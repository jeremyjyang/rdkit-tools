#!/bin/sh
#############################################################################
# See:
#	$RDBASE/Code/PgSQL/rdkit/README
#	$RDBASE/Code/PgSQL/rdkit/sql/
#	https://rdkit.readthedocs.org/en/latest/Cartridge.html#reference-guide
#	https://github.com/rdkit/rdkit/issues/1762
#############################################################################
# Old way; build from source:
# cd $RDBASE/Code/PgSQL/rdkit
# make
# make install
# make installcheck
#############################################################################
# New way; Ubuntu 20.04
# apt install postgresql-12
# apt-cache policy postgresql-12
# apt install postgresql-server-dev-12
# apt install postgresql-12-rdkit
#############################################################################
#
DBNAME="scratch"
DBSCHEMA="pdb"
#
sudo -u postgres psql -d $DBNAME -c 'CREATE EXTENSION rdkit'
#
psql -d $DBNAME <<__EOF__
SELECT
	*
	INTO ${DBSCHEMA}.mols
FROM
	(SELECT
		structureid AS id,
		mol_from_smiles(smiles::cstring) AS mol
	FROM
		${DBSCHEMA}.pdb_ligands_druglike
	) tmp
WHERE
	mol IS NOT NULL
	;
--
__EOF__
#
psql -d $DBNAME -c "CREATE INDEX molidx ON ${DBSCHEMA}.mols USING gist(mol)"
#
#
### One way is to create separate table for FPs.
#psql -d $DBNAME -c "SELECT id,torsionbv_fp(mol) AS torsionbv,morganbv_fp(mol) AS mfp,featmorganbv_fp(mol) AS ffp INTO ${DBSCHEMA}.fps FROM ${DBSCHEMA}.mols"
#
### Another way is to add FPs to mols table.
psql -d $DBNAME -c "ALTER TABLE ${DBSCHEMA}.mols ADD COLUMN torsionbv BFP"
psql -d $DBNAME -c "ALTER TABLE ${DBSCHEMA}.mols ADD COLUMN mfp BFP"
psql -d $DBNAME -c "ALTER TABLE ${DBSCHEMA}.mols ADD COLUMN ffp BFP"
psql -d $DBNAME -c "UPDATE ${DBSCHEMA}.mols SET torsionbv = torsionbv_fp(mol)"
psql -d $DBNAME -c "UPDATE ${DBSCHEMA}.mols SET mfp = morganbv_fp(mol)"
psql -d $DBNAME -c "UPDATE ${DBSCHEMA}.mols SET ffp = featmorganbv_fp(mol)"
#
psql -d $DBNAME -c "CREATE INDEX fps_ttbv_idx ON ${DBSCHEMA}.mols USING gist(torsionbv)"
psql -d $DBNAME -c "CREATE INDEX fps_mfp_idx ON ${DBSCHEMA}.mols USING gist(mfp)"
psql -d $DBNAME -c "CREATE INDEX fps_ffp_idx ON ${DBSCHEMA}.mols USING gist(ffp)"
#
#
### Convenience function:
sudo -u postgres psql -d $DBNAME <<__EOF__
CREATE OR REPLACE FUNCTION
	${DBSCHEMA}_simsearch_mfp(smiles text)
RETURNS TABLE(id VARCHAR, mol mol, similarity double precision) AS
	$$
	SELECT
		id,mol,tanimoto_sml(morganbv_fp(mol_from_smiles($1::cstring)),mfp) AS similarity
	FROM
		${DBSCHEMA}.mols
	WHERE
		morganbv_fp(mol_from_smiles($1::cstring))%mfp
	ORDER BY
		morganbv_fp(mol_from_smiles($1::cstring))<%>mfp
		;
	$$
LANGUAGE SQL STABLE
	;
__EOF__
#
#
