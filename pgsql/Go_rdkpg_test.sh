#!/bin/sh
#
# https://rdkit.readthedocs.org/en/latest/Cartridge.html#reference-guide
#
DB="scratch"
DBSCHEMA="pdb"
#
#
set -x
#
psql -d $DB -c "SELECT id,mol FROM ${DBSCHEMA}.mols WHERE mol@>'c1cccc2c1nncc2' LIMIT 10"
psql -d $DB -c "SELECT id,mol FROM ${DBSCHEMA}.mols WHERE mol@>'c1[o,s]ncn1'::qmol LIMIT 10"
psql -d $DB -c "SELECT * FROM pdb_simsearch_mfp('c1ccc(-c2nsc(N3CCNCC3)n2)cc1')"
#
psql -d $DB <<__EOF__
SET rdkit.tanimoto_threshold=0.0;
SELECT
	id,
	mol,
	TO_CHAR(similarity,'0.99') AS sim_tanimoto
FROM
	pdb_simsearch_mfp('c1ccc(-c2nsc(N3CCNCC3)n2)cc1')
WHERE
	similarity > 0.4
	;
__EOF__
#
