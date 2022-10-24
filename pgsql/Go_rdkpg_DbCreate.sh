#!/bin/bash
###
# https://www.rdkit.org/docs/Cartridge.html
###
#
T0=$(date +%s)
#
DBNAME="REPLACE_WITH_DB_NAME"
DBSCHEMA="public"
DBHOST="localhost"
#
cwd=$(pwd)
#
dropdb $DBNAME
createdb $DBNAME
psql -d $DBNAME -c "COMMENT ON DATABASE $DBNAME IS '${REPLACE_WITH_DB_NAME} Db'";
###
# Create mols table for RDKit structural searching.
#sudo -u postgres psql -d $DBNAME -c 'CREATE EXTENSION rdkit'
psql -d $DBNAME -c 'CREATE EXTENSION rdkit'
#
# Molecules table:
psql -d $DBNAME <<__EOF__
CREATE TABLE mols (
	id SERIAL PRIMARY KEY,
	name VARCHAR(100),
	smiles VARCHAR(2000),
	cansmi VARCHAR(2000),
	molecule MOL
	);
__EOF__
#
# Metadata table:
psql -d $DBNAME <<__EOF__
CREATE TABLE meta (
	field VARCHAR(18) PRIMARY KEY,
	value VARCHAR(100),
	description VARCHAR(2000)
	);
__EOF__
###
# Load mols table from TSV file.
cat $DBDIR/REPLACE_WITH_FILE_NAME.tsv \
	|awk -F '\t' '{print INSERT INTO mols (id, name, smiles) VALUES (\'$1\', \'$2\', \'$3\');' \
	|psql -q $DBNAME
###
#
psql -d $DBNAME -c "UPDATE mols SET mol = mol_from_smiles(smiles::cstring)"
psql -d $DBNAME -c "UPDATE mols SET cansmi = mol_to_smiles(mol)"
psql -d $DBNAME -c "CREATE INDEX molidx ON mols USING gist(molecule)"
#
### Add FPs to mols table.
# https://www.rdkit.org/docs/GettingStartedInPython.html
# Path-based, Daylight-like.
psql -d $DBNAME -c "SET rdkit.rdkit_fp_size=2048"
psql -d $DBNAME -c "ALTER TABLE mols DROP COLUMN IF EXISTS fp"
psql -d $DBNAME -c "ALTER TABLE mols ADD COLUMN fp BFP"
psql -d $DBNAME -c "UPDATE mols SET fp = rdkit_fp(molecule)"
psql -d $DBNAME -c "CREATE INDEX fps_fp_idx ON mols USING gist(fp)"
#
# Morgan (Circular) Fingerprints (with radius=2 ECFP4-like).
psql -d $DBNAME -c "SET rdkit.morgan_fp_size=2048"
psql -d $DBNAME -c "ALTER TABLE mols DROP COLUMN IF EXISTS ecfp"
psql -d $DBNAME -c "ALTER TABLE mols ADD COLUMN ecfp BFP"
psql -d $DBNAME -c "UPDATE mols SET ecfp = morganbv_fp(molecule, 2)"
psql -d $DBNAME -c "CREATE INDEX fps_ecfp_idx ON mols USING gist(ecfp)"
#
###
# Updated by postprocess script.
psql -d $DBNAME -c "INSERT INTO meta (field, value, description) VALUES ('timestamp', CURRENT_DATE::VARCHAR, 'Date of build.')"
###
#
printf "Elapsed time: %ds\n" "$[$(date +%s) - ${T0}]"
#
