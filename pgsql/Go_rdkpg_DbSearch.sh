#!/bin/bash
###
# https://rdkit.readthedocs.org/en/latest/Cartridge.html#reference-guide
# This script inputs a SMILES file and for each SMILES, runs a similarity
# search on a PostgreSQL Db with the RDKit Cartridge, as built with
# the accompanying script Go_rdkpg_DbCreate.sh.
###
#
set -e
#
function SimilaritySearch {
	DBNAME=$1
	QSMILES=$2
	OUTPUT_FILE=$3
	#
	(psql -qAF $'\t' -d $DBNAME <<__EOF__
SET rdkit.morgan_fp_size=2048;
SET rdkit.tanimoto_threshold=0.0;
SELECT
	name,
	smiles,
	ROUND(tanimoto_sml(morganbv_fp(mol_from_smiles('${QSMILES}'), 2), ecfp)::NUMERIC, 3) AS similarity
FROM
	mols
WHERE
	morganbv_fp(mol_from_smiles('${QSMILES}'), 2)%ecfp
ORDER BY
	morganbv_fp(mol_from_smiles('${QSMILES}'), 2)<%>ecfp,
	LENGTH(smiles)
LIMIT 1000
	;
__EOF__
		) |sed -e '$d' \
	>$OUTPUT_FILE
}
#
T0=$(date +%s)
#
#
if [ $# -ne 2 ]; then
	printf "Syntax: %s DBNAME SMIFILE (Format: SMILES<tab>NAME)\n" $0
	exit
fi
#
PROG=$(basename $0 |sed -e 's/\..*$//')
DBNAME=$1
SMIFILE=$2
N=$(cat $SMIFILE |wc -l)
#
I="0"
while [ "$I" -lt "$N" ]; do
	I=$[$I + 1]
	smi=$(cat $SMIFILE |sed "${I}q;d" |awk -F '\t' '{print $1}')
	name=$(cat $SMIFILE |sed "${I}q;d" |awk -F '\t' '{print $2}')
	if [ ! "${smi}" ]; then
		printf "ERROR: Failed to parse SMILES from line %d\n" "$I"
		continue
	fi
	if [ ! "${name}" ]; then
		printf "ERROR: Failed to parse NAME from line %d\n" "$I"
		continue
	fi
	ofile="${PROG}_output_${name}.tsv"
	#
	printf "${I}.  NAME:${name}; SMILES: ${smi}; OFILE:${ofile}\n"
	SimilaritySearch "${DBNAME}" "${smi}" "${ofile}"
done
#
printf "Elapsed time: %ds\n" "$[$(date +%s) - ${T0}]"
#
