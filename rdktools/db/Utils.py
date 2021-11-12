#!/usr/bin/env python3
#############################################################################
import os,sys,re,time,logging
import pandas as pd
import psycopg2,psycopg2.extras

import rdkit
from rdkit import Chem

from ..util import db as util_db

#############################################################################
def TestRDKit(tables, dbcon, dbschema=util_db.DBSCHEMA, fout=None):
  """Test for RDKit extension, and objects in specified tables."""
  df = pd.read_sql(f"SELECT * FROM pg_extension WHERE extname='rdkit'", dbcon)
  logging.info("RDKit extension {0}FOUND.".format( "NOT " if df.shape[0]==0 else ""))
  if fout: df.to_csv(fout, "\t", index=False)

  sql=f"""\
SELECT
	n.nspname||'.'||p.proname as "Function",
	pg_catalog.pg_get_function_arguments(p.oid) as "Argument datatypes",
	pg_catalog.pg_get_function_result(p.oid) as "Result datatype"
FROM pg_catalog.pg_proc p
LEFT JOIN pg_catalog.pg_namespace n ON n.oid = p.pronamespace
WHERE pg_catalog.pg_function_is_visible(p.oid) AND n.nspname = 'public'
ORDER BY 1
"""
  df = pd.read_sql(sql, dbcon)
  if fout: df.to_csv(fout, "\t", index=False)
  return df

#############################################################################
def MolSubstructureSearch(smiles, moltable, dbcon, dbschema=util_db.DBSCHEMA, fout=None):
  sql=f"""\
SELECT
	cansmi,
	name
FROM
	{dbschema}.{moltable}
WHERE
	mol@>'{smiles}'
ORDER BY
	LENGTH(cansmi) ASC
"""
  logging.debug(sql)
  df = pd.read_sql(sql, dbcon)
  if fout: df.to_csv(fout, "\t", index=False)
  logging.info(f"n_out: {df.shape[0]}")
  return df

#############################################################################
def MolSimilaritySearch(smiles, moltable, dbcon, dbschema=util_db.DBSCHEMA, fout=None):
  sql=f"""\
SELECT
	cansmi,
	name,
	tanimoto_sml(rdkit_fp(mol_from_smiles('{smiles}'::cstring)), {moltable}.fp)::NUMERIC AS tanimoto
FROM
	{dbschema}.{moltable}
ORDER BY
	tanimoto DESC
"""
  logging.debug(sql)
  df = pd.read_sql(sql, dbcon)
  if fout: df.to_csv(fout, "\t", index=False)
  logging.info(f"n_out: {df.shape[0]}")
  return df
