#!/usr/bin/env python3
"""
Generic db utilities, mostly focused on PostgreSql, for use with
the RDKit chemical cartridge.
"""
#############################################################################
import os,sys,re,logging,yaml
import pandas as pd

import psycopg2,psycopg2.extras

DBSCHEMA="public"

#############################################################################
def Connect(dbhost, dbport, dbname, dbusr, dbpw):
  """Connect to db; specify default cursor type DictCursor."""
  dsn = (f"host='{dbhost}' port='{dbport}' dbname='{dbname}' user='{dbusr}' password='{dbpw}'")
  dbcon = psycopg2.connect(dsn)
  dbcon.cursor_factory = psycopg2.extras.DictCursor
  return dbcon

#############################################################################
def ListTables(dbcon, dbschema=DBSCHEMA, fout=None):
  '''Listing the tables.'''
  sql = (f"SELECT table_name FROM information_schema.tables WHERE table_schema = '{dbschema}'")
  df = pd.read_sql(sql, dbcon)
  if fout: df.to_csv(fout, "\t", index=False)
  logging.info(f"n_out: {df.shape[0]}")
  return df

#############################################################################
def ListTablesRowCounts(dbcon, dbschema=DBSCHEMA, fout=None):
  '''Listing the table rowcounts.'''
  df=None;
  sql1 = (f"SELECT table_name FROM information_schema.tables WHERE table_schema = '{dbschema}'")
  df1 = pd.read_sql(sql1, dbcon)
  for tname in df1.table_name:
    sql2 = (f"SELECT COUNT(*) AS rowcount FROM {dbschema}.{tname}")
    df_this = pd.read_sql(sql2, dbcon)
    df_this["schema"] = dbschema
    df_this["table"] = tname
    df = df_this if df is None else pd.concat([df, df_this])
  df = df[["schema", "table", "rowcount"]]
  if fout: df.to_csv(fout, "\t", index=False)
  logging.info(f"n_out: {df.shape[0]}")
  return df

#############################################################################
def ListColumns(tables, dbcon, dbschema=DBSCHEMA, fout=None):
  df=None;
  for tname in tables:
    sql = (f"SELECT column_name,data_type FROM information_schema.columns WHERE table_schema = '{dbschema}' AND table_name = '{tname}'")
    df_this = pd.read_sql(sql, dbcon)
    df_this["schema"] = dbschema
    df_this["table"] = tname
    df = df_this if df is None else pd.concat([df, df_this])
  df = df[["schema", "table", "column_name", "data_type"]]
  if fout: df.to_csv(fout, "\t", index=False)
  logging.info(f"n_out: {df.shape[0]}")
  return df

