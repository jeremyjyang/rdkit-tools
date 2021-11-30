#!/usr/bin/env python3
"""
DrugCentral PostgreSql db client.
"""
import os,sys,argparse,re,time,logging
import pandas as pd

from .. import db
from .. import util
from ..util import yaml as util_yaml
from ..util import db as util_db

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description="RDKit-tools Pg db client utility", epilog="")
  ops = [
	"list_tables",
	"list_tables_rowCounts",
	"list_columns",
	"test_rdkit",
	"show_params",
	"molSubstructureSearch",
	"molSimilaritySearch",
	]
  parser.add_argument("op", choices=ops, help="OPERATION (select one)")
  parser.add_argument("--i", dest="ifile", help="input ID file")
  parser.add_argument("--ids", help="input IDs (comma-separated)")
  parser.add_argument("--o", dest="ofile", help="output (TSV)")
  parser.add_argument("--moltable", default="mols", help="table with RDKit molecules")
  parser.add_argument("--smiles", help="SMILES query")
  parser.add_argument("--dbhost")
  parser.add_argument("--dbport")
  parser.add_argument("--dbname")
  parser.add_argument("--dbusr")
  parser.add_argument("--dbpw")
  parser.add_argument("--param_file", default=os.environ['HOME']+"/.rdktools.yaml")
  parser.add_argument("--dbschema", default="public")
  parser.add_argument("-v", "--verbose", dest="verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  if os.path.isfile(args.param_file):
    params = util_yaml.ReadParamFile(args.param_file)
  if args.dbhost: params['DBHOST'] = args.dbhost 
  if args.dbport: params['DBPORT'] = args.dbport 
  if args.dbname: params['DBNAME'] = args.dbname 
  if args.dbusr: params['DBUSR'] = args.dbusr 
  if args.dbpw: params['DBPW'] = args.dbpw 

  if args.op=='show_params':
    pd.DataFrame({key:[params[key]] for key in params.keys()}).transpose().to_csv(sys.stdout, "\t", index=True, header=False)
    sys.exit()

  fout = open(args.ofile, "w+") if args.ofile else sys.stdout

  if args.ifile:
    fin = open(args.ifile)
    while True:
      line = fin.readline()
      if not line: break
      ids.append(line.rstrip())
    logging.info(f"Input IDs: {len(ids)}")
    fin.close()
  elif args.ids:
    ids = re.split(r'[,\s]+', args.ids)

  t0 = time.time()

  dbcon = util_db.Connect(params['DBHOST'], params['DBPORT'], params['DBNAME'], params['DBUSR'], params['DBPW'])
  if dbcon is None: sys.exit()

  if args.op=='list_tables':
    util_db.ListTables(dbcon, args.dbschema, fout)

  elif args.op=='list_tables_rowCounts':
    util_db.ListTablesRowCounts(dbcon, args.dbschema, fout)

  elif args.op=='list_columns':
    util_db.ListColumns(dbcon, args.dbschema, fout)

  elif args.op=='test_rdkit':
    db.TestRDKit(args.moltable, dbcon, args.dbschema, fout)

  elif args.op=='molSubstructureSearch':
    db.MolSubstructureSearch(args.smiles, args.moltable, dbcon, args.dbschema, fout)

  elif args.op=='molSimilaritySearch':
    db.MolSimilaritySearch(args.smiles, args.moltable, dbcon, args.dbschema, fout)

  else:
    parser.error(f"Invalid operation: {args.op}")

  dbcon.close()
  logging.info('Elapsed time: %s'%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))))


