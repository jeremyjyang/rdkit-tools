#!/usr/bin/env python3
#############################################################################
import os,sys,re,time,argparse,logging

import rdkit
#import rdkit.Chem.AllChem
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles

from .. import standard

#############################################################################
def GetNorms(args):
  if args.norms=="unm":
    norms = standard.MyNorms()
  elif args.ifile_norms:
    fin = open(args.ifile_norms)
    norms = standard.ReadNormsFile(fin)
  else:
    norms = MolStandardize.normalize.NORMALIZATIONS
  return norms

#############################################################################
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="RDKit chemical standardizer", epilog="")
  OPS = ["standardize", "list_norms", "show_params", "demo"]
  parser.add_argument("op", choices=OPS, help="OPERATION")
  parser.add_argument("--i", dest="ifile", help="input file, SMI or SDF")
  parser.add_argument("--o", dest="ofile", help="output file, SMI or SDF")
  parser.add_argument("--delim", default=" \t", help="SMILES/TSV delimiter")
  parser.add_argument("--smilesColumn", type=int, default=0, help="")
  parser.add_argument("--nameColumn", type=int, default=1, help="")
  parser.add_argument("--header", action="store_true", help="SMILES/TSV has header line")
  parser.add_argument("--norms", choices=["default", "unm"], default="default", help="normalizations")
  parser.add_argument("--i_norms", dest="ifile_norms", help="input normalizations file, format: SMIRKS<space>NAME")
  parser.add_argument("--remove_isomerism", action="store_true", help="if true, output SMILES isomerism removed")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f"RDKit version: {rdkit.__version__}")

  t0=time.time()

  if args.op=="show_params":
    standard.ShowParameters()

  elif args.op=="list_norms":
    norms = GetNorms(args)
    fout = open(args.ofile, "w") if args.ofile else sys.stdout
    standard.ListNormalizations(norms, fout)

  elif args.op=="standardize":
    if not (args.ifile and args.ofile): parser.error('--i and --o required.')
    molReader = util.File2Molreader(ifile, args.delim, args.smilesColumn, args.nameColumn, args.header)
    molWriter = util.File2Molwriter(args.ofile, args.delim, args.header, isomericSmiles=(not args.remove_isomerism), kekuleSmiles=True)
    norms = GetNorms(args)
    stdzr = standard.MyStandardizer(norms)
    standard.Standardize(stdzr, args.remove_isomerism, molReader, molWriter)

  elif args.op=="demo":
    standard.Demo()

  else:
    parser.error(f"Unsupported operation: {args.op}")

  logging.info('Elapsed: %s'%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))))

