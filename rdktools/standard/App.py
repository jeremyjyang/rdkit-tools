#!/usr/bin/env python3
#############################################################################
import os,sys,re,time,argparse,logging

import rdkit
from rdkit.Chem import MolStandardize

from .. import standard as standard_utils
from .. import util

#############################################################################
def GetNorms(args):
  if args.normset=="UNM":
    norms = standard_utils.MyNorms()
  elif args.ifile_normset:
    fin = open(args.ifile_normset)
    norms = standard_utils.ReadNormsFile(fin)
  else: #DEFAULT
    norms = MolStandardize.normalize.NORMALIZATIONS
  return norms

#############################################################################
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="RDKit chemical standardizer", epilog="")
  NORMSETS = ["DEFAULT", "UNM"]
  OPS = ["standardize", "canonicalize", "list_norms", "show_params", "demo"]
  parser.add_argument("op", choices=OPS, help="OPERATION")
  parser.add_argument("--i", dest="ifile", required=True, help="input file, SMI or SDF")
  parser.add_argument("--o", dest="ofile", default="-", help="output file, SMI or SDF")
  parser.add_argument("--delim", default="\t", help="SMILES/TSV delimiter")
  parser.add_argument("--smilesColumn", type=int, default=0, help="")
  parser.add_argument("--nameColumn", type=int, default=1, help="")
  parser.add_argument("--header", action="store_true", help="SMILES/TSV has header line")
  parser.add_argument("--sanitize", action="store_true", help="Sanitize molecules as read.")
  parser.add_argument("--kekuleSmiles", action="store_true", help="Kekule SMILES output.")
  parser.add_argument("--normset", choices=NORMSETS, default="DEFAULT", help="normalization sets")
  parser.add_argument("--i_normset", dest="ifile_normset", help="input normalizations file, format: SMIRKS<space>NAME")
  parser.add_argument("--isomericSmiles", action="store_true", help="If false, output SMILES isomerism removed")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f"RDKit version: {rdkit.__version__}")

  t0=time.time()

  if args.op=="show_params":
    standard_utils.ShowParameters()

  elif args.op=="list_norms":
    norms = GetNorms(args)
    fout = open(args.ofile, "w") if args.ofile else sys.stdout
    standard_utils.ListNormalizations(norms, fout)

  elif args.op=="standardize":
    molReader = util.File2Molreader(args.ifile, args.delim, args.smilesColumn, args.nameColumn, args.header, False) #Sanitize separately.
    molWriter = util.File2Molwriter(args.ofile, args.delim, args.header, isomericSmiles=args.isomericSmiles, kekuleSmiles=args.kekuleSmiles)
    norms = GetNorms(args)
    stdzr = standard_utils.MyStandardizer(norms)
    standard_utils.Standardize(stdzr, args.sanitize, args.isomericSmiles, molReader, molWriter)

  elif args.op=="canonicalize":
    molReader = util.File2Molreader(args.ifile, args.delim, args.smilesColumn, args.nameColumn, args.header, False) #Sanitize separately.
    molWriter = util.File2Molwriter(args.ofile, args.delim, args.header, isomericSmiles=args.isomericSmiles, kekuleSmiles=args.kekuleSmiles)
    standard_utils.Canonicalize(args.isomericSmiles, args.kekuleSmiles, molReader, molWriter)

  elif args.op=="demo":
    standard_utils.Demo()

  else:
    parser.error(f"Unsupported operation: {args.op}")

  logging.info('Elapsed: %s'%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))))

