#!/usr/bin/env python3
#
"""
https://www.rdkit.org/docs/source/rdkit.Chem.rdChemReactions.html
"""
import sys,os,re,argparse,logging

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from .. import util
from .. import reactions

#############################################################################
if __name__=='__main__':
  EPILOG = "Reactants specified as disconnected components of single molecule, or from separate input files."
  parser = argparse.ArgumentParser(description="RDKit chemical reactions utility", epilog=EPILOG)
  OPS = [ "enumerateLibrary", "react", "demo", "demo2", "demo3"]
  parser.add_argument("op", choices=OPS, help="OPERATION")
  parser.add_argument("--i", dest="ifiles", help="input file[s] (SMILES/TSV or SDF)")
  parser.add_argument("--o", dest="ofile", default="-", help="output file (specify '-' for stdout)")
  parser.add_argument("--smirks", help="SMIRKS reaction transform")
  parser.add_argument("--kekulize", action="store_true", help="Kekulize")
  parser.add_argument("--sanitize", action="store_true", help="Sanitize")
  parser.add_argument("--header", action="store_true", help="input SMILES/TSV file has header line")
  parser.add_argument("--delim", default="\t", help="delimiter for SMILES/TSV")
  parser.add_argument("--smilesColumn", type=int, default=0, help="input SMILES column")
  parser.add_argument("--nameColumn", type=int, default=1, help="input name column")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  molReaders=[];
  if args.ifiles is not None:
    for ifile in re.split(r'[\s*,\s*]', args.ifiles):
      molReader = util.File2Molreader(ifile, args.delim, args.smilesColumn, args.nameColumn, args.header)
      molReaders.append(molReader)
  molWriter = util.File2Molwriter(args.ofile, args.delim, args.header)


  if args.op == "react":
    reactions.React(args.smirks, molReaders, molWriter)

  elif args.op == "enumerateLibrary":
    reactions.EnumerateLibrary(args.smirks, molReaders, molWriter)

  elif args.op == "demo": reactions.Demo()
  elif args.op == "demo2": reactions.Demo2()
  elif args.op == "demo3": reactions.Demo3()
  else:
    parser.error(f"Unsupported operation: {args.op}")

