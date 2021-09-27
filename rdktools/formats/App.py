#!/usr/bin/env python3
#
import sys,os,re,argparse,logging

import rdkit
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolToSmiles, MolFromSmiles
from rdkit.Chem import AllChem

from .. import formats

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='RDKit chemical format utility')
  OPS = ["mdl2smi", "mdl2tsv", "smi2mdl", "smiclean", "mdlclean", "mol2inchi", "mol2inchikey", "demo"]
  parser.add_argument("op", choices=OPS, help="operation")
  parser.add_argument("--i", dest="ifile", help="input file (SMILES/TSV or SDF)")
  parser.add_argument("--o", dest="ofile", help="output file (specify '-' for stdout)")
  parser.add_argument("--kekulize", action="store_true", help="Kekulize")
  parser.add_argument("--sanitize", action="store_true", help="Sanitize")
  parser.add_argument("--header", action="store_true", help="input SMILES/TSV file has header line")
  parser.add_argument("--delim", default="\t", help="delimiter for SMILES/TSV")
  parser.add_argument("--smicol", type=int, default=0, help="input SMILES column")
  parser.add_argument("--namcol", type=int, default=1, help="input name column")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  if args.op in ("mdl2smi", "mdlclean"):
    molReader = SDMolSupplier(args.ifile, sanitize=True, removeHs=False)
  elif args.op in ("smi2mdl", "smiclean"):
    molReader = SmilesMolSupplier(args.ifile, delimiter=args.delim, smilesColumn=args.smicol, nameColumn=args.namcol, titleLine=True, sanitize=True)
  else:
    if re.sub(r'.*\.', '', args.ifile).lower()in ('smi', 'smiles'):
      molReader = SmilesMolSupplier(args.ifile, delimiter=args.delim, smilesColumn=args.smicol, nameColumn=args.namcol, titleLine=True, sanitize=True)
    elif re.sub(r'.*\.', '', args.ifile).lower() in ('sdf','sd','mdl','mol'):
      molReader = SDMolSupplier(args.ifile, sanitize=True, removeHs=True)
    else:
      logging.error(f'Invalid file extension: {args.ifile}')

  if args.op in ("mdl2smi", "smiclean"):
    molWriter = SmilesWriter(args.ofile if args.ofile else "-", delimiter=args.delim, includeHeader=False, isomericSmiles=True, kekuleSmiles=True)
  elif args.op in ("mdl2tsv", ):
    molWriter = SmilesWriter(args.ofile if args.ofile else "-", delimiter=args.delim, nameHeader='Name', includeHeader=True, isomericSmiles=True, kekuleSmiles=True)
  elif args.op in ("smi2mdl", "mdlclean"):
    molWriter = SDWriter(args.ofile if args.ofile else "-")
  elif args.op in ("mol2inchikey", "mol2inchi"):
    fout = open(args.ofile, "w+") if args.ofile else sys.stdout
  else:
    molWriter = SmilesWriter("-", delimiter='\t', nameHeader='Name', includeHeader=True, isomericSmiles=True, kekuleSmiles=True)

  if args.op == "mdl2smi":
    formats.Utils.Mdl2Smi(molReader, molWriter)

  elif args.op == "mdl2tsv":
    formats.Utils.Mdl2Tsv(molReader, molWriter)

  elif args.op == "smi2mdl":
    formats.Utils.Smi2Mdl(molReader, molWriter)

  elif args.op == "mol2inchi":
    formats.Utils.Mol2Inchi(molReader, fout)

  elif args.op == "mol2inchikey":
    formats.Utils.Mol2Inchikey(molReader, fout)

  else:
    parser.error(f"Unsupported operation: {args.op}")


