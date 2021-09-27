#!/usr/bin/env python3
#
import sys,os,argparse,logging

import rdkit
from rdkit.Chem import AllChem

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='RDKit chemical format utility')
  parser.add_argument("--i", dest="ifile", help="input file (SMILES/TSV)")
  parser.add_argument("--kekulize", action="store_true", help="Kekulize")
  parser.add_argument("--sanitize", action="store_true", help="Sanitize")
  parser.add_argument("--header", action="store_true", help="input file has header line")
  parser.add_argument("--o", dest="ofile", help="output file")
  parser.add_argument("--delim", default="\t", help="input molecule file")
  parser.add_argument("--smicol", type=int, default=0, help="input SMILES column")
  parser.add_argument("--namcol", type=int, default=1, help="input name column")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  molsuppl = rdkit.Chem.SmilesMolSupplier(args.ifile, delimiter=args.delim, smilesColumn=args.smicol, nameColumn=args.namcol, titleLine=args.header, sanitize=args.sanitize)

  sdwriter = rdkit.Chem.SDWriter(args.ofile if args.ofile else "-")

  i_mol=0
  for mol in molsuppl:
    i_mol+=1
    logging.info(f"{i_mol}. {mol.GetProp('_Name')}") 
    rdkit.Chem.SanitizeMol(mol)
    AllChem.Compute2DCoords(mol)
    sdwriter.write(mol)

  logging.info(f"n_out: {i_mol}")
