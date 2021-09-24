#!/usr/bin/env python3
"""
Calculate Kier-Hall electrotopological descriptors
  [x] smiles or molfile input and output with data
  [ ] stdout and stdin
  [ ] sumDeltaI -- can it be done with RDKit?
 
References:
  1. L.H. Hall, B. Mohney and L.B. Kier., "The electrotopological state:
     structure information at the atomic level for molecular graphs",
     JCICS 31 (1) 76-82 (1991)
  2. L.B. Kier and L.H. Hall _Molecular Structure Description:
     The Electrotopological State"_  Academic Press (1999)
"""
###       
import sys,os,argparse,logging

from rdkit import Chem
from rdkit.Chem import AllChem

from .. import properties

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='RDKit molecular properties utility')
  OPS = ["logp", "estate", "demo"]
  parser.add_argument("op", choices=OPS, help="OPERATION")
  parser.add_argument("--i", required=True, dest="ifile", help="input molecule file")
  parser.add_argument("--o", dest="ofile", help="output file with data (TSV)")
  parser.add_argument("--inheader", action="store_true", help="input file has header line")
  parser.add_argument("--outheader", action="store_true", help="include TSV header line with smiles output")
  parser.add_argument("--kekulize", action="store_true", help="Kekulize")
  parser.add_argument("--sanitize", action="store_true", help="Sanitize")
  parser.add_argument("--delim", default="\t", help="input molecule file")
  parser.add_argument("--smicol", type=int, default=0, help="input SMILES column")
  parser.add_argument("--namcol", type=int, default=1, help="input name column")

  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  if args.ifile[-4:]=='.smi':
    molreader = Chem.SmilesMolSupplier(args.ifile, delimiter=args.delim, smilesColumn=args.smicol, nameColumn=args.namcol, titleLine=args.inheader, sanitize=args.sanitize)
  elif args.ifile[-4:] in ('.sdf','.sd','.mdl','.mol'):
    molreader = Chem.SDMolSupplier(args.ifile, sanitize=args.sanitize, removeHs=True)
  else:
    logging.error(f"Invalid file extension: {args.ifile}")

  if args.ofile is None:
    molwriter = Chem.SDWriter("-")
  elif args.ofile[-4:]=='.smi':
    molwriter = Chem.SmilesWriter(args.ofile, delimiter='\t', nameHeader='Name', includeHeader=args.outheader, isomericSmiles=True, kekuleSmiles=args.kekulize)
  elif args.ofile[-4:] in ('.sdf','.sd','.mdl','.mol'):
    molwriter = Chem.SDWriter(args.ofile)
  else:
    logging.error(f"Invalid file extension: {args.ofile}")

  if args.op == "estate":
    i_mol=0
    for mol in molreader:
      i_mol+=1
      properties.Utils.CalcEStateIndices(mol)
      if i_mol==1:
        molwriter.SetProps(mol.GetPropNames())
      molwriter.write(mol)
    logging.info(f"n_out: {i_mol}")

  elif args.op == "logp":
    i_mol=0
    for mol in molreader:
      i_mol+=1
      name = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
      logging.info(f"{i_mol}. {name}")
      #AllChem.Compute2DCoords(mol)
      #Chem.SanitizeMol(mol)
      properties.Utils.CalcCrippenLogP(mol)
      if i_mol==1:
        molwriter.SetProps(mol.GetPropNames())
      molwriter.write(mol)
    logging.info(f"n_out: {i_mol}")

  else:
    logging.error(f"Invalid operation: {args.op}")


