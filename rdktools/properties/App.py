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

from .. import util
from .. import properties

#############################################################################
if __name__=='__main__':
  EPILOG = """estate: Kier-Hall electrotopological descriptors;
lipinski: rdkit.Chem.Lipinski descriptors;
logp: Wildman-Crippen LogP;
descriptors: Chem.Descriptors (various); descriptors3d: Chem.Descriptors3D (various) (requires 3D coordinates);
freesasa: Accessible Surface Area (SASA) and related descriptors (requires 3D coordinates)
"""
  parser = argparse.ArgumentParser(description='RDKit molecular properties utility')
  OPS = [ "descriptors", "descriptors3d", "lipinski", "logp", "estate", "freesasa", "demo"]
  parser.add_argument("op", choices=OPS, help="OPERATION")
  parser.add_argument("--i", dest="ifile", help="input molecule file")
  parser.add_argument("--o", dest="ofile", help="output file with data (TSV)")
  parser.add_argument("--iheader", action="store_true", help="input file has header line")
  parser.add_argument("--oheader", action="store_true", help="include TSV header line with smiles output")
  parser.add_argument("--kekulize", action="store_true", help="Kekulize")
  parser.add_argument("--sanitize", action="store_true", help="Sanitize")
  parser.add_argument("--delim", default=" \t", help="SMILES/TSV delimiter")
  parser.add_argument("--smilesColumn", type=int, default=0, help="input SMILES column")
  parser.add_argument("--nameColumn", type=int, default=1, help="input name column")

  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  molReader = util.File2Molreader(args.ifile, args.delim, args.smilesColumn, args.nameColumn, args.iheader) if args.ifile is not None else None

  molWriter = util.File2Molwriter(args.ofile, args.delim, args.oheader) if args.ofile is not None else Chem.SDWriter("-")

  if args.op == "descriptors":
    properties.Run_CalcDescriptors(molReader, molWriter)

  elif args.op == "descriptors3d":
    properties.Run_CalcDescriptors3D(molReader, molWriter)

  elif args.op == "estate":
    properties.Run_CalcEStateIndices(molReader, molWriter)

  elif args.op == "logp":
    properties.Run_CalcCrippenLogP(molReader, molWriter)

  elif args.op == "lipinski":
    properties.Run_CalcLipinski(molReader, molWriter)

  elif args.op == "freesasa":
    properties.Run_CalcFreeSASA(molReader, molWriter)

  elif args.op == "demo":
    properties.Demo()

  else:
    logging.error(f"Invalid operation: {args.op}")

