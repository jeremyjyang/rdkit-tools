#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/
https://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints

"""
###
## RDKFingerprint # returns ExplicitBitVect
## args: minPath=1,maxPath=7,fpSize=2048,nBitsPerHash=4,useHs=True, tgtDensity=0.0,minSize=128
##
## GenMACCSKeys (no args) # returns SparseBitVect
##
## AllChem.GetMorganFingerprint(mol,2) # returns SparseBitVect
## AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
##
## When comparing the ECFP/FCFP fingerprints and Morgan fingerprints by
## RDKit, remember that 4 in ECFP4 corresponds to diameter of atom environments,
## while Morgan fingerprints take a radius parameter. So radius=2 roughly
## equivalent to ECFP4 and FCFP4.
#rdkit.DataStructs.BitVectToText
#rdkit.DataStructs.BitVectToBinaryText
#rdkit.DataStructs.BitVectToFPSText
#rdkit.DataStructs.CreateFromBinaryText
#rdkit.DataStructs.CreateFromBitString
#rdkit.DataStructs.CreateFromFPSText
#############################################################################
import os,sys,re,json,time,inspect,argparse,logging
import matplotlib as mpl
#from matplotlib import pyplot as plt

import rdkit
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles, Draw, rdDepictor

from .. import fp
from .. import util

#############################################################################
def Demo():
  mols=[];
  for smi in util.Utils.DEMOSMIS:
    mol = MolFromSmiles(re.sub(r'\s.*$', '', smi))
    mols.append(mol)
  molWriter = SmilesWriter("-", delimiter="\t", nameHeader='Name', includeHeader=True, isomericSmiles=True, kekuleSmiles=False)
  fp.Utils.Mols2FPs_RDK(mols, molWriter)
 
#############################################################################
def DemoMorgan():
  mols=[];
  for smi in util.Utils.DEMOSMIS:
    mol = MolFromSmiles(re.sub(r'\s.*$', '', smi))
    mols.append(mol)
  molWriter = SmilesWriter("-", delimiter="\t", nameHeader='Name', includeHeader=True, isomericSmiles=True, kekuleSmiles=False)
  fp.Utils.Mols2FPs_Morgan(mols, 2, 1024, molWriter)
 
#############################################################################
def DemoMACCSKeys():
  mols=[];
  for smi in util.Utils.DEMOSMIS:
    mol = MolFromSmiles(re.sub(r'\s.*$', '', smi))
    mols.append(mol)
  molWriter = SmilesWriter("-", delimiter="\t", nameHeader='Name', includeHeader=True, isomericSmiles=True, kekuleSmiles=False)
  fp.Utils.Mols2FPs_MACCSKeys(mols, molWriter)
 
#############################################################################
if __name__ == "__main__":
  MORGAN_NBITS=1024; MORGAN_RADIUS=2;
  parser = argparse.ArgumentParser(description="RDKit fingerprint generator", epilog="")
  OPS = ["demo", "demo_maccs", "demo_morgan", "fpgen", "fpgen_morgan", "fpgen_maccs"]
  parser.add_argument("op", choices=OPS, help="OPERATION")
  parser.add_argument("--i", dest="ifile", help="input file, TSV or SDF")
  parser.add_argument("--o", dest="ofile", help="output file, TSV|SDF")
  parser.add_argument("--scratchdir", default="/tmp")
  parser.add_argument("--smicol", type=int, default=1, help="SMILES column from TSV (counting from 1)")
  parser.add_argument("--namcol", type=int, default=2, help="name column from TSV (counting from 1)")
  parser.add_argument("--idelim", default="\t", help="delim for input TSV")
  parser.add_argument("--odelim", default="\t", help="delim for output TSV")
  parser.add_argument("--iheader", action="store_true", help="input TSV has header")
  parser.add_argument("--oheader", action="store_true", help="output TSV has header")
  parser.add_argument("--morgan_nbits", type=int, default=MORGAN_NBITS)
  parser.add_argument("--morgan_radius", type=int, default=MORGAN_RADIUS)
  parser.add_argument("-v", "--verbose", action="count", default=0)

  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f"RDKit version: {rdkit.rdBase.rdkitVersion}")
  #logging.info(f"MatplotLib version: {mpl.__version__}")

  t0=time.time()

  if args.op=="demo":
    Demo()
    sys.exit()

  elif args.op=="demo_morgan":
    DemoMorgan()
    sys.exit()

  elif args.op=="demo_maccs":
    DemoMACCSKeys()
    sys.exit()

  if not (args.ifile): parser.error('--i required.')

  molReader = util.Utils.File2Molreader(args.ifile, args.idelim, args.smicol-1, args.namcol-1, args.iheader)
  molWriter = util.Utils.File2Molwriter(args.ofile, args.odelim, args.oheader)
  mols = util.Utils.ReadMols(molReader)

  if args.op=="fpgen":
    fp.Utils.Mols2FPs_RDK(mols, molWriter)

  elif args.op=="fpgen_morgan":
    fp.Utils.Mols2FPs_Morgan(mols, args.morgan_radius, args.morgan_nbits, molWriter)

  elif args.op=="fpgen_maccs":
    fp.Utils.Mols2FPs_MACCS(mols, molWriter)

  else:
    parser.error(f"Unsupported operation: {args.op}")

  logging.info('Elapsed: %s'%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))))

