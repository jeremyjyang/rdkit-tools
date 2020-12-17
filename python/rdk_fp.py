#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/
https://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints

"""
#############################################################################
import os,sys,re,json,time,inspect,argparse,logging
import matplotlib as mpl
#from matplotlib import pyplot as plt

import rdkit
#from rdkit import DataStructs
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles, Draw, rdDepictor
from rdkit.Chem import RDKFingerprint, PatternFingerprint, LayeredFingerprint, LayeredFingerprint_substructLayers
from rdkit.Chem.MACCSkeys import GenMACCSKeys
from rdkit.Chem.AllChem import GetMorganFingerprint

from rdkit.Chem.Draw import rdMolDraw2D


DEMOSMIS = [
'C[C@H]1CN(CCCN1[S+]([O-])(=O)C2=CC=CC3=C2C(=CN=C3)C)C(=O)CN Rho Kinase Inhibitor IV',
'N[S+]([O-])(=O)C1=C(Cl)C=C2NC(N[S+]([O-])(=O)C2=C1)C(Cl)Cl trichlormethiazide',
'C[S+]([O-])C1=CC=C(C=C1)\C=C2\C(=C(\CC(O)=O)C3=C2C=CC(=C3)F)C Sulindac',
'O=[N+]([O-])c1cc(S(=O)(=O)N2CCCC2)ccc1NN=Cc1ccccn1 PUBCHEM_CID:4036736',
'Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12 flucloxacillin',
'CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O ampicillin',
'CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O.[Na] penicillin',
'Cc1onc(-c2ccccc2)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12 oxacillin'
	]

#############################################################################
def moltosvg(mol, molSize=(450,250), kekulize=True):
    mc = rdMolDraw2D.PrepareMolForDrawing(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
    opts = drawer.drawOptions()
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg.replace('svg:', '')

#############################################################################
def ReadMols(molReader):
  mols=[]; 
  for mol in molReader:
    if mol is not None: mols.append(mol)
  logging.info(f'{len(mols)} mols read from {molReader}')
  return mols

#############################################################################
def File2Molreader(ifile, idelim, smicol, namcol, iheader):
  if re.sub(r'.*\.', '', ifile).lower()in ('smi', 'smiles', 'csv', 'tsv'):
    molReader = SmilesMolSupplier(ifile, delimiter=idelim, smilesColumn=smicol, nameColumn=namcol, titleLine=iheader, sanitize=True)
  elif re.sub(r'.*\.', '', ifile).lower() in ('sdf','sd','mdl','mol'):
    molReader = SDMolSupplier(ifile, sanitize=True, removeHs=True)
  else:
    molReader = None
    logging.error(f'Invalid file extension: {ifile}')
  return molReader

#############################################################################
def File2Molwriter(ofile, odelim, oheader):
  if not ofile:
    molWriter = SmilesWriter("-", delimiter=odelim, nameHeader='Name', includeHeader=oheader, isomericSmiles=True, kekuleSmiles=False)
  elif re.sub(r'.*\.', '', ofile).lower() in ('sdf','sd','mdl','mol'):
    molWriter = SDWriter(ofile)
  elif re.sub(r'.*\.', '', ofile).lower() in ('smi', 'smiles', 'csv', 'tsv'):
    molWriter = SmilesWriter(ofile, delimiter=odelim, nameHeader='Name', includeHeader=oheader, isomericSmiles=True, kekuleSmiles=False)
  else:
    logging.error(f'Invalid file extension: {ofile}')
    molWriter = None
  return molWriter

#############################################################################
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
#############################################################################
def Mols2FPs(mols, fper, molWriter):
  #rdkit.DataStructs.CreateFromBinaryText
  #rdkit.DataStructs.CreateFromFPSText
  #rdkit.DataStructs.CreateFromBitString

  for molA in mols:
    fpA = fper(molA)
    #logging.debug(rdkit.DataStructs.BitVectToText(fpA))
    #logging.debug(rdkit.DataStructs.BitVectToBinaryText(fpA))
    logging.debug(rdkit.DataStructs.BitVectToFPSText(fpA))
    
    count=0;
    for i in range(fpA.GetNumBits()):
      if fpA.GetBit(i): count+=1
    logging.debug(f"count={count}; size={fpA.GetNumBits()}")
    if fper==GenMACCSKeys:
      onbits = tuple(fpA.GetOnBits())
      logging.debug(str(onbits))
      for i in onbits:
        logging.debug(rdkit.Chem.MACCSkeys.smartsPatts[i][0])

#############################################################################
def Demo(fper):
  logging.debug(str(fper))
  mols=[];
  for smi in DEMOSMIS:
    mol = MolFromSmiles(re.sub(r'\s.*$', '', smi))
    mols.append(mol)
  molWriter = SmilesWriter("-", delimiter="\t", nameHeader='Name', includeHeader=True, isomericSmiles=True, kekuleSmiles=False)
  Mols2FPs(mols, fper, molWriter)
 
#############################################################################
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="RDKit fingerprint generator", epilog="")
  OPS = ["demo", "fpgen"]
  parser.add_argument("op", choices=OPS, help="OPERATION")
  parser.add_argument("--i", dest="ifile", help="input file, TSV or SDF")
  parser.add_argument("--o", dest="ofile", help="output file, TSV|SDF")
  parser.add_argument("--scratchdir", default="/tmp")
  parser.add_argument("--smicol", type=int, default=0, help="SMILES column from TSV (counting from 0)")
  parser.add_argument("--namcol", type=int, default=1, help="name column from TSV (counting from 0)")
  parser.add_argument("--idelim", default="\t", help="delim for input TSV")
  parser.add_argument("--odelim", default="\t", help="delim for output TSV")
  parser.add_argument("--iheader", action="store_true", help="input TSV has header")
  parser.add_argument("--oheader", action="store_true", help="output TSV has header")
  parser.add_argument("-v", "--verbose", action="count", default=0)

  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f"RDKit version: {rdkit.rdBase.rdkitVersion}")
  logging.info(f"MatplotLib version: {mpl.__version__}")

  t0=time.time()

  fpers = [RDKFingerprint, GenMACCSKeys, GetMorganFingerprint]
  #fper = RDKFingerprint
  fper = GenMACCSKeys
  #fper = GetMorganFingerprint

  if args.op=="demo":
    Demo(fper)
    sys.exit()

  if not (args.ifile): parser.error('--i required.')

  if args.op=="fpgen":
    molReader = File2Molreader(args.ifile, args.idelim, args.smicol, args.namcol, args.iheader)
    molWriter = File2Molwriter(args.ofile, args.odelim, args.oheader)
    mols = ReadMols(molReader)
    Mols2FPs(mols, fper, molWriter)

  else:
    parser.error(f"Unsupported operation: {args.op}")

  logging.info('Elapsed: %s'%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))))

