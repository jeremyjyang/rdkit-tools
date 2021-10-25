#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.MurckoScaffold.html
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.rdScaffoldNetwork.html
rdScaffoldNetwork available RDKit 2020.03.1+.
"""
#############################################################################
import os,sys,re,logging,argparse,json,time,inspect
import pandas as pd

import matplotlib as mpl
#from matplotlib import pyplot as plt

import rdkit
import rdkit.Chem
import rdkit.Chem.AllChem
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles, Draw, rdDepictor

# fingerprints:
from rdkit.Chem import RDKFingerprint, PatternFingerprint, LayeredFingerprint, LayeredFingerprint_substructLayers
from rdkit.Chem.MACCSkeys import GenMACCSKeys
from rdkit.Chem.AllChem import GetMorganFingerprint,GetMorganFingerprintAsBitVect

#############################################################################
def ShowDetails(details):
  for key,val in inspect.getmembers(details):
    if not key.startswith('_'): # ignore private and protected functions
      if not inspect.ismethod(val): # ignore other methods 
        logging.debug(f"{key:>16}: {val}")

#############################################################################
def DemoPath():
  mols=[];
  for smi in rdktools_util.Utils.DEMOSMIS:
    mol = MolFromSmiles(re.sub(r'\s.*$', '', smi))
    mols.append(mol)
  molWriter = SmilesWriter("-", delimiter="\t", nameHeader='Name', includeHeader=True, isomericSmiles=True, kekuleSmiles=False)
  Mols2FPs_RDK(mols, molWriter)

#############################################################################
def DemoMorgan():
  mols=[];
  for smi in rdktools_util.Utils.DEMOSMIS:
    mol = MolFromSmiles(re.sub(r'\s.*$', '', smi))
    mols.append(mol)
  molWriter = SmilesWriter("-", delimiter="\t", nameHeader='Name', includeHeader=True, isomericSmiles=True, kekuleSmiles=False)
  Mols2FPs_Morgan(mols, 2, 1024, molWriter)

#############################################################################
def DemoMACCSKeys():
  mols=[];
  for smi in rdktools_util.Utils.DEMOSMIS:
    mol = MolFromSmiles(re.sub(r'\s.*$', '', smi))
    mols.append(mol)
  molWriter = SmilesWriter("-", delimiter="\t", nameHeader='Name', includeHeader=True, isomericSmiles=True, kekuleSmiles=False)
  Mols2FPs_MACCSKeys(mols, molWriter)

#############################################################################
def ListMACCSkeys():
  data=[];
  for i,val in rdkit.Chem.MACCSkeys.smartsPatts.items():
    smarts,n = val
    logging.debug(f"{i}\t{n}\t{smarts}")
    data.append([i, n, smarts])
  df = pd.DataFrame(data, columns=['I', 'N', 'SMARTS'])
  df.to_csv(sys.stdout, "\t", index=False)

#############################################################################
def Mols2FPs_RDK(mols, molWriter=None):
  for molA in mols:
    fpA = RDKFingerprint(molA)
    #fpstxt = rdkit.DataStructs.BitVectToFPSText(fpA)
    #logging.debug(f"({len(fpstxt)} bits) {fpstxt}")
    logging.debug(f"{fpA.ToBitString()} ({fpA.GetNumBits()} bits)")
    molA.SetProp("_Name", (molA.GetProp("_Name") if molA.HasProp("_Name") else "")+"\t"+fpstxt)
    count=0;
    for i in range(fpA.GetNumBits()):
      if fpA.GetBit(i): count+=1
    logging.debug(f"count={count}; size={fpA.GetNumBits()}")
    if molWriter is not None:
      molWriter.write(molA)

#############################################################################
def Mols2FPs_Morgan(mols, radius=2, nbits=1024, molWriter=None):
  for molA in mols:
    fpA = GetMorganFingerprintAsBitVect(molA, radius, nbits)
    #fpstxt = rdkit.DataStructs.BitVectToFPSText(fpA)
    #logging.debug(f"({len(fpstxt)} bits) {fpstxt}")
    logging.debug(f"{fpA.ToBitString()} ({fpA.GetNumBits()} bits)")
    molA.SetProp("_Name", (molA.GetProp("_Name") if molA.HasProp("_Name") else "")+"\t"+fpstxt)
    count=0;
    for i in range(fpA.GetNumBits()):
      if fpA.GetBit(i): count+=1
    logging.debug(f"count={count}; size={fpA.GetNumBits()}")
    if molWriter is not None:
      molWriter.write(molA)

#############################################################################
def Mols2FPs_MACCSKeys(mols, molWriter=None):
  for molA in mols:
    fpA = GenMACCSKeys(molA)
    #fpstxt = rdkit.DataStructs.BitVectToFPSText(fpA)
    #logging.debug(f"({len(fpstxt)} bits) {fpstxt}")
    logging.debug(f"{fpA.ToBitString()} ({fpA.GetNumBits()} bits)")
    molA.SetProp("_Name", (molA.GetProp("_Name") if molA.HasProp("_Name") else "")+"\t"+fpstxt)
    count=0;
    for i in range(fpA.GetNumBits()):
      if fpA.GetBit(i): count+=1
    logging.debug(f"count={count}; size={fpA.GetNumBits()}")
    onbits = tuple(fpA.GetOnBits())
    smarts = {i:rdkit.Chem.MACCSkeys.smartsPatts[i][0] for i in onbits}
    logging.debug(str(smarts))
    if molWriter is not None:
      molWriter.write(molA)

#############################################################################
def Demo2():
  fpers = [ RDKFingerprint, GenMACCSKeys ]
  #fpers = [ GenMACCSKeys ]

  ##
  ## RDKFingerprint args:
  ## minPath=1, maxPath=7, fpSize=2048, nBitsPerHash=4, useHs=True,
  ## tgtDensity=0.0, minSize=128
  ## returns ExplicitBitVect
  ##
  ## GenMACCSKeys (no args)
  ## returns SparseBitVect
  ##
  smiA='NCCc1cc(O)c(O)cc1'
  smiB='CNC(C)Cc1cc(O)c(O)cc1'
  molA = rdkit.Chem.MolFromSmiles(smiA)
  molB = rdkit.Chem.MolFromSmiles(smiB)

  for fper in fpers:
    fpA = fper(molA)
    fptxtA = rdkit.DataStructs.BitVectToText(fpA)
    print(fptxtA)
    size=fpA.GetNumBits()
    count=0
    for i in range(fpA.GetNumBits()):
      if fpA.GetBit(i): count+=1
    print(count, size)
    if fper==GenMACCSKeys:
      onbits = tuple(fpA.GetOnBits())
      print(str(onbits))
      for i in onbits:
        print(rdkit.Chem.MACCSkeys.smartsPatts[i][0])

    fpB = fper(molB)
    fptxtB = rdkit.DataStructs.BitVectToText(fpB)
    print(fptxtB)
    size=fpB.GetNumBits()
    count=0
    for i in range(fpB.GetNumBits()):
      if fpB.GetBit(i): count+=1
    print(count, size)
    if fper==GenMACCSKeys:
      onbits = tuple(fpB.GetOnBits())
      print(str(onbits))
      for i in onbits:
        print(rdkit.Chem.MACCSkeys.smartsPatts[i][0])

    sim = rdkit.DataStructs.TanimotoSimilarity(fpA, fpB)
    print(sim)
    sim = rdkit.DataStructs.TverskySimilarity(fpA, fpB, 0.9, 0.1)
    print(sim)
    bcom = rdkit.DataStructs.NumBitsInCommon(fpA, fpB)
    print(bcom)

#############################################################################
