#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.MurckoScaffold.html
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.rdScaffoldNetwork.html
rdScaffoldNetwork available RDKit 2020.03.1+.
"""
#############################################################################
import os,sys,re,logging,json,time,inspect

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
def Mols2FPs_RDK(mols, molWriter=None):
  for molA in mols:
    fpA = RDKFingerprint(molA)
    fpstxt = rdkit.DataStructs.BitVectToFPSText(fpA)
    logging.debug(fpstxt)
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
    fpstxt = rdkit.DataStructs.BitVectToFPSText(fpA)
    logging.debug(fpstxt)
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
    fpstxt = rdkit.DataStructs.BitVectToFPSText(fpA)
    logging.debug(fpstxt)
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
