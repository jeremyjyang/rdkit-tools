#!/usr/bin/env python3
#
import sys,os
from rdkit import Chem
from rdkit.Chem.MACCSkeys import GenMACCSKeys

from rdkit.Chem import RDKFingerprint
from rdkit import DataStructs

fpers = [ RDKFingerprint,GenMACCSKeys ]
#fpers = [ GenMACCSKeys ]

##
## RDKFingerprint args:
## minPath=1,maxPath=7,fpSize=2048,nBitsPerHash=4,useHs=True,
## tgtDensity=0.0,minSize=128
## returns ExplicitBitVect
##
## GenMACCSKeys (no args)
## returns SparseBitVect
##

smiA='NCCc1cc(O)c(O)cc1'
smiB='CNC(C)Cc1cc(O)c(O)cc1'
molA = Chem.MolFromSmiles(smiA)
molB = Chem.MolFromSmiles(smiB)


for fper in fpers:
  fpA = fper(molA)
  fptxtA = DataStructs.BitVectToText(fpA)
  print(fptxtA)
  size=fpA.GetNumBits()
  count=0
  for i in range(fpA.GetNumBits()):
    if fpA.GetBit(i): count+=1
  print(count,size)
  if fper==GenMACCSKeys:
    onbits = tuple(fpA.GetOnBits())
    print(str(onbits))
    for i in onbits:
      print(Chem.MACCSkeys.smartsPatts[i][0])

  fpB = fper(molB)
  fptxtB = DataStructs.BitVectToText(fpB)
  print(fptxtB)
  size=fpB.GetNumBits()
  count=0
  for i in range(fpB.GetNumBits()):
    if fpB.GetBit(i): count+=1
  print(count,size)
  if fper==GenMACCSKeys:
    onbits = tuple(fpB.GetOnBits())
    print(str(onbits))
    for i in onbits:
      print(Chem.MACCSkeys.smartsPatts[i][0])


  sim = DataStructs.TanimotoSimilarity(fpA,fpB)
  print(sim)
  sim = DataStructs.TverskySimilarity(fpA,fpB,0.9,0.1)
  print(sim)
  bcom = DataStructs.NumBitsInCommon(fpA,fpB)
  print(bcom)
