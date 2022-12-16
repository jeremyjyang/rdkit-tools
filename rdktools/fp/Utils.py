#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.MurckoScaffold.html
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.rdScaffoldNetwork.html
rdScaffoldNetwork available RDKit 2020.03.1+.

https://www.rdkit.org/docs/source/rdkit.Chem.AtomPairs.Pairs.html
"""
#############################################################################
import os,sys,re,logging,argparse,json,time,inspect,tqdm,pickle
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
import rdkit.Chem.AtomPairs.Pairs

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
def ListMACCSkeys(fout):
  data=[];
  for i,val in rdkit.Chem.MACCSkeys.smartsPatts.items():
    smarts,n = val
    logging.debug(f"{i}\t{n}\t{smarts}")
    data.append([i, n, smarts])
  df = pd.DataFrame(data, columns=['I', 'N', 'SMARTS'])
  if fout is not None: df.to_csv(fout, "\t", index=False)
  logging.info(f"n_out: {df.shape[0]}")
  return df

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
def AtomPairFingerprints(details, reportFreq):
  # AtomPairs.Pairs.GetAtomPairFingerprintAsBitVect(mol)
  # AtomPairs.Pairs.pyScorePair(at1, at2, dist, atomCodes=None, includeChirality=False)
  # AtomPairs.Pairs.ExplainPairScore(score, includeChirality=False)
  # AtomPairs.Pairs.GetAtomPairFingerprintAsBitVect()
  # AtomPairs.Pairs.GetAtomPairFingerprintAsIntVect()
  # AtomPairs.Pairs.GetHashedAtomPairFingerprint()
  if details.inFileName and details.useSmiles:
    from rdkit.ML.Data import DataUtils
    conn = None
    if not details.idName:
      details.idName = 'ID'
    try:
      dataSet = DataUtils.TextFileToData(details.inFileName, onlyCols=[details.idName, details.smilesName])
    except IOError:
      import traceback
      logging.error(f"Problems reading from file {details.inFileName}")
      traceback.print_exc()
    idCol = 0
    smiCol = 1
  elif details.inFileName and details.useSD:
    conn = None
    if not details.idName:
      details.idName = 'ID'
    dataSet = []
    molReader = Chem.SDMolSupplier(details.inFileName)
    for mol in molReader:
      dataSet.append(mol)
    for i, mol in enumerate(dataSet):
      if mol.HasProp(details.idName):
        nm = mol.GetProp(details.idName)
      else:
        nm = mol.GetProp('_Name')
      dataSet[i] = (nm, mol)
  else:
    dataSet = None

  fps = None
  if dataSet :
    if details.useSD:
      fps = AtomPairFingerprintsFromMols(dataSet, reportFreq=reportFreq, **details.__dict__)
    else:
      data = dataSet.GetNamedData()
      if details.molPklName:
        fps = AtomPairFingerprintsFromPickles(data, idCol, smiCol, reportFreq=reportFreq, **details.__dict__)
      else:
        fps = AtomPairFingerprintsFromSmiles(data, idCol, smiCol, reportFreq=reportFreq, **details.__dict__)

  if fps:
    if details.outFileName:
      logging.info(f"Writing pickled FPs to {details.outFileName}")
      outF = open(details.outFileName, 'wb+')
      for i in range(len(fps)):
        pickle.dump(fps[i], outF)
      outF.close()
    dbName = details.outDbName or details.dbName
    if details.outTableName and dbName:
      from rdkit.Dbase.DbConnection import DbConnect
      from rdkit.Dbase import DbUtils, DbModule
      conn = DbConnect(dbName)
      #
      #  Figure out column-types.
      colTypes = DbUtils.TypeFinder(data, len(data), len(data[0]))
      typeStrs = DbUtils.GetTypeStrings([details.idName, details.smilesName], colTypes, keyCol=details.idName)
      cols = f"{typeStrs[0]}, {details.fpColName} {DbModule.binaryTypeName}"

      # create the new table
      if details.replaceTable or \
         details.outTableName.upper() not in [x.upper() for x in conn.GetTableNames()]:
        conn.AddTable(details.outTableName, cols)

      # And add the data
      for ID, fp in fps:
        tpl = ID, DbModule.binaryHolder(fp.ToBinary())
        conn.InsertData(details.outTableName, tpl)
      conn.Commit()
  return fps

#############################################################################
def AtomPairFingerprintsFromSmiles(dataSource, idCol, smiCol, reportFreq=10, maxMols=-1, **fpArgs):
  res=[]; nDone=0; tq=None;
  for entry in dataSource:
    ID, smi = str(entry[idCol]), str(entry[smiCol])
    mol = rdkit.Chem.MolFromSmiles(smi)
    if mol is not None:
      fp = rdkit.Chem.AtomPairs.Pairs.GetAtomPairFingerprintAsBitVect(mol)
      logging.debug(f"fp.GetNumOnBits()/fp.GetNumBits():{fp.GetNumOnBits()}/{fp.GetNumBits()}")
      #logging.debug(f"{fp.ToBitString()}")
      res.append((ID, fp))
      if not tq: tq = tqdm.tqdm(total=len(dataSource), unit="molecules")
      tq.update()
      nDone+=1
      if maxMols>0 and nDone>=maxMols:
        break
    else:
      logging.error(f"Problems parsing SMILES: {smi}")
  if tq is not None: tq.close()
  return res

#############################################################################
def AtomPairFingerprintsFromPickles(dataSource, idCol, pklCol, reportFreq=10, maxMols=-1,**fpArgs):
  res=[]; nDone=0; tq=None;
  for entry in dataSource:
    ID, pkl = str(entry[idCol]), str(entry[pklCol])
    mol = rdkit.Chem.Mol(pkl)
    if mol is not None:
      fp = rdkit.Chem.AtomPairs.Pairs.GetAtomPairFingerprintAsBitVect(mol)
      res.append((ID, fp))
      if not tq: tq = tqdm.tqdm(total=len(dataSource), unit="molecules")
      tq.update()
      nDone += 1
      if maxMols>0 and nDone>=maxMols:
        break
    else:
      logging.error(f"Problems parsing pickle for ID: {ID}")
  if tq is not None: tq.close()
  return res

#############################################################################
def AtomPairFingerprintsFromMols(dataSet, reportFreq=10, maxMols=-1, **fpArgs):
  res=[]; nDone=0; tq=None;
  for ID, mol in dataSet:
    if mol:
      fp = rdkit.Chem.AtomPairs.Pairs.GetAtomPairFingerprintAsBitVect(mol)
      res.append((ID, fp))
      if not tq: tq = tqdm.tqdm(total=len(dataSet), unit="molecules")
      tq.update()
      nDone+=1
      if maxMols>0 and nDone>=maxMols:
        break
    else:
      logging.error(f"Problems parsing SMILES: {smi}")
  if tq is not None: tq.close()
  return res

#############################################################################

