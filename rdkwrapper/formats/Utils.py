#!/usr/bin/env python3
#############################################################################
import sys,os,re,logging

import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.Chem.inchi

#############################################################################
def Mol2Inchi():
  if not rdkit.Chem.inchi.INCHI_AVAILABLE:
    logging.error(f"INCHI_AVAILABLE={rdkit.Chem.inchi.INCHI_AVAILABLE}")
    exit(1)
  smi = "NCCc1cc(O)c(O)cc1"
  mol = rdkit.Chem.MolFromSmiles(smi)
  #inchi,auxinfo = rdkit.Chem.inchi.MolToInchiAndAuxInfo(mol, options='', logLevel= None, treatWarningAsError=False)
  inchi = rdkit.Chem.inchi.MolToInchi(mol, options='', logLevel=None, treatWarningAsError=False)
  inchikey = rdkit.Chem.inchi.MolToInchiKey(mol, options='')
  logging.debug(f"SMILES: \"{smi}\"; INCHI: \"{inchi}\"; INCHIKEY: \"{inchikey}\"")

#############################################################################
if __name__=='__main__':
  logging.basicConfig(level=logging.DEBUG)
  logging.info("Kekulization test...")
  #sma='[R]-O'
  #sma='*-[#6]'
  #smi='NCCc1ccc(O)c(O)c1'
  #smi="NC(C=O)CC1=CC=C(O)C=C1"
  #smi="C(C(C=O)N)C1=CN=CN1"
  smi="c1c([nH]cn1)CC(C=O)N"
  sma='[R]'

  pat = rdkit.Chem.MolFromSmarts(sma)
  if not pat:
    logging.error(f'Bad smarts: {sma}')
  mol = rdkit.Chem.MolFromSmiles(smi)
  matches = mol.GetSubstructMatches(pat, uniquify=True, useChirality=False)
  logging.info(matches)
  alist=[]
  for match in matches:
    for a in match:
      if a not in alist:
        alist.append(a)
  alist.sort()
  logging.info(alist)
  rdkit.Chem.AllChem.Kekulize(mol, clearAromaticFlags=True)
  keksmi = rdkit.Chem.MolToSmiles(mol, isomericSmiles=False, kekuleSmiles=True)
  keksmi += (' |ha:{0}|'.format(','.join(map(lambda x:str(x), alist))))
  logging.info(keksmi)
