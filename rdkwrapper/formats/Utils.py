#!/usr/bin/env python3
#############################################################################
import sys,os,re,logging

import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.Chem.inchi

#############################################################################
def Mol2Inchi():
  if not rdkit.Chem.inchi.INCHI_AVAILABLE:
    logging.error("INCHI_AVAILABLE={}".format(rdkit.Chem.inchi.INCHI_AVAILABLE))
    exit(1)

  smi = "NCCc1cc(O)c(O)cc1"
  mol = rdkit.Chem.MolFromSmiles(smi)

  #inchi,auxinfo = rdkit.Chem.inchi.MolToInchiAndAuxInfo(mol, options='', logLevel= None, treatWarningAsError=False)
  inchi = rdkit.Chem.inchi.MolToInchi(mol, options='', logLevel=None, treatWarningAsError=False)
  #inchikey = rdkit.Chem.inchi.InchiToInchiKey(inchi)
  inchikey = rdkit.Chem.inchi.MolToInchiKey(mol, options='')

  logging.debug("SMILES: \"{}\"; INCHI: \"{}\"; INCHIKEY: \"{}\"".format(smi, inchi, inchikey))

#############################################################################
if __name__=='__main__':

  #sma='[R]-O'
  #sma='*-[#6]'
  #smi='NCCc1ccc(O)c(O)c1'
  #smi="NC(C=O)CC1=CC=C(O)C=C1"
  #smi="C(C(C=O)N)C1=CN=CN1"
  smi="c1c([nH]cn1)CC(C=O)N"
  sma='[R]'

  pat=rdkit.Chem.MolFromSmarts(sma)
  if not pat:
    print('Bad smarts: %s'%(sma), file=sys.stderr)

  mol=rdkit.Chem.MolFromSmiles(smi)

  matches=mol.GetSubstructMatches(pat,uniquify=True,useChirality=False)

  print(matches)
  alist=[]
  for match in matches:
    for a in match:
      if a not in alist:
        alist.append(a)
  alist.sort()
  print(alist)

  rdkit.Chem.AllChem.Kekulize(mol,clearAromaticFlags=True)

  keksmi=rdkit.Chem.MolToSmiles(mol,isomericSmiles=False,kekuleSmiles=True)

  print(keksmi)
  keksmi+=(' |ha:%s|'%(','.join(map(lambda x:str(x),alist))))
  print(keksmi)
