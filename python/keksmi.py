#!/usr/bin/env python3
#############################################################################
import rdkit.Chem
import rdkit.Chem.AllChem
import sys,os,getopt,re

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
