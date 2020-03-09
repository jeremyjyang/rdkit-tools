#!/usr/bin/env python3
#############################################################################
import rdkit.rdBase
import rdkit.Chem
import sys,os

#############################################################################
### DeduplicateMatches() - remove matches not USA (Unique Sets of Atoms)
### uumatches - unique set of matches, sorted atom order
### umatches - unique set of matches, original atom order
#############################################################################
def DeduplicateMatches(matches):
  uumatches=[]
  umatches=[]
  for match in matches:
    uumatch=list(match)
    uumatch.sort()
    if uumatch not in uumatches:
      uumatches.append(uumatch)
      umatches.append(match)
  return tuple(umatches)

#############################################################################
if __name__=='__main__':

  print('RDKit version: %s'%rdkit.rdBase.rdkitVersion, file=sys.stderr) 

  smas=[
	'*~1~*~*~*~*~*1',
	'*~1~*~*~*~*~*~*~*~*~*1',
	'*~1~*~*~*~*~*~*~*~*~*~*~*~*~*1',
  ]
  smis=[
	'C1CCCCC1',
	'C1CCCC2C1CCCC2',
	'N1CCCC2C1CCCC2',
	'C1CCC3C4CCCCC4CCC3C1',
	'CC1CCC3C4CCCCC4CCC3C1',
  ]

  for smi in smis:
    print('smi: %s'%smi)
    for sma in smas:
      print('\tsma: %s'%sma)
      pat=rdkit.Chem.MolFromSmarts(sma)
      mol=rdkit.Chem.MolFromSmiles(smi)
      matches=mol.GetSubstructMatches(pat,uniquify=True,useChirality=False)
      for match in matches: print('\t',match)
      nmatches=len(matches),
      print('\tmatches: %d'%nmatches,)
      umatches=DeduplicateMatches(matches)
      nmatches_usa=len(umatches),
      print('\tusa matches: %d'%nmatches_usa,)
      if nmatches==nmatches_usa:
        print()
      else:
        print('<-- LOOK !!!')

