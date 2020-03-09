#!/usr/bin/env python
#############################################################################
### mol2img_rdk_utils.py
###
### Jeremy Yang
### 16 Oct 2009
#############################################################################
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import sys,os


#############################################################################
if __name__=='__main__':
  import warnings
  #import tempfile
  #fpath=tempfile.mktemp(PROG+'.jpg')
  fpath='z.png'
  print >>sys.stderr, fpath

  mdlstr='''\
DOPAMINE
  -OEChem-04230909332D

 11 11  0     0  0  0  0  0  0999 V2000
   -3.4670    1.9950    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5995    1.4976    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7350    2.0001    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8675    1.5027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8675    0.4975    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8675    0.4975    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7328   -0.0038    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8675    1.5027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7350    2.0001    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000    2.0104    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4 11  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  7  1  0  0  0  0
  7  8  1  0  0  0  0
  7  9  2  0  0  0  0
  9 10  1  0  0  0  0
  9 11  1  0  0  0  0
M  END
'''

  mol=Chem.MolFromMolBlock(mdlstr)

  #print >>sys.stderr, 'NumConfs: %d'%mol.GetNumConformers()
  AllChem.Compute2DCoords(mol,clearConfs=False)
  #print >>sys.stderr, 'NumConfs: %d'%mol.GetNumConformers()

  sma='a[R](a)-[!R]'
  pat=Chem.MolFromSmarts(sma)
  matches=mol.GetSubstructMatches(pat,uniquify=True,useChirality=False)

  hitatoms=[]
  for match in matches:
    print >>sys.stderr, match
    hitatoms.extend(match)
    hitatoms.sort()  ##repeats ok?

  print Chem.MolToSmiles(mol)
