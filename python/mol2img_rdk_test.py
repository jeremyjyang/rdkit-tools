#!/usr/bin/env python
#############################################################################
### mol2img_rdk_test.py
###
### Jeremy Yang
###  4 Jun 2014
#############################################################################
import sys,os
import StringIO

#import warnings

import rdkit.Chem.Draw
import rdkit.Chem.AllChem
import rdkit.Chem.rdmolfiles

#import mol2img_rdk_utils


#############################################################################
if __name__=='__main__':

  fpath='z.png'

  mdlstr='''\
[His] histidine (H)
  -OEChem-11160913072D

 10 10  0     1  0  0  0  0  0999 V2000
   -3.6947   -0.5273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0915    0.2684    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   -1.4054    0.2766    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4082    1.9369    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1221    0.2784    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3031   -1.4474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3288   -1.6726    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2358   -2.6686    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1596   -3.0629    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8168   -2.3070    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  2  5  1  0  0  0  0
  1  6  1  0  0  0  0
  6 10  1  0  0  0  0
  6  7  2  0  0  0  0
  7  8  1  0  0  0  0
  8  9  2  0  0  0  0
  9 10  1  0  0  0  0
M  END
> <COMPOUND_CID>
999999

> <LogP>
9.999999

$$$$
'''

  #Works, but ignores properties.
  #mol=rdkit.Chem.MolFromMolBlock(mdlstr, sanitize=True, removeHs=True, strictParsing=True)

  #Works, and reads properties.
  sdms = rdkit.Chem.AllChem.SDMolSupplier()
  sdms.SetData(mdlstr, sanitize=True, removeHs=True)
  mol = sdms[0]

  print >>sys.stderr, '%s'%(mol.GetProp('_Name'))
  mol_kek = rdkit.Chem.AllChem.Mol(mol)
  rdkit.Chem.Kekulize(mol_kek)
  print >>sys.stderr, '%s'%(rdkit.Chem.MolToSmiles(mol_kek,kekuleSmiles=True))
  print >>sys.stderr, '%s'%(rdkit.Chem.MolToSmiles(mol,isomericSmiles=True))
  print >>sys.stderr, '\t%d atoms, %d bonds'%(mol.GetNumAtoms(),mol.GetNumBonds())
  for propname in mol.GetPropNames():
    print >>sys.stderr, '\tProperty ["%s"]: %s'%(propname,mol.GetProp(propname))

  #smi='c1c(CN)cccc1C(=O)O'

  #mol=rdkit.Chem.MolFromSmiles(smi,sanitize=True,replacements={})
  #mol=rdkit.Chem.MolFromMol2Block(mdlstr, sanitize=True, removeHs=True)
  #mol=rdkit.Chem.MolFromPDB2Block(mdlstr, sanitize=True, removeHs=True, flavor=None)
  #mol=rdkit.Chem.MolFromTPL2Block(mdlstr, sanitize=True, skipFirstConf=False)

  print >>sys.stderr, 'NumConfs: %d'%mol.GetNumConformers()
  rdkit.Chem.AllChem.Compute2DCoords(mol,clearConfs=False)
  print >>sys.stderr, 'NumConfs: %d'%mol.GetNumConformers()

  smarts='a:a'
  pat=rdkit.Chem.MolFromSmarts(smarts)
  matches=mol.GetSubstructMatches(pat,uniquify=True,useChirality=False)

  hitatoms=[]
  for match in matches:
    #print >>sys.stderr, match
    hitatoms.extend(match)
    hitatoms.sort()  ##repeats ok?

  width = 400
  height = 300
  #warnings.filterwarnings("ignore")

  img = rdkit.Chem.Draw.MolToImage(mol,size=(width,height),kekulize=True,highlightAtoms=hitatoms,wedgeBonds=True)
  print >>sys.stderr, 'NumConfs: %d'%mol.GetNumConformers()

  strio = StringIO.StringIO()
  img.save(strio,format='PNG')
  imgbytes = strio.getvalue()
  strio.close()
  f = open(fpath,'wb')
  f.write(imgbytes)
  f.close()

  print >>sys.stderr, 'Image written to: %s'%fpath

  #img.save(fpath) ##also works ok
  #img.save(sys.stdout,'PNG') ##also works ok

