#!/usr/bin/env python
#############################################################################
### mdl2img_rdk.py
###
### to do:
###   [ ] kekule
###
### Jeremy J Yang
### 21 Apr 2011
#############################################################################
import os,sys,cgi,re,tempfile,base64,urllib
#import cgitb; cgitb.enable()   ### debugger
import warnings,exceptions

import htm_utils

from rdkit import Chem
from rdkit.Chem import AllChem

import mol2img_rdk_utils

#############################################################################
def Mdl2Img(form):
  mdl=urllib.unquote(form.getvalue('mdlcode'))
  mdl=base64.decodestring(mdl)
  #mdl=zlib.decompress(mdl)
  mdl=htm_utils.GunzipBytes(mdl)
  if not mdl: sys.exit()
  try: w=int(form.getvalue('w'))
  except: w=200
  try: h=int(form.getvalue('h'))
  except: h=160
  smarts=None
  for tag in ('smarts','smartscode'):
    if form.getvalue(tag,''):
      smarts=form.getvalue(tag)
  kekule=bool(form.getvalue('kekule',''))
  mol=Chem.MolFromMolBlock(mdl)
  hitatoms=[]
  if smarts:
    pat=Chem.MolFromSmarts(smarts)
    matches=mol.GetSubstructMatches(pat,uniquify=True,useChirality=False)
    for match in matches:
      hitatoms.extend(match)
      hitatoms.sort()  ##repeats ok?
  warnings.filterwarnings('ignore')	##for some deprecated PIL code
  img=mol2img_rdk_utils.MyMolToImage(mol,size=(w,h),
  	highlightAtoms=hitatoms,kekulize=kekule,wedgeBonds=True)
  sys.stdout.write('Content-type: image/png\n\n')
  sys.stdout.flush()
  img.save(sys.stdout,'png')
  sys.stdout.flush()

#############################################################################
if __name__=='__main__':
  form=cgi.FieldStorage(keep_blank_values=1)
  Mdl2Img(form)
