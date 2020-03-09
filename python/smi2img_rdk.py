#!/usr/bin/env python3
#############################################################################
### smi2img_rdk.py
###
### Jeremy Yang
#############################################################################
import os,sys,cgi,re
#import cgitb; cgitb.enable()   ### debugger

import rdkit.Chem.AllChem
import rdkit.Chem.Draw

import mol2img_rdk_utils

#############################################################################
def Smi2Img(form):
  smi=None
  for tag in ('smi','smiles','smicode'):
    if form.getvalue(tag,''):
      smi=form.getvalue(tag)
      break
  if not smi: sys.exit()

  try: width=int(form.getvalue('w'))
  except: width=200
  try: height=int(form.getvalue('h'))
  except: height=160

  smarts=None
  for tag in ('smarts','smartscode'):
    if form.getvalue(tag,''):
      smarts=form.getvalue(tag)

  kekule=bool(form.getvalue('kekule',''))

  smi=re.sub(r'\s.*$','',smi) #remove name if present

  mol=rdkit.Chem.MolFromSmiles(smi)
  rdkit.Chem.AllChem.Compute2DCoords(mol)

  hitatoms=[]
  if smarts:
    pat=rdkit.Chem.MolFromSmarts(smarts)
    matches=mol.GetSubstructMatches(pat,uniquify=True,useChirality=False)
    for match in matches:
      hitatoms.extend(match)
      hitatoms.sort()  ##repeats ok?

  #img=rdkit.Chem.Draw.MolToImage(mol, size=(width,height), kekulize=kekule, highlightAtoms=hitatoms, wedgeBonds=True)

  img=mol2img_rdk_utils.Mol2Image(mol,width=width,height=height,
	highlightAtoms=hitatoms,kekulize=True,wedgeBonds=True,verbose=True)

  sys.stdout.write('Content-type: image/png\n\n')
  img.save(sys.stdout,format='PNG')
  sys.stdout.flush()

#############################################################################
if __name__=='__main__':
  form=cgi.FieldStorage(keep_blank_values=1)
  Smi2Img(form)
