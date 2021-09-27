#!/usr/bin/env python3
"""
Molecule formats: SMI, MDL, MOL2, PDB
Note that RDKit requires PIL, Python Imaging Library.
"""
###
import os,sys,cgi,logging
import re,base64

import rdkit.Chem.AllChem
import rdkit.Chem.Draw

from .. import util
from .. import depict

#############################################################################
def Mol2Img(form):
  mol=None; smi=None; molname=None;
  for tag in ('smi', 'smiles', 'smicode'):
    if form.getvalue(tag, ''):
      smi = form.getvalue(tag)
      if re.search(r'\s', smi):
        smi,molname = re.split(r'\s', smi, 1)
      mol = rdkit.Chem.MolFromSmiles(smi)

  molfmt = form.getvalue('molfmt', 'MDL')

  if form.getvalue('fcode', ''):
    fdata = base64.decodestring(form.getvalue('fcode'))
    if form.getvalue('gz', 'false')=='true':
      fdata = util.http.Utils.GunzipBytes(fdata)
    if molfmt.upper() in ('MDL', 'SDF', 'MOL'): 
      sdms = rdkit.Chem.AllChem.SDMolSupplier()
      sdms.SetData(fdata, sanitize=True, removeHs=True)
      mol = sdms[0]
    elif molfmt.upper() == 'PDB':
      mol = rdkit.Chem.MolFromPDB2Block(fdata, sanitize=True, removeHs=True, flavor=None)
    elif molfmt.upper() == 'MOL2':
      mol = rdkit.Chem.MolFromMol2Block(fdata, sanitize=True, removeHs=True)

  if not mol: sys.exit()

  if not molname:
    molname = form.getvalue('title', '')

  try: width = int(form.getvalue('w'))
  except: width=300
  try: height = int(form.getvalue('h'))
  except: height=180

  smarts=None
  for tag in ('smarts', 'smartscode'):
    if form.getvalue(tag, ''):
      smarts = form.getvalue(tag)

  kekule = bool(form.getvalue('kekule', ''))

  gen2d = bool(form.getvalue('gen2d') or form.getvalue('ignore2d'))

  nconf = mol.GetNumConformers()

  if gen2d or nconf==0:
    rdkit.Chem.AllChem.Compute2DCoords(mol, clearConfs=gen2d)

  hitatoms=[]
  #highmap = {} ##dict of (atom,color) pairs
  if smarts:
    pat = rdkit.Chem.MolFromSmarts(smarts)
    matches = mol.GetSubstructMatches(pat, uniquify=True, useChirality=False)
    for match in matches:
      hitatoms.extend(match)
      hitatoms.sort()  ##repeats ok?

    ###  Does not work, would be nice for arbitrary atom coloring.
    ###  How are colors specified?
    #for atom in mol.GetAtoms():
    #  highmap[atom]=rdkit.sping.pid.darkblue
    #  highmap[atom]=(0,0,0)
    #for aidx in hitatoms:
    #  highmap[mol.GetAtomWithIdx(aidx)]=rdkit.sping.pid.red
    #  highmap[mol.GetAtomWithIdx(aidx)]=(1,0,0)

    # By default the RDKit colors atoms by element in depictions.
    # We can turn this off by replacing the element dictionary:
    import collections
    from rdkit.Chem.Draw import MolDrawing
    from rdkit.Chem.Draw.MolDrawing import DrawingOptions

    DrawingOptions.elemDict = collections.defaultdict(lambda : (0,0,0))

  #img = rdkit.Chem.Draw.MolToImage(mol, size=(width,height), kekulize=kekule, highlightAtoms=hitatoms, wedgeBonds=True, legend=molname)

  # PIL.Image.Image
  img = depict.Utils.Mol2Image(mol, width=width, height=height, highlightAtoms=hitatoms, kekulize=True, wedgeBonds=True)

  sys.stdout.buffer.write(b'Content-type: image/png\n\n')
  img.save(sys.stdout.buffer, format='PNG')
  sys.stdout.flush()

#############################################################################
if __name__=='__main__':
  logging.basicConfig(level=logging.DEBUG)
  form = cgi.FieldStorage(keep_blank_values=1)
  Mol2Img(form)
