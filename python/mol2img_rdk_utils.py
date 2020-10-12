#!/usr/bin/env python3
#############################################################################
### mol2img_rdk_utils.py
#############################################################################
import sys,os,logging
from PIL import Image

import rdkit.Chem
import rdkit.Chem.Draw
#import rdkit.Chem.AllChem
#import rdkit.rdBase

#############################################################################
def InitializeCanvas(verbose=0):
  '''Initialize Canvas, depending on library availability.'''
  global RDKIT_useAGG, RDKIT_useCAIRO, RDKIT_useSPING
  RDKIT_useAGG=False
  RDKIT_useCAIRO=False
  RDKIT_useSPING=False

  try:
    import cairo
    from rdkit.Chem.Draw.cairoCanvas import Canvas
    RDKIT_useCAIRO=True
    logging.debug('RDKIT_useCAIRO=True')
  except:
    logging.debug('No cairo.')
    try:
      import aggdraw
      from rdkit.Chem.Draw.aggCanvas import Canvas
      RDKIT_useAGG=True
      logging.debug('RDKIT_useAGG=True')
    except:
      logging.debug('No aggdraw.')
      from rdkit.Chem.Draw.spingCanvas import Canvas
      RDKIT_useSPING=True
      logging.debug('RDKIT_useSPING=True')
  return Canvas

#############################################################################
def Mol2Image(mol, width=300, height=300, kekulize=True, wedgeBonds=True, highlightAtoms=None):
  """	We need this custom function because the RDKit function
	rdkit.Chem.Draw.MolToImage() does not allow atom highlighting.
	Aha!  Now (2014) this is possible.
  """
  logging.debug('In mol2img_rdk_utils.')
  global Canvas
  try:
    x = Canvas #Check that Canvas defined.
  except:
    Canvas = InitializeCanvas(False)

  if not mol:
    raise ValueError('NULL molecule.') #python3

  img = Image.new("RGBA",(width,height),"white")
  canvas = Canvas(img)
  if RDKIT_useAGG:
    canvas.setantialias(True)
  drawing = rdkit.Chem.Draw.MolDrawing(canvas)

  if kekulize:
    mol = rdkit.Chem.Mol(mol.ToBinary())
    rdkit.Chem.Kekulize(mol)

  drawing.wedgeDashedBonds=wedgeBonds

  #Scaling needed or larger image just has more border.
  #...but fails with error.
  #drawing.scaleAndCenter(mol,mol.GetConformer(0),coordCenter=True)

  if highlightAtoms:
    drawing.AddMol(mol, highlightAtoms=highlightAtoms)
  else:
    drawing.AddMol(mol)
  canvas.flush()

  return img

#############################################################################
