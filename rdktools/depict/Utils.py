#!/usr/bin/env python3
#############################################################################
import sys,os,platform,re,logging
import StringIO
import PIL
from PIL import Image

import rdkit.Chem
import rdkit.Chem.Draw
import rdkit.rdBase
import rdkit.Chem.AllChem

from .. import util

#############################################################################
def SelectMolsupplier(ifmt, ifile):
  ifmt = SelectFormat(ifmt, ifile)
  molSupplier=None
  if ifmt=='MDL':
    molSupplier = rdkit.Chem.rdmolfiles.SDMolSupplier(ifile)
  elif ifmt=='SMI':
    fin = open(ifile)
    fdata = fin.read()
    fin.close()
    molSupplier = rdkit.Chem.rdmolfiles.SmilesMolSupplier()
    ### Fails to read 1st SMILES without "title line".
    ### Also, spaces in names not handled.
    molSupplier.SetData('smiles name\n'+fdata) ##Kludge
  else:
    logging.info(f'Unknown format: {ifmt}')
  return molSupplier

#############################################################################
def SelectFormat(ifmt, ifile):
  fext = os.path.splitext(ifile)[1]
  if ifmt.upper() == 'AUTO' and fext:
    if fext.upper() in ('.MDL','.SDF','.MOL'):
      ifmt='MDL'
    elif fext.upper() in ('.SMI','.SMILES','.ISM'):
      ifmt='SMI'
  return ifmt

#############################################################################
def ReadMols(ifile,ifmt):
  molSupplier = SelectMolsupplier(ifmt,ifile)
  mols=[];
  i_mol=0;
  for mol in molSupplier:
    i_mol+=1
    molname = mol.GetProp('_Name') if mol else ""
    logging.debug(f'{i_mol}. "{molname}"')
    mols.append(mol)
  logging.debug(f'mols read: {i_mol}')
  return mols

#############################################################################
def Mols2Images(mols, width, height, kekulize, wedgeBonds):
  if not mols:
    logging.info('ERROR: no mols.')
    return []
  imgs=[]
  img_blank = PIL.Image.new('RGB', (width, height))
  img_blank.info['name'] = ""
  for mol in mols:
    if not mol:
      imgs.append(img_blank)
      continue
    rdkit.Chem.AllChem.Compute2DCoords(mol, clearConfs=True)
    #Fails with RDKit 2017.03.1
    #img = rdkit.Chem.Draw.MolToImage(mol, size=(width, height), kekulize=kekulize, wedgeBonds=wedgeBonds, fitImage=True)
    img = Mol2Image(mol, width, height, kekulize, wedgeBonds, None)
    img.info['name']=mol.GetProp('_Name')
    imgs.append(img)
  return imgs

#############################################################################
def ReadMols2Images(ifile, ifmt, width, height, kekulize, wedgebonds):
  mols = ReadMols(ifile, ifmt)
  return Mols2Images(mols,  width, height, kekulize, wedgebonds)

#############################################################################
def WriteImage2ImageFile(img, ofmt, ofile):
  logging.debug(f"Image name: {img.info['name']}")
  logging.debug(f"Image mode: {img.mode}; format: {img.format}; size: {img.width}x{img.height}")
  logging.debug(f"Output file: {ofile} ; format: {ofmt}")
  img.save(ofile, format=ofmt)

#############################################################################
def WriteImages2ImageFiles(imgs, ofmt, batch_dir, batch_prefix):
  i_img=0;
  for img in imgs:
    i_img+=1
    ofile = (f"{batch_dir}/{batch_prefix}_{i_img:02d}.{ofmt.lower()}")
    WriteImage2ImageFile(img, ofmt, ofile)

#############################################################################
def WriteImages2PDFFile(imgs, img_width, img_height, title, ofile):
  '''Converts PIL Image objects to pylatex StandAloneGraphics objects, then
  write LaTeX tex and pdf with grid of depictions.'''
  sagps = [util.latex.PILImage2SAGPlus(img) for img in imgs]
  npr = int(440/img_width)
  npc = int(560/img_height)
  logging.debug(f"w={img_width} ; h={img_height} ; npr={npr} ; npc={npc}")
  util.latex.WriteImageGrid(sagps, img_width, img_height, npr, npc, title, ofile)

#############################################################################
def InitializeCanvas():
  '''Initialize Canvas, depending on library availability.'''
  global RDKIT_useAGG, RDKIT_useCAIRO, RDKIT_useSPING
  RDKIT_useAGG=False
  RDKIT_useCAIRO=False
  RDKIT_useSPING=False

  try:
    import cairo
    from rdkit.Chem.Draw.cairoCanvas import Canvas
    RDKIT_useCAIRO=True
    logging.info('RDKIT_useCAIRO=True')
  except:
    logging.info('No cairo.')
    try:
      import aggdraw
      from rdkit.Chem.Draw.aggCanvas import Canvas
      RDKIT_useAGG=True
      logging.info('RDKIT_useAGG=True')
    except:
      logging.info('No aggdraw.')
      from rdkit.Chem.Draw.spingCanvas import Canvas
      RDKIT_useSPING=True
      logging.info('RDKIT_useSPING=True')
  return Canvas

#############################################################################
def Mol2Image(mol, width=300, height=300, kekulize=True, wedgeBonds=True, highlightAtoms=None):
  """	We need this custom function because the RDKit function
	rdkit.Chem.Draw.MolToImage() does not allow atom highlighting.
	Aha!  Now (2014) this is possible.
  """
  global Canvas
  try:
    x = Canvas #Check that Canvas defined.
  except:
    Canvas = InitializeCanvas()

  if not mol:
    raise ValueError('NULL molecule.') #python3

  img = Image.new("RGBA", (width, height), "white")
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

  sdms = rdkit.Chem.AllChem.SDMolSupplier()
  sdms.SetData(mdlstr, sanitize=True, removeHs=True)
  mol = sdms[0]

  logging.info(mol.GetProp('_Name'))
  mol_kek = rdkit.Chem.AllChem.Mol(mol)
  rdkit.Chem.Kekulize(mol_kek)
  logging.info(rdkit.Chem.MolToSmiles(mol_kek, kekuleSmiles=True))
  logging.info(rdkit.Chem.MolToSmiles(mol, isomericSmiles=True))
  logging.info(f'\t{mol.GetNumAtoms()} atoms, {mol.GetNumBonds()} bonds')
  for propname in mol.GetPropNames():
    logging.info(f'\tProperty ["{propname}"]: {mol.GetProp(propname)}')

  logging.info(f'NumConfs: {mol.GetNumConformers()}')
  rdkit.Chem.AllChem.Compute2DCoords(mol, clearConfs=False)
  logging.info(f'NumConfs: {mol.GetNumConformers()}')

  smarts='a:a'
  pat=rdkit.Chem.MolFromSmarts(smarts)
  matches=mol.GetSubstructMatches(pat, uniquify=True, useChirality=False)

  hitatoms=[]
  for match in matches:
    hitatoms.extend(match)
    hitatoms.sort()  ##repeats ok?

  width = 400
  height = 300

  img = rdkit.Chem.Draw.MolToImage(mol, size=(width, height), kekulize=True, highlightAtoms=hitatoms, wedgeBonds=True)
  logging.info('NumConfs: %d'%mol.GetNumConformers())

  strio = StringIO.StringIO()
  img.save(strio, format='PNG')
  imgbytes = strio.getvalue()
  strio.close()
  f = open(fpath, 'wb')
  f.write(imgbytes)
  f.close()

  logging.info(f'Image written to: {fpath}')

  #img.save(fpath) ##also works ok
  #img.save(sys.stdout, 'PNG') ##also works ok

