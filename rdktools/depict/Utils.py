#!/usr/bin/env python3
#############################################################################
import sys,os,platform,re,logging
import pandas as pd
import PIL
from PIL import Image

import rdkit.Chem
import rdkit.Chem.Draw
import rdkit.rdBase
import rdkit.Chem.AllChem

from ..util import latex as util_latex

PDF_DOCTYPES = ["usletter", "letterpaper", "legalpaper", "a4paper", "a5paper", "b5paper"]

#############################################################################
def SelectMolsupplier(ifile, ifmt, smiCol, namCol, delim, header):
  ifmt = SelectFormat(ifmt, ifile)
  molSupplier=None
  if ifmt=='MDL':
    molSupplier = rdkit.Chem.rdmolfiles.SDMolSupplier(ifile)
  elif ifmt=='SMI':
    molSupplier = rdkit.Chem.rdmolfiles.SmilesMolSupplier(ifile, delimiter=delim, smilesColumn=smiCol, nameColumn=namCol, titleLine=header, sanitize=True)
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
def ReadMols(ifile, ifmt, smiCol, namCol, delim, header, parse_as_smarts):
  #logging.debug(f'ReadMols parse_as_smarts={parse_as_smarts}')
  if parse_as_smarts:
    return ReadMolsFromSmarts(ifile, smiCol, namCol, delim, header)
  molSupplier = SelectMolsupplier(ifile, ifmt, smiCol, namCol, delim, header)
  mols=[]; i_mol=0;
  for mol in molSupplier:
    i_mol+=1
    molname = mol.GetProp('_Name') if mol else ""
    logging.debug(f'{i_mol}. "{molname}"')
    mols.append(mol)
  logging.debug(f'mols read: {i_mol}')
  return mols

#############################################################################
def ReadMolsFromSmarts(ifile, smiCol, namCol, delim, header):
  mols=[]; i_mol=0;
  with pd.read_csv(ifile, sep=delim, header=(0 if header else None), iterator=True, chunksize=10, engine="python") as reader:
    for df_this in reader:
      for i in range(df_this.shape[0]):
        i_mol+=1
        smarts = df_this.iloc[i,smiCol]
        molname = df_this.iloc[i,namCol]
        try:
          mol = rdkit.Chem.rdmolfiles.MolFromSmarts(smarts)
        except Exception as e:
          logging.error(f'{i_mol}. "{molname}"; {e}')
          mol = rdkit.Chem.Mol() #empty mol
        if mol is None: mol = rdkit.Chem.Mol() #empty mol
        logging.debug(f"type(mol): {type(mol)}")
        logging.debug(f"type(molname): {type(molname)}")
        mol.SetProp('_Name', str(molname) if bool(molname) else f"mol_{i_mol}")
        logging.debug(f'{i_mol}. "{molname}"')
        mols.append(mol)
  logging.debug(f'mols read (from SMARTS): {i_mol}')
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
def ReadMols2Images(ifile, ifmt, smiCol, namCol, delim, header, kekulize, wedgebonds, parse_as_smarts, width, height):
  mols = ReadMols(ifile, ifmt, smiCol, namCol, delim, header, parse_as_smarts)
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
def WriteImages2PDFFile(ifile, ifmt, smilesColumn, nameColumn, delim, header, kekulize, wedgebonds, parse_as_smarts, grid_width, grid_height, nPerRow, nPerCol, doctype, landscape, title, ofile):
  '''Converts PIL Image objects to pylatex StandAloneGraphics objects, then
  write LaTeX tex and pdf with grid of depictions.'''
  img_width = int(grid_width/nPerRow)
  img_height = int(grid_height/nPerCol)
  logging.debug(f"grid_width={grid_width}; grid_height={grid_height}; w={img_width}; h={img_height}; nPerRow={nPerRow}; nPerCol={nPerCol}")
  imgs = ReadMols2Images(ifile, ifmt, smilesColumn, nameColumn, delim, header, kekulize, wedgebonds, parse_as_smarts, img_width, img_height) 
  sagps = [util_latex.Utils.PILImage2SAGPlus(img) for img in imgs]
  util_latex.WriteImageGrid(sagps, 
	img_width, img_height, 
	nPerRow, nPerCol, 
	doctype, landscape,
	title, ofile)

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
