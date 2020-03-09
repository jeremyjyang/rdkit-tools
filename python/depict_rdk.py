#!/usr/bin/env python3
#############################################################################
### depict_rdk.py
###
### Jeremy J Yang
### 10 Jan 2018
#############################################################################
import sys,os,platform,re,argparse
import PIL

import rdkit.rdBase
import rdkit.Chem
import rdkit.Chem.Draw
import rdkit.Chem.AllChem

import mol2img_rdk_utils
import pylatex_utils

#############################################################################
def SelectMolsupplier(ifmt,ifile):
  ifmt = SelectFormat(ifmt,ifile)
  molSupplier=None
  if ifmt=='MDL':
    molSupplier=rdkit.Chem.rdmolfiles.SDMolSupplier(ifile)
  elif ifmt=='SMI':
    fin = open(ifile)
    fdata = fin.read()
    fin.close()
    molSupplier=rdkit.Chem.rdmolfiles.SmilesMolSupplier()
    ### Fails to read 1st SMILES without "title line".
    ### Also, spaces in names not handled.
    molSupplier.SetData('smiles name\n'+fdata) ##Kludge
  else:
    print >>sys.stderr, ('Unknown format: %s'%(ifmt))
  return molSupplier

#############################################################################
def SelectFormat(ifmt,ifile):
  fext=os.path.splitext(ifile)[1]
  if ifmt.upper() == 'AUTO' and fext:
    if fext.upper() in ('.MDL','.SDF','.MOL'):
      ifmt='MDL'
    elif fext.upper() in ('.SMI','.SMILES','.ISM'):
      ifmt='SMI'
  return ifmt

#############################################################################
def ReadMols(ifile,ifmt,verbose):
  molSupplier = SelectMolsupplier(ifmt,ifile)
  mols=[];
  i_mol=0;
  for mol in molSupplier:
    i_mol+=1
    molname = mol.GetProp('_Name') if mol else ""
    if verbose>1:
      print >>sys.stderr, '%d. "%s"'%(i_mol,molname)
    mols.append(mol)
  if verbose:
    print >>sys.stderr, 'mols read: %d'%(i_mol)
  return mols

#############################################################################
def Mols2Images(mols,width,height,kekulize,wedgeBonds,verbose):
  if not mols:
    print >>sys.stderr, 'ERROR: no mols.'
    return []
  imgs=[]
  img_blank = PIL.Image.new('RGB', (width, height))
  img_blank.info['name'] = ""
  for mol in mols:
    if not mol:
      imgs.append(img_blank)
      continue
    rdkit.Chem.AllChem.Compute2DCoords(mol,clearConfs=True)
    #Fails with RDKit 2017.03.1
    #img = rdkit.Chem.Draw.MolToImage(mol,size=(width,height),kekulize=kekulize,wedgeBonds=wedgeBonds,fitImage=True)
    img = mol2img_rdk_utils.Mol2Image(mol,width,height,kekulize,wedgeBonds,None,bool(verbose>0))
    img.info['name']=mol.GetProp('_Name')
    imgs.append(img)
  return imgs

#############################################################################
def ReadMols2Images(ifile,ifmt,width,height,kekulize,wedgebonds,verbose):
  mols = ReadMols(ifile,ifmt,verbose)
  return Mols2Images(mols, width,height,kekulize,wedgebonds,verbose)

#############################################################################
def WriteImage2ImageFile(img, ofmt, ofile, verbose):
  if verbose:
    print >>sys.stderr, 'Image name: %s'%(img.info['name'])
    print >>sys.stderr, 'Image mode: %s ; format: %s ; size: %dx%d'%(img.mode, img.format, img.width,img.height)
    print >>sys.stderr, 'Output file: %s ; format: %s'%(ofile, ofmt)
  img.save(ofile, format=ofmt)

#############################################################################
def WriteImages2ImageFiles(imgs, ofmt, batch_dir, batch_prefix, verbose):
  i_img=0;
  for img in imgs:
    i_img+=1
    ofile = '%s/%s_%02d.%s'%(batch_dir,batch_prefix,i_img,ofmt.lower())
    WriteImage2ImageFile(img, ofmt, ofile, verbose)

#############################################################################
def WriteImages2PDFFile(imgs, img_width, img_height, title, ofile, verbose):
  '''Converts PIL Image objects to pylatex StandAloneGraphics objects, then
  write LaTeX tex and pdf with grid of depictions.'''
  sagps = [pylatex_utils.PILImage2SAGPlus(img) for img in imgs]
  npr = int(440/img_width)
  npc = int(560/img_height)
  print >>sys.stderr, 'DEBUG: w=%d ; h=%d ; npr=%d ; npc=%d'%(img_width,img_height,npr,npc)
  pylatex_utils.WriteImageGrid(sagps,img_width,img_height,npr,npc,title,ofile,verbose)

#############################################################################
if __name__=='__main__':

  PROG=os.path.basename(sys.argv[0])

  parser = argparse.ArgumentParser(
        description='RDKit molecule depiction utility', epilog='Modes: single = one image ; batch = multiple images ; pdf = multi-page')
  ifmts = ['AUTO','SMI','MDL']
  ofmts = ['PNG','JPEG','PDF']
  ops = ['single', 'batch', 'pdf']
  parser.add_argument("op",choices=ops,default='single',help='operation')
  parser.add_argument("--i",dest="ifile",help="input molecule file")
  parser.add_argument("--ifmt",choices=ifmts,help='input file format',default='AUTO')
  parser.add_argument("--ofmt",choices=ofmts,help='output file format',default='PNG')
  parser.add_argument("--height",type=int,help='height of image',default=120)
  parser.add_argument("--width",type=int,help='width of image',default=140)
  parser.add_argument("--kekulize",action="store_true",help="display Kekule form")
  parser.add_argument("--wedgebonds",action="store_true",help="stereo wedge bonds")
  parser.add_argument("--pdf_title",help="PDF doc title")
  parser.add_argument("--batch_dir",help="destination for batch files",default='/tmp')
  parser.add_argument("--batch_prefix",help="prefix for batch files",default=re.sub('\.py$','',PROG))
  parser.add_argument("--o",dest="ofile",help="output file")
  parser.add_argument("-v","--verbose",action="count")

  args = parser.parse_args()
  if not args.verbose: args.verbose=0

  if args.verbose:
    print >>sys.stderr, 'RDKit version: %s'%rdkit.rdBase.rdkitVersion
    print >>sys.stderr, 'Boost version: %s'%rdkit.rdBase.boostVersion
    print >>sys.stderr, 'Python version: %s'%platform.python_version()
    print >>sys.stderr, 'PIL version: %s'%PIL.VERSION
    print >>sys.stderr, 'PIL_PILLOW version: %s'%PIL.PILLOW_VERSION

  if args.op == 'batch':
    if not args.ifile: parser.error('Input file required.')
    imgs = ReadMols2Images(args.ifile,args.ifmt,args.width,args.height,args.kekulize,args.wedgebonds,args.verbose) 
    if not args.batch_dir:
      parser.error('--batch_dir required for --batch_mode.')
      parser.print_help()
    if not os.access(args.batch_dir, os.W_OK):
      parser.error('batch_dir %s non-existent or non-writeable.'%args.batch_dir)
      parser.print_help()
    WriteImages2ImageFiles(imgs, args.ofmt, args.batch_dir, args.batch_prefix, args.verbose)

  elif args.op == 'pdf':
    if not args.ifile: parser.error('Input file required.')
    imgs = ReadMols2Images(args.ifile,args.ifmt,args.width,args.height,args.kekulize,args.wedgebonds,args.verbose) 
    title = args.pdf_title if args.pdf_title else os.path.basename(args.ifile)
    WriteImages2PDFFile(imgs, args.width, args.height, title, args.ofile, args.verbose)

  elif args.op == 'single':
    if not args.ifile: parser.error('Input file required.')
    imgs = ReadMols2Images(args.ifile,args.ifmt,args.width,args.height,args.kekulize,args.wedgebonds,args.verbose) 
    WriteImage2ImageFile(imgs[0], args.ofmt, args.ofile, args.verbose)

  else:
    parser.print_help()

