#!/usr/bin/env python3
#############################################################################
import sys,os,platform,re,argparse,logging
import tempfile
import rdkit
import PIL

from .. import depict

#############################################################################
def Demo():
  with tempfile.NamedTemporaryFile("w+b", suffix=".smi", delete=False) as ftmp:
    ftmp.write(b"NCCc1ccc(O)c(O)c1 dopamine\n")
    fname = ftmp.name
  imgs = depict.ReadMols2Images(fname, "SMI", 400, 300, True, False) 
  os.remove(fname)
  imgs[0].show()

#############################################################################
def Demo2():
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

  smarts = 'a:a'
  pat = rdkit.Chem.MolFromSmarts(smarts)
  matches = mol.GetSubstructMatches(pat, uniquify=True, useChirality=False)

  hitatoms=[]
  for match in matches:
    hitatoms.extend(match)
    hitatoms.sort()  ##repeats ok?

  width = 400
  height = 300

  img = rdkit.Chem.Draw.MolToImage(mol, size=(width, height), kekulize=True, highlightAtoms=hitatoms, wedgeBonds=True)
  img.show()

  #strio = io.StringIO()
  #img.save(strio, format='PNG')
  #imgbytes = strio.getvalue()
  #strio.close()

  #fpath = "/tmp/z.png"
  #f = open(fpath, 'wb')
  #f.write(imgbytes)
  #f.close()
  #logging.info(f'Image written to: {fpath}')
  #img.save(fpath) ##also works ok
  #img.save(sys.stdout, 'PNG') ##also works ok

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='RDKit molecule depiction utility', epilog='Modes: single = one image; batch = multiple images; pdf = multi-page')
  ifmts = ['AUTO','SMI','MDL']
  ofmts = ['PNG','JPEG','PDF']
  ops = ['single', 'batch', 'pdf', 'demo', 'demo2']
  parser.add_argument("op", choices=ops, default='single', help='OPERATION')
  parser.add_argument("--i", dest="ifile", help="input molecule file")
  parser.add_argument("--ifmt", choices=ifmts, help='input file format', default='AUTO')
  parser.add_argument("--ofmt", choices=ofmts, help='output file format', default='PNG')
  parser.add_argument("--smilesColumn", type=int, default=0,  help='')
  parser.add_argument("--nameColumn", type=int, default=1,  help='')
  parser.add_argument("--header", action="store_true", help="SMILES/TSV file has header")
  parser.add_argument("--delim", default=" \t", help="SMILES/TSV field delimiter")
  parser.add_argument("--height", type=int, help='height of image', default=120)
  parser.add_argument("--width", type=int, help='width of image', default=140)
  parser.add_argument("--kekulize", action="store_true", help="display Kekule form")
  parser.add_argument("--wedgebonds", action="store_true", help="stereo wedge bonds")
  parser.add_argument("--pdf_title", help="PDF doc title")
  parser.add_argument("--batch_dir", help="destination for batch files", default='/tmp')
  parser.add_argument("--batch_prefix", help="prefix for batch files", default="RDKDEPICT")
  parser.add_argument("--o", dest="ofile", help="output file")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f"Python version: {platform.python_version()}")
  logging.info(f"RDKit version: {rdkit.rdBase.rdkitVersion}")
  logging.info(f"Boost version: {rdkit.rdBase.boostVersion}")
  logging.info(f"PIL_PILLOW version: {PIL.__version__}")

  if args.op == 'batch':
    if not args.ifile: parser.error('Input file required.')
    imgs = depict.ReadMols2Images(args.ifile, args.ifmt, args.smilesColumn, args.nameColumn, args.delim, args.header, args.width, args.height, args.kekulize, args.wedgebonds) 
    if not args.batch_dir:
      parser.error('--batch_dir required for --batch_mode.')
      parser.print_help()
    if not os.access(args.batch_dir, os.W_OK):
      parser.error(f"batch_dir {args.batch_dir} non-existent or non-writeable.")
      parser.print_help()
    depict.WriteImages2ImageFiles(imgs, args.ofmt, args.batch_dir, args.batch_prefix)

  elif args.op == 'pdf':
    if not args.ifile: parser.error('Input file required.')
    imgs = depict.ReadMols2Images(args.ifile, args.ifmt, args.smilesColumn, args.nameColumn, args.delim, args.header, args.width, args.height, args.kekulize, args.wedgebonds) 
    title = args.pdf_title if args.pdf_title else os.path.basename(args.ifile)
    depict.WriteImages2PDFFile(imgs, args.width, args.height, title, args.ofile)

  elif args.op == 'single':
    if not args.ifile: parser.error('Input file required.')
    imgs = depict.ReadMols2Images(args.ifile, args.ifmt, args.smilesColumn, args.nameColumn, args.delim, args.header, args.width, args.height, args.kekulize, args.wedgebonds) 
    depict.WriteImage2ImageFile(imgs[0], args.ofmt, args.ofile)

  elif args.op == 'demo':
    Demo()

  elif args.op == 'demo2':
    Demo2()

  else:
    parser.print_help()

