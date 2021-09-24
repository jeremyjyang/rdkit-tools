#!/usr/bin/env python3
#############################################################################
import sys,os,platform,re,argparse,logging
import rdkit
import PIL

from .. import depict

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='RDKit molecule depiction utility', epilog='Modes: single = one image; batch = multiple images; pdf = multi-page')
  ifmts = ['AUTO','SMI','MDL']
  ofmts = ['PNG','JPEG','PDF']
  ops = ['single', 'batch', 'pdf', 'demo']
  parser.add_argument("op", choices=ops, default='single', help='OPERATION')
  parser.add_argument("--i", dest="ifile", help="input molecule file")
  parser.add_argument("--ifmt", choices=ifmts, help='input file format', default='AUTO')
  parser.add_argument("--ofmt", choices=ofmts, help='output file format', default='PNG')
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
    imgs = depict.Utils.ReadMols2Images(args.ifile, args.ifmt, args.width, args.height, args.kekulize, args.wedgebonds) 
    if not args.batch_dir:
      parser.error('--batch_dir required for --batch_mode.')
      parser.print_help()
    if not os.access(args.batch_dir, os.W_OK):
      parser.error(f"batch_dir {args.batch_dir} non-existent or non-writeable.")
      parser.print_help()
    depict.Utils.WriteImages2ImageFiles(imgs, args.ofmt, args.batch_dir, args.batch_prefix)

  elif args.op == 'pdf':
    if not args.ifile: parser.error('Input file required.')
    imgs = depict.Utils.ReadMols2Images(args.ifile, args.ifmt, args.width, args.height, args.kekulize, args.wedgebonds) 
    title = args.pdf_title if args.pdf_title else os.path.basename(args.ifile)
    depict.Utils.WriteImages2PDFFile(imgs, args.width, args.height, title, args.ofile)

  elif args.op == 'single':
    if not args.ifile: parser.error('Input file required.')
    imgs = depict.Utils.ReadMols2Images(args.ifile, args.ifmt, args.width, args.height, args.kekulize, args.wedgebonds) 
    depict.Utils.WriteImage2ImageFile(imgs[0], args.ofmt, args.ofile)

  elif args.op == 'demo':
    import tempfile
    with tempfile.NamedTemporaryFile("w+b", suffix=".smi", delete=False) as ftmp:
      ftmp.write(b"NCCc1ccc(O)c(O)c1 dopamine\n")
      fname = ftmp.name
    imgs = depict.Utils.ReadMols2Images(fname, "SMI", 400, 300, True, False) 
    os.remove(fname)
    imgs[0].show()

  else:
    parser.print_help()

