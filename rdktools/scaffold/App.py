#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.MurckoScaffold.html
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.rdScaffoldNetwork.html
rdScaffoldNetwork available RDKit 2020.03.1+.
"""
#############################################################################
import os,sys,stat,re,json,time,inspect,argparse,logging
import matplotlib as mpl
import pyvis
from pyvis.network import Network

import rdkit
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles, Draw
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Scaffolds import rdScaffoldNetwork

from .. import scaffold
from .. import util

SCRATCHDIR = f"{os.environ['HOME']}/tmp/rdktools"

#############################################################################
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="RDKit scaffold analysis", epilog="")
  OPS = ["bmscaf", "scafnet", "demobm", "demonet_img", "demonet_html"]
  parser.add_argument("op", choices=OPS, default="mol2scaf", help="OPERATION")
  parser.add_argument("--i", dest="ifile", help="input file, TSV or SDF")
  parser.add_argument("--o", dest="ofile", help="output file, TSV|SDF")
  parser.add_argument("--scratchdir", default=SCRATCHDIR)
  parser.add_argument("--o_png", dest="ofile_png", default=f"{SCRATCHDIR}/rdktools_scafnet.png", help="visualization output file, PNG")
  parser.add_argument("--o_html", dest="ofile_html", default=f"{SCRATCHDIR}/rdktools_scafnet.html", help="visualization output file, HTML")
  parser.add_argument("--smilesColumn", type=int, default=0, help="SMILES column from TSV (counting from 0)")
  parser.add_argument("--nameColumn", type=int, default=1, help="name column from TSV (counting from 0)")
  parser.add_argument("--idelim", default="\t", help="delim for input TSV")
  parser.add_argument("--odelim", default="\t", help="delim for output TSV")
  parser.add_argument("--iheader", action="store_true", help="input TSV has header")
  parser.add_argument("--oheader", action="store_true", help="output TSV has header")
  parser.add_argument("--brics", action="store_true", help="BRICS fragmentation rules (Degen, 2008)")
  parser.add_argument("--scafname", default="RDKit Scaffold Analysis", help="title for output")
  parser.add_argument("--display", action="store_true", help="Display scafnet interactively.")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f"RDKit version: {rdkit.rdBase.rdkitVersion}")
  logging.info(f"MatplotLib version: {mpl.__version__}")
  logging.info(f"Pyvis version: {pyvis.__version__}")

  t0=time.time()

  if not os.path.isdir(args.scratchdir): os.mkdir(args.scratchdir)

  if args.op=="demobm":
    scaffold.DemoBM()
    sys.exit()

  elif args.op=="demonet_img":
    scaffold.DemoNetImg(args.scratchdir)
    sys.exit()

  elif args.op=="demonet_html":
    scaffold.DemoNetHtml(args.scratchdir)
    sys.exit()

  if not (args.ifile): parser.error('--i required.')

  if args.op=="bmscaf":
    molReader = util.File2Molreader(args.ifile, args.idelim, args.smilesColumn, args.nameColumn, args.iheader)
    molWriter = util.File2Molwriter(args.ofile, args.odelim, args.oheader)
    mols = util.ReadMols(molReader)
    scafmols = scaffold.Mols2BMScaffolds(mols, molWriter)
    if args.ofile_png:
      img = rdkit.Chem.Draw.MolsToGridImage(scafmols, molsPerRow=8)
      img.save(args.ofile_png, format="PNG")

  elif args.op=="scafnet":
    molReader = util.File2Molreader(args.ifile, args.idelim, args.smilesColumn, args.nameColumn, args.iheader)
    mols = util.ReadMols(molReader)
    scafnet = scaffold.Mols2ScafNet(mols, args.brics, args.ofile)
    if args.ofile_png:
      scaffold.Scafnet2Img(scafnet, args.ofile_png)
    if args.ofile_html:
      scaffold.Scafnet2Html(scafnet, args.scafname, args.scratchdir, args.ofile_html)

  else:
    parser.error(f"Unsupported operation: {args.op}")

  logging.info(f"""Elapsed: {time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))}""")

