#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.MurckoScaffold.html
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.rdScaffoldNetwork.html
rdScaffoldNetwork available RDKit 2020.03.1+.
"""
#############################################################################
import os,sys,re,json,time,inspect,argparse,logging
import matplotlib as mpl
import pyvis
from pyvis.network import Network

import rdkit
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles, Draw
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Scaffolds import rdScaffoldNetwork

from .. import scaffold
from .. import util

#############################################################################
def DemoBM():
  scafmols=[];
  for smi in util.Utils.DEMOSMIS:
    mol = MolFromSmiles(re.sub(r'\s.*$', '', smi))
    scafmol = MurckoScaffold.GetScaffoldForMol(mol) if mol else None
    scafmols.append(scafmol)
    smi_std = MolToSmiles(scafmol, isomericSmiles=False) if scafmol else None
    logging.info(f"{smi:>28s} >> {smi_std}")
  img = rdkit.Chem.Draw.MolsToGridImage(scafmols, molsPerRow=4)
  img.show()

#############################################################################
def DemoNet(brics):
  smi = ('Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12 flucloxacillin')
  mols = [MolFromSmiles(re.sub(r'\s.*$', '', smi))]
  scafnet = scaffold.Utils.Mols2ScafNet(mols, brics)
  scafmols = [MolFromSmiles(m) for m in scafnet.nodes]
  #title="RDKit_ScafNet:"+re.sub(r'^[^\s]*\s+(.*)$', r'\1', smi)) #How to add title?
  img = rdkit.Chem.Draw.MolsToGridImage(scafmols, legends=[f'{i}, counts: {c}' for i,c in enumerate(scafnet.counts)], molsPerRow=4)
  img.show()

#############################################################################
def DemoNetVis(brics, scratchdir, ofile):
  DEMOSMI = ('Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12 flucloxacillin')
  mols = [MolFromSmiles(re.sub(r'\s.*$', '', DEMOSMI))]
  scafnet = scaffold.Utils.Mols2ScafNet(mols, brics)
  if not os.path.isdir(scratchdir): os.mkdir(scratchdir)
    
  g = pyvis.network.Network(notebook=False, height='800px', width='1000px', heading="RDKit_ScafNet:"+re.sub(r'^[^\s]*\s+(.*)$', r'\1', DEMOSMI))

  for i,n in enumerate(scafnet.nodes):
    svg = util.Utils.moltosvg(rdkit.Chem.MolFromSmiles(n))
    with open(f'{scratchdir}/{i}.svg', 'w') as outf:
      outf.write(svg)
    g.add_node(i, shape="image", label=' ', image=f'{scratchdir}/{i}.svg', title=svg, size=60)
  for e in scafnet.edges:
    g.add_edge(e.beginIdx, e.endIdx, label=str(e.type))

  VisJS_options = {
    "edges": {
     "font":{
     "size":20
     }
    },
    "nodes": {
      "font": {
        "color": "rgba(214,47,66,1)",
        "size": 16,
        "face": "tahoma"
      }
    },
    "physics": {
      "forceAtlas2Based": {
        "gravitationalConstant": -120,
        "springLength": 200,
        "avoidOverlap": 0.42
      },
      "minVelocity": 0.75,
      "solver": "forceAtlas2Based"
    }
  }
  g.set_options(options=json.dumps(VisJS_options))
  #g.show_buttons()
  logging.info(f"Writing SCAFNET to: {ofile}")
  g.show(ofile)

#############################################################################
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="RDKit scaffold analysis", epilog="")
  OPS = ["bmscaf", "scafnet", "demobm", "demonet", "demonetvis"]
  parser.add_argument("op", choices=OPS, default="mol2scaf", help="OPERATION")
  parser.add_argument("--i", dest="ifile", help="input file, TSV or SDF")
  parser.add_argument("--o", dest="ofile", help="output file, TSV|SDF")
  parser.add_argument("--o_html", dest="ofile_html", default="/tmp/rdk_scafnet.html", help="output file, HTML")
  parser.add_argument("--scratchdir", default="/tmp")
  parser.add_argument("--smicol", type=int, default=0, help="SMILES column from TSV (counting from 0)")
  parser.add_argument("--namcol", type=int, default=1, help="name column from TSV (counting from 0)")
  parser.add_argument("--idelim", default="\t", help="delim for input TSV")
  parser.add_argument("--odelim", default="\t", help="delim for output TSV")
  parser.add_argument("--iheader", action="store_true", help="input TSV has header")
  parser.add_argument("--oheader", action="store_true", help="output TSV has header")
  parser.add_argument("--brics", action="store_true", help="BRICS fragmentation rules (Degen, 2008)")
  parser.add_argument("-v", "--verbose", action="count", default=0)

  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f"RDKit version: {rdkit.rdBase.rdkitVersion}")
  logging.info(f"MatplotLib version: {mpl.__version__}")
  logging.info(f"Pyvis version: {pyvis.__version__}")

  t0=time.time()

  if args.op=="demobm":
    DemoBM()
    sys.exit()

  elif args.op=="demonet":
    DemoNet(args.brics)
    sys.exit()

  elif args.op=="demonetvis":
    DemoNetVis(args.brics, args.scratchdir, args.ofile_html)
    sys.exit()

  if not (args.ifile): parser.error('--i required.')

  if args.op=="bmscaf":
    molReader = util.Utils.File2Molreader(args.ifile, args.idelim, args.smicol, args.namcol, args.iheader)
    molWriter = util.Utils.File2Molwriter(args.ofile, args.odelim, args.oheader)
    mols = util.Utils.ReadMols(molReader)
    scaffold.Utils.Mols2BMScaffolds(mols, molWriter)

  elif args.op=="scafnet":
    molReader = util.Utils.File2Molreader(args.ifile, args.idelim, args.smicol, args.namcol, args.iheader)
    mols = util.Utils.ReadMols(molReader)
    scaffold.Utils.Mols2ScafNet(mols, args.brics, args.ofile)

  else:
    parser.error(f"Unsupported operation: {args.op}")

  logging.info(f"""Elapsed: {time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))}""")

