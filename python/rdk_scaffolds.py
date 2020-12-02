#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.MurckoScaffold.html
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.rdScaffoldNetwork.html
rdScaffoldNetwork available RDKit 2020.03.1+.
"""
#############################################################################
import os,sys,re,json,time,inspect,argparse,logging
import matplotlib as mpl
#from matplotlib import pyplot as plt
import pyvis
from pyvis.network import Network

import rdkit
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles, Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Scaffolds import rdScaffoldNetwork

DEMOSMIS = [
'C[C@H]1CN(CCCN1[S+]([O-])(=O)C2=CC=CC3=C2C(=CN=C3)C)C(=O)CN Rho Kinase Inhibitor IV',
'N[S+]([O-])(=O)C1=C(Cl)C=C2NC(N[S+]([O-])(=O)C2=C1)C(Cl)Cl trichlormethiazide',
'C[S+]([O-])C1=CC=C(C=C1)\C=C2\C(=C(\CC(O)=O)C3=C2C=CC(=C3)F)C Sulindac',
'O=[N+]([O-])c1cc(S(=O)(=O)N2CCCC2)ccc1NN=Cc1ccccn1 PUBCHEM_CID:4036736',
'Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12 flucloxacillin',
'CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O ampicillin',
'CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O.[Na] penicillin',
'Cc1onc(-c2ccccc2)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12 oxacillin'
	]

#############################################################################
def moltosvg(mol, molSize=(450,250), kekulize=True):
    mc = rdMolDraw2D.PrepareMolForDrawing(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
    opts = drawer.drawOptions()
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg.replace('svg:', '')

#############################################################################
def ReadMols(molReader):
  mols=[]; 
  for mol in molReader:
    if mol is not None: mols.append(mol)
  logging.info(f'{len(mols)} mols read from {molReader}')
  return mols

#############################################################################
def Mols2BMScaffolds(mols, molWriter):
  for i,mol in enumerate(mols):
    molname = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    logging.debug(f'{i+1}. {molname}')
    scafmol = MurckoScaffold.GetScaffoldForMol(mol)
    molWriter.write(scafmol)
  logging.info(f'{len(mols)} mols written to {molWriter}')

#############################################################################
def Mols2ScafNet(mols, brics=False, ofile=None):
  if brics:
    params = rdScaffoldNetwork.BRICSScaffoldParams()
  else:
    params = rdScaffoldNetwork.ScaffoldNetworkParams()
    params.flattenChirality = True
    params.flattenIsotopes = True
    params.flattenKeepLargest = True
    params.includeGenericBondScaffolds = False
    params.includeGenericScaffolds = False
    params.includeScaffoldsWithAttachments = False
    params.includeScaffoldsWithoutAttachments = True
    params.keepOnlyFirstFragment = False
    params.pruneBeforeFragmenting = True

  attrs = [a for a in inspect.getmembers(params) if not(a[0].startswith('__'))]
  for a in attrs:
    logging.info(f"{a[0]}: {a[1]}")

  for i,mol in enumerate(mols):
    molname = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    logging.debug(f'{i+1}. {molname}:')

  scafnet = rdScaffoldNetwork.CreateScaffoldNetwork(mols, params)
  fout = open(ofile, "w") if ofile else sys.stdout
  for i in range(len(scafnet.nodes)):
    fout.write(f"node\t{i}\t{scafnet.nodes[i]}\t{scafnet.counts[i]}\n")
  for i in range(len(scafnet.edges)):
    fout.write(f"edge\t{i}\t{scafnet.edges[i].beginIdx}\t{scafnet.edges[i].endIdx}\t{scafnet.edges[i].type}\n")
  fout.flush()
  logging.info(f"nodes: {len(scafnet.nodes)}; edges:{len(scafnet.edges)}")
  return scafnet

#############################################################################
def DemoBM():
  scafmols=[];
  for smi in DEMOSMIS:
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
  scafnet = Mols2ScafNet(mols, brics)
  scafmols = [MolFromSmiles(m) for m in scafnet.nodes]
  #title="RDKit_ScafNet:"+re.sub(r'^[^\s]*\s+(.*)$', r'\1', smi)) #How to add title?
  img = rdkit.Chem.Draw.MolsToGridImage(scafmols, legends=[f'{i}, counts: {c}' for i,c in enumerate(scafnet.counts)], molsPerRow=4)
  img.show()

#############################################################################
def DemoNetVis(brics, scratchdir, ofile):
  DEMOSMI = ('Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12 flucloxacillin')
  mols = [MolFromSmiles(re.sub(r'\s.*$', '', DEMOSMI))]
  scafnet = Mols2ScafNet(mols, brics)
  if not os.path.isdir(scratchdir): os.mkdir(scratchdir)
    
  g = pyvis.network.Network(notebook=False, height='800px', width='1000px', heading="RDKit_ScafNet:"+re.sub(r'^[^\s]*\s+(.*)$', r'\1', DEMOSMI))

  for i,n in enumerate(scafnet.nodes):
    svg = moltosvg(rdkit.Chem.MolFromSmiles(n))
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
def File2Molreader(ifile, idelim, smicol, namcol, iheader):
  if re.sub(r'.*\.', '', ifile).lower()in ('smi', 'smiles', 'csv', 'tsv'):
    molReader = SmilesMolSupplier(ifile, delimiter=idelim, smilesColumn=smicol, nameColumn=namcol, titleLine=iheader, sanitize=True)
  elif re.sub(r'.*\.', '', ifile).lower() in ('sdf','sd','mdl','mol'):
    molReader = SDMolSupplier(ifile, sanitize=True, removeHs=True)
  else:
    molReader = None
    logging.error(f'Invalid file extension: {ifile}')
  return molReader

#############################################################################
def File2Molwriter(ofile, odelim, oheader):
  if not ofile:
    molWriter = SmilesWriter("-", delimiter=odelim, nameHeader='Name', includeHeader=oheader, isomericSmiles=True, kekuleSmiles=False)
  elif re.sub(r'.*\.', '', ofile).lower() in ('sdf','sd','mdl','mol'):
    molWriter = SDWriter(ofile)
  elif re.sub(r'.*\.', '', ofile).lower() in ('smi', 'smiles', 'csv', 'tsv'):
    molWriter = SmilesWriter(ofile, delimiter=odelim, nameHeader='Name', includeHeader=oheader, isomericSmiles=True, kekuleSmiles=False)
  else:
    logging.error(f'Invalid file extension: {ofile}')
    molWriter = None
  return molWriter

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
    molReader = File2Molreader(args.ifile, args.idelim, args.smicol, args.namcol, args.iheader)
    molWriter = File2Molwriter(args.ofile, args.odelim, args.oheader)
    mols = ReadMols(molReader)
    Mols2BMScaffolds(mols, molWriter)

  elif args.op=="scafnet":
    molReader = File2Molreader(args.ifile, args.idelim, args.smicol, args.namcol, args.iheader)
    mols = ReadMols(molReader)
    Mols2ScafNet(mols, args.brics, args.ofile)

  else:
    parser.error(f"Unsupported operation: {args.op}")

  logging.info('Elapsed: %s'%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))))

