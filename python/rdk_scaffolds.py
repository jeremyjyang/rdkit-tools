#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.MurckoScaffold.html
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.rdScaffoldNetwork.html
rdScaffoldNetwork available RDKit 2020.03.1+.
"""
#############################################################################
import os,sys,re,time,argparse,logging

import rdkit
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Scaffolds import rdScaffoldNetwork

DEMOSMIS = [
	'c1ccc(cc1)N(=O)=O',
	'c1ccc(cc1)[N+](=O)[O-]',
	'C[C@H]1CN(CCCN1[S+]([O-])(=O)C2=CC=CC3=C2C(=CN=C3)C)C(=O)CN',
	'N[S+]([O-])(=O)C1=C(Cl)C=C2NC(N[S+]([O-])(=O)C2=C1)C(Cl)Cl',
	'C[S+]([O-])C1=CC=C(C=C1)\C=C2\C(=C(\CC(O)=O)C3=C2C=CC(=C3)F)C',
	'O=[N+]([O-])c1cc(S(=O)(=O)N2CCCC2)ccc1NN=Cc1ccccn1'
	]

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
  for key in ('flattenChirality', 'flattenIsotopes', 'flattenKeepLargest', 'includeGenericBondScaffolds', 'includeGenericScaffolds', 'includeScaffoldsWithAttachments', 'includeScaffoldsWithoutAttachments', 'keepOnlyFirstFragment', 'pruneBeforeFragmenting'):
    logging.info(f"{key}: {params.__getattribute__(key)}") 
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
  for smi in DEMOSMIS:
    mol1 = MolFromSmiles(smi)
    mol2 = MurckoScaffold.GetScaffoldForMol(mol1) if mol1 else None
    smi_std = MolToSmiles(mol2, isomericSmiles=False) if mol2 else None
    logging.info(f"{smi:>28s} >> {smi_std}")

#############################################################################
def DemoNet(brics):
  mols = [MolFromSmiles(smi) for smi in DEMOSMIS]
  Mols2ScafNet(mols, brics)

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
  OPS = ["bmscaf", "scafnet", "demobm", "demonet"]
  parser.add_argument("op", choices=OPS, default="mol2scaf", help="OPERATION")
  parser.add_argument("--i", dest="ifile", help="input file, TSV or SDF")
  parser.add_argument("--o", dest="ofile", help="output file, TSV or SDF")
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

  logging.info(f"RDKit version: {rdkit.__version__}")

  t0=time.time()

  if args.op=="demobm":
    DemoBM()
    sys.exit()

  elif args.op=="demonet":
    DemoNet(args.brics)
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

