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

#############################################################################
def Mol2BMScaffold(mol):
  scafmol = MurckoScaffold.GetScaffoldForMol(mol)
  return scafmol

#############################################################################
def Mols2BMScaffolds(molReader, molWriter):
  n_mol=0; 
  for mol in molReader:
    n_mol+=1
    molname = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    logging.debug('%d. %s:'%(n_mol, molname))
    scafmol = Mol2BMScaffold(mol)
    molWriter.write(scafmol)
  #logging.info(f'{n_mol} mols written to {args.ofile}')

#############################################################################
def Mols2ScafNet(molReader, molWriter):
  n_mol=0; 
  params = rdScaffoldNetwork.BRICSScaffoldParams()
  for key in ('flattenChirality', 'flattenIsotopes', 'flattenKeepLargest', 'includeGenericBondScaffolds', 'includeGenericScaffolds', 'includeScaffoldsWithAttachments', 'includeScaffoldsWithoutAttachments', 'keepOnlyFirstFragment', 'pruneBeforeFragmenting'):
    logging.info("{}: {}".format(key, params.__getattribute__(key)))
  mols = []
  for mol in molReader:
    n_mol+=1
    molname = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    logging.debug('%d. %s:'%(n_mol, molname))
    mols.append(mol)

  scafnet = rdScaffoldNetwork.CreateScaffoldNetwork(mols, params)
  logging.info("nodes: {}; edges:{}".format(len(scafnet.nodes), len(scafnet.edges)))

#############################################################################
def Demo():
  smis = [
	'c1ccc(cc1)N(=O)=O',
	'c1ccc(cc1)[N+](=O)[O-]',
	'C[C@H]1CN(CCCN1[S+]([O-])(=O)C2=CC=CC3=C2C(=CN=C3)C)C(=O)CN',
	'N[S+]([O-])(=O)C1=C(Cl)C=C2NC(N[S+]([O-])(=O)C2=C1)C(Cl)Cl',
	'C[S+]([O-])C1=CC=C(C=C1)\C=C2\C(=C(\CC(O)=O)C3=C2C=CC(=C3)F)C',
	'O=[N+]([O-])c1cc(S(=O)(=O)N2CCCC2)ccc1NN=Cc1ccccn1'
	]
  for smi in smis:
    mol1 = MolFromSmiles(smi)
    mol2 = Mol2BMScaffold(mol1) if mol1 else None
    smi_std = MolToSmiles(mol2, isomericSmiles=False) if mol2 else None
    logging.info(f"{smi:>28s} >> {smi_std}")

#############################################################################
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="RDKit scaffold analysis", epilog="")
  OPS = ["bmscaf", "scafnet", "demo"]
  parser.add_argument("op", choices=OPS, default="mol2scaf", help="OPERATION")
  parser.add_argument("--i", dest="ifile", help="input file, TSV or SDF")
  parser.add_argument("--o", dest="ofile", help="output file, TSV or SDF")
  parser.add_argument("--smicol", type=int, default=0, help="SMILES column from TSV (counting from 0)")
  parser.add_argument("--namcol", type=int, default=1, help="name column from TSV (counting from 0)")
  parser.add_argument("--idelim", default="\t", help="delim for input TSV")
  parser.add_argument("--odelim", default="\t", help="delim for output TSV")
  parser.add_argument("--iheader", action="store_true", help="input TSV has header")
  parser.add_argument("--oheader", action="store_true", help="output TSV has header")
  parser.add_argument("-v", "--verbose", action="count", default=0)

  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f"RDKit version: {rdkit.__version__}")

  t0=time.time()

  if args.op=="demo":
    Demo()
    sys.exit()

  if not (args.ifile and args.ofile): parser.error('--i and --o required.')
  if re.sub(r'.*\.', '', args.ifile).lower()in ('smi', 'smiles', 'csv', 'tsv'):
    molReader = SmilesMolSupplier(args.ifile, delimiter=args.idelim, smilesColumn=args.smicol, nameColumn=args.namcol, titleLine=args.iheader, sanitize=True)
  elif re.sub(r'.*\.', '', args.ifile).lower() in ('sdf','sd','mdl','mol'):
    molReader = SDMolSupplier(args.ifile, sanitize=True, removeHs=True)
  else:
    parser.error('Invalid file extension: {}'.format(args.ifile))

  if re.sub(r'.*\.', '', args.ofile).lower() in ('sdf','sd','mdl','mol'):
    molWriter = SDWriter(args.ofile)
  elif re.sub(r'.*\.', '', args.ofile).lower() in ('smi', 'smiles', 'csv', 'tsv'):
    molWriter = SmilesWriter(args.ofile, delimiter=args.odelim, nameHeader='Name', includeHeader=args.oheader, isomericSmiles=True, kekuleSmiles=False)
  else:
    logging.error(f'Invalid file extension: {args.ofile}')

  if args.op=="bmscaf":
    Mols2BMScaffolds(molReader, molWriter)

  elif args.op=="scafnet":
    Mols2ScafNet(molReader, molWriter)

  else:
    parser.error(f"Unsupported operation: {args.op}")

  logging.info('Elapsed: %s'%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))))

