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
def Mols2BMScaffolds(molReader, molWriter):
  n_mol=0; 
  for mol in molReader:
    n_mol+=1
    molname = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    logging.debug('%d. %s:'%(n_mol, molname))
    scafmol = MurckoScaffold.GetScaffoldForMol(mol)
    molWriter.write(scafmol)
  logging.info(f'{n_mol} mols written to {}'.format(str(molWriter))

#############################################################################
def Mols2ScafNet(molReader, ofile):
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
  fout = open(ofile, "w") if ofile else sys.stdout
  for i in range(len(scafnet.nodes)):
    fout.write("node\t{}\t{}\t{}\n".format(i, scafnet.nodes[i], scafnet.counts[i]))
  for i in range(len(scafnet.edges)):
    fout.write("edge\t{}\t{}\t{}\t{}\n".format(i, scafnet.edges[i].beginIdx, scafnet.edges[i].endIdx, scafnet.edges[i].type))
  fout.flush()
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

  if not (args.ifile): parser.error('--i required.')

  if args.op=="bmscaf":
    molReader = File2Molreader(args.ifile, args.idelim, args.smicol, args.namcol, args.iheader)
    molWriter = File2Molwriter(args.ofile, args.odelim, args.oheader)
    Mols2BMScaffolds(molReader, molWriter)

  elif args.op=="scafnet":
    molReader = File2Molreader(args.ifile, args.idelim, args.smicol, args.namcol, args.iheader)
    Mols2ScafNet(molReader, args.ofile)

  else:
    parser.error(f"Unsupported operation: {args.op}")

  logging.info('Elapsed: %s'%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))))

