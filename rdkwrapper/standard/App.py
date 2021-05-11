#!/usr/bin/env python3
#############################################################################
import os,sys,re,time,argparse,logging

import rdkit
#import rdkit.Chem.AllChem
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles

from .. import standard

#############################################################################
def Demo(norms):
  stdzr = standard.MyStandardizer(norms)
  smis = [
	'CCC(=O)O',
	'c1ccc(cc1)N(=O)=O',
	'c1ccc(cc1)[N+](=O)[O-]',
	'CCC([O-])[OH+]',
	'CCN(=O)([O-])[OH+]',
	'C[C@H]1CN(CCCN1[S+]([O-])(=O)C2=CC=CC3=C2C(=CN=C3)C)C(=O)CN',
	'N[S+]([O-])(=O)C1=C(Cl)C=C2NC(N[S+]([O-])(=O)C2=C1)C(Cl)Cl',
	'C[S+]([O-])C1=CC=C(C=C1)\C=C2\C(=C(\CC(O)=O)C3=C2C=CC(=C3)F)C',
	'O=[N+]([O-])c1cc(S(=O)(=O)N2CCCC2)ccc1NN=Cc1ccccn1',
	'CCC(CC)COC(=O)[C@H](C)N[P@](=O)(OC[C@H]1O[C@](C#N)([C@H](O)[C@@H]1O)C1=CC=C2N1N=CN=C2N)OC1=CC=CC=C1',
	]
  for smi in smis:
    mol1 = MolFromSmiles(smi)
    mol2 = stdzr.standardize(mol1) if mol1 else None
    smi_std = MolToSmiles(mol2, isomericSmiles=True) if mol2 else None
    logging.info("{:>28s} >> {}".format(smi, smi_std))

#############################################################################
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="RDKit chemical standardizer", epilog="")
  OPS = ["standardize", "list_norms", "show_params", "demo"]
  parser.add_argument("op", choices=OPS, default="standardize", help="operation")
  parser.add_argument("--i", dest="ifile", help="input file, SMI or SDF")
  parser.add_argument("--o", dest="ofile", help="output file, SMI or SDF")
  parser.add_argument("--norms", choices=["default", "unm"], default="default", help="normalizations")
  parser.add_argument("--i_norms", dest="ifile_norms", help="input normalizations file, format: SMIRKS<space>NAME")
  parser.add_argument("--remove_isomerism", action="store_true", help="if true, output SMILES isomerism removed")
  parser.add_argument("-v", "--verbose", action="count", default=0)

  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f"RDKit version: {rdkit.__version__}")

  t0=time.time()

  if args.norms=="unm":
    norms = standard.Utils.MyNorms()
  elif args.ifile_norms:
    fin = open(args.ifile_norms)
    norms = standard.Utils.ReadNormsFile(fin)

  else:
    norms = MolStandardize.normalize.NORMALIZATIONS

  if args.op=="list_norms":
    fout = open(args.ofile, "w") if args.ofile else sys.stdout
    standard.Utils.ListNormalizations(norms, fout)

  elif args.op=="standardize":
    if not (args.ifile and args.ofile): parser.error('--i and --o required.')
    if re.sub(r'.*\.', '', args.ifile).lower()in ('smi', 'smiles'):
      molReader = SmilesMolSupplier(args.ifile, delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=True, sanitize=True)
    elif re.sub(r'.*\.', '', args.ifile).lower() in ('sdf','sd','mdl','mol'):
      molReader = SDMolSupplier(args.ifile, sanitize=True, removeHs=True)
    else:
      logging.error(f'Invalid file extension: {args.ifile}')

    if re.sub(r'.*\.', '', args.ofile).lower() in ('sdf','sd','mdl','mol'):
      molWriter = SDWriter(args.ofile)
    elif re.sub(r'.*\.', '', args.ofile).lower()in ('smi', 'smiles'):
      molWriter = SmilesWriter(args.ofile, delimiter='\t', nameHeader='Name',
        includeHeader=True, isomericSmiles=(not args.remove_isomerism), kekuleSmiles=False)
    else:
      logging.error(f'Invalid file extension: {args.ofile}')

    stdzr = standard.Utils.MyStandardizer(norms)
    standard.Utils.Standardize(stdzr, args.remove_isomerism, molReader, molWriter)

  elif args.op=="show_params":
    standard.Utils.ShowParameters()

  elif args.op=="demo":
    Demo(norms)

  else:
    parser.error(f"Unsupported operation: {args.op}")

  logging.info('Elapsed: %s'%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))))

