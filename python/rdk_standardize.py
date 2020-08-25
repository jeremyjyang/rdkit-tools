#!/usr/bin/env python3
#############################################################################
import os,sys,re,time,argparse,logging

from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles
#import rdkit.Chem.AllChem

#############################################################################
def Standardize(stdzr, molReader, molWriter):
  n_mol=0; 
  for mol in molReader:
    n_mol+=1
    molname = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    logging.debug('%d. %s:'%(n_mol, molname))
    mol2 = StdMol(stdzr, mol)
    molWriter.write(mol2)
  logging.info('%d mols written to %s' %(n_mol, args.ofile))

#############################################################################
def MyNorms():
  norms = list(MolStandardize.normalize.NORMALIZATIONS)
  for i in range(len(norms)-1, 0, -1):
   norm = norms[i]
   if norm.name == "Sulfoxide to -S+(O-)-":
     del(norms[i])
  norms.append(MolStandardize.normalize.Normalization("[S+]-[O-] to S=O",
	"[S+:1]([O-:2])>>[S+0:1](=[O-0:2])"))
  logging.info("Normalizations: {}".format(len(norms)))
  return(norms)

#############################################################################
def ReadNormsFile(fin):
  norms=[];
  while True:
    line = fin.readline()
    if not line: break
    smirks, name = re.split(r'[\s]+', line.rstrip(), 1)
    norms.append(MolStandardize.normalize.Normalization(name, smirks))
  logging.info("Normalizations: {}".format(len(norms)))
  return(norms)

#############################################################################
def ShowParameters():
  logging.info("PREFER_ORGANIC: {}".format(MolStandardize.fragment.PREFER_ORGANIC))
  logging.info("MAX_RESTARTS: {}".format(MolStandardize.normalize.MAX_RESTARTS))
  logging.info("MAX_TAUTOMERS: {}".format(MolStandardize.tautomer.MAX_TAUTOMERS))
  logging.info("ACID_BASE_PAIRS:\n\t{}".format(
	"\n".join(map(lambda ab: ("\t{}\t{}\t{}".format(ab.name, MolToSmiles(ab.acid), MolToSmiles(ab.base))), MolStandardize.charge.ACID_BASE_PAIRS))))
  logging.info("CHARGE_CORRECTIONS:\n\t{}".format(
	"\n".join(map(lambda cc: ("\t{}\t{}\t{}".format(cc.name, MolToSmiles(cc.smarts), cc.charge)), MolStandardize.charge.CHARGE_CORRECTIONS))))
  logging.info("TAUTOMER_TRANSFORMS:\n\t{}".format(
	"\n".join(map(lambda tt: ("\t{}\t{}\t{}\t{}".format(tt.name, tt.tautomer_str, tt.bonds, tt.charges)), MolStandardize.tautomer.TAUTOMER_TRANSFORMS))))
  logging.info("TAUTOMER_SCORES:\n\t{}".format(
	"\n".join(map(lambda ts: ("\t{}\t{}\t{}".format(ts.name, ts.smarts_str, ts.score)), MolStandardize.tautomer.TAUTOMER_SCORES))))

#############################################################################
def MyStandardizer(norms):
  stdzr = MolStandardize.Standardizer(
	normalizations = norms,
	max_restarts = MolStandardize.normalize.MAX_RESTARTS,
	prefer_organic = MolStandardize.fragment.PREFER_ORGANIC,
	acid_base_pairs = MolStandardize.charge.ACID_BASE_PAIRS,
	charge_corrections = MolStandardize.charge.CHARGE_CORRECTIONS,
	tautomer_transforms = MolStandardize.tautomer.TAUTOMER_TRANSFORMS,
	tautomer_scores = MolStandardize.tautomer.TAUTOMER_SCORES,
	max_tautomers = MolStandardize.tautomer.MAX_TAUTOMERS
	)
  return(stdzr)

#############################################################################
def ListNormalizations(norms, fout):
  for norm in norms:
    fout.write("{}\t{}\n".format(norm.transform_str, norm.name))
  logging.info("Normalizations: {}".format(len(norms)))

#############################################################################
def StdMol(stdzr, mol):
  smi = MolToSmiles(mol, isomericSmiles=False) if mol else None
  mol_std = stdzr.standardize(mol) if mol else None
  smi_std = MolToSmiles(mol_std, isomericSmiles=False) if mol_std else None
  logging.debug("{:>28s} >> {}".format(smi, smi_std))
  return(mol_std)

#############################################################################
def Demo(norms):
  stdzr = MyStandardizer(norms)
  smis = [
	'CCC(=O)O',
	'c1ccc(cc1)N(=O)=O',
	'c1ccc(cc1)[N+](=O)[O-]',
	'CCC([O-])[OH+]',
	'CCN(=O)([O-])[OH+]',
	'C[C@H]1CN(CCCN1[S+]([O-])(=O)C2=CC=CC3=C2C(=CN=C3)C)C(=O)CN',
	'N[S+]([O-])(=O)C1=C(Cl)C=C2NC(N[S+]([O-])(=O)C2=C1)C(Cl)Cl',
	'C[S+]([O-])C1=CC=C(C=C1)\C=C2\C(=C(\CC(O)=O)C3=C2C=CC(=C3)F)C',
	'O=[N+]([O-])c1cc(S(=O)(=O)N2CCCC2)ccc1NN=Cc1ccccn1'
	]
  for smi in smis:
    mol1 = MolFromSmiles(smi)
    mol2 = stdzr.standardize(mol1) if mol1 else None
    smi_std = MolToSmiles(mol2, isomericSmiles=False) if mol2 else None
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
  parser.add_argument("-v", "--verbose", action="count", default=0)

  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  t0=time.time()

  if args.norms=="unm":
    norms = MyNorms()
  elif args.ifile_norms:
    fin = open(args.ifile_norms)
    norms = ReadNormsFile(fin)

  else:
    norms = MolStandardize.normalize.NORMALIZATIONS

  if args.op=="list_norms":
    fout = open(args.ofile, "w") if args.ofile else sys.stdout
    ListNormalizations(norms, fout)

  elif args.op=="standardize":
    if not (args.ifile and args.ofile): parser.error('--i and --o required.')
    if re.sub(r'.*\.', '', args.ifile).lower()in ('smi', 'smiles'):
      molReader = SmilesMolSupplier(args.ifile, delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=True, sanitize=True)
    elif re.sub(r'.*\.', '', args.ifile).lower() in ('sdf','sd','mdl','mol'):
      molReader = SDMolSupplier(args.ifile, sanitize=True, removeHs=True)
    else:
      parser.error('Invalid file extension: {}'.format(args.ifile))

    if re.sub(r'.*\.', '', ofile).lower() in ('sdf','sd','mdl','mol'):
      molWriter = SDWriter(ofile)
    elif re.sub(r'.*\.', '', ofile).lower()in ('smi', 'smiles'):
      molWriter = SmilesWriter(ofile, delimiter='\t', nameHeader='Name',
        includeHeader=True, isomericSmiles=True, kekuleSmiles=False)
    else:
      logging.error('Invalid file extension: {}'.format(ofile))

    stdzr = MyStandardizer(norms)
    Standardize(stdzr, molReader, molWriter)

  elif args.op=="show_params":
    ShowParameters()

  elif args.op=="demo":
    Demo(norms)

  else:
    parser.error("Unsupported operation: {}".format(args.op))

  logging.info('Elapsed: %s'%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))))

