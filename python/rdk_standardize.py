#!/usr/bin/env python3
#############################################################################
import os,sys,re,time,argparse,logging

import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.Chem.MolStandardize

#############################################################################
def MyStandardizer():
  norms = list(rdkit.Chem.MolStandardize.normalize.NORMALIZATIONS)

  for i in range(len(norms)-1, 0, -1):
   norm = norms[i]
   if norm.name == "Sulfoxide to -S+(O-)-":
     del(norms[i])

  norms.append(rdkit.Chem.MolStandardize.normalize.Normalization("[S+]-[O-] to S=O",
	"[S+:1]([O-:2])>>[S+0:1](=[O-0:2])"))

  logging.info("Normalizations: {}".format(len(norms)))
  for i,n in enumerate(norms):
    logging.info("{}. {} : {}".format(i, n.name, n.transform_str))

  s = rdkit.Chem.MolStandardize.Standardizer(
	normalizations = norms,
	max_restarts = rdkit.Chem.MolStandardize.normalize.MAX_RESTARTS,
	prefer_organic = rdkit.Chem.MolStandardize.fragment.PREFER_ORGANIC,
	acid_base_pairs = rdkit.Chem.MolStandardize.charge.ACID_BASE_PAIRS,
	charge_corrections = rdkit.Chem.MolStandardize.charge.CHARGE_CORRECTIONS,
	tautomer_transforms = rdkit.Chem.MolStandardize.tautomer.TAUTOMER_TRANSFORMS,
	tautomer_scores = rdkit.Chem.MolStandardize.tautomer.TAUTOMER_SCORES,
	max_tautomers = rdkit.Chem.MolStandardize.tautomer.MAX_TAUTOMERS
	)
  return(s)

#############################################################################
def StdMol(s, mol):
  smi = rdkit.Chem.MolToSmiles(mol, isomericSmiles=False) if mol else None
  mol_std = s.standardize(mol) if mol else None
  smi_std = rdkit.Chem.MolToSmiles(mol_std, isomericSmiles=False) if mol_std else None
  logging.debug("{:>28s}\t{}".format(smi, smi_std))
  return(mol_std)

#############################################################################
def Demo():
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
  s = MyStandardizer()
  for smi in smis:
    mol1 = rdkit.Chem.MolFromSmiles(smi)
    mol2 = s.standardize(mol1) if mol1 else None
    smi_std = rdkit.Chem.MolToSmiles(mol2, isomericSmiles=False) if mol2 else None
    logging.info("{:>28s}\t{}\n".format(smi, smi_std))

#############################################################################
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="RDKit chemical standardizer", epilog="")
  parser.add_argument("--i", dest="ifile", help="input file, SMI or SDF")
  parser.add_argument("--o", dest="ofile")
  parser.add_argument("-v", "--verbose", action="count", default=0)

  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  if not (args.ifile and args.ofile): parser.error('--i and --o required.')

  if re.sub(r'.*\.', '', args.ifile).lower()in ('smi', 'smiles'):
    molreader = rdkit.Chem.SmilesMolSupplier(args.ifile, delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=True, sanitize=True)
  elif re.sub(r'.*\.', '', args.ifile).lower() in ('sdf','sd','mdl','mol'):
    molreader = rdkit.Chem.SDMolSupplier(args.ifile, sanitize=True, removeHs=True)
  else:
    parser.error('Invalid file extension: %s'%args.ifile)

  if re.sub(r'.*\.', '', args.ofile).lower() in ('sdf','sd','mdl','mol'):
    molwriter = rdkit.Chem.SDWriter(args.ofile)
  elif re.sub(r'.*\.', '', args.ofile).lower()in ('smi', 'smiles'):
    molwriter = rdkit.Chem.SmilesWriter(args.ofile, delimiter='\t', nameHeader='Name',
        includeHeader=True, isomericSmiles=True, kekuleSmiles=False)
  else:
    parser.error('Invalid file extension: %s'%args.ofile)

  t0=time.time()
  n_mol=0; 
  s = MyStandardizer()
  for mol in molreader:
    n_mol+=1
    molname = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    logging.debug('%d. %s:'%(n_mol, molname))
    mol2 = StdMol(s, mol)
    molwriter.write(mol2)

  logging.info('%d mols written to %s' %(n_mol, args.ofile))
  logging.info('Elapsed: %s'%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))))

