#!/usr/bin/env python3
#############################################################################
import os,sys,re,argparse,logging

import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.Chem.MolStandardize

def StdSmiles(s, smi, fout):
  mol1 = rdkit.Chem.MolFromSmiles(smi)
  mol2 = s.standardize(mol1) if mol1 else None
  smi_std = rdkit.Chem.MolToSmiles(mol2, isomericSmiles=False) if mol2 else None
  fout.write("{:>28s}\t{}\n".format(smi, smi_std))

#############################################################################
if __name__ == "__main__":
  #logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
  logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

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
    StdSmiles(s, smi, sys.stdout)

  if len(sys.argv)>1:
    with open(sys.argv[1]) as fin:
      while True:
        line = fin.readline()
        if not line: break
        smi = line.rstrip()
        StdSmiles(s, smi, sys.stdout)

