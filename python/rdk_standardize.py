#!/usr/bin/env python3
#############################################################################
import os,sys,re,argparse,logging

import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.Chem.MolStandardize

s = rdkit.Chem.MolStandardize.Standardizer()

smis = [
	'C1=CC=CC=C1',
	'CCC(=O)O',
	'c1ccc(cc1)N(=O)=O',
	'c1ccc(cc1)[N+](=O)[O-]',
	'CCC([O-])[OH+]',
	'CCN(=O)([O-])[OH+]',
	'O=[N+]([O-])c1cc(S(=O)(=O)N2CCCC2)ccc1NN=Cc1ccccn1'
	]

for smi in smis:
  mol1 = rdkit.Chem.MolFromSmiles(smi)
  mol2 = s.standardize(mol1) if mol1 else None
  smi_std = rdkit.Chem.MolToSmiles(mol2, isomericSmiles=False) if mol2 else None
  print("{:>28s} => {}".format(smi, smi_std))
