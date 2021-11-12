#!/usr/bin/env python3
#############################################################################
import os,sys,re,time,argparse,logging
import pandas as pd

import rdkit
from rdkit import Chem
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles

from .. import util

#############################################################################
def StdMol(stdzr, mol, sanitize, isomeric=True):
  smi = MolToSmiles(mol, isomericSmiles=isomeric) if mol else None
  if mol and sanitize:
    flag = Chem.SanitizeMol(mol)
    logging.debug(f"Sanitize_flag: {flag}")
  mol_std = stdzr.standardize(mol) if mol else None
  if mol_std: mol_std.SetProp("_Name", util.MolName(mol)) #Sometimes name is lost?
  smi_std = MolToSmiles(mol_std, isomericSmiles=isomeric) if mol_std else None
  logging.debug(f"{smi:>28s} >> {smi_std}")
  return mol_std

#############################################################################
def Standardize(stdzr, sanitize, isomeric, molReader, molWriter):
  n_mol=0; n_out=0; n_empty_in=0; n_empty_out=0; n_err=0;
  for mol in molReader:
    n_mol+=1
    if mol is None:
      n_err+=1
      logging.error(f"[N={n_mol}] Failed to read mol.")
      mol_out = Chem.Mol() #empty mol
    elif mol.GetNumAtoms()==0:
      logging.info(f"[N={n_mol}] {util.MolName(mol)}: Empty molecule -- no atoms.")
      n_empty_in+=1
      mol_out = mol
    else:
      logging.debug(f"[N={n_mol}] {util.MolName(mol)}: {MolToSmiles(mol, isomericSmiles=isomeric)}")
      try:
        mol_out = StdMol(stdzr, mol, sanitize, isomeric)
      except Exception as e:
        logging.error(f"[N={n_mol}]: standardize failed: {e}")
        n_err+=1
        mol_out = mol
    if mol_out.GetNumAtoms()==0: n_empty_out+=1
    molWriter.write(mol_out)
    n_out+=1
  logging.info(f"Mols in: {n_mol}; empty mols in: {n_empty_in}; mols out: {n_out}; empty mols out: {n_empty_out}; errors: {n_err}")

#############################################################################
def MyNorms():
  norms = list(MolStandardize.normalize.NORMALIZATIONS)
  for i in range(len(norms)-1, 0, -1):
   norm = norms[i]
   if norm.name == "Sulfoxide to -S+(O-)-":
     del(norms[i])
  norms.append(MolStandardize.normalize.Normalization("[S+]-[O-] to S=O", "[S+:1]([O-:2])>>[S+0:1](=[O-0:2])"))
  logging.info(f"Normalizations: {len(norms)}")
  return(norms)

#############################################################################
def ReadNormsFile(fin):
  norms=[];
  while True:
    line = fin.readline()
    if not line: break
    smirks, name = re.split(r'[\s]+', line.rstrip(), 1)
    norms.append(MolStandardize.normalize.Normalization(name, smirks))
  logging.info(f"Normalizations: {len(norms)}")
  return(norms)

#############################################################################
def ShowParameters():
  logging.info(f"PREFER_ORGANIC: {MolStandardize.fragment.PREFER_ORGANIC}")
  logging.info(f"MAX_RESTARTS: {MolStandardize.normalize.MAX_RESTARTS}")
  logging.info(f"MAX_TAUTOMERS: {MolStandardize.tautomer.MAX_TAUTOMERS}")
  logging.info("ACID_BASE_PAIRS:\n"+("\n".join(map(lambda ab: (f"{ab.name}\t{MolToSmiles(ab.acid)}\t{MolToSmiles(ab.base)}"), MolStandardize.charge.ACID_BASE_PAIRS))))
  logging.info("CHARGE_CORRECTIONS:\n"+("\n".join(map(lambda cc: (f"{cc.name}\t{MolToSmiles(cc.smarts)}\t{cc.charge}"), MolStandardize.charge.CHARGE_CORRECTIONS))))
  logging.info("TAUTOMER_TRANSFORMS:\n"+("\n".join(map(lambda tt: (f"{tt.name}\t{tt.tautomer_str}\t{tt.bonds}\t{tt.charges}"), MolStandardize.tautomer.TAUTOMER_TRANSFORMS))))
  logging.info("TAUTOMER_SCORES:\n"+("\n".join(map(lambda ts: (f"{ts.name}\t{ts.smarts_str}\t{ts.score}"), MolStandardize.tautomer.TAUTOMER_SCORES))))

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
  return stdzr

#############################################################################
def ListNormalizations(norms, fout):
  df = pd.DataFrame([[norm.name, norm.transform_str] for norm in norms])
  if fout is not None: df.to_csv(fout, "\t", index=True)
  logging.info(f"Normalizations: {len(norms)}")
  return df

#############################################################################
def Demo():
  norms = MolStandardize.normalize.NORMALIZATIONS
  ListNormalizations(norms, sys.stdout)
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
        'O=[N+]([O-])c1cc(S(=O)(=O)N2CCCC2)ccc1NN=Cc1ccccn1',
        'CCC(CC)COC(=O)[C@H](C)N[P@](=O)(OC[C@H]1O[C@](C#N)([C@H](O)[C@@H]1O)C1=CC=C2N1N=CN=C2N)OC1=CC=CC=C1',
        ]
  for smi in smis:
    mol = MolFromSmiles(smi)
    mol_out = stdzr.standardize(mol) if mol else None
    smi_std = MolToSmiles(mol_out, isomericSmiles=True) if mol_out else None
    logging.info(f"{smi:>28s} >> {smi_std}")

#############################################################################
