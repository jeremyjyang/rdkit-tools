#!/usr/bin/env python3
#############################################################################
import os,sys,re,time,argparse,logging

import rdkit
#import rdkit.Chem.AllChem
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles

#############################################################################
def Standardize(stdzr, remove_isomerism, molReader, molWriter):
  n_mol=0; 
  for mol in molReader:
    n_mol+=1
    molname = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    logging.debug('%d. %s:'%(n_mol, molname))
    mol2 = StdMol(stdzr, mol, remove_isomerism)
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
    fout.write(f"{norm.transform_str}\t{norm.name}\n")
  logging.info(f"Normalizations: {len(norms)}")

#############################################################################
def StdMol(stdzr, mol, remove_isomerism=False):
  smi = MolToSmiles(mol, isomericSmiles=(not remove_isomerism)) if mol else None
  mol_std = stdzr.standardize(mol) if mol else None
  smi_std = MolToSmiles(mol_std, isomericSmiles=(not remove_isomerism)) if mol_std else None
  logging.debug(f"{smi:>28s} >> {smi_std}")
  return(mol_std)

#############################################################################