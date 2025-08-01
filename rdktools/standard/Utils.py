#!/usr/bin/env python3
#############################################################################
import os,sys,re,time,argparse,logging
import pandas as pd

import rdkit
from rdkit import Chem
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles, SaltRemover

from .. import util

#############################################################################
def StdMol(stdzr, mol, sanitize, isomeric=True, nameSDField=None):
  if not mol: return None
  if isomeric:
    Set_Chiral_Flags(mol)
    CheckStereo(mol)
  smi = MolToSmiles(mol, isomericSmiles=isomeric)
  if mol and sanitize:
    flag = Chem.SanitizeMol(mol, Chem.rdmolops.SanitizeFlags.SANITIZE_ALL)
    logging.debug(f"Sanitize_flag: {flag}")
  mol_std = stdzr.standardize(mol)
  if mol_std:
    if nameSDField is not None: mol_std.SetProp("_Name", mol.GetProp(nameSDField))
    else: mol_std.SetProp("_Name", util.MolName(mol))
  smi_std = MolToSmiles(mol_std, isomericSmiles=isomeric)
  logging.debug(f"{smi:>28s} >> {smi_std}")
  return mol_std

#############################################################################
def Canonicalize(isomeric, kekulize, nameSDField, molReader, molWriter):
  """Canonicalize only."""
  n_mol=0; n_out=0; n_empty_in=0; n_empty_out=0; n_err=0;
  for mol in molReader:
    n_mol+=1
    if mol is None:
      n_err+=1
      logging.error(f"[N={n_mol}] Failed to read mol.")
      mol = Chem.Mol() #empty mol
    elif mol.GetNumAtoms()==0:
      logging.info(f"[N={n_mol}] {util.MolName(mol)}: Empty molecule -- no atoms.")
      n_empty_in+=1
      n_empty_out+=1
    else:
      try:
        if kekulize:
          Chem.Kekulize(mol, clearAromaticFlags=True)
        else:
          Chem.SanitizeMol(mol) # setAromaticity=True (default)
      except Exception as e:
        logging.error(f"[N={n_mol}]: canonicalize failed: {e}")
        n_err+=1
      logging.debug(f"[N={n_mol}] {util.MolName(mol)}: {MolToSmiles(mol, isomericSmiles=isomeric)}")
      if nameSDField is not None: mol.SetProp("_Name", mol.GetProp(nameSDField))
      else: mol.SetProp("_Name", util.MolName(mol))
    molWriter.write(mol)
    n_out+=1
  logging.info(f"Mols in: {n_mol}; empty mols in: {n_empty_in}; mols out: {n_out}; empty mols out: {n_empty_out}; errors: {n_err}")

#############################################################################
def NeutralizeAtoms(mol):
  """From https://www.rdkit.org/docs/Cookbook.html"""
  pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
  at_matches = mol.GetSubstructMatches(pattern)
  at_matches_list = [y[0] for y in at_matches]
  for at_idx in at_matches_list:
    atom = mol.GetAtomWithIdx(at_idx)
    chg = atom.GetFormalCharge()
    hcount = atom.GetTotalNumHs()
    atom.SetFormalCharge(0)
    atom.SetNumExplicitHs(hcount - chg)
    logging.debug(f"NeutralizeAtoms pattern matched: [{atom.GetSymbol()}{chg}H{hcount}] >> [{atom.GetSymbol()}{atom.GetFormalCharge()}H{atom.GetTotalNumHs()}]")
    atom.UpdatePropertyCache()
  return mol

#############################################################################
def SaltRemove(metalRemove, largestFragment, neutralize, nameSDField, molReader, molWriter):
  """Salt-remove plus optionally, metals and non-largest fragments."""
  salt_smarts = "[Cl,Br,Fe,Cu,Ni,Zn]" if metalRemove else "[Cl,Br]"
  saltRemover = SaltRemover.SaltRemover(defnData=salt_smarts)
  LargestFragmentChooser = MolStandardize.rdMolStandardize.LargestFragmentChooser()
  n_mol=0; n_out=0; n_mod=0; n_empty_in=0; n_empty_out=0; n_err=0;
  for mol in molReader:
    n_mol+=1
    if mol is None:
      n_err+=1
      logging.error(f"[N={n_mol}] Failed to read mol.")
      mol_out = Chem.Mol() #empty mol
    elif mol.GetNumAtoms()==0:
      logging.info(f"[N={n_mol}] {util.MolName(mol)}: Empty molecule -- no atoms.")
      n_empty_in+=1
      n_empty_out+=1
    else:
      try:
        mol_out = saltRemover.StripMol(mol)
        if largestFragment: mol_out = LargestFragmentChooser.choose(mol_out)
        if neutralize: mol_out = NeutralizeAtoms(mol_out)
      except Exception as e:
        logging.error(f"[N={n_mol}]: salt-remove failed: {e}")
        n_err+=1
      if mol.GetNumAtoms()==0:
        logging.warning(f"[N={n_mol}]: salt-remove removed all atoms!")
      logging.debug(f"[N={n_mol}] {util.MolName(mol)}: {MolToSmiles(mol)}")
      if nameSDField is not None: mol_out.SetProp("_Name", mol.GetProp(nameSDField))
      else: mol_out.SetProp("_Name", util.MolName(mol_out))
      if mol_out.GetNumAtoms() < mol.GetNumAtoms(): n_mod+=1
    molWriter.write(mol_out)
    n_out+=1
  logging.info(f"Mols in: {n_mol}; empty mols in: {n_empty_in}; mols out: {n_out}; empty mols out: {n_empty_out}; mols modified: {n_mod}; errors: {n_err}")

#############################################################################
def Standardize(stdzr, sanitize, isomeric, nameSDField, molReader, molWriter):
  """Sanitize, standardize, and canonicalize."""
  n_mol=0; n_out=0; n_empty_in=0; n_empty_out=0; n_err=0;
  for mol in molReader:
    mol_in = Chem.Mol(mol)
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
      try:
        logging.debug(f"[N={n_mol}] {util.MolName(mol)}: {MolToSmiles(mol, isomericSmiles=isomeric)}")
        mol_out = StdMol(stdzr, mol, sanitize, isomeric, nameSDField)
      except Exception as e:
        logging.error(f"[N={n_mol}]: standardize (StdMol) failed: {e}")
        n_err+=1
        mol_out = Chem.Mol(mol_in)
    if mol_out.GetNumAtoms()==0: n_empty_out+=1
    try:
      molWriter.write(mol_out)
      n_out+=1
    except Exception as e:
      logging.error(f"[N={n_mol}]: standardize (write) failed: {e}")
      n_err+=1
  logging.info(f"Mols in: {n_mol}; empty mols in: {n_empty_in}; mols out: {n_out}; empty mols out: {n_empty_out}; errors: {n_err}")

#############################################################################
def ValidateMolecules(molReader):
  """Validate and report any problems."""
  n_mol=0; n_empty_in=0; n_invalid=0; n_err=0;
  for mol in molReader:
    n_mol+=1
    if mol is None:
      n_err+=1
      logging.error(f"[N={n_mol}] Failed to read mol.")
    elif mol.GetNumAtoms()==0:
      logging.info(f"[N={n_mol}] {util.MolName(mol)}: Empty molecule -- no atoms.")
      n_empty_in+=1
    else:
      try:
        smi = MolToSmiles(mol)
        logging.debug(f"[N={n_mol}] {util.MolName(mol)}: {smi}")
      except Exception as e:
        logging.error(f"[N={n_mol}]: validation failed: {e}")
        n_err+=1
        n_invalid+=1
  logging.info(f"Mols in: {n_mol}; empty mols in: {n_empty_in}; errors: {n_err}; invalid/non-empty: {n_invalid}")

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
  if fout is not None: df.to_csv(fout, "\t", index=False)
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
def CheckStereo(mol_in):
  if not mol_in: return None
  mol = Chem.Mol(mol_in) #copy mol
  smi_pre = MolToSmiles(mol, isomericSmiles=True)
  error_code = Chem.SanitizeMol(mol, Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY, catchErrors=True)
  if error_code != 0:
    logging.error(f"Sanitize_error: {error_code}")
  smi_san = MolToSmiles(mol, isomericSmiles=True)
  if smi_san != smi_pre:
    logging.info(f"Stereo fixed: {smi_pre} >> {smi_san}")
  else:
    logging.info(f"Stereo ok: {smi_san}")
  return mol

#############################################################################
def Set_Chiral_Flags(mol):
  for a in mol.GetAtoms():
    if a.HasProp("molParity"):
      try:
        parity = int(a.GetProp("molParity"))
      except ValueError:
        parity=None
      if parity: logging.debug(f"{a.GetSymbol()} [idx={a.GetIdx()}] (parity={parity})")
      if parity==1:
        a.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
      elif parity==2:
        a.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
      elif parity==3:
        a.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
  return mol
