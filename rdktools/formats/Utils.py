#!/usr/bin/env python3
#############################################################################
import sys,os,re,logging

import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.Chem.inchi

#############################################################################
def Mdl2Smi(molReader, molWriter, nameField):
  n_mol=0; n_err=0; n_out=0;
  for mol in molReader:
    n_mol+=1
    try:
      if nameField is not None: mol.SetProp('_Name', mol.GetProp(nameField))
      logging.debug(f"{n_mol}. {mol.GetProp('_Name')}")
      mol_out = rdkit.Chem.RemoveHs(mol, implicitOnly=False)
      rdkit.Chem.SanitizeMol(mol_out)
      rdkit.Chem.AllChem.Kekulize(mol_out, clearAromaticFlags=True)
      molWriter.write(mol_out)
      n_out+=1
    except Exception as e:
      logging.error(f"{n_mol}. {e}")
      n_err+=1
  logging.info(f"n_mol: {n_mol}; n_out: {n_out}; n_err: {n_err}")

#############################################################################
def Mdl2Tsv(molReader, molWriter):
  n_mol=0;
  for mol in molReader:
    logging.debug(f"{n_mol+1}. {mol.GetProp('_Name')}")
    mol_out = rdkit.Chem.RemoveHs(mol, implicitOnly=False)
    rdkit.Chem.SanitizeMol(mol_out)
    rdkit.Chem.AllChem.Kekulize(mol_out, clearAromaticFlags=True)
    if n_mol==0: molWriter.SetProps(mol_out.GetPropNames())
    molWriter.write(mol_out)
    n_mol+=1
  logging.info(f"n_out: {n_mol}")

#############################################################################
def Smi2Mdl(molReader, molWriter):
  n_mol=0;
  for mol in molReader:
    logging.debug(f"{n_mol+1}. {mol.GetProp('_Name')}")
    rdkit.Chem.SanitizeMol(mol)
    rdkit.Chem.AllChem.Compute2DCoords(mol)
    molWriter.write(mol)
    n_mol+=1
  logging.info(f"n_out: {n_mol}")

#############################################################################
def Mol2Inchi(molReader, fout=None):
  n_mol=0; n_out=0; n_err=0;
  for mol in molReader:
    n_mol+=1
    try:
      molId = mol.GetProp('_Name')
      logging.debug(f"{n_mol}. {molId}")
      inchi = rdkit.Chem.inchi.MolToInchi(mol, options='', logLevel=None, treatWarningAsError=False)
      fout.write(f"{molId}\t{inchi}\n")
      n_out+=1
    except Exception as e:
      fout.write("\n")
      logging.error(f"{e}")
      n_err+=1
  logging.info(f"n_mol: {n_mol}; n_out: {n_out}; n_err: {n_err} (blank lines output on errors)")

#############################################################################
def Mol2Inchikey(molReader, fout=None):
  n_mol=0; n_out=0; n_err=0;
  for mol in molReader:
    n_mol+=1
    try:
      molId = mol.GetProp('_Name')
      logging.debug(f"{n_mol}. {molId}")
      inchikey = rdkit.Chem.inchi.MolToInchiKey(mol, options='')
      fout.write(f"{molId}\t{inchikey}\n")
      n_out+=1
    except Exception as e:
      fout.write("\n")
      logging.error(f"{e}")
      n_err+=1
  logging.info(f"n_mol: {n_mol}; n_out: {n_out}; n_err: {n_err} (blank lines output on errors)")

#############################################################################
if __name__=='__main__':
  logging.basicConfig(level=logging.DEBUG)
  logging.info("Kekulization test...")
  #sma='[R]-O'
  #sma='*-[#6]'
  #smi='NCCc1ccc(O)c(O)c1'
  #smi="NC(C=O)CC1=CC=C(O)C=C1"
  #smi="C(C(C=O)N)C1=CN=CN1"
  smi="c1c([nH]cn1)CC(C=O)N"
  sma='[R]'

  pat = rdkit.Chem.MolFromSmarts(sma)
  if not pat:
    logging.error(f'Bad smarts: {sma}')
  mol = rdkit.Chem.MolFromSmiles(smi)
  matches = mol.GetSubstructMatches(pat, uniquify=True, useChirality=False)
  logging.info(matches)
  alist=[]
  for match in matches:
    for a in match:
      if a not in alist:
        alist.append(a)
  alist.sort()
  logging.info(alist)
  rdkit.Chem.AllChem.Kekulize(mol, clearAromaticFlags=True)
  keksmi = rdkit.Chem.MolToSmiles(mol, isomericSmiles=False, kekuleSmiles=True)
  keksmi += (' |ha:{0}|'.format(','.join(map(lambda x:str(x), alist))))
  logging.info(keksmi)
