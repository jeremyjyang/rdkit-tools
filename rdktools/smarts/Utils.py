#!/usr/bin/env python3
#############################################################################
import sys,os,argparse,re,logging
import rdkit.Chem
from rdkit.Chem import MolToSmiles, MolFromSmiles
import rdkit.Chem.AllChem
import rdkit.rdBase

#############################################################################
### DeduplicateMatches() - remove matches not USA (Unique Sets of Atoms)
### uumatches - unique set of matches, sorted atom order
### umatches - unique set of matches, original atom order
#############################################################################
def DeduplicateMatches(matches):
  uumatches=[]
  umatches=[]
  for match in matches:
    uumatch=list(match)
    uumatch.sort()
    if uumatch not in uumatches:
      uumatches.append(uumatch)
      umatches.append(match)
  return tuple(umatches)

#############################################################################
class Options(object):
  def __init__(self):
    self.ifile=None;
    self.smarts=None;
    self.ofile=None;
    self.usa=False;
    self.verbose=0;

#############################################################################
def MatchFilter(pat, molReader, molWriter):
  n_mol=0; n_mol_matched=0;
  for mol in molReader:
    matches = mol.GetSubstructMatches(pat, uniquify=True, useChirality=False)
    if len(matches)>0:
      n_mol_matched+=1
      molWriter.write(mol)
    n_mol+=1
  logging.info(f"n_mol: {n_mol}; n_mol_matched: {n_mol_matched}")

#############################################################################
def MatchCounts(pat, molReader, molWriter):
  n_mol=0; n_mol_matched=0;
  for mol in molReader:
    #name = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    #smi = MolToSmiles(mol, isomericSmiles=True)
    matches = mol.GetSubstructMatches(pat, uniquify=True, useChirality=False)
    n_matches = len(matches)
    n_mol_matched += (1 if n_matches>0 else 0)
    matches_usa = DeduplicateMatches(matches)
    n_matches_usa = len(matches_usa)
    mol.SetProp("n_matches", f"{n_matches}")
    mol.SetProp("n_matches_usa", f"{n_matches_usa}")
    if n_mol==0: molWriter.SetProps(mol.GetPropNames())
    molWriter.write(mol)
    n_mol+=1
  logging.info(f"n_mol: {n_mol}; n_mol_matched: {n_mol_matched}")

#############################################################################
