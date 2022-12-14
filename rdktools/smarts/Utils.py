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
def MatchFilter(smarts, molReader, molWriter):
  pat = rdkit.Chem.MolFromSmarts(smarts)
  n_mol=0; n_mol_matched=0;
  for mol in molReader:
    matches = mol.GetSubstructMatches(pat, uniquify=True, useChirality=False)
    if len(matches)>0:
      n_mol_matched+=1
      molWriter.write(mol)
    n_mol+=1
  logging.info(f"n_mol: {n_mol}; n_mol_matched: {n_mol_matched}")

#############################################################################
def MatchFilterMulti(smartsfile, molReader, molWriter):
  smartses=[];
  with open(smartsfile, "r") as fin:
    while True:
      line = fin.readline()
      if not line: break
      smartses.append(line.rstrip())
  logging.info(f"SMARTS read from file {smartsfile}: {len(smartses)}")
  querys = [{'smarts':smarts, 'pat':rdkit.Chem.MolFromSmarts(smarts), 'n_mol_matched':0
	} for smarts in smartses]

  n_mol=0; n_mol_matched=0;
  for mol in molReader:
    for j,query in enumerate(querys):
      matches = mol.GetSubstructMatches(query['pat'], uniquify=True, useChirality=False)
      if len(matches)>0: query['n_mol_matched']+=1
    matchcounts = [q['n_mol_matched'] for q in querys]
    if min(matchcounts)>0:
      molWriter.write(mol)
      n_mol_matched+=1
    n_mol+=1
  logging.info(f"n_mol: {n_mol}; n_mol_matched: {n_mol_matched}")

#############################################################################
def MatchCounts(smarts, usa, molReader, molWriter):
  pat = rdkit.Chem.MolFromSmarts(smarts)
  n_mol=0; n_mol_matched=0;
  for mol in molReader:
    #name = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    #smi = MolToSmiles(mol, isomericSmiles=True)
    matches = mol.GetSubstructMatches(pat, uniquify=True, useChirality=False)
    if usa: matches = DeduplicateMatches(matches)
    n_matches = len(matches)
    n_mol_matched += (1 if n_matches>0 else 0)
    mol.SetProp("n_matches", f"{n_matches}")
    if n_mol==0: molWriter.SetProps(mol.GetPropNames())
    molWriter.write(mol)
    n_mol+=1
  logging.info(f"n_mol: {n_mol}; n_mol_matched: {n_mol_matched}")

#############################################################################
def MatchCountsMulti(smartsfile, usa, molReader, molWriter):
  smartses=[];
  with open(smartsfile, "r") as fin:
    while True:
      line = fin.readline()
      if not line: break
      smartses.append(line.rstrip())
  logging.info(f"SMARTS read from file {smartsfile}: {len(smartses)}")
  querys = [{'smarts':smarts, 'pat':rdkit.Chem.MolFromSmarts(smarts), 'n_mol_matched':0
	} for smarts in smartses]
  n_mol=0;
  for mol in molReader:
    for j,query in enumerate(querys):
      matches = mol.GetSubstructMatches(query['pat'], uniquify=True, useChirality=False)
      if usa: matches = DeduplicateMatches(matches)
      n_matches = len(matches)
      if n_matches>0: query['n_mol_matched']+=1
      mol.SetProp(f"""n_matches(query_{j+1:02d} = "{query['smarts']}")""", f"{n_matches}")
    if n_mol==0: molWriter.SetProps(mol.GetPropNames())
    molWriter.write(mol)
    n_mol+=1
  logging.info(f"n_mol: {n_mol}; mols matched: "+(','.join([f"query_{j+1:02d}:{q['n_mol_matched']}" for j,q in enumerate(querys)])))

#############################################################################
