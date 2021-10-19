#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/rdkit.Chem.rdChemReactions.html
"""
#############################################################################
import sys,os,re,logging

from rdkit import Chem
from rdkit.Chem import MolFromSmiles, AllChem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.rdChemReactions import PreprocessReaction

#############################################################################
def React(smirks, molReader, molWriter):
  rxn = rdChemReactions.ReactionFromSmarts(smirks)
  for m1 in molReader:
    frags = Chem.rdmolops.GetMolFrags(m1, asMols=True)
    r = AllChem.EnumerateLibraryFromReaction(rxn, frags)
    prods = [x[0] for x in r]
    for prod in prods:
      logging.debug(f"{Chem.MolToSmiles(prod)})")
      molWriter.write(prod)

#############################################################################
def Demo():
  rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')
  reactants = (Chem.MolFromSmiles('C(=O)O'),Chem.MolFromSmiles('CNC'))
  products = rxn.RunReactants(reactants) # tuple of tuples
  logging.info(f"{Chem.MolToSmiles(products[0][0])}")

#############################################################################
def Demo2():
  testFile = '/home/app/src/rdkit/rdkit/Chem/SimpleEnum/test_data/boronic1.rxn'
  with open(testFile) as f:
    txt = f.read()
  logging.info(f"{testFile}: {txt}")
  rxn = AllChem.ReactionFromRxnFile(testFile)
  rxn.Initialize()

  r1 = rxn.GetReactantTemplate(0)
  m1 = Chem.MolFromSmiles('CCBr')
  m2 = Chem.MolFromSmiles('c1ccccc1Br')

  # These both match because the reaction file itself just has R1-Br:
  logging.info(f"{m1.HasSubstructMatch(r1)}")
  logging.info(f"{m2.HasSubstructMatch(r1)}")

  d = PreprocessReaction(rxn)
  reactantLabels = d[-1]
  logging.info(f"reactantLabels: {reactantLabels}")

  # After preprocessing, we only match the aromatic Br:
  logging.info(f"{m1.HasSubstructMatch(r1)}")
  logging.info(f"{m2.HasSubstructMatch(r1)}")

  #testFile = '/home/app/src/rdkit/rdkit/Chem/SimpleEnum/test_data/azide_reaction.rxn'
  #testFile = '/home/app/src/rdkit/Code/GraphMol/ChemReactions/testData/AmideBond.rxn'
