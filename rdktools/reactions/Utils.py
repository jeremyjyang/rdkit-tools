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

from .. import util

#############################################################################
def React(smirks, molReader, molWriter):
  rxn = rdChemReactions.ReactionFromSmarts(smirks)

#############################################################################
def EnumerateLibrary(smirks, molReaders, molWriter):
  """
s1=[Chem.MolFromSmiles(x) for x in ('NC','NCC')]
s2=[Chem.MolFromSmiles(x) for x in ('OC=O','OC(=O)C')]
rxn = AllChem.ReactionFromSmarts('[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]')
r = AllChem.EnumerateLibraryFromReaction(rxn,[s2,s1])
"""
  logging.debug(f"{smirks}")
  rxn = rdChemReactions.ReactionFromSmarts(smirks)
  molsets=[];
  for molReader in molReaders:
    mols = util.ReadMols(molReader)
    molsets.append(mols)
  logging.debug(f"len(molsets): {len(molsets)}")

  lib = AllChem.EnumerateLibraryFromReaction(rxn, molsets)
  prods = [x[0] for x in list(lib)]
  logging.debug(f"len(prods): {len(prods)}")
  for prod in prods:
    logging.debug(f"{Chem.MolToSmiles(prod)})")
    molWriter.write(prod)

#  for m1 in molReader:
#    frags = Chem.rdmolops.GetMolFrags(m1, asMols=True)
#    logging.debug(f"len(frags): {len(frags)}")
#    for f in frags:
#      logging.debug(f"{Chem.MolToSmiles(f)}")
#    #lib = rdChemReactions.EnumerateLibrary(rxn, [(f,) for f in frags])
#    lib = AllChem.EnumerateLibraryFromReaction(rxn, [(f,) for f in frags])
#    prods = [x[0] for x in list(lib)]
#    logging.debug(f"len(prods): {len(prods)}")
#    for prod in prods:
#      logging.debug(f"{Chem.MolToSmiles(prod)})")
#      molWriter.write(prod)

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

#############################################################################
def Demo3():
  s1=[Chem.MolFromSmiles(x) for x in ('NC','NCC')]
  s2=[Chem.MolFromSmiles(x) for x in ('OC=O','OC(=O)C')]
  rxn = rdChemReactions.ReactionFromSmarts('[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]')
  #lib = rdChemReactions.EnumerateLibrary(rxn,[s2,s1])
  lib = AllChem.EnumerateLibraryFromReaction(rxn, [s2,s1])
  prods = [x[0] for x in list(lib)]
  logging.debug(f"len(prods): {len(prods)}")
  for prod in prods:
    logging.debug(f"{Chem.MolToSmiles(prod)}")
