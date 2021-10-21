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
  for molmix in molReader:
    reactants = Chem.rdmolops.GetMolFrags(molmix, asMols=True)
    logging.debug(f"REACTANTS: {len(reactants)}")
    for r in reactants:
      logging.debug(f"REACTANT: {Chem.MolToSmiles(r)}")
    products = rxn.RunReactants(reactants)
    productmols = [p[0] for p in products]
    logging.info(f"products: {len(productmols)}")
    for m in productmols:
      logging.debug(f"{Chem.MolToSmiles(m)})")
      molWriter.write(m)

#############################################################################
def Demo3():
  smirks = '[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]'
  molReaders = [
	Chem.SmilesMolSupplierFromText("""\
OC=O
OC(=O)C
OC(=O)CCc1ccccc1
""", delimiter="\t", smilesColumn=0, nameColumn=0, titleLine=False),
	Chem.SmilesMolSupplierFromText("""\
NC
NCC
NCc1ccccc1
""", delimiter="\t", smilesColumn=0, nameColumn=0, titleLine=False)
	]
  molWriter = Chem.SmilesWriter("-", delimiter='\t', includeHeader=False, isomericSmiles=True, kekuleSmiles=True)
  EnumerateLibrary(smirks, molReaders, molWriter)

#############################################################################
def EnumerateLibrary(smirks, molReaders, molWriter):
  """Molset order must agree with SMIRKS reactant order."""
  logging.info(f"molReaders: {len(molReaders)}")
  rxn = rdChemReactions.ReactionFromSmarts(smirks)
  logging.info(f"{smirks}")
  molsets=[];
  for molReader in molReaders:
    molsets.append([m for m in molReader])
  for i,molset in enumerate(molsets):
    logging.debug(f"MOLSET {i+1}: "+(",".join([Chem.MolToSmiles(m) for m in molset])))
  products = AllChem.EnumerateLibraryFromReaction(rxn, molsets)
  #products = rdChemReactions.EnumerateLibrary(rxn, molsets)
  productmols = [p[0] for p in products]
  for m in productmols:
    logging.debug(f"PRODUCT: {Chem.MolToSmiles(m)}")
    molWriter.write(m)
  logging.info(f"products: {len(productmols)}")

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
