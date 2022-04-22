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
  logging.info(f"SMIRKS: {smirks}")
  rxn = rdChemReactions.ReactionFromSmarts(smirks)
  for molmix in molReader:
    reactants = Chem.rdmolops.GetMolFrags(molmix, asMols=True)
    logging.debug(f"REACTANTS: {len(reactants)}")
    for k,r in enumerate(reactants):
      logging.debug(f"REACTANT {k+1}: {Chem.MolToSmiles(r)}")
    product_sets = rxn.RunReactants(reactants)
    for i,product_set in enumerate(product_sets):
      for j,product in enumerate(product_set):
        logging.info(f"PRODUCT {i+1}.{j+1}: {Chem.MolToSmiles(product)}")
        molWriter.write(product)

#############################################################################
def EnumerateLibrary(smirks, molReaders, molWriter):
  """Reactant-set order must agree with SMIRKS reactant order."""
  logging.info(f"SMIRKS: {smirks}")
  logging.info(f"molReaders: {len(molReaders)}")
  rxn = rdChemReactions.ReactionFromSmarts(smirks)
  reactant_sets=[];
  for molReader in molReaders:
    reactant_sets.append([m for m in molReader])
  for i,reactant_set in enumerate(reactant_sets):
    logging.debug(f"REACTANT_SET {i+1}: "+(",".join([Chem.MolToSmiles(m) for m in reactant_set])))
  product_sets = AllChem.EnumerateLibraryFromReaction(rxn, reactant_sets) # tuple of tuples
  for i,product_set in enumerate(product_sets):
    for j,product in enumerate(product_set):
      logging.info(f"PRODUCT {i+1}.{j+1}: {Chem.MolToSmiles(product)}")
      molWriter.write(product)

#############################################################################
def Demo():
  smirks = '[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]'
  logging.info(f"SMIRKS: {smirks}")
  reactant_smis = ['C(=O)O', 'CNC']
  for i,smi in enumerate(reactant_smis):
    logging.info(f"REACTANT {i+1}: {smi}")
  rxn = rdChemReactions.ReactionFromSmarts(smirks)
  reactants = [Chem.MolFromSmiles(smi) for smi in reactant_smis]
  product_sets = rxn.RunReactants(reactants) # tuple of tuples
  for i,product_set in enumerate(product_sets):
    for j,product in enumerate(product_set):
      logging.info(f"PRODUCT {i+1}.{j+1}: {Chem.MolToSmiles(product)}")

#############################################################################
def Demo2():
  #rxnFile = '/home/app/src/rdkit/rdkit/Chem/SimpleEnum/test_data/azide_reaction.rxn'
  #rxnFile = '/home/app/src/rdkit/Code/GraphMol/ChemReactions/testData/AmideBond.rxn'
  rxnFile = f"{os.environ['HOME']}/src/rdkit-tools/data/boronic1.rxn"
  with open(rxnFile) as f:
    txt = f.read()
  logging.debug(f"{rxnFile}: {txt}")
  rxn = AllChem.ReactionFromRxnFile(rxnFile)
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


#############################################################################
def Demo3():
  smirks = '[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]'
  molReader1 = Chem.SmilesMolSupplierFromText("""\
OC=O
OC(=O)C
OC(=O)CCc1ccccc1
""", delimiter="\t", smilesColumn=0, nameColumn=0, titleLine=False)
  molReader2 = Chem.SmilesMolSupplierFromText("""\
NC
NCC
NCc1ccccc1
""", delimiter="\t", smilesColumn=0, nameColumn=0, titleLine=False)
  molReaders = [ molReader1, molReader2 ]
  molWriter = Chem.SmilesWriter("-", delimiter='\t', includeHeader=False, isomericSmiles=True, kekuleSmiles=True)
  EnumerateLibrary(smirks, molReaders, molWriter)

