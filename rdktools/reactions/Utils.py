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
def React(smirks, molReader, output_mode, fout):
  """For this function, each molecule record from the input stream must be a disconnected mixture of reactants."""
  osmis=[];
  logging.debug(f"SMIRKS: {smirks}")
  rxn = rdChemReactions.ReactionFromSmarts(smirks)
  for i,molmix in enumerate(molReader):
    reactants = Chem.rdmolops.GetMolFrags(molmix, asMols=True)
    rsmis = [Chem.MolToSmiles(r) for r in reactants]
    logging.debug("REACTANTS: "+(".".join([f"{rsmi}" for rsmi in rsmis])))
    product_sets = rxn.RunReactants(reactants)
    for j,product_set in enumerate(product_sets):
      psmis = [Chem.MolToSmiles(p) for p in product_set]
      logging.debug("PRODUCTS: "+(".".join([f"{psmi}" for psmi in psmis])))
      rxnsmi = f"{'.'.join(rsmis)}>>{'.'.join(psmis)}"
      logging.debug(f"REACTION {i+1}.{j+1}: {rxnsmi}")
      fout.write(f"{psmi if output_mode=='products' else rxnsmi}\n")
      osmis.append(psmi if output_mode=='products' else rxnsmi)
  return osmis

#############################################################################
def EnumerateLibrary(smirks, molReaders, output_mode, fout):
  """Each input stream corresponds with a reactant-set, with the same reaction role. Reactant-set order must agree with SMIRKS reactant order."""
  logging.info(f"SMIRKS: {smirks}")
  logging.info(f"molReaders: {len(molReaders)}")
  rxn = rdChemReactions.ReactionFromSmarts(smirks)
  reactant_sets=[];
  for molReader in molReaders:
    reactant_sets.append([m for m in molReader])
  for k,reactant_set in enumerate(reactant_sets):
    rsmis = [Chem.MolToSmiles(r) for r in reactant_set]
    logging.debug(f"REACTANT_SET {k+1}: {','.join(rsmis)}")
  product_sets = AllChem.EnumerateLibraryFromReaction(rxn, reactant_sets) # tuple of tuples
  for i,product_set in enumerate(product_sets):
    for j,product in enumerate(product_set):
      psmi = Chem.MolToSmiles(product)
      logging.info(f"PRODUCT {i+1}.{j+1}: {psmi}")
      fout.write(psmi)

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
  EnumerateLibrary(smirks, molReaders, sys.stdout)

#############################################################################
def Demo4():
  smirks = '[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]'
  molReader = Chem.SmilesMolSupplierFromText("""\
OC=O.NC
OC(=O)C.NCC
OC(=O)CCc1ccccc1.NCc1ccccc1
""", delimiter="\t", smilesColumn=0, nameColumn=0, titleLine=False)
  React(smirks, molReader, sys.stdout)
