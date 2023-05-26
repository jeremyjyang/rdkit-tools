#!/usr/bin/env python3
"""
http://rdkit.org/docs/source/rdkit.Chem.fmcs.fmcs.html
https://www.rdkit.org/docs/source/rdkit.Chem.rdFMCS.html
See example in rdkit source: Code/GraphMol/FMCS/Wrap/testFMCS.py
"""
#############################################################################
import os,sys,re,logging,json,time,inspect,tempfile,stat

import rdkit
import rdkit.Chem
#import rdkit.Chem.AllChem
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles, Draw, rdDepictor
from rdkit.Chem.AllChem import Compute2DCoords

from rdkit.Chem import rdFMCS

from .. import util

#############################################################################
def Mols2FMCS(mols, atom_compare, bond_compare, fout=None):

  params = rdFMCS.MCSParameters()
  setattr(params, "AtomTyper", atom_compare)
  setattr(params, "BondTyper", bond_compare)

  params.MaximizeBonds = True

  params.AtomCompareParameters.MatchFormalCharge = False
  params.AtomCompareParameters.RingMatchesRingOnly = True

  params.BondCompareParameters.RingMatchesRingOnly = True
  params.BondCompareParameters.CompleteRingsOnly = True
  params.BondCompareParameters.MatchStereo = False

  mcs = rdFMCS.FindMCS(mols, params)

  logging.info(f"MCS atoms: {mcs.numAtoms}; bonds: {mcs.numBonds}")
  logging.info(f"MCS SMARTS: {mcs.smartsString}")

  qm = rdkit.Chem.MolFromSmarts(mcs.smartsString)

  return qm

#############################################################################
