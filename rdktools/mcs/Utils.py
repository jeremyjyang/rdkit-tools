#!/usr/bin/env python3
"""
http://rdkit.org/docs/source/rdkit.Chem.fmcs.fmcs.html
https://www.rdkit.org/docs/source/rdkit.Chem.rdFMCS.html
"""
#############################################################################
import os,sys,re,logging,json,time,inspect,tempfile,stat

import rdkit
import rdkit.Chem
import rdkit.Chem.AllChem
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles, Draw, rdDepictor
from rdkit.Chem.AllChem import Compute2DCoords

from rdkit.Chem.MCS import rdFMCS

from .. import util

#############################################################################
def Mols2FMCS(mols, fout=None):

  params = rdFMCS.MCSParameters()
  setattr(params, "AtomTyper", rdFMCS.AtomCompare.CompareElements)
  setattr(params, "BondTyper", rdFMCS.BondCompare.CompareOrder)

  #params_atom = rdFMCS.MCSAtomCompareParameters()
  #params_bond = rdFMCS.MCSBondCompareParameters()

  mcs_result = rdFMCS.FindMCS(mols, params)

  logging.info(f"MCS atoms: {mcs.numAtoms}; bonds: {mcs.numBonds}")
  logging.info(f"MCS SMARTS: {mcs.smartsString}")

  qm = Chem.MolFromSmarts(mcs.smartsString)

  return qm

#############################################################################
