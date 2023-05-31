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

MCS_ATOM_FLAGS = {
	"CompleteRingsOnly": "results cannot include lone ring atoms",
	"MatchChiralTag": "include atom chirality in the match",
	"MatchFormalCharge": "include formal charge in the match",
	"MatchIsotope": "use isotope atom queries in MCSResults",
	"MatchValences": "include atom valences in the match",
	"RingMatchesRingOnly": "ring atoms are only allowed to match other ring atoms",
	}
MCS_BOND_FLAGS = {
	"CompleteRingsOnly": "results cannot include partial rings",
	"MatchFusedRings": "enforce check on ring fusion",
	"MatchFusedRingsStrict": "only enforced if MatchFusedRings is True; the ring fusion must be the same in both query and target",
	"MatchStereo": "include bond stereo in the comparison",
	"RingMatchesRingOnly": "ring bonds are only allowed to match other ring bonds",
	}

#############################################################################
def Mols2FMCS(mols, mces, atom_compare, bond_compare, atom_flags, bond_flags, fout=None):

  for flag in atom_flags:
    if flag not in MCS_ATOM_FLAGS.keys():
      logging.error(f"Atom flag invalid: {flag}")
  for flag in bond_flags:
    if flag not in MCS_BOND_FLAGS.keys():
      logging.error(f"Bond flag invalid: {flag}")

  params = rdFMCS.MCSParameters()
  setattr(params, "AtomTyper", atom_compare)
  setattr(params, "BondTyper", bond_compare)

  params.MaximizeBonds = mces

  params.AtomCompareParameters.MatchFormalCharge = bool("MatchFormalCharge" in atom_flags)
  params.AtomCompareParameters.RingMatchesRingOnly = bool("RingMatchesRingOnly" in atom_flags)

  params.BondCompareParameters.RingMatchesRingOnly = bool("RingMatchesRingOnly" in bond_flags)
  params.BondCompareParameters.CompleteRingsOnly = bool("CompleteRingsOnly" in bond_flags)
  params.BondCompareParameters.MatchStereo = bool("MatchStereo" in bond_flags)

  mcs = rdFMCS.FindMCS(mols, params)

  logging.info(f"MCS atoms: {mcs.numAtoms}; bonds: {mcs.numBonds}")
  logging.info(f"MCS SMARTS: {mcs.smartsString}")

  qm = rdkit.Chem.MolFromSmarts(mcs.smartsString)

  return qm

#############################################################################
