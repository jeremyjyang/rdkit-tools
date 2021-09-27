#!/usr/bin/env python3
#
import sys,os,logging

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.EState.EState import EStateIndices

#############################################################################
def CalcCrippenLogP(mol):
  logp = MolLogP(mol)
  mol.SetProp('WildmanCrippenLogP', f"{logp:.3f}")

#############################################################################
def CalcEStateIndices(mol):
  """
Calculate Kier-Hall electrotopological descriptors
  [x] smiles or molfile input and output with data
  [ ] stdout and stdin
  [ ] sumDeltaI -- can it be done with RDKit?
 
References:
  1. L.H. Hall, B. Mohney and L.B. Kier., "The electrotopological state:
     structure information at the atomic level for molecular graphs",
     JCICS 31 (1) 76-82 (1991)
  2. L.B. Kier and L.H. Hall _Molecular Structure Description:
     The Electrotopological State"_  Academic Press (1999)
  """
  esvals = EStateIndices(mol, force=1)
  mol.SetProp('RDKit_EStates', (','.join(map(lambda x:(f"{x:.3f}"), esvals))))
  mol.SetProp('RDKit_MaxEState', (f"{max(esvals):.3f}"))

#############################################################################
