#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.MurckoScaffold.html
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.rdScaffoldNetwork.html
rdScaffoldNetwork available RDKit 2020.03.1+.
"""
#############################################################################
import os,sys,re,logging,json,time,inspect

import matplotlib as mpl
#from matplotlib import pyplot as plt

import rdkit
import rdkit.Chem
import rdkit.Chem.AllChem
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles, Draw, rdDepictor

# fingerprints:
from rdkit.Chem import RDKFingerprint, PatternFingerprint, LayeredFingerprint, LayeredFingerprint_substructLayers
from rdkit.Chem.MACCSkeys import GenMACCSKeys
from rdkit.Chem.AllChem import GetMorganFingerprint,GetMorganFingerprintAsBitVect

# scaffolds:
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Scaffolds import rdScaffoldNetwork

import pyvis
from pyvis.network import Network

DEMOSMIS = [
'C[C@H]1CN(CCCN1[S+]([O-])(=O)C2=CC=CC3=C2C(=CN=C3)C)C(=O)CN Rho Kinase Inhibitor IV',
'N[S+]([O-])(=O)C1=C(Cl)C=C2NC(N[S+]([O-])(=O)C2=C1)C(Cl)Cl trichlormethiazide',
'C[S+]([O-])C1=CC=C(C=C1)\C=C2\C(=C(\CC(O)=O)C3=C2C=CC(=C3)F)C Sulindac',
'O=[N+]([O-])c1cc(S(=O)(=O)N2CCCC2)ccc1NN=Cc1ccccn1 PUBCHEM_CID:4036736',
'Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12 flucloxacillin',
'CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O ampicillin',
'CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O.[Na] penicillin',
'Cc1onc(-c2ccccc2)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12 oxacillin'
	]

#############################################################################
def Mols2BMScaffolds(mols, molWriter):
  for i,mol in enumerate(mols):
    molname = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    logging.debug(f'{i+1}. {molname}')
    scafmol = MurckoScaffold.GetScaffoldForMol(mol)
    molWriter.write(scafmol)
  logging.info(f'{len(mols)} mols written to {molWriter}')

#############################################################################
def Mols2ScafNet(mols, brics=False, ofile=None):
  if brics:
    params = rdScaffoldNetwork.BRICSScaffoldParams()
  else:
    params = rdScaffoldNetwork.ScaffoldNetworkParams()
    params.flattenChirality = True
    params.flattenIsotopes = True
    params.flattenKeepLargest = True
    params.includeGenericBondScaffolds = False
    params.includeGenericScaffolds = False
    params.includeScaffoldsWithAttachments = False
    params.includeScaffoldsWithoutAttachments = True
    params.keepOnlyFirstFragment = False
    params.pruneBeforeFragmenting = True

  attrs = [a for a in inspect.getmembers(params) if not(a[0].startswith('__'))]
  for a in attrs:
    logging.info(f"{a[0]}: {a[1]}")

  for i,mol in enumerate(mols):
    molname = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    logging.debug(f'{i+1}. {molname}:')

  scafnet = rdScaffoldNetwork.CreateScaffoldNetwork(mols, params)
  fout = open(ofile, "w") if ofile else sys.stdout
  for i in range(len(scafnet.nodes)):
    fout.write(f"node\t{i}\t{scafnet.nodes[i]}\t{scafnet.counts[i]}\n")
  for i in range(len(scafnet.edges)):
    fout.write(f"edge\t{i}\t{scafnet.edges[i].beginIdx}\t{scafnet.edges[i].endIdx}\t{scafnet.edges[i].type}\n")
  fout.flush()
  logging.info(f"nodes: {len(scafnet.nodes)}; edges:{len(scafnet.edges)}")
  return scafnet

#############################################################################
