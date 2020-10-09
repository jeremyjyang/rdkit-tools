#!/usr/bin/env python3
"""
"""
import os,sys,re,logging

import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.Chem.inchi


#############################################################################
if __name__=='__main__':
  if not rdkit.Chem.inchi.INCHI_AVAILABLE:
    logging.error("INCHI_AVAILABLE={}".format(rdkit.Chem.inchi.INCHI_AVAILABLE))
    exit(1)

  smi = "NCCc1cc(O)c(O)cc1"
  mol = rdkit.Chem.MolFromSmiles(smi)

  #inchi,auxinfo = rdkit.Chem.inchi.MolToInchiAndAuxInfo(mol, options='', logLevel=None, treatWarningAsError=False)
  inchi = rdkit.Chem.inchi.MolToInchi(mol, options='', logLevel=None, treatWarningAsError=False)
  #inchikey = rdkit.Chem.inchi.InchiToInchiKey(inchi)
  inchikey = rdkit.Chem.inchi.MolToInchiKey(mol, options='')

  logging.debug("SMILES: \"{}\"; INCHI: \"{}\"; INCHIKEY: \"{}\"".format(smi, inchi, inchikey))
