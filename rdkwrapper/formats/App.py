#!/usr/bin/env python3
#
import sys,os,logging

from rdkit import Chem
from rdkit.Chem import AllChem


if __name__=='__main__':
  fpath_in=sys.argv[1]
  fpath_out=sys.argv[2]
  #sniffer=csv.Sniffer()
  #d=file(fpath_in,'r').read(1000)
  #delim=sniffer.sniff(d).delimiter
  molsuppl = Chem.SmilesMolSupplier(fpath_in,delimiter=' ', smilesColumn=0,nameColumn=1,titleLine=False, sanitize=False)

  sdwriter=Chem.SDWriter(fpath_out)

  i_mol=0
  for mol in molsuppl:
    i_mol+=1
    logging.info('%d. %s' %(i_mol,mol.GetProp('_Name'))) 
    Chem.SanitizeMol(mol)
    AllChem.Compute2DCoords(mol)
    sdwriter.write(mol)

  logging.info('%d mols written to %s' %(i_mol,fpath_out)) 
