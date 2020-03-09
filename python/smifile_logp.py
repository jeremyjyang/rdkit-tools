#!/usr/bin/env python3
#
import sys,os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Crippen import MolLogP

#import csv

if __name__=='__main__':
  fpath_in=sys.argv[1]
  fpath_out=sys.argv[2]
  #sniffer=csv.Sniffer()
  #d=file(fpath_in,'r').read(1000)
  #delim=sniffer.sniff(d).delimiter
  molsuppl=Chem.SmilesMolSupplier(fpath_in,delimiter=' ',
	smilesColumn=0,nameColumn=1,titleLine=False,
	sanitize=False)

  sdwriter=Chem.SDWriter(fpath_out)

  i_mol=0
  for mol in molsuppl:
    i_mol+=1
    if mol.HasProp('_Name'):
      name=mol.GetProp('_Name')
    else:
      name=''
    print('%d. %s' %(i_mol,name), file=sys.stderr)
    Chem.SanitizeMol(mol)
    AllChem.Compute2DCoords(mol)
    logp=MolLogP(mol)
    mol.SetProp('WildmanCrippenLogP','%.2f'%logp)
    sdwriter.write(mol)

  print('%d mols written to %s' %(i_mol,fpath_out), file=sys.stderr)
