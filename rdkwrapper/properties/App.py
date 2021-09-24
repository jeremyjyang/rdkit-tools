#!/usr/bin/env python3
#############################################################################
###  estatecalc_rdk.py - calculate Kier-Hall electrotopological descriptors
###  [x] smiles or molfile input and output with data
###  [ ] stdout and stdin
###  [ ] sumDeltaI -- can it be done with RDKit?
### 
### References:
###    1. L.H. Hall, B. Mohney and L.B. Kier., "The electrotopological state:
###       structure information at the atomic level for molecular graphs",
###       JCICS 31 (1) 76-82 (1991)
###    2. L.B. Kier and L.H. Hall _Molecular Structure Description:
###       The Electrotopological State"_  Academic Press (1999)
###       
#############################################################################
import sys,os,argparse,logging

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EState.EState import EStateIndices

from .. import properties

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='RDKit molecular property utility')
  parser.add_argument("--i", required=True, dest="ifile", help="input molecule file")
  parser.add_argument("--o", dest="ofile", help="output file with data (TSV)")
  parser.add_argument("--inheader", help="header line with smiles input")
  parser.add_argument("--outheader", help="include TSV header line with smiles output")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))


  if args.ifile[-4:]=='.smi':
    molreader=Chem.SmilesMolSupplier(args.ifile,delimiter=' ',
	smilesColumn=0,nameColumn=1,titleLine=inheader,
	sanitize=True)
  elif args.ifile[-4:] in ('.sdf','.sd','.mdl','.mol'):
    molreader=Chem.SDMolSupplier(args.ifile,sanitize=True,removeHs=True)
  else:
    logging.error('Unrecognized file extension: %s'%args.ifile)

  if args.ofile[-4:]=='.smi':
    molwriter=Chem.SmilesWriter(args.ofile,delimiter='\t',nameHeader='Name',
	includeHeader=outheader,isomericSmiles=True,kekuleSmiles=False)
  elif args.ofile[-4:] in ('.sdf','.sd','.mdl','.mol'):
    molwriter=Chem.SDWriter(args.ofile)
  else:
    logging.error('Unrecognized file extension: %s'%args.ofile)

  i_mol=0
  for mol in molreader:
    i_mol+=1
    esvals = EStateIndices(mol,force=1)
    mol.SetProp('RDKit_EStates',(','.join(map(lambda x:('%.2f'%x),esvals))))
    mol.SetProp('RDKit_MaxEState',('%.2f'%max(esvals)))
    if i_mol==1:
      molwriter.SetProps(mol.GetPropNames())
    molwriter.write(mol)

  logging.info('%d mols written to %s' %(i_mol,args.ofile))
