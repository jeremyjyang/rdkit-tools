#!/usr/bin/env python
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
### Jeremy Yang
###  22 Oct 2009
#############################################################################
import sys,os,getopt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EState.EState import EStateIndices

PROG=os.path.basename(sys.argv[0])

#############################################################################
def ErrorExit(msg):
  print >>sys.stderr,msg
  sys.exit(1)

#############################################################################
if __name__=='__main__':
  usage='''
%(PROG)s [options]

  required:
    --i=<INFILE>  ... .smi or .sdf
    --o=<OUTFILE> ... .smi or .sdf with data

  options:
    --inheader ... header line with smiles input
    --outheader ... include csv header line with smiles output
    --v     ... verbose
    --h     ... help

'''%{'PROG':PROG}

  ifile=None; ofile=None; verbose=0; inheader=False; outheader=False;
  opts,pargs=getopt.getopt(sys.argv[1:],'',['h','v','vv','i=','o=',
	'inheader','outheader'])
  if not opts: ErrorExit(usage)
  for (opt,val) in opts:
    if opt=='--h': ErrorExit(usage)
    elif opt=='--i': ifile=val
    elif opt=='--o': ofile=val
    elif opt=='--inheader': inheader=True
    elif opt=='--outheader': outheader=True
    elif opt=='--v': verbose=1
    elif opt=='--vv': verbose=2
    else: ErrorExit('Illegal option: %s'%val)

  if not ifile and ofile: ErrorExit('-ifile and -ofile required.')

  if ifile[-4:]=='.smi':
    molreader=Chem.SmilesMolSupplier(ifile,delimiter=' ',
	smilesColumn=0,nameColumn=1,titleLine=inheader,
	sanitize=True)
  elif ifile[-4:] in ('.sdf','.sd','.mdl','.mol'):
    molreader=Chem.SDMolSupplier(ifile,sanitize=True,removeHs=True)
  else:
    ErrorExit('ERROR: unrecognized file extension: %s'%ifile)

  if ofile[-4:]=='.smi':
    molwriter=Chem.SmilesWriter(ofile,delimiter='\t',nameHeader='Name',
	includeHeader=outheader,isomericSmiles=True,kekuleSmiles=False)
  elif ofile[-4:] in ('.sdf','.sd','.mdl','.mol'):
    molwriter=Chem.SDWriter(ofile)
  else:
    ErrorExit('ERROR: unrecognized file extension: %s'%ofile)

  i_mol=0
  for mol in molreader:
    i_mol+=1
    
    esvals=EStateIndices(mol,force=1)

    mol.SetProp('RDKit_EStates',(','.join(map(lambda x:('%.2f'%x),esvals))))
    mol.SetProp('RDKit_MaxEState',('%.2f'%max(esvals)))

    if i_mol==1:
      molwriter.SetProps(mol.GetPropNames())

    molwriter.write(mol)

  print >>sys.stderr, '%d mols written to %s' %(i_mol,ofile)
