#!/usr/bin/env python3
"""
Conformer generation vi RDKit distance geometry method.

Problem: RDKit hard coded to write to both stderr and stdout.  Logging a
known issue.

Jeremy Yang
"""
#############################################################################
import os,sys,re,time,getopt
import rdkit.rdBase
import rdkit.Chem
#import rdkit.Chem.AllChem

import time_utils
import rdk_utils

PROG=os.path.basename(sys.argv[0])

#############################################################################
def ErrorExit(msg):
  print(msg, file=sys.stderr)
  sys.exit(1)

FFS=['UFF','MMFF']

#############################################################################
if __name__=='__main__':
  nconf=1;
  ff='MMFF';
  optiters=200;
  etol=1e-6;
  USAGE='''\
%(PROG)s - RDKit Conformer Generation
Based on distance geometry method by Blaney et al.

required:
	--i INFILE .......... SMI or SDF
	--o OUTFILE ......... SDF with 3D 

options:
	--ff FORCEFIELD ..... %(FFS)s [%(FF)s]
	--nconf NCONF ....... # confs per mol [%(NCONF)d]
	--optiters N  ....... optimizer iterations per conf [%(OPTITERS)d]
	--etol ETOL ......... energy tolerance [%(ETOL)g]
	--v[v] .............. verbose [very]
	--h ................. help

RDKit: %(RDK_VERSION)s
'''%{	'PROG':PROG,
	'NCONF':nconf,
	'OPTITERS':optiters,
	'ETOL':etol,
	'FF':ff,
	'FFS':('|'.join(FFS)),
	'RDK_VERSION':rdkit.rdBase.rdkitVersion
	}

  ifile=None; ofile=None; verbose=0;
  inheader=False;
  opts,pargs=getopt.getopt(sys.argv[1:],'',['h','v','vv','vvv',
	'i=','o=','nconf=','ff=','optiters=','etol=',
        ])
  if not opts: ErrorExit(USAGE)
  for (opt,val) in opts:
    if opt=='--h': ErrorExit(USAGE)
    elif opt=='--i': ifile=val
    elif opt=='--o': ofile=val
    elif opt=='--nconf': nconf=int(val)
    elif opt=='--optiters': optiters=int(val)
    elif opt=='--etol': etol=float(val)
    elif opt=='--ff': ff=val
    elif opt=='--v': verbose=1
    elif opt=='--vv': verbose=2
    elif opt=='--vvv': verbose=3
    else: ErrorExit('Illegal option: %s'%val)

  if not ifile and ofile: ErrorExit('--i and --o required.')

  if ifile[-4:].lower()=='.smi':
    molreader=rdkit.Chem.SmilesMolSupplier(ifile,delimiter=' ',
        smilesColumn=0,nameColumn=1,titleLine=inheader,
        sanitize=True)
  elif ifile[-4:].lower() in ('.sdf','.sd','.mdl','.mol'):
    molreader=rdkit.Chem.SDMolSupplier(ifile,sanitize=True,removeHs=True)
  else:
    ErrorExit('ERROR: unrecognized file extension: %s'%ifile)

  if ofile[-4:].lower() in ('.sdf','.sd','.mdl','.mol'):
    molwriter=rdkit.Chem.SDWriter(ofile)
  else:
    ErrorExit('ERROR: unrecognized file extension: %s'%ofile)

  if ff.upper() not in FFS:
    ErrorExit('ERROR: unrecognized force field: "%s"'%ff)

  t0=time.time()
  n_mol=0; n_conf=0;
  for mol in molreader:
    n_mol+=1
    molname = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    if verbose>1:
      print('%d. %s:'%(n_mol,molname), file=sys.stderr)
    ###
    ### redirect sys.stderr
    #fmsg = open('/tmp/z.err','w')
    #old_target, sys.stderr = sys.stderr, fmsg
    ###

    mol, confIds = rdk_utils.GenerateConformations(mol,nconf,ff,optiters,etol,verbose)

    ###
    #print >>sys.stderr, 'DEBUG: fmsg = "%s"'%(open('/tmp/z.err').read())
    #os.remove('/tmp/z.err')
    ### restore sys.stderr
    #sys.stderr = old_target
    ###

    for confId in confIds:
      molwriter.write(mol, confId = confId)
    n_conf+=len(confIds)

  print('%s: %d mols, %d confs written to %s' %(PROG,n_mol,n_conf,ofile), file=sys.stderr)
  print('%s: total elapsed time: %s'%(PROG,time_utils.NiceTime(time.time()-t0)), file=sys.stderr)
