#!/usr/bin/env python3
"""
Conformer generation vi RDKit distance geometry method.

Problem: RDKit hard coded to write to both stderr and stdout.  Logging a
known issue.
"""
#############################################################################
import os,sys,re,time,argparse,logging

import rdkit.rdBase
import rdkit.Chem

import rdk_utils

#############################################################################
if __name__=='__main__':
  FFS=['UFF', 'MMFF'];
  NCONF=1; OPTITERS=200; ETOL=1e-6;

  parser = argparse.ArgumentParser(description="RDKit Conformer Generation", epilog="Based on distance geometry method by Blaney et al.")
  parser.add_argument("--i", dest="ifile", help="input file, SMI or SDF")
  parser.add_argument("--o", dest="ofile", help="output SDF with 3D")
  parser.add_argument("--ff", choices=FFS, default="MMFF", help="force-field")
  parser.add_argument("--optiters", type=int, default=200, help="optimizer iterations per conf")
  parser.add_argument("--nconf", type=int, default=1, help="# confs per mol")
  parser.add_argument("--etol", type=float, default=ETOL, help="energy tolerance")
  parser.add_argument("--title_in_header", action="store_true", help="title line in header")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f"RDK_VERSION: {rdkit.rdBase.rdkitVersion}")

  if not (args.ifile and args.ofile): parser.error('--i and --o required.')

  if re.sub(r'.*\.', '', args.ifile).lower()=='smi':
    molreader = rdkit.Chem.SmilesMolSupplier(args.ifile, delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=args.title_in_header, sanitize=True)
  
  elif re.sub(r'.*\.', '', args.ifile).lower() in ('sdf','sd','mdl','mol'):
    molreader = rdkit.Chem.SDMolSupplier(args.ifile, sanitize=True, removeHs=True)
  else:
    parser.error(f"Unrecognized file extension: {args.ifile}")

  if re.sub(r'.*\.', '', args.ofile).lower() in ('sdf', 'sd', 'mdl', 'mol'):
    molwriter = rdkit.Chem.SDWriter(args.ofile)
  else:
    parser.error(f"Unrecognized file extension: {args.ofile}")

  if args.ff.upper() not in FFS:
    parser.error(f'''Unrecognized force field: "{args.ff}"''')

  t0=time.time()
  n_mol=0; n_conf=0;
  for mol in molreader:
    n_mol+=1
    molname = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    logging.debug(f"{n_mol}. {molname}:")
    ###
    ### redirect sys.stderr
    #fmsg = open('/tmp/z.err','w')
    #old_target, sys.stderr = sys.stderr, fmsg
    ###

    mol, confIds = rdk_utils.GenerateConformations(mol, args.nconf, args.ff, args.optiters, args.etol)

    ###
    #logging.debug(f'''fmsg = "{open('/tmp/z.err').read()}"''')
    #os.remove('/tmp/z.err')
    ### restore sys.stderr
    #sys.stderr = old_target
    ###

    for confId in confIds:
      molwriter.write(mol, confId = confId)
    n_conf+=len(confIds)

  logging.info(f"{n_mol} mols, {n_conf} confs written to {args.ofile}")
  logging.info(f"total elapsed time: {time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))}")
