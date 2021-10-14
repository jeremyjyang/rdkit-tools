#!/usr/bin/env python3
"""
Conformer generation vi RDKit distance geometry method.

Problem: RDKit hard coded to write to both stderr and stdout.  Logging a known issue.
"""
#############################################################################
import os,sys,io,re,time,argparse,logging,tempfile

import rdkit.rdBase
import rdkit.Chem

from .. import util
from .. import conform

#############################################################################
def GenerateConformers(molreader, molwriter, ff, nconf, optiters, etol):
  n_mol=0; n_conf=0;
  for mol in molreader:
    n_mol+=1
    molname = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    logging.info(f"{n_mol}. {molname}: {rdkit.Chem.MolToSmiles(mol, isomericSmiles=False)}")
    ###
    ### redirect sys.stderr
    #fmsg = open("/tmp/z.err", "w")
    #old_target, sys.stderr = sys.stderr, fmsg
    ###

    mol, confIds = conform.GenerateConformations(mol, nconf, ff, optiters, etol)

    ###
    #logging.debug('''fmsg = "{0}"'''.format(open("/tmp/z.err").read()))
    #os.remove("/tmp/z.err")
    ### restore sys.stderr
    #sys.stderr = old_target
    ###

    for confId in confIds:
      molwriter.write(mol, confId = confId)
    n_conf+=len(confIds)
  logging.info(f"mols: {n_mol}; confs: {n_conf}")

#############################################################################
if __name__=="__main__":
  FFS = ["UFF", "MMFF"];
  MDLEXTS = ["sdf", "sd", "mdl", "mol"]
  NCONF=1; OPTITERS=200; ETOL=1e-6;
  OPS = ["generate", "demo"]
  parser = argparse.ArgumentParser(description="RDKit Conformer Generation", epilog="Based on distance geometry method by Blaney et al.")
  parser.add_argument("op", choices=OPS, help="OPERATION")
  parser.add_argument("--i", dest="ifile", help="input file, SMI or SDF")
  parser.add_argument("--o", dest="ofile", help="output SDF with 3D")
  parser.add_argument("--ff", choices=FFS, default="MMFF", help="force-field")
  parser.add_argument("--optiters", type=int, default=200, help="optimizer iterations per conf")
  parser.add_argument("--nconf", type=int, default=1, help="# confs per mol")
  parser.add_argument("--etol", type=float, default=ETOL, help="energy tolerance")
  parser.add_argument("--delim", default=" \t", help="SMILES/TSV delimiter")
  parser.add_argument("--smilesColumn", type=int, default=0, help="")
  parser.add_argument("--nameColumn", type=int, default=1, help="")
  parser.add_argument("--header", action="store_true", help="SMILES/TSV has header line")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f"RDK_VERSION: {rdkit.rdBase.rdkitVersion}")

  t0=time.time()

  if args.op == "demo":
    fin = tempfile.NamedTemporaryFile(mode="w+b", suffix=".smi", delete=False)
    fin.write(b"NCCc1ccc(O)c(O)c1\tdopamine\n")
    fin.close()
    molreader = rdkit.Chem.SmilesMolSupplier(fin.name, delimiter=args.delim, smilesColumn=args.smilesColumn, nameColumn=args.nameColumn, titleLine=args.header, sanitize=True)
    fout = tempfile.NamedTemporaryFile(delete=False)
    molwriter = rdkit.Chem.SDWriter(fout.name)
    GenerateConformers(molreader, molwriter, args.ff, args.nconf, args.optiters, args.etol)
    print(fout.read().decode('utf8'))
    fout.close()
    os.remove(fin.name)
    os.remove(fout.name)
  else:
    if not (args.ifile and args.ofile): parser.error('--i and --o required.')
    if args.ff.upper() not in FFS: parser.error(f"Invalid force field: {args.ff}; allowed values: {','.join(FFS)}")

    if re.sub(r'.*\.', '', args.ifile).lower()=='smi':
      molreader = rdkit.Chem.SmilesMolSupplier(args.ifile, delimiter=args.delim, smilesColumn=args.smilesColumn, nameColumn=args.nameColumn, titleLine=args.header, sanitize=True)
  
    elif re.sub(r'.*\.', '', args.ifile).lower() in MDLEXTS:
      molreader = rdkit.Chem.SDMolSupplier(args.ifile, sanitize=True, removeHs=True)
    else:
      parser.error(f"Invalid file extension: {args.ifile}")

    if re.sub(r'.*\.', '', args.ofile).lower() in MDLEXTS:
      molwriter = rdkit.Chem.SDWriter(args.ofile)
    else:
      parser.error(f"Invalid file extension: {args.ofile}")

    GenerateConformers(molreader, molwriter, args.ff, args.nconf, args.optiters, args.etol)

  logging.info(f"Total elapsed time: {time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))}")
