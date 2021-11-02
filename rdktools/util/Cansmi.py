#!/usr/bin/env python3
#############################################################################
import os,sys,re,time,argparse,logging

import rdkit
from rdkit import Chem

from .. import util

#############################################################################
def CanonicalizeSmiles(isomeric, molReader, molWriter):
  n_mol=0; n_out=0; n_err=0; n_empty=0;
  cansmis = set();
  for mol in molReader:
    n_mol+=1
    if mol is None:
      n_err+=1
      logging.error(f"[N={n_mol}] Failed to read mol.")
      continue
    molname = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    if mol.GetNumAtoms()==0:
      logging.info(f"[N={n_mol}] {molname}: Empty molecule -- no atoms.")
      n_empty+=1
    else:
      logging.debug(f"[N={n_mol}] {molname}: {Chem.MolToSmiles(mol, isomericSmiles=isomeric)}")
      cansmi = Chem.MolToSmiles(mol, isomericSmiles=isomeric)
      cansmis.add(cansmi)
    molWriter.write(mol)
    n_out+=1
  logging.info(f"Mols in: {n_mol}; mols out: {n_out}")
  logging.info(f"Empty mols: {n_empty}; non-empty mols: {n_mol-n_empty}")
  logging.info(f"Unique canonical {'ISOMERIC ' if isomeric else ''}SMILES: {len(cansmis)}")
  logging.info(f"Errors: {n_err}")

#############################################################################
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="RDKit SMILES canonicalization", epilog="")
  OPS = ["canonicalize", "demo"]
  parser.add_argument("op", choices=OPS, help="OPERATION")
  parser.add_argument("--i", dest="ifile", required=True, help="input file, SMI or SDF")
  parser.add_argument("--o", dest="ofile", default="-", help="output file, SMI")
  parser.add_argument("--delim", default="\t", help="SMILES/TSV delimiter")
  parser.add_argument("--smilesColumn", type=int, default=0, help="")
  parser.add_argument("--nameColumn", type=int, default=1, help="")
  parser.add_argument("--header", action="store_true", help="SMILES/TSV has header line")
  parser.add_argument("--sanitize", action="store_true", help="Sanitize molecules.")
  parser.add_argument("--kekule", action="store_true", help="Kekulize molecules, Kekule SMILES.")
  parser.add_argument("--isomeric", action="store_true", help="If false, output SMILES isomerism removed")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f"RDKit version: {rdkit.__version__}")

  t0=time.time()

  if args.op=="canonicalize":
    molReader = util.File2Molreader(args.ifile, args.delim, args.smilesColumn, args.nameColumn, args.header, args.sanitize)
    molWriter = util.File2Molwriter(args.ofile, args.delim, args.header, isomericSmiles=args.isomeric, kekuleSmiles=args.kekule)
    CanonicalizeSmiles(args.isomeric, molReader, molWriter)

  elif args.op=="demo":
    Demo()

  else:
    parser.error(f"Unsupported operation: {args.op}")

  logging.info('Elapsed: %s'%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))))

