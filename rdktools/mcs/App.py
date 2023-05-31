#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/rdkit.Chem.rdFMCS.html
http://rdkit.org/docs/source/rdkit.Chem.fmcs.fmcs.html
https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-5-S1-O6
"""
#############################################################################
import os,sys,stat,re,json,time,inspect,argparse,logging

import rdkit
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles, Draw
from rdkit.Chem import rdFMCS

from .. import mcs
from .. import util

SCRATCHDIR = f"{os.environ['HOME']}/tmp/rdktools"

#############################################################################
if __name__ == "__main__":
  epilog = """atom_flags: {}; bond_flags: {}""".format(
	", ".join([f"{key} ({val})" for key,val in mcs.MCS_ATOM_FLAGS.items()]),
	", ".join([f"{key} ({val})" for key,val in mcs.MCS_BOND_FLAGS.items()]))
  parser = argparse.ArgumentParser(description="RDKit FMCS analysis", epilog=epilog)
  OPS = ["fmcs", "demo", ]
  parser.add_argument("op", choices=OPS, help="OPERATION")
  parser.add_argument("--i", dest="ifile", help="input file, TSV or SDF")
  parser.add_argument("--o", dest="ofile", help="output file, TSV|SDF")
  parser.add_argument("--atom_compare", choices=rdFMCS.AtomCompare.names.keys(), default="CompareElements", help="atom-compare criterion")
  parser.add_argument("--bond_compare", choices=rdFMCS.BondCompare.names.keys(), default="CompareOrder", help="bond-compare criterion")
  parser.add_argument("--mces", action="store_true", help="maximize common edges/bonds")

  parser.add_argument("--atom_flags", default="", help="comma-separated MCS atom flags")
  parser.add_argument("--bond_flags", default="", help="comma-separated MCS bond flags")

  parser.add_argument("--scratchdir", default=SCRATCHDIR)
  parser.add_argument("--smilesColumn", type=int, default=0, help="SMILES column from TSV (counting from 0)")
  parser.add_argument("--nameColumn", type=int, default=1, help="name column from TSV (counting from 0)")
  parser.add_argument("--idelim", default="\t", help="delim for input TSV")
  parser.add_argument("--odelim", default="\t", help="delim for output TSV")
  parser.add_argument("--iheader", action="store_true", help="input TSV has header")
  parser.add_argument("--oheader", action="store_true", help="output TSV has header")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f"RDKit version: {rdkit.rdBase.rdkitVersion}")

  t0=time.time()

  if not os.path.isdir(args.scratchdir): os.mkdir(args.scratchdir)

  fout = open(args.ofile, "w+") if args.ofile is not None else sys.stdout

  if args.op=="demo":
    mcs.Demo()
    sys.exit()

  if not (args.ifile): parser.error('--i required.')

  if args.op=="fmcs":
    molReader = util.File2Molreader(args.ifile, args.idelim, args.smilesColumn, args.nameColumn, args.iheader)
    mols = util.ReadMols(molReader)
    fmcsmols = mcs.Mols2FMCS(mols,
	args.mces,
	rdFMCS.AtomCompare.names[args.atom_compare],
	rdFMCS.BondCompare.names[args.bond_compare],
	re.split(f'\s*,\s*', args.atom_flags.strip()),
	re.split(f'\s*,\s*', args.bond_flags.strip()),
	fout)
  else:
    parser.error(f"Unsupported operation: {args.op}")

  logging.info(f"""Elapsed: {time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))}""")

