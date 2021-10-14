#!/usr/bin/env python3
#############################################################################
import os,sys,re,time,argparse,logging

import rdkit
import rdkit.Chem
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter
import rdkit.Chem.AllChem
from rdkit.Chem.rdmolfiles import SmilesMolSupplierFromText

from .. import smarts

#############################################################################
def Demo():
  molReader = SmilesMolSupplierFromText("""\
SMILES	Name
CN1C(=O)N(C)C(=O)C(N(C)C=N2)=C12	caffeine
c1ncccc1C1CCCN1C	nicotine
NC(C)Cc1ccccc1	adderall
c1ccccc1C(=O)OC2CC(N3C)CCC3C2C(=O)OC	cocaine
COc1cc2c(ccnc2cc1)C(O)C4CC(CC3)C(C=C)CN34	quinine
C123C5C(O)C=CC2C(N(C)CC1)Cc(ccc4O)c3c4O5	morphine
C1C(C)=C(C=CC(C)=CC=CC(C)=CCO)C(C)(C)C1	Vitamin A
CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O	benzylpenicillin
Oc1cc(ccc1)/C=C/c2cc(O)cc(O)c2	resveratrol
CC/C(=C(/c1ccc(OCCN(C)C)cc1)c1ccccc1)c1ccccc1	Tamoxifen
CN/C(=C\[N+](=O)[O-])/NCCSCC1=CC=C(O1)CN(C)C.Cl	Zantac
CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)OCC	Tamiflu
CNCCC(c1ccccc1)Oc2ccc(cc2)C(F)(F)F.Cl	Prozac
CN1c2ccc(cc2C(=NCC1=O)c3ccccc3)Cl	Valium
""", delimiter="\t", smilesColumn=0, nameColumn=1, titleLine=True)
  sma = "[$([N;!H0]-[#6]);!$(N-C=[O,N,S])]" #amine
  pat = rdkit.Chem.MolFromSmarts(sma)
  molWriter = SmilesWriter("-", delimiter='\t', nameHeader='Name', includeHeader=True, isomericSmiles=True, kekuleSmiles=False)
  smarts.MatchCounts(pat, molReader, molWriter)

#############################################################################
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="RDKit SMARTS utility", epilog="")
  OPS = ["matchCounts", "matchFilter", "demo"]
  parser.add_argument("op", choices=OPS, help="OPERATION")
  parser.add_argument("--i", dest="ifile", help="input file, SMI or SDF")
  parser.add_argument("--o", dest="ofile", help="output file, TSV")
  parser.add_argument("--smarts", help="query SMARTS")
  parser.add_argument("--usa", action="store_true", help="unique set-of-atoms match counts")
  parser.add_argument("--delim", default=" \t", help="delimiter for SMILES/TSV")
  parser.add_argument("--smilesColumn", type=int, default=0, help="")
  parser.add_argument("--nameColumn", type=int, default=1, help="")
  parser.add_argument("--header", action="store_true", help="SMILES/TSV has header line")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f"RDKit version: {rdkit.__version__}")

  t0=time.time()
  
  if args.op=="demo":
    Demo()
    sys.exit()

  if re.sub(r'.*\.', '', args.ifile).lower()in ('smi', 'smiles'):
    molReader = SmilesMolSupplier(args.ifile, delimiter=args.delim, smilesColumn=args.smilesColumn, nameColumn=args.nameColumn, titleLine=args.header, sanitize=True)
  elif re.sub(r'.*\.', '', args.ifile).lower() in ('sdf','sd','mdl','mol'):
    molReader = SDMolSupplier(args.ifile, sanitize=True, removeHs=True)
  else:
    logging.error(f'Invalid file extension: {args.ifile}')

  if args.ofile is None:
    molWriter = SmilesWriter("-", delimiter=args.delim, nameHeader='Name', includeHeader=True, isomericSmiles=True, kekuleSmiles=False)
  elif re.sub(r'.*\.', '', args.ofile).lower() in ('sdf','sd','mdl','mol'):
    molWriter = SDWriter(args.ofile)
  elif re.sub(r'.*\.', '', args.ofile).lower()in ('smi', 'smiles'):
    molWriter = SmilesWriter(args.ofile, delimiter=args.delim, nameHeader='Name', includeHeader=True, isomericSmiles=True, kekuleSmiles=False)
  else:
      logging.error(f'Invalid file extension: {args.ofile}')

  if args.op=="matchCounts":
    pat = rdkit.Chem.MolFromSmarts(args.smarts)
    smarts.MatchCounts(pat, molReader, molWriter)

  elif args.op=="matchFilter":
    pat = rdkit.Chem.MolFromSmarts(args.smarts)
    smarts.MatchFilter(pat, molReader, molWriter)

  else:
    parser.error(f"Unsupported operation: {args.op}")

  logging.info('Elapsed: %s'%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))))

