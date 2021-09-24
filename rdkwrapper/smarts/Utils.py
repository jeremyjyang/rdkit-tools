#!/usr/bin/env python3
#############################################################################
import sys,os,argparse,re,logging
import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.rdBase

#############################################################################
### DeduplicateMatches() - remove matches not USA (Unique Sets of Atoms)
### uumatches - unique set of matches, sorted atom order
### umatches - unique set of matches, original atom order
#############################################################################
def DeduplicateMatches(matches):
  uumatches=[]
  umatches=[]
  for match in matches:
    uumatch=list(match)
    uumatch.sort()
    if uumatch not in uumatches:
      uumatches.append(uumatch)
      umatches.append(match)
  return tuple(umatches)

#############################################################################
class Options(object):
  def __init__(self):
    self.ifile=None;
    self.smarts=None;
    self.ofile=None;
    self.usa=False;
    self.verbose=0;

#############################################################################
def Demo():
  smas=[
	'*~1~*~*~*~*~*1',
	'*~1~*~*~*~*~*~*~*~*~*1',
	'*~1~*~*~*~*~*~*~*~*~*~*~*~*~*1',
  ]
  smis=[
	'C1CCCCC1',
	'C1CCCC2C1CCCC2',
	'N1CCCC2C1CCCC2',
	'C1CCC3C4CCCCC4CCC3C1',
	'CC1CCC3C4CCCCC4CCC3C1',
  ]

  for smi in smis:
    logging.info('smi: %s'%smi)
    for sma in smas:
      logging.info('\tsma: %s'%sma)
      pat = rdkit.Chem.MolFromSmarts(sma)
      mol = rdkit.Chem.MolFromSmiles(smi)
      matches = mol.GetSubstructMatches(pat,uniquify=True,useChirality=False)
      for match in matches: logging.info('\t',match)
      nmatches=len(matches),
      logging.info('\tmatches: %d'%nmatches,)
      umatches = DeduplicateMatches(matches)
      nmatches_usa = len(umatches),
      logging.info(f"""\tusa matches: {nmatches_usa}{'<-- LOOK !!!' if nmatches!=nmatches_usa else ''}""")

#############################################################################
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='RDKit SMARTS utility')
  parser.add_argument("--i", required=True, dest="ifile", help="input molecule file")
  parser.add_argument("--smarts", required=True, help="input SMARTS")
  parser.add_argument("--usa", action="store_true", help="unique set of atoms match counts")

  parser.add_argument("--o", dest="ofile", help="output file (TSV)")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f'RDKit version: {rdkit.rdBase.rdkitVersion}')

  fout= open(ofile,'w+') if args.ofile else sys.stdout
  fout.write('smiles\tname\tn_match\n')

  pat = rdkit.Chem.MolFromSmarts(args.smarts)
  if not pat:
    logging.error('Bad smarts: %s'%(args.smarts))

  n_mol=0; n_mol_matched=0;
  fin = open(args.ifile, 'r')
  while True:
    line=fin.readline()
    if not line: break
    n_mol+=1
    line=line.rstrip()
    m = re.match('(\S*)\s+(\S.*)$',line)
    if m:
      smi,name = m.group(1),m.group(2)
    else:
      smi,name = line,''
    mol = rdkit.Chem.MolFromSmiles(smi)

    matches = mol.GetSubstructMatches(pat,uniquify=True,useChirality=False)
    if args.usa:
      matches = DeduplicateMatches(matches)
    n_matches = len(matches)
    n_mol_matched += (1 if n_matches>0 else 0)
    fout.write(f"{smi}\t{name}\t{n_matches}\n")
  fin.close()
  fout.close()

  logging.info('mols: %d'%(n_mol))
  logging.info('matched: %d'%(n_mol_matched))
