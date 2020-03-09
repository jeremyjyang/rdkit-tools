#!/usr/bin/env python3
#############################################################################
### smartstest_rdk.py
###
### Jeremy Yang
#############################################################################
import sys,os,getopt,re
import rdkit.Chem
import rdkit.Chem.AllChem

PROG=os.path.basename(sys.argv[0])

def ErrorExit(msg):
  print(msg, file=sys.stderr)
  sys.exit(1)

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
if __name__=='__main__':
  usage='''
%(PROG)s [options]

  required:
    --i IFILE ...... input molecule file (SMILES) 
    --smarts SMA ... input SMARTS

  options:
    --o OUTFILE .... csv showing smiles and match counts
    --usa .......... unique set of atoms match counts
    --v ............ verbose
    --h ............ help

'''%{'PROG':PROG}
  options=Options()
  opts,pargs=getopt.getopt(sys.argv[1:],'',['h','v','vv','i=',
  'smarts=','o=','usa'])
  if not opts: ErrorExit(usage)
  for (opt,val) in opts:
    if opt=='--h': ErrorExit(usage)
    elif opt=='--i': options.ifile=val
    elif opt=='--smarts': options.smarts=val
    elif opt=='--o': options.ofile=val
    elif opt=='--usa': options.usa=True
    elif opt=='--v': options.verbose=1
    elif opt=='--vv': options.verbose=2
    else: ErrorExit('Illegal option: %s'%val)

  if not options.smarts and options.ifile: ErrorExit('-smarts and -i required.')

  if options.ofile: fout=file(ofile,'w+')
  else: fout=sys.stdout
  fout.write('smiles,name,n_match\n')

  pat=rdkit.Chem.MolFromSmarts(options.smarts)
  if not pat:
    ErrorExit('Bad smarts: %s'%(options.smarts))

  n_mol=0; n_mol_matched=0;
  fin=file(options.ifile,'r')
  while True:
    line=fin.readline()
    if not line: break
    n_mol+=1
    line=line.rstrip()
    m=re.match('(\S*)\s+(\S.*)$',line)
    if m:
      smi,name=m.group(1),m.group(2)
    else:
      smi,name=line,''
    mol=rdkit.Chem.MolFromSmiles(smi)

    matches=mol.GetSubstructMatches(pat,uniquify=True,useChirality=False)
    if options.usa:
      matches=DeduplicateMatches(matches)
    n_matches=len(matches)
    n_mol_matched += (1 if n_matches>0 else 0)
    fout.write('"%s","%s",%d\n'%(smi,name,n_matches))
  fin.close()
  fout.close()

  print('%s: mols: %d'%(PROG,n_mol), file=sys.stderr)
  print('%s: matched: %d'%(PROG,n_mol_matched), file=sys.stderr)
