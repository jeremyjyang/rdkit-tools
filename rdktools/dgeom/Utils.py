#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/
"""
#############################################################################
import os,sys,re,logging,json,time,inspect
from threading import Thread

import matplotlib as mpl
#from matplotlib import pyplot as plt

import rdkit
import rdkit.Chem
import rdkit.Chem.AllChem
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles, Draw, rdDepictor

from .. import util

#############################################################################
def GenerateConformations(mol, nconf, ffname, optiters, etol):
  '''verbose=2 means MMFF generates much data, unfortunately to stdout.'''
  verbose = 2 if logging.getLogger().getEffectiveLevel()==logging.DEBUG else 1 if logging.getLogger().getEffectiveLevel()==logging.INFO else 0
  molh = rdkit.Chem.AddHs(mol)
  confIds=rdkit.Chem.AllChem.EmbedMultipleConfs(molh, numConfs=nconf)
  ok_opts = {confId:False for confId in confIds} #optimization status
  molh_0 = rdkit.Chem.Mol(molh)
  i_conf=0;
  for confId in confIds:
    i_conf+=1
    if ffname.upper()=='MMFF':
      ff = rdkit.Chem.AllChem.MMFFGetMoleculeForceField(molh, rdkit.Chem.AllChem.MMFFGetMoleculeProperties(molh, "MMFF94", verbose), confId=confId)
    else:
      ff = rdkit.Chem.AllChem.UFFGetMoleculeForceField(molh, confId=confId)
    try:
      e_i = ff.CalcEnergy()
      rval = ff.Minimize(maxIts=optiters, energyTol=etol)
      ok_opts[confId] = bool(rval==0)
      e_f = ff.CalcEnergy()
    except Exception as e:
      logging.info('%s'%str(e))
      e_i,e_f,ok_opts[confId] = None,None,False
    try:
      rmsd = rdkit.Chem.AllChem.AlignMol(molh, molh_0, prbCid=confId, refCid=confId)
    except Exception as e:
      rmsd = 0.0
      logging.info(e)
    logging.debug("\t{i_conf:3d}. RMSD = {rmsd:.2f} ; Ei = {e_i:.2f}, Ef = {e_f:.2f}, converged = {ok_opts[confId]}")
  if not confIds:
    logging.error('EmbedMultipleConfs() failed.') 

  return molh, list(confIds)

#############################################################################
class dgeom_thread(Thread):
  """
Used by RDK DGeom CGI, to run DGeom in spawned thread.

From threading docs: (https://docs.python.org/3/library/threading.html)
No other methods (except for the constructor) should be overridden in a subclass. In
other words, only override the __init__() and run() methods of this class.
Once a thread object is created, its activity must be started by calling the thread's
start() method. This invokes the run() method in a separate thread of control.
Once the thread's activity is started, the thread is considered alive. It
stops being alive when its run() method terminates - either normally, or by raising
an unhandled exception. The is_alive() method tests whether the thread is alive.
Thread methods: start(), join(), setName(), getName(), is_alive()
  """
  def __init__ (self, ifile, ofile, ff, nconf):
    Thread.__init__(self)
    self.name='dgeom'
    self.ifile=ifile
    self.ofile=ofile
    self.nconf=nconf
    self.ff=ff
    self.n_in=0
    self.n_fail=0
    self.log=''
  def run(self):
    if self.ifile[-4:].lower() in ('.sdf','.sd','.mdl','.mol'):
      molreader = rdkit.Chem.SDMolSupplier(self.ifile)
    else:
      molreader = rdkit.Chem.SmilesMolSupplier(self.ifile, delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=False)
    molwriter = rdkit.Chem.SDWriter(self.ofile)

    for mol in molreader:
      self.n_in+=1
      mol_out = rdkit.Chem.AddHs(mol)
      confIds = rdkit.Chem.AllChem.EmbedMultipleConfs(mol_out, numConfs=self.nconf)
      mol_out_0 = rdkit.Chem.Mol(mol_out)
      i_conf=0;
      for confId in confIds:
        i_conf+=1
        ff = rdkit.Chem.AllChem.UFFGetMoleculeForceField(mol_out, confId=confId)
        e_i = ff.CalcEnergy()
        rval = ff.Minimize()
        e_f = ff.CalcEnergy()
        rmsd = rdkit.Chem.AllChem.AlignMol(mol_out, mol_out_0, prbCid=confId, refCid=confId)
        self.log+=(f"{i_conf:3d}. Ei = {e_i:.2f}, Ef = {e_f:.2f}, RMSD = {rmsd:.2f}, opt_status = {rval}\n")

      if confIds:
        for confId in list(confIds):
          molwriter.write(mol_out,confId=confId)
      else:
        self.n_fail+=1
        molwriter.write(mol)

#############################################################################
class countmols_thread(Thread):
  '''Expect mdlfile or smiles.'''
  def __init__ (self,ifile):
    Thread.__init__(self)
    self.ifile=ifile
    self.n_mol=0
  def run(self):
    if self.ifile[-4:].lower() in ('.sdf','.sd','.mdl','.mol'):
      molreader = rdkit.Chem.SDMolSupplier(self.ifile)
    else:
      molreader = rdkit.Chem.SmilesMolSupplier(self.ifile,delimiter=' ',smilesColumn=0,nameColumn=1,titleLine=False)
    for mol in molreader:
      self.n_mol+=1

#############################################################################
def DgeomProgress(th, t0, tpoll, progresswin, job='Dgeom'):
  sys.stdout.write(f"""<SCRIPT>
var pwin=window.open('','{progresswin}');
pwin.document.writeln('{job}...<BR>');
</SCRIPT>
""")
  sys.stdout.flush()
  while th.is_alive():
    sys.stdout.write(f"""<SCRIPT>
pwin.document.writeln('{job} [{time.strftime('%Hh:%Mm:%Ss',time.gmtime(time.time()-t0))}]: {th.n_in} in, {th.n_in-th.n_fail} out, {th.n_fail} failed...<BR>');
if (navigator.appName.match('Explorer')) pwin.scrollTo(0,99999);
else pwin.scrollTo(0,pwin.document.body.offsetHeight);
</SCRIPT>
""")
    sys.stdout.flush()
    th.join(tpoll)
  sys.stdout.write(f"""<SCRIPT>
pwin.document.writeln('{job} [{time.strftime('%Hh:%Mm:%Ss',time.gmtime(time.time()-t0))}]: {th.n_in} in, {th.n_in-th.n_fail} out, {th.n_fail} failed (done).<BR>');
if (navigator.appName.match('Explorer')) pwin.scrollTo(0,99999);
else pwin.scrollTo(0,pwin.document.body.offsetHeight);
</SCRIPT>
""")
  sys.stdout.flush()
  return th.n_in,th.n_in-th.n_fail,th.n_fail

#############################################################################
if __name__=='__main__':
  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG))
  logging.debug('TESTING:')
  if len(sys.argv) != 5:
    logging.debug(f'Syntax: {__name__} INFILE OUTFILE FORCEFIELD MAXCONF')
    sys.exit(1)
  infile = sys.argv[1]
  outfile = sys.argv[2]
  ff = sys.argv[3]
  maxconf = int(sys.argv[4])
  th = dgeom_thread(infile, outfile, ff, maxconf)
  t0=time.time()
  logging.debug('DEBUG: start()...')
  th.start()
  logging.debug('DEBUG: back from start()...')
  n_mol,n_ok,n_fail = DgeomProgress(th, t0, 1, 'TESTING')
  logging.debug(f"execution time: {time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))}")
  logging.debug(th.log)

