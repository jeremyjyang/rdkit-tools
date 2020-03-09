#!/usr/bin/env python3
#############################################################################
### rdk_threads.py 
###
### Used by RDK DGeom CGI, to run DGeom in spawned thread.
#############################################################################
### From threading docs: (https://docs.python.org/3/library/threading.html)
### No other methods (except for the constructor) should be overridden in a subclass. In
### other words, only override the __init__() and run() methods of this class.
### Once a thread object is created, its activity must be started by calling the thread's
### start() method. This invokes the run() method in a separate thread of control.
### Once the thread's activity is started, the thread is considered alive. It
### stops being alive when its run() method terminates - either normally, or by raising
### an unhandled exception. The is_alive() method tests whether the thread is alive.
### Thread methods: start(), join(), setName(), getName(), is_alive()
#############################################################################
import os,re,time,sys,logging
from threading import Thread

import rdkit.Chem
import rdkit.Chem.AllChem

import rdk_utils

#############################################################################
class dgeom_thread(Thread):
  def __init__ (self,ifile,ofile,ff,nconf):
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
      molreader=rdkit.Chem.SDMolSupplier(self.ifile)
    else:
      molreader=rdkit.Chem.SmilesMolSupplier(self.ifile,delimiter=' ',smilesColumn=0,nameColumn=1,titleLine=False)
    molwriter=rdkit.Chem.SDWriter(self.ofile)

    for mol in molreader:
      self.n_in+=1
      #mol_out, confIds = rdk_utils.GenerateConformations(mol,self.ff,self.nconf)

      mol_out = rdkit.Chem.AddHs(mol)
      confIds=rdkit.Chem.AllChem.EmbedMultipleConfs(mol_out,numConfs=self.nconf)
      mol_out_0 = rdkit.Chem.Mol(mol_out)
      i_conf=0;
      for confId in confIds:
        i_conf+=1
        ff = rdkit.Chem.AllChem.UFFGetMoleculeForceField(mol_out,confId=confId)
        e_i = ff.CalcEnergy()
        rval = ff.Minimize()
        e_f = ff.CalcEnergy()
        rmsd = rdkit.Chem.AllChem.AlignMol(mol_out,mol_out_0,prbCid=confId,refCid=confId)
        self.log+=('%3d. Ei = %.2f, Ef = %.2f, RMSD = %.2f, opt_status = %s\n'%(i_conf,e_i,e_f,rmsd,str(rval)))

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
      molreader=rdkit.Chem.SDMolSupplier(self.ifile)
    else:
      molreader=rdkit.Chem.SmilesMolSupplier(self.ifile,delimiter=' ',smilesColumn=0,nameColumn=1,titleLine=False)
    for mol in molreader:
      self.n_mol+=1

#############################################################################
def DgeomProgress(th,t0,tpoll,progresswin,job='Dgeom'):
  sys.stdout.write('''<SCRIPT>
var pwin=window.open('','%(PROGRESSWIN)s');
pwin.document.writeln('%(JOB)s...<BR>');
</SCRIPT>
'''%{	'PROGRESSWIN':progresswin,
        'JOB':job
        })
  sys.stdout.flush()
  while th.isAlive():
    sys.stdout.write('''<SCRIPT>
pwin.document.writeln('%(JOB)s [%(T)s]: %(NIN)d in, %(NOUT)d out, %(NFAIL)d failed...<BR>');
if (navigator.appName.match('Explorer')) pwin.scrollTo(0,99999);
else pwin.scrollTo(0,pwin.document.body.offsetHeight);
</SCRIPT>
'''%{   'T':time.strftime('%Hh:%Mm:%Ss',time.gmtime(time.time()-t0)),
        'NIN':th.n_in,
        'NOUT':th.n_in-th.n_fail,
        'NFAIL':th.n_fail,
        'JOB':job
        })
    sys.stdout.flush()
    th.join(tpoll)
  sys.stdout.write('''<SCRIPT>
pwin.document.writeln('%(JOB)s [%(T)s]: %(NIN)d in, %(NOUT)d out, %(NFAIL)d failed (done).<BR>');
if (navigator.appName.match('Explorer')) pwin.scrollTo(0,99999);
else pwin.scrollTo(0,pwin.document.body.offsetHeight);
</SCRIPT>
'''%{   'T':time.strftime('%Hh:%Mm:%Ss',time.gmtime(time.time()-t0)),
        'NIN':th.n_in,
        'NOUT':th.n_in-th.n_fail,
        'NFAIL':th.n_fail,
        'JOB':job
        })
  sys.stdout.flush()
  return th.n_in,th.n_in-th.n_fail,th.n_fail

#############################################################################

if __name__=='__main__':
  PROG=os.path.basename(sys.argv[0])
  logging.info('TESTING:')
  if len(sys.argv) != 5:
    logging.info('Syntax: %s INFILE OUTFILE FORCEFIELD MAXCONF'%PROG)
    sys.exit(1)
  infile = sys.argv[1]
  outfile = sys.argv[2]
  ff = sys.argv[3]
  maxconf = int(sys.argv[4])

  th = dgeom_thread(infile,outfile,ff,maxconf)

  t0=time.time()
  logging.info('DEBUG: start()...')
  th.start()
  logging.info('DEBUG: back from start()...')
  n_mol,n_ok,n_fail = DgeomProgress(th,t0,1,'TESTING',PROG)
  logging.info("%s: execution time: %s"%(PROG, time.strftime('%Hh:%Mm:%Ss',time.gmtime(time.time()-t0))))
  logging.info("%s"%(th.log))

