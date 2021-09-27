#!/usr/bin/env python3
#############################################################################
### rdk_mp.py  - multiprocessing
###
### Used by RDK DGeom CGI, to run DGeom in spawned processes.
#############################################################################
### Minimize parameters:
### maxIts [200]
### forceTol [1e-4]
### energyTol [1e-6]
#############################################################################
import os,re,time,sys,json,logging
from multiprocessing import Process,Manager,Pool

import rdkit.Chem
import rdkit.Chem.AllChem

import rdk_utils

#############################################################################
def dgeom_process(ifile, smifile_header, ofile, logfile, statusfile, nconf, ffname, iters, etol):
  status_tags = ['n_in','n_mol_out','n_conf_out','n_err','n_conf_converged']
  status = { tag:0 for tag in status_tags }
  if ifile[-4:].lower() in ('.sdf','.sd','.mdl','.mol'):
    molreader=rdkit.Chem.SDMolSupplier(ifile)
  else:
    molreader=rdkit.Chem.SmilesMolSupplier(ifile, delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=smifile_header)
  molwriter=rdkit.Chem.SDWriter(ofile)
  flog = open(logfile, 'w+')
  for mol in molreader:
    status['n_in']+=1
    molname = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    mol_out = rdkit.Chem.AddHs(mol)
    confIds=rdkit.Chem.AllChem.EmbedMultipleConfs(mol_out, numConfs=nconf)
    ok_opts = {confId:False for confId in confIds} #optimization status
    mol_out_0 = rdkit.Chem.Mol(mol_out)
    flog.write('%d. "%s"\n'%(status['n_in'], molname))
    i_conf=0;
    for confId in confIds:
      i_conf+=1
      ff = rdkit.Chem.AllChem.MMFFGetMoleculeForceField(mol_out, rdkit.Chem.AllChem.MMFFGetMoleculeProperties(mol_out, 'MMFF94', 0), confId=confId) if ffname.upper()=='MMFF' else rdkit.Chem.AllChem.UFFGetMoleculeForceField(mol_out, confId=confId)

      try:
        e_i = ff.CalcEnergy()
        rval = ff.Minimize(maxIts=iters, energyTol=etol)
        ok_opts[confId] = bool(rval==0)
        if rval==0: status['n_conf_converged']+=1
        e_f = ff.CalcEnergy()
      except Exception as e:
        flog.write('\tforce field exception: %s\n'%str(e))
        e_i,e_f,ok_opts[confId] = None,None,False

      try:
        rmsd = rdkit.Chem.AllChem.AlignMol(mol_out, mol_out_0, prbCid=confId, refCid=confId)
      except Exception as e:
        rmsd = 0.0
        flog.write('\talignment exception: %s\n'%str(e))

      flog.write('\t%3d. RMSD = %.2f ; Ei = %s, Ef = %s, converged = %s\n'%(i_conf, rmsd, ('%.2f'%e_i if e_i else 'None'), ('%.2f'%e_f if e_f else 'None'), str(ok_opts[confId])))

    if confIds:
      status['n_mol_out']+=1
      for confId in list(confIds):
        molwriter.write(mol_out, confId=confId)
        status['n_conf_out']+=1
    else:
      status['n_err']+=1
      status['n_mol_out']+=1
      molwriter.write(mol)
    fstat = open(statusfile, 'w+')
    fstat.write(json.dumps(status))
    fstat.close()
  flog.close()

#############################################################################
def DgeomProgress(proc, statusfile, t0, tpoll, progresswin, job='Dgeom'):
  sys.stdout.write('''<SCRIPT>
var pwin=window.open('', '%(PROGRESSWIN)s');
pwin.document.writeln('%(JOB)s...<BR>');
</SCRIPT>
'''%{   'PROGRESSWIN':progresswin,
        'JOB':job
        })
  sys.stdout.flush()
  while proc.is_alive():
    try:
      status = json.load(open(statusfile))
    except Exception as e:
      status = {}
    sys.stdout.write('''<SCRIPT>
pwin.document.writeln('%(JOB)s [%(T)s]: %(NIN)d in, %(NOUT)d out, %(NERR)d errors ...<BR>');
if (navigator.appName.match('Explorer')) pwin.scrollTo(0,99999);
else pwin.scrollTo(0,pwin.document.body.offsetHeight);
</SCRIPT>
'''%{   'T':time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0)),
	'NIN':(status['n_in'] if 'n_in' in status else 0),
	'NOUT':(status['n_mol_out'] if 'n_mol_out' in status else 0),
	'NERR':(status['n_err'] if 'n_err' in status else 0),
        'JOB':job
        })
    sys.stdout.flush()
    proc.join(tpoll)

  sys.stdout.write('''<SCRIPT>
pwin.document.writeln('%(JOB)s [%(T)s] (done).<BR>');
if (navigator.appName.match('Explorer')) pwin.scrollTo(0,99999);
else pwin.scrollTo(0,pwin.document.body.offsetHeight);
</SCRIPT>
'''%{   'T':time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0)),
        'JOB':job
        })
  sys.stdout.flush()
  return

#############################################################################
if __name__=='__main__':
  PROG=os.path.basename(sys.argv[0])
  if len(sys.argv) != 6:
    logging.info('Syntax: %s INFILE OUTFILE LOGFILE MAXCONF FORCEFIELD'%PROG)
    sys.exit(1)
  infile = sys.argv[1]
  smifile_header = False
  outfile = sys.argv[2]
  logfile = sys.argv[3]
  maxconf = int(sys.argv[4])
  ff = sys.argv[5]
  proc = Process(target=dgeom_process, args=(infile,smifile_header,outfile,logfile,maxconf,ff,200,1e-6))
  logging.info('TESTING:')
  t0=time.time()
  logging.info('DEBUG: start()...', file=sys.stderr) 
  proc.start()
  logging.info('DEBUG: back from start()...', file=sys.stderr) 
  DgeomProgress(proc,t0,1,'TESTING', PROG)
  logging.info("execution time: %s"%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))))
  logging.info("%s"%(open(logfile).read()))

