#!/usr/bin/env python3
#############################################################################
### purgescratchdirs.py -- for purging temporary files created by
### webapps, where the names indicate creation time.
#############################################################################
import sys,os,re,time
#
PROG = os.path.basename(sys.argv[0])
#
#############################################################################
def PurgeScratchDirs(dirlist,retire_secs,verbose):
  ok=True
  exp_date=time.strftime('%Y%m%d%H%M%S',time.localtime(time.time()-retire_secs))
  #
  ndeleted=0
  ntotal=0
  for dir in dirlist:
    filenames=os.listdir(dir)
    for filename in filenames:
      ntotal+=1
      path=dir+'/'+filename
      ### should be YYYYMMDDHHMMSS
      date = filename.split('.')[1] if '.' in filename else filename
      if not len(date)==14 or not re.match('\d\d\d\d\d\d\d\d\d\d\d\d\d\d$',date):
        if verbose:
          sys.stderr.write('%s: ERROR: bad date format: "%s"\n'%(PROG,date))
          ok=False
        date=time.strftime('%Y%m%d%H%M%S',time.localtime(os.stat(path).st_ctime))
        if date<exp_date:
          if verbose:
            sys.stderr.write('%s: WARNING: stale ctime: "%s"\n'%(PROG,path))
          ok=False
        continue
      if date<exp_date:	### ascii compare is fine
        if not os.access(path,os.F_OK):
          if verbose:
            sys.stderr.write('%s: ERROR: cannot find "%s"\n'%(PROG,path))
          ok=False
        if os.access(path,os.W_OK):
          os.unlink(path)
          if verbose:
            sys.stdout.write('%s: deleted "%s"\n'%(PROG,path))
          ndeleted+=1
        else:
          if verbose:
            sys.stderr.write('%s: ERROR: cannot delete "%s"\n'%(PROG,path))
          ok=False
      else:
        if verbose:
          sys.stdout.write('%s: not deleting "%s"\n'%(PROG,path))
  
  if verbose:
    sys.stdout.write('%s: %d / %d deleted\n'%(PROG,ndeleted,ntotal))
  return ok

#############################################################################
if __name__=='__main__':
  SCRATCHDIRS = [ os.environ['HOME']+'/web/cgi-bin/omegascratch',
                  os.environ['HOME']+'/web/cgi-bin/rocsscratch',
                  os.environ['HOME']+'/web/cgi-bin/fredscratch']
  RETIRE_HRS=1	### Beyond this age files deleted.
  retire_secs=RETIRE_HRS*3600
  verbose=1
  ok=PurgeScratchDirs(SCRATCHDIRS,retire_secs,verbose)
  if not ok: sys.stdout.write('%s: ERROR detected\n'%PROG)
