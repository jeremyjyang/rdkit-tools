#! /usr/bin/env python3
"""
Front end for RDKit conformer generation via distance geometry.
Ref: http://www.rdkit.org/Python_Docs/
"""
import os,sys,stat,cgi,re,time,tempfile,random,urllib,base64,json
from multiprocessing import Process,Manager,Pool

import rdkit.Chem
import rdkit.rdBase
import rdkit.Chem.AllChem
  
import rdk_utils
import rdk_mp

import env_cgi
import mol2imghtm
import htm_utils
import purgescratchdirs

#############################################################################
def JavaScript(progname):
  return '''\
function go_run_dgeom(form)
{
  if (!checkform(form)) return;
  var x,y;
  if (!navigator.appName.match('Explorer'))
  {
    x=window.screenX+300;
    y=window.screenY+300;
  }
  else
  {
    x=window.screenLeft+300;
    y=window.screenTop+300;
  }
  if (form.runmode.value=='multi')
  {
    var pwin=window.open('','%(PROGRESS_WIN_NAME)s',
	'width=400,height=100,left='+x+',top='+y+',scrollbars=1,resizable=1,location=0,status=0,toolbar=0');
    if (!pwin) {
      alert('ERROR: popup windows must be enabled for progress indicator.');
      return false;
    }
    pwin.focus();
    pwin.document.close(); //if window exists, clear
    pwin.document.open('text/html');
    pwin.document.writeln('<HTML><HEAD>');
    pwin.document.writeln('<LINK REL=\"stylesheet\" type=\"text/css\" HREF=\"%(HTML_SUBDIR)s/css/biocomp.css\" />');
    pwin.document.writeln('</HEAD><BODY BGCOLOR=\"#DDDDDD\">');
    if (navigator.appName.match('Explorer'))
      pwin.document.title='%(PROG)s progress'; //not-ok for IE
  }
  form.run_dgeom.value='TRUE';
  form.submit();
}
function checkform(form)
{
  if (!form.intxt.value && !form.infile.value) {
    alert('ERROR: No input specified');
    return false;
  }
  if (!form.maxconf.value) {
    alert('ERROR: maxconf must be specified.');
    return false;
  }
  return true;
}
function go_init(form)
{
  var i;
  for (i=0;i<form.ifmt.length;++i)
    if (form.ifmt.options[i].value==1) // smi
      form.ifmt.options[i].selected=true;
  form.maxconf.value='1';
  form.optiters.value='400';
  form.etol.value='1e-06';
  form.ff.value='mmff';
  form.smifile_header.checked=false;
  form.verbose.checked=false;
}
/// JSME stuff:
function StartJSME()
{
  window.open('/jsme_win.html','JSME','width=500,height=450,scrollbars=0,resizable=1,location=0');
}
function fromJSME(smiles) // function called from JSME window
{
  if (smiles=='')
  {
    alert("ERROR: no molecule submitted");
    return;
  }
  var form=document.mainform;
  form.intxt.value=smiles;
  for (i=0;i<form.ifmt.length;++i)
  {
    if (form.ifmt.options[i].value==1) form.ifmt.options[i].selected=true;
  }
}
'''%{'PROG':progname,'PROGRESS_WIN_NAME':PROGRESS_WIN_NAME,'HTML_SUBDIR':env_cgi.HTML_SUBDIR}

#############################################################################
def PrintForm():
  ifmtmenu='<SELECT NAME="ifmt">\n'
  for fmt in ('smiles','sdf'):
    if fmt==IFMT: s=' SELECTED'
    else: s=''
    ifmtmenu+='<OPTION VALUE="%s"%s>%s'%(fmt,s,FMTS[fmt])
  ifmtmenu+='</SELECT>'

  RUNMODE_SINGLE=''; RUNMODE_MULTI='';
  if RUNMODE=='single': RUNMODE_SINGLE='CHECKED';
  elif RUNMODE=='multi': RUNMODE_MULTI='CHECKED';

  FF_MMFF=''; FF_UFF='';
  if FF=='mmff': FF_MMFF='CHECKED';
  elif FF=='uff': FF_UFF='CHECKED';

  print('<FORM NAME="mainform" ACTION="%s" METHOD="POST" ENCTYPE="multipart/form-data">'%(CGIURL))
  print('<INPUT TYPE=HIDDEN NAME="run_dgeom">')
  print('<TABLE WIDTH="100%%"><TR><TD><H2>%s</H2></TD>'%APPNAME)
  print('<TD><TD>- RDKit distance-geometry conformer generation</TD>')
  print('<TD ALIGN="right">')
  print('<BUTTON TYPE=BUTTON onClick="void window.open(\'%s?help=TRUE\',\'helpwin\',\'width=600,height=400,scrollbars=1,resizable=1\')">'%CGIURL)
  print('<B>Help</B></BUTTON>')
  print('<BUTTON TYPE=BUTTON onClick="window.location.replace(\'%s\')"><B>Reset</B></BUTTON>&nbsp;'%CGIURL)
  print('</TD></TR></TABLE>')
  print('<HR>')
  print('<TABLE WIDTH=100% CELLPADDING=5>')
  print('<TR>')
  print('<TD VALIGN=TOP>')
  print('input format:%s \n'%ifmtmenu)
  print('file2txt:<INPUT TYPE=CHECKBOX NAME="infile2txt" VALUE="CHECKED"')
  print(' %s><BR>'%INFILE2TXT)
  print('upload:<INPUT TYPE="FILE" NAME="infile">&nbsp;')
  print('or paste...<BR>')
  print('<TEXTAREA NAME="intxt" ROWS=12 COLS=50 WRAP=OFF>%s</TEXTAREA><BR>'%INTXT)
  print('or draw:')
  print('<BUTTON TYPE=BUTTON onClick="StartJSME()">JSME</BUTTON>')
  print('</TD>')
  print('<TD VALIGN=TOP>')
  print('<TABLE WIDTH="100%" BGCOLOR="#EEEEEE" CELLPADDING=2>')
  print('<TR><TD COLSPAN=2><B>options:</B></TD></TR>')
  print('<TR><TD ALIGN="RIGHT" VALIGN="TOP">mode:</TD><TD><INPUT TYPE=RADIO NAME="runmode" VALUE="single" %s>single <I>(interactive)</I><BR>'%(RUNMODE_SINGLE))
  print('<INPUT TYPE=RADIO NAME="runmode" VALUE="multi" %s>multi <I>(batch)</I></TD></TR>'%(RUNMODE_MULTI))
  print('<TR><TD ALIGN=RIGHT>maxconf:</TD>')
  print('<TD><INPUT TYPE=TEXT SIZE=3 NAME="maxconf" VALUE="%s"></TD></TR>'%MAXCONF)
  print('<TR><TD ALIGN="RIGHT">forcefield:</TD><TD><INPUT TYPE=RADIO NAME="ff" VALUE="mmff" %s>MMFF</I>'%(FF_MMFF))
  print('<INPUT TYPE=RADIO NAME="ff" VALUE="uff" %s>UFF</TD></TR>'%(FF_UFF))
  print('<TR><TD ALIGN=RIGHT>optimizer iters:</TD>')
  print('<TD><INPUT TYPE=TEXT SIZE=3 NAME="optiters" VALUE="%s"></TD></TR>'%OPTITERS)
  print('<TR><TD ALIGN=RIGHT>optimizer E-tolerance:</TD>')
  print('<TD><INPUT TYPE=TEXT SIZE=5 NAME="etol" VALUE="%s"></TD></TR>'%ETOL)
  print('<TR><TD COLSPAN=2><B>misc:</B></TD></TR>')
  print('<TR><TD ALIGN=RIGHT>smifile header:</TD>')
  print('<TD ALIGN=LEFT><INPUT TYPE=CHECKBOX NAME="smifile_header" VALUE="CHECKED" %s></TD></TR>'%SMIFILE_HEADER)
  print('<TR><TD ALIGN=RIGHT>verbose:</TD>')
  print('<TD ALIGN=LEFT><INPUT TYPE=CHECKBOX NAME="verbose" VALUE="CHECKED" %s></TD></TR>'%VERBOSE)
  print('</TABLE>')
  print('</TD></TR>')
  print('<TR><TD COLSPAN=2 ALIGN=CENTER>')
  print('<BUTTON TYPE=BUTTON onClick="go_run_dgeom(this.form)"')
  print('><B>Go %s</B></BUTTON>'%APPNAME)
  print('</TD></TR>')
  print('</TABLE>')
  print('</FORM>')

#############################################################################
def DgeomSingle(imol,outfile,maxconf,ff,optiters,etol):
  if not imol:
    OUTPUTS.append("<H2>ERROR: no input molecule</H2>")
    return
  OUTPUTS.append("<H2>input:</H2>")
  molname = imol.GetProp('_Name') if imol.HasProp('_Name') else ''
  thtm=('<TABLE><TR><TD VALIGN=TOP><TT>%s</TT></TD></TR>'%molname)
  ismi=rdkit.Chem.MolToSmiles(imol,isomericSmiles=True)
  thtm+=("<TR><TD>%s</TD></TR>"%(ismi))
  imghtm=(mol2imghtm.Smi2ImgHtm((ismi+' '+molname if molname else ismi),'',180,220,SMI2IMG))
  thtm+=('<TR><TD BGCOLOR="white">%s</TD></TR></TABLE>'%(imghtm))
  OUTPUTS.append("<BLOCKQUOTE>%s</BLOCKQUOTE>"%thtm)

  OUTPUTS.append('<H2>output:</H2>')

  molwriter=rdkit.Chem.SDWriter(outfile)
  outmol, confIds = rdk_utils.GenerateConformations(imol,maxconf,ff,optiters,etol)
  for confId in confIds:
    molwriter.write(outmol,confId=confId)
  if not confIds:
    OUTPUTS.append('<B>Dgeom failed.</B>')
    return

  filecode=urllib.parse.quote(SCRATCHDIRURL+'/'+os.path.basename(outfile), safe='')
  bhtm=('<FORM><BUTTON TYPE=BUTTON')
  bhtm+=(" onClick=\"go_view3d('%s','%s',%d,'view3d_dgeom','multiconf')\">"%(VIEW3D,filecode,600))
  bhtm+=('<IMG SRC="/images/Jmol_icon_128.png" HEIGHT="60" VALIGN="middle">JSmol; ')
  bhtm+=('view output (%d confs)</BUTTON></FORM>'%(len(confIds)))
  OUTPUTS.append('<BLOCKQUOTE>%s</BLOCKQUOTE>'%bhtm)

  fname='dgeom_out.sdf'
  bhtm=('<FORM ACTION="%s/%s" METHOD="POST">'%(CGIURL,fname))
  bhtm+=('<INPUT TYPE=HIDDEN NAME="downloadfile" VALUE="%s">'%outfile)
  bhtm+=('<BUTTON TYPE=BUTTON onClick="this.form.submit()">')
  bhtm+=('%s (%s)</BUTTON></FORM>'%(fname,htm_utils.NiceBytes(os.stat(outfile).st_size)))
  OUTPUTS.append('<BLOCKQUOTE><b>Download:</b> %s</BLOCKQUOTE>'%bhtm)

  DATE=time.strftime('%Y%m%d%H%M',time.localtime())
  LOG.write('%s\t%s\t%d\n' % (DATE,os.environ['REMOTE_ADDR'],1))
  LOG.close()

#############################################################################
def DgeomLaunchProcess(infile,smifile_header,outfile,logfile,statusfile,maxconf,ff,optiters,etol,job):
  proc = Process(target=rdk_mp.dgeom_process, args=(infile,smifile_header,outfile,logfile,statusfile,maxconf,ff,optiters,etol))
  t0=time.time()
  proc.start()
  rdk_mp.DgeomProgress(proc,statusfile,t0,1,PROGRESS_WIN_NAME,job)
  ERRORS.append("%s execution time: %s"%(job, time.strftime('%Hh:%Mm:%Ss',time.gmtime(time.time()-t0))
))
  if VERBOSE:
    ERRORS.append("<PRE>%s</PRE>"%(open(logfile).read()))
  return

#############################################################################
def ParseStatusFile(statusfile):
  try:
    status = json.load(open(statusfile))
  except Exception as e:
    status = {}
  return status

#############################################################################
def DgeomResultsMulti(status,logfile):
  if OUTPUT_TRUNCATED:
    OUTPUTS.append('Warning: output truncated; limit reached NMAX=%d.'%NMAX)

  n_in = status['n_in'] if 'n_in' in status else 0
  n_mol_out = status['n_mol_out'] if 'n_mol_out' in status else 0
  n_conf_out = status['n_conf_out'] if 'n_conf_out' in status else 0
  n_conf_converged = status['n_conf_converged'] if 'n_conf_converged' in status else 0
  n_err = status['n_err'] if 'n_err' in status else 0
  OUTPUTS.append('<H2>output:</H2>')
  OUTPUTS.append('mols in: %d'%n_in)
  OUTPUTS.append('mols out: %d'%n_mol_out)
  OUTPUTS.append('confs out: %d'%n_conf_out)
  OUTPUTS.append('confs converged: %d ; NOT converged: %d'%(n_conf_converged,n_conf_out-n_conf_converged))
  OUTPUTS.append('errors: %d'%n_err)

  OUTPUTS.append('<H2>download:</H2>')

  fname=('dgeom_out.sdf')
  bhtm=('<FORM ACTION="%s/%s" METHOD="POST">'%(CGIURL,fname))
  bhtm+=('<INPUT TYPE=HIDDEN NAME="downloadfile" VALUE="%s">'%OUTFILE)
  bhtm+=('<BUTTON TYPE=BUTTON onClick="this.form.submit()">')
  bhtm+=('%s (%s)</BUTTON></FORM>'%(fname,htm_utils.NiceBytes(os.stat(OUTFILE).st_size)))
  OUTPUTS.append('<BLOCKQUOTE><b>%s</b> - output conformations (SDF)</BLOCKQUOTE>'%bhtm)

  fname=('dgeom_log.txt')
  bhtm=('<FORM ACTION="%s/%s" METHOD="POST">'%(CGIURL,fname))
  bhtm+=('<INPUT TYPE=HIDDEN NAME="downloadfile" VALUE="%s">'%logfile)
  bhtm+=('<BUTTON TYPE=BUTTON onClick="this.form.submit()">')
  bhtm+=('%s (%s)</BUTTON></FORM>'%(fname,htm_utils.NiceBytes(os.stat(logfile).st_size)))
  OUTPUTS.append('<BLOCKQUOTE><b>%s</b> - log file</BLOCKQUOTE>'%bhtm)

  DATE=time.strftime('%Y%m%d%H%M',time.localtime())
  LOG.write('%s\t%s\t%d\n' % (DATE,os.environ['REMOTE_ADDR'],n_in))
  LOG.close()

#############################################################################
def Initialize():
  global FORM,VERBOSE,INTXT,FIXTXT,LOG,INFILE2TXT,APPNAME
  global ERRORS,OUTPUTS,CGIURL,PROG,FMTS,FMTS3D,IFMT
  global SMI2IMG,VIEW3D,FIXFMT,SMIFILE_HEADER
  global RUNMODE
  ERRORS=[]; OUTPUTS=[];
  FORM=cgi.FieldStorage(keep_blank_values=1)
  CGIURL=os.environ['SCRIPT_NAME']
  PROG=os.path.basename(sys.argv[0])
  APPNAME='DGeom'

  logohtm="<TABLE CELLSPACING=5 CELLPADDING=5>\n<TR><TD>"
  href="http://medicine.unm.edu/informatics/"
  imghtm=('<IMG SRC="'+env_cgi.HTML_SUBDIR+'/images/biocomp_logo_only.gif">')
  tiphtm=('%s web app from UNM Translational Informatics Divsion'%APPNAME)
  logohtm+=(htm_utils.HtmTipper(imghtm,'',tiphtm,href,200,"white"))

  href="http://www.rdkit.org"
  imghtm=('<IMG HEIGHT="70" SRC="'+env_cgi.HTML_SUBDIR+'/images/rdkit_logo_only.png">')
  tiphtm=('RDKit: Cheminformatics and Machine Learning Software')
  logohtm+=(htm_utils.HtmTipper(imghtm,'',tiphtm,href,200,"white"))

  href="http://www.jmol.org"
  imghtm=('<IMG HEIGHT="60" SRC="'+env_cgi.HTML_SUBDIR+'/images/Jmol_icon_128.png">')
  tiphtm=('JSmol 3D viewer from the Jmol project')
  logohtm+=(htm_utils.HtmTipper(imghtm,'',tiphtm,href,200,"white"))

  logohtm+="</TD></TR></TABLE>"
  ERRORS.append("<CENTER>"+logohtm+"</CENTER>")

  VERBOSE=FORM.getvalue('verbose')
  INTXT=FORM.getvalue('intxt','')
  FIXTXT=FORM.getvalue('fixtxt','')
  INFILE2TXT=FORM.getvalue('infile2txt')
  URLHOST=os.environ['SERVER_NAME']
  urldir=os.path.dirname(os.environ['REQUEST_URI'])
  SMI2IMG=urldir+'/mol2img.cgi'
  VIEW3D=urldir+'/view3d_jsmol.cgi'
  SMIFILE_HEADER=FORM.getvalue('smifile_header')

  global NMAX,NOLIMIT,OUTPUT_TRUNCATED
  NOLIMIT=FORM.getvalue('nolimit')
  NMAX=2000
  if NOLIMIT and env_cgi.ENABLE_NOLIMIT:
    NMAX=0
    MAXFILESIZE=0
  OUTPUT_TRUNCATED=False
  RUNMODE=FORM.getvalue('runmode','single')

  global FF,MAXCONF,OPTITERS,ETOL
  FF=FORM.getvalue('ff')
  MAXCONF=int(FORM.getvalue('maxconf','1'))
  OPTITERS=int(FORM.getvalue('optiters','400'))
  ETOL=float(FORM.getvalue('etol','1e-6'))

  IFMT=FORM.getvalue('ifmt','smiles')
  FMTS={
	'smi':'SMILES',
	'smiles':'SMILES',
	'mdl':'MDL Molfile',
	'mol':'MDL Molfile',
	'sdf':'MDL SDfile'
	}

  global PROGRESS_WIN_NAME
  PROGRESS_WIN_NAME='dgeom_progress_win'

  logfile = './logs/'+PROG+'.log'	# web-app logfile
  logfields=['ip','N']
  if os.access(logfile,os.R_OK):
    tmp=open(logfile)
    lines=tmp.readlines()
    tmp.close()
    if len(lines)>1:
      date=lines[1].split('\t')[0][:8]
      date='%s/%s/%s'%(date[4:6],date[6:8],date[:4])
      ERRORS.append('%s has been used %d times since %s.'%(PROG,len(lines)-1,date))
  if os.access(logfile,os.W_OK):
    LOG=open(logfile,'a')
    if os.path.getsize(logfile)==0:
      LOG.write('date\t' + '\t'.join(logfields) + '\n')
  else:
    LOG=open(logfile,'w+')
    LOG.write('date\t' + '\t'.join(logfields) + '\n')
  if not LOG:
    ERRORS.append('LOGFILE error.')
    return False

  ERRORS.append('RDKit version: %s'%rdkit.rdBase.rdkitVersion)
  ERRORS.append('Boost version: %s'%rdkit.rdBase.boostVersion)

  global SCRATCHDIR,SCRATCHDIRURL,SCRATCHDIR_LIFETIME
  SCRATCHDIR=env_cgi.scratchdir
  SCRATCHDIRURL=env_cgi.scratchdirurl
  SCRATCHDIR_LIFETIME=env_cgi.scratchdir_lifetime

  global TS,PREFIX
  TS=time.strftime('%Y%m%d%H%M%S',time.localtime())
  PREFIX=SCRATCHDIR+'/'+re.sub('\..*$','',PROG)+'.'+TS+'.'+str(random.randint(100,999))

  ### Assure that writeable tmp directory exists.
  if not os.access(SCRATCHDIR,os.F_OK):
    ERRORS.append('ERROR: missing scratch dir "%s"'%SCRATCHDIR)
    return False
  elif not os.access(SCRATCHDIR,os.W_OK):
    ERRORS.append('ERROR: non-writable scratch dir "%s"'%SCRATCHDIR)
    return False

  if not FORM.getvalue('run_dgeom'):
    return True

  ### NOW STUFF FOR A RUN ############################################

  global INFILE,OUTFILE,LOGFILE_DGEOM,STATUSFILE_DGEOM
  os.chdir(SCRATCHDIR)
  tempfile.tempdir=os.getcwd()
  fd,INFILE=tempfile.mkstemp('.'+IFMT,PREFIX)
  os.close(fd)
  fd,OUTFILE=tempfile.mkstemp('.sdf',PREFIX)
  os.close(fd)
  fd,LOGFILE_DGEOM=tempfile.mkstemp('.log',PREFIX)
  os.close(fd)
  fd,STATUSFILE_DGEOM=tempfile.mkstemp('.status',PREFIX)
  os.close(fd)

  global DELFILES
  DELFILES=[INFILE,LOGFILE_DGEOM,STATUSFILE_DGEOM]

  f = open(INFILE,'w+')
  intxt=''
  if FORM['infile'].file and FORM['infile'].filename:
    intxt = FORM['infile'].value.decode("utf-8") #bytes
    if FORM.getvalue('infile2txt'): INTXT=intxt
  elif FORM['intxt'].value:
    intxt = FORM['intxt'].value
  else:
    ERRORS.append('ERROR: No input data.')
    return False
  intxt = re.sub('\t',' ',intxt) ##Fix TABs
  intxt = re.sub('\r\n*','\n',intxt) ##Fix CRs
  f.write(intxt)
  f.close()

  if IFMT in ('smi','smiles'):
    molreader=rdkit.Chem.SmilesMolSupplier(INFILE,delimiter=' ',smilesColumn=0,nameColumn=1,titleLine=bool(SMIFILE_HEADER),sanitize=True)
  else:
    molreader=rdkit.Chem.SDMolSupplier(INFILE,sanitize=True,removeHs=True)

  global IMOL
  for mol in molreader:
    try:
      IMOL = rdkit.Chem.Mol(mol)
    except Exception as e:
      ERRORS.append(str(e))
      return False
    break
  molreader.reset()

  if not IMOL:
    ERRORS.append('ERROR: No input data.')
    return False

  if VERBOSE:
    ERRORS.append('server: '+(', '.join(os.uname())))

  return True

#############################################################################
def Help():
  print('''\
<HR>
<H3>%(APPNAME)s help</H3>
<P>
RDKit distance-geometry conformer generation.
</P><P>
<B>User controls:</B>
<UL>
<LI><B>Single mode</B> processes one molecule, with detailed output and JSmol interactive visualization and
modeling.  <B>Multi mode</B> processes multiple molecules in batch mode, with summary output and downloadable
output files.
<LI><B>Maxconfs</B> specifies the attempted and maximum number of conformations per molecule.
<LI><B>Force field</B> specifies algorithm used for optimization (energy minimization).
<LI><B>Optimization iters</B> specifies iterations of the FF algorithm, and may or may not result in convergence to
<LI><B>Energy tolerance</B> defines convergence as an iteration where energy change within tolerance.
</UL>
</P><P>

<B>Force fields:</B>
<UL>
<LI>MMFF: Merck Molecular Force Field (MMFF94)<BR>
Usually more accurate than UFF, but fails if atom types undefined (e.g. boron).
<LI>UFF: Universal Force Field<BR>
Usually faster than MMFF, applicable to more atom types.
</UL>
</P>

<B>References:</B>
<UL>
<LI> <A HREF="http://www.rdkit.org/">RDKit</A>.
<LI> <A HREF="http://www.rdkit.org/Python_Docs/">RDKit Python API</A>.
<LI> <A HREF="http://jmol.sourceforge.net">JSmol</A>.

<LI> Blaney, JM, Dixon, JS, "Distance Geometry in Molecular Modeling", Reviews in
Computational Chemistry; VCH: New York, 1994.

<LI> Blaney, JM, "DGEOM", Program 159 of the Quantum Chemical Program, Exchange, University of Indiana,
Bloomington, IN (1990).

<LI> Rubicon Manual, Daylight Chemical Information Systems, 
http://www.daylight.com/dayhtml/doc/rubicon/.

<LI>Wenger, JC, Smith, DH, J. Chem. Inf. Comp. Sci., 1982, 22, 29.

<LI> Crippen, GM, Distance Geometry and Conformational Calculations; Bawden, D., Ed.;
Research Studies Press (Wiley): New York, 1981.

<LI> Crippen, GM, Havel, TF, Distance Geometry and Molecular Conformation; Bawden, D,
Ed.; Research Studies Press (Wiley): New York, 1988.

<LI>"Molecular modeling software and methods for medicinal chemistry", Cohen, NC,
Blaney, JM, Humblet, C, Gund, P,  Barry, DC, J. Med. Chem., 1990, 33(3), pp883-894.

<LI>Halgren, TA, "Merck molecular force field. I. Basis, form, scope, parameterization, and performance of
MMFF94", J Comp Chem., 17:490-19 (1996).

<LI>Rapp&eacute;, AK, Casewit, CJ, Colwell, KS, Goddard III, WA, Skiff,
WM, "UFF, a full periodic table force field for molecular mechanics and molecular
dynamics simulations", J. Am. Chem. Soc. 114:10024-35 (1992).

<LI> "Freely available conformer generation methods: how good are they?",
Ebejer JP1, Morris GM, Deane CM, J Chem Inf Model. 2012 May 25;52(5):1146-58, doi: 10.1021/ci2004658,
Epub 2012 Apr 19 (http://www.ncbi.nlm.nih.gov/pubmed/22482737).

</UL>
<P>
NMAX = %(NMAX)d<BR>
Depictions by %(DEPICT_TOOL)s.<BR>
''' % {	'APPNAME':APPNAME,
	'NMAX':NMAX,
	'DEPICT_TOOL':env_cgi.DEPICT_TOOL
	})

#############################################################################
if __name__=='__main__':
  ok=Initialize()
  if not ok:
    htm_utils.PrintHeader(APPNAME,JavaScript(PROG))
    htm_utils.PrintFooter(ERRORS)
  elif FORM.getvalue('run_dgeom'):
    htm_utils.PrintHeader(APPNAME,JavaScript(PROG))
    PrintForm()
    if RUNMODE=='single':
      DgeomSingle(IMOL,OUTFILE,MAXCONF,FF,OPTITERS,ETOL)
    else:
      DgeomLaunchProcess(INFILE,bool(SMIFILE_HEADER),OUTFILE,LOGFILE_DGEOM,STATUSFILE_DGEOM,MAXCONF,FF,OPTITERS,ETOL,'Dgeom')
      status = ParseStatusFile(STATUSFILE_DGEOM)
      DgeomResultsMulti(status,LOGFILE_DGEOM)
      print('<SCRIPT>pwin.parent.focus(); pwin.focus(); pwin.close();</SCRIPT>')
    htm_utils.PrintOutput(OUTPUTS)
    htm_utils.PrintFooter(ERRORS)
    htm_utils.Cleanup(DELFILES)
    purgescratchdirs.PurgeScratchDirs([SCRATCHDIR],SCRATCHDIR_LIFETIME,0)
  elif FORM.getvalue('downloadtxt'):
    htm_utils.DownloadString(FORM.getvalue('downloadtxt'))
  elif FORM.getvalue('downloadfile'):
    htm_utils.DownloadFile(FORM.getvalue('downloadfile'))
  elif FORM.getvalue('help'):
    htm_utils.PrintHeader(APPNAME,JavaScript(PROG))
    Help()
    htm_utils.PrintFooter(ERRORS)
  elif FORM.getvalue('test'):
    htm_utils.PrintTestText(APPNAME,{})
  else:
    htm_utils.PrintHeader(APPNAME,JavaScript(PROG))
    PrintForm()
    print('<SCRIPT>go_init(window.document.mainform);</SCRIPT>')
    htm_utils.PrintFooter(ERRORS)
