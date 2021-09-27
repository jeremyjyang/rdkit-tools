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
  
from .. import util

#############################################################################
def JavaScript(progname):
  return """\
var DEMO_SMI='c1(ccc(cc1)OC(c2ccccc2)CCNC)C(F)(F)F';
function go_dgeom(form)
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
    var pwin=window.open('','"""+PROGRESS_WIN_NAME+"""',
	'width=400,height=100,left='+x+',top='+y+',scrollbars=1,resizable=1,location=0,status=0,toolbar=0');
    if (!pwin) {
      alert('ERROR: popup windows must be enabled for progress indicator.');
      return false;
    }
    pwin.focus();
    pwin.document.close(); //if window exists, clear
    pwin.document.open('text/html');
    pwin.document.writeln('<HTML><HEAD>');
    pwin.document.writeln('<LINK REL=\"stylesheet\" type=\"text/css\" HREF=\""""+util.http.env_cgi.HTML_SUBDIR+"""/css/biocomp.css\" />');
    pwin.document.writeln('</HEAD><BODY BGCOLOR=\"#DDDDDD\">');
    if (navigator.appName.match('Explorer'))
      pwin.document.title='"""+progname+""" progress'; //not-ok for IE
  }
  form.run_dgeom.value='TRUE';
  form.submit();
}
function go_demo(form)
{
  go_init(form)
  form.intxt.value=DEMO_SMI;
  go_dgeom(form)
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
  window.open('"""+util.http.env_cgi.HTML_SUBDIR+"""/jsme_win.html','JSME','width=500,height=450,scrollbars=0,resizable=1,location=0');
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
"""

#############################################################################
def PrintForm():
  ifmtmenu='<SELECT NAME="ifmt">\n'
  for fmt in ('smiles', 'sdf'):
    s = ' SELECTED' if fmt==IFMT else ''
    ifmtmenu+=(f'<OPTION VALUE="{fmt}"{s}>{FMTS[fmt]}')
  ifmtmenu+='</SELECT>'

  RUNMODE_SINGLE=''; RUNMODE_MULTI='';
  if RUNMODE=='single': RUNMODE_SINGLE='CHECKED';
  elif RUNMODE=='multi': RUNMODE_MULTI='CHECKED';

  FF_MMFF=''; FF_UFF='';
  if FF=='mmff': FF_MMFF='CHECKED';
  elif FF=='uff': FF_UFF='CHECKED';

  print(f'<FORM NAME="mainform" ACTION="{CGIURL}" METHOD="POST" ENCTYPE="multipart/form-data">')
  print('<INPUT TYPE=HIDDEN NAME="run_dgeom">')
  print(f'<TABLE WIDTH="100%"><TR><TD><H2>{APPNAME}</H2></TD>')
  print('<TD><TD>- RDKit distance-geometry conformer generation</TD>')
  print('<TD ALIGN="right">')
  print(f"""<BUTTON TYPE=BUTTON onClick="void window.open('{CGIURL}?help=TRUE','helpwin','width=600,height=400,scrollbars=1,resizable=1')"><B>Help</B></BUTTON>""")
  print('<BUTTON TYPE="button" onClick="go_demo(this.form)"><B>Demo</B></BUTTON>')
  print(f"""<BUTTON TYPE=BUTTON onClick="window.location.replace('{CGIURL}')"><B>Reset</B></BUTTON>&nbsp;""")
  print('</TD></TR></TABLE>')
  print('<HR>')
  print('<TABLE WIDTH=100% CELLPADDING=5>')
  print('<TR>')
  print('<TD VALIGN=TOP>')
  print(f'input format:{ifmtmenu}\n')
  print(f'file2txt:<INPUT TYPE=CHECKBOX NAME="infile2txt" VALUE="CHECKED" {INFILE2TXT}><BR>')
  print('upload:<INPUT TYPE="FILE" NAME="infile">&nbsp;')
  print('or paste...<BR>')
  print(f'<TEXTAREA NAME="intxt" ROWS=12 COLS=50 WRAP=OFF>{INTXT}</TEXTAREA><BR>')
  print('or draw:')
  print('<BUTTON TYPE=BUTTON onClick="StartJSME()">JSME</BUTTON>')
  print('</TD>')
  print('<TD VALIGN=TOP>')
  print('<TABLE WIDTH="100%" BGCOLOR="#EEEEEE" CELLPADDING=2>')
  print('<TR><TD COLSPAN=2><B>options:</B></TD></TR>')
  print(f'<TR><TD ALIGN="RIGHT" VALIGN="TOP">mode:</TD><TD><INPUT TYPE=RADIO NAME="runmode" VALUE="single" {RUNMODE_SINGLE}>single <I>(interactive)</I><BR>')
  print(f'<INPUT TYPE=RADIO NAME="runmode" VALUE="multi" {RUNMODE_MULTI}>multi <I>(batch)</I></TD></TR>')
  print('<TR><TD ALIGN=RIGHT>maxconf:</TD>')
  print(f'<TD><INPUT TYPE=TEXT SIZE=3 NAME="maxconf" VALUE="{MAXCONF}"></TD></TR>')
  print(f'<TR><TD ALIGN="RIGHT">forcefield:</TD><TD><INPUT TYPE=RADIO NAME="ff" VALUE="mmff" {FF_MMFF}>MMFF</I>')
  print(f'<INPUT TYPE=RADIO NAME="ff" VALUE="uff" {FF_UFF}>UFF</TD></TR>')
  print('<TR><TD ALIGN=RIGHT>optimizer iters:</TD>')
  print(f'<TD><INPUT TYPE=TEXT SIZE=3 NAME="optiters" VALUE="{OPTITERS}"></TD></TR>')
  print('<TR><TD ALIGN=RIGHT>optimizer E-tolerance:</TD>')
  print(f'<TD><INPUT TYPE=TEXT SIZE=5 NAME="etol" VALUE="{ETOL}"></TD></TR>')
  print('<TR><TD COLSPAN=2><B>misc:</B></TD></TR>')
  print('<TR><TD ALIGN=RIGHT>smifile header:</TD>')
  print(f'<TD ALIGN=LEFT><INPUT TYPE=CHECKBOX NAME="smifile_header" VALUE="CHECKED" {SMIFILE_HEADER}></TD></TR>')
  print('<TR><TD ALIGN=RIGHT>verbose:</TD>')
  print(f'<TD ALIGN=LEFT><INPUT TYPE=CHECKBOX NAME="verbose" VALUE="CHECKED" {VERBOSE}></TD></TR>')
  print('</TABLE>')
  print('</TD></TR>')
  print('<TR><TD COLSPAN=2 ALIGN=CENTER>')
  print(f'<BUTTON TYPE=BUTTON onClick="go_dgeom(this.form)"><B>Go {APPNAME}</B></BUTTON>')
  print('</TD></TR>')
  print('</TABLE>')
  print('</FORM>')

#############################################################################
def DgeomSingle(imol, outfile, maxconf, ff, optiters, etol):
  if not imol:
    OUTPUTS.append("<H2>ERROR: no input molecule</H2>")
    return
  molname = imol.GetProp('_Name') if imol.HasProp('_Name') else ''
  ismi = rdkit.Chem.MolToSmiles(imol, isomericSmiles=True)
  imghtm = (util.http.mol2imghtm.Smi2ImgHtm((ismi+' '+molname if molname else ismi), '', 180, 220, SMI2IMG))
  thtm = ('<TABLE><TR><TD VALIGN=TOP><TT>'+molname+'</TT></TD></TR>'
 	+ "<TR><TD>"+ismi+"</TD></TR>"
 	+ '<TR><TD BGCOLOR="white">'+imghtm+'</TD></TR></TABLE>')
  OUTPUTS.append("<H2>input:</H2>")
  OUTPUTS.append("<BLOCKQUOTE>"+thtm+"</BLOCKQUOTE>")
  OUTPUTS.append('<H2>output:</H2>')

  molwriter = rdkit.Chem.SDWriter(outfile)
  outmol, confIds = util.Utils.GenerateConformations(imol,maxconf,ff,optiters,etol)
  for confId in confIds:
    molwriter.write(outmol,confId=confId)
  if not confIds:
    OUTPUTS.append('<B>Dgeom failed.</B>')
    return

  # furl (file URL) must refer to URL, since JSMOL (JavaScript) runs in client!
  furl = 'https://'+URLHOST+SCRATCHDIRURL+'/'+os.path.basename(outfile)
  #ERRORS.append("DEBUG: furl: "+furl)
  filecode = urllib.parse.quote(furl, safe='')
  bhtm = (f"""<FORM><BUTTON TYPE=BUTTON onClick="go_view3d('{VIEW3D}','{filecode}',{600},'view3d_dgeom','multiconf')">"""
	+ f'''<IMG SRC="'''+util.http.env_cgi.HTML_SUBDIR+f'''/images/Jmol_icon_128.png" HEIGHT="60" VALIGN="middle">JSmol; view output ({len(confIds)} confs)</BUTTON></FORM>''')
  OUTPUTS.append('<BLOCKQUOTE>'+bhtm+'</BLOCKQUOTE>')

  fname = 'dgeom_out.sdf'
  bhtm = (f'<FORM ACTION="{CGIURL}/{fname}" METHOD="POST">'
	+ '<INPUT TYPE=HIDDEN NAME="downloadfile" VALUE="'+outfile+'">'
	+ f"""<BUTTON TYPE=BUTTON onClick="this.form.submit()"> {fname} ({util.http.Utils.NiceBytes(os.stat(outfile).st_size)})</BUTTON></FORM>""")
  OUTPUTS.append('<BLOCKQUOTE><b>Download:</b> '+bhtm+'</BLOCKQUOTE>')

#############################################################################
def DgeomLaunchProcess(infile, smifile_header, outfile, logfile, statusfile, maxconf, ff, optiters, etol, job):
  proc = Process(target=util.mp.Utils.dgeom_process, args=(infile, smifile_header, outfile, logfile, statusfile, maxconf, ff, optiters, etol))
  t0 = time.time()
  proc.start()
  util.mp.Utils.DgeomProgress(proc, statusfile, t0, 1, PROGRESS_WIN_NAME, job)
  ERRORS.append(f"{job} execution time: "+time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0)))
  if VERBOSE:
    ERRORS.append("<PRE>"+open(logfile).read()+"</PRE>")
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
    OUTPUTS.append(f'Warning: output truncated; limit reached NMAX={NMAX}.')

  n_in = status['n_in'] if 'n_in' in status else 0
  n_mol_out = status['n_mol_out'] if 'n_mol_out' in status else 0
  n_conf_out = status['n_conf_out'] if 'n_conf_out' in status else 0
  n_conf_converged = status['n_conf_converged'] if 'n_conf_converged' in status else 0
  n_err = status['n_err'] if 'n_err' in status else 0
  OUTPUTS.append('<H2>output:</H2>')
  OUTPUTS.append('mols in: '+n_in)
  OUTPUTS.append('mols out: '+n_mol_out)
  OUTPUTS.append('confs out: '+n_conf_out)
  OUTPUTS.append(f'confs converged: {n_conf_converged} ; NOT converged: {n_conf_out-n_conf_converged}')
  OUTPUTS.append('errors: '+n_err)
  OUTPUTS.append('<H2>download:</H2>')

  fname = 'dgeom_out.sdf'
  bhtm = (f'<FORM ACTION="{CGIURL}/{fname}" METHOD="POST">'
	+ f'<INPUT TYPE=HIDDEN NAME="downloadfile" VALUE="{OUTFILE}">'
	+ f'<BUTTON TYPE=BUTTON onClick="this.form.submit()"><B>{fname} ({util.http.Utils.NiceBytes(os.stat(OUTFILE).st_size)})</B></BUTTON></FORM>')
  OUTPUTS.append('<BLOCKQUOTE>'+bhtm+' - output conformations (SDF)</BLOCKQUOTE>')

  fname = 'dgeom_log.txt'
  bhtm = (f'<FORM ACTION="{CGIURL}/{fname}" METHOD="POST">'
	+ f'<INPUT TYPE=HIDDEN NAME="downloadfile" VALUE="{logfile}">'
	+ f'<BUTTON TYPE=BUTTON onClick="this.form.submit()"><B>{fname} ({util.http.Utils.NiceBytes(os.stat(logfile).st_size)})</B></BUTTON></FORM>')
  OUTPUTS.append('<BLOCKQUOTE><b>'+bhtm+' - log file</BLOCKQUOTE>')

#############################################################################
def Initialize():
  global FORM,VERBOSE,INTXT,FIXTXT,INFILE2TXT,APPNAME
  global ERRORS,OUTPUTS,CGIURL,PROG,FMTS,FMTS3D,IFMT
  global SMI2IMG,VIEW3D,FIXFMT,SMIFILE_HEADER
  global RUNMODE,URLHOST
  ERRORS=[]; OUTPUTS=[];
  FORM = cgi.FieldStorage(keep_blank_values=1)
  CGIURL = os.environ['SCRIPT_NAME']
  PROG = os.path.basename(sys.argv[0])
  APPNAME = 'DGeom'

  logohtm = "<TABLE CELLSPACING=5 CELLPADDING=5>\n<TR><TD>"
  href = "http://medicine.unm.edu/informatics/"
  imghtm = ('<IMG SRC="'+util.http.env_cgi.HTML_SUBDIR+'/images/biocomp_logo_only.gif">')
  tiphtm = (f'{APPNAME} web app from UNM Translational Informatics Divsion')
  logohtm += (util.http.Utils.HtmTipper(imghtm, '', tiphtm, href, 200, "white"))

  href = "http://www.rdkit.org"
  imghtm = ('<IMG HEIGHT="70" SRC="'+util.http.env_cgi.HTML_SUBDIR+'/images/rdkit_logo_only.png">')
  tiphtm = ('RDKit: Cheminformatics and Machine Learning Software')
  logohtm += (util.http.Utils.HtmTipper(imghtm, '', tiphtm, href, 200, "white"))

  href = "http://www.jmol.org"
  imghtm = ('<IMG HEIGHT="60" SRC="'+util.http.env_cgi.HTML_SUBDIR+'/images/Jmol_icon_128.png">')
  tiphtm = ('JSmol 3D viewer from the Jmol project')
  logohtm += (util.http.Utils.HtmTipper(imghtm, '', tiphtm, href, 200, "white"))

  logohtm += "</TD></TR></TABLE>"
  ERRORS.append("<CENTER>"+logohtm+"</CENTER>")

  VERBOSE = FORM.getvalue('verbose')
  INTXT = FORM.getvalue('intxt','')
  FIXTXT = FORM.getvalue('fixtxt','')
  INFILE2TXT = FORM.getvalue('infile2txt')
  URLHOST = os.environ['HTTP_HOST'] if 'HTTP_HOST' in os.environ else os.environ['SERVER_NAME']
  urldir = os.path.dirname(os.environ['REQUEST_URI'])
  SMI2IMG = urldir+'/mol2img.cgi'
  VIEW3D = urldir+'/view3d_jsmol.cgi'
  SMIFILE_HEADER = FORM.getvalue('smifile_header')

  global NMAX,OUTPUT_TRUNCATED
  NMAX=2000
  OUTPUT_TRUNCATED = False
  RUNMODE = FORM.getvalue('runmode','single')

  global FF,MAXCONF,OPTITERS,ETOL
  FF = FORM.getvalue('ff')
  MAXCONF = int(FORM.getvalue('maxconf','1'))
  OPTITERS = int(FORM.getvalue('optiters','400'))
  ETOL = float(FORM.getvalue('etol','1e-6'))

  IFMT = FORM.getvalue('ifmt','smiles')
  FMTS={
	'smi':'SMILES',
	'smiles':'SMILES',
	'mdl':'MDL Molfile',
	'mol':'MDL Molfile',
	'sdf':'MDL SDfile'
	}

  global PROGRESS_WIN_NAME
  PROGRESS_WIN_NAME = 'dgeom_progress_win'

  ERRORS.append('RDKit version: '+rdkit.rdBase.rdkitVersion)
  ERRORS.append('Boost version: '+rdkit.rdBase.boostVersion)

  global SCRATCHDIR,SCRATCHDIRURL,SCRATCHDIR_LIFETIME
  SCRATCHDIR = util.http.env_cgi.scratchdir
  SCRATCHDIRURL = util.http.env_cgi.scratchdirurl
  SCRATCHDIR_LIFETIME = util.http.env_cgi.scratchdir_lifetime

  global TS,PREFIX
  TS = time.strftime('%Y%m%d%H%M%S',time.localtime())
  PREFIX = SCRATCHDIR+'/'+re.sub('\..*$','',PROG)+'.'+TS+'.'+str(random.randint(100,999))

  ### Assure that writeable tmp directory exists.
  if not os.access(SCRATCHDIR,os.F_OK):
    ERRORS.append("ERROR: missing SCRATCHDIR: "+SCRATCHDIR)
    return False
  elif not os.access(SCRATCHDIR,os.W_OK):
    ERRORS.append("ERROR: non-writable SCRATCHDIR: "+SCRATCHDIR)
    return False
  #ERRORS.append("DEBUG: SCRATCHDIR: "+SCRATCHDIR)
  #ERRORS.append("DEBUG: SCRATCHDIRURL: "+SCRATCHDIRURL)

  if not FORM.getvalue('run_dgeom'):
    return True

  ### NOW STUFF FOR A RUN ############################################

  global INFILE,OUTFILE,LOGFILE_DGEOM,STATUSFILE_DGEOM
  os.chdir(SCRATCHDIR)
  tempfile.tempdir = os.getcwd()
  fd,INFILE = tempfile.mkstemp('.'+IFMT, PREFIX); os.close(fd)
  fd,OUTFILE = tempfile.mkstemp('.sdf', PREFIX); os.close(fd)
  fd,LOGFILE_DGEOM = tempfile.mkstemp('.log', PREFIX); os.close(fd)
  fd,STATUSFILE_DGEOM = tempfile.mkstemp('.status', PREFIX); os.close(fd)

  global DELFILES
  DELFILES = [INFILE, LOGFILE_DGEOM, STATUSFILE_DGEOM]

  intxt=''
  if FORM['infile'].file and FORM['infile'].filename:
    intxt = FORM['infile'].value.decode("utf-8") #bytes
    if FORM.getvalue('infile2txt'): INTXT=intxt
  elif FORM['intxt'].value:
    intxt = FORM['intxt'].value
  else:
    ERRORS.append('ERROR: No input data.')
    return False
  intxt = re.sub('\t', ' ', intxt) ##Fix TABs
  intxt = re.sub('\r\n*', '\n', intxt) ##Fix CRs
  with open(INFILE,'w+') as f:
    f.write(intxt)

  if IFMT in ('smi','smiles'):
    molreader = rdkit.Chem.SmilesMolSupplier(INFILE, delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=bool(SMIFILE_HEADER), sanitize=True)
  else:
    molreader = rdkit.Chem.SDMolSupplier(INFILE, sanitize=True, removeHs=True)

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
  print(f"""\
<HR>
<H3>{APPNAME} help</H3>
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
NMAX = {NMAX}<BR>
Depictions by {util.http.env_cgi.DEPICT_TOOL}.<BR>
""")

#############################################################################
if __name__=='__main__':
  ok=Initialize()
  if not ok:
    util.http.Utils.PrintHeader(APPNAME, JavaScript(PROG))
    util.http.Utils.PrintFooter(ERRORS)
  elif FORM.getvalue('run_dgeom'):
    util.http.Utils.PrintHeader(APPNAME, JavaScript(PROG))
    PrintForm()
    if RUNMODE=='single':
      DgeomSingle(IMOL, OUTFILE, MAXCONF, FF, OPTITERS, ETOL)
    else:
      DgeomLaunchProcess(INFILE, bool(SMIFILE_HEADER), OUTFILE, LOGFILE_DGEOM, STATUSFILE_DGEOM, MAXCONF, FF, OPTITERS, ETOL, 'Dgeom')
      status = ParseStatusFile(STATUSFILE_DGEOM)
      DgeomResultsMulti(status, LOGFILE_DGEOM)
      print('<SCRIPT>pwin.parent.focus(); pwin.focus(); pwin.close();</SCRIPT>')
    util.http.Utils.PrintOutput(OUTPUTS)
    util.http.Utils.PrintFooter(ERRORS)
    util.http.Utils.Cleanup(DELFILES)
    util.http.purgescratchdirs.PurgeScratchDirs([SCRATCHDIR], SCRATCHDIR_LIFETIME, 0)
  elif FORM.getvalue('downloadtxt'):
    util.http.Utils.DownloadString(FORM.getvalue('downloadtxt'))
  elif FORM.getvalue('downloadfile'):
    util.http.Utils.DownloadFile(FORM.getvalue('downloadfile'))
  elif FORM.getvalue('help'):
    util.http.Utils.PrintHeader(APPNAME, JavaScript(PROG))
    Help()
    util.http.Utils.PrintFooter(ERRORS)
  elif FORM.getvalue('test'):
    util.http.Utils.PrintTestText(APPNAME, {})
  else:
    util.http.Utils.PrintHeader(APPNAME, JavaScript(PROG))
    PrintForm()
    print('<SCRIPT>go_init(window.document.mainform);</SCRIPT>')
    util.http.Utils.PrintFooter(ERRORS)
