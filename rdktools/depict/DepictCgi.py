#! /usr/bin/env python3
#############################################################################
import os,sys,cgi,platform
import re,time,urllib,tempfile,base64,zlib,random
#import cgitb; cgitb.enable()   ### debugger

import rdkit.rdBase
import rdkit.Chem.AllChem
import rdkit.Chem.Draw

from .. import util

#############################################################################
def JavaScript():
  return '''\
function go_init(form)
{
  var i;
  for (i=0;i<form.ifmt.length;++i)
    if (form.ifmt.options[i].value=='auto')
        form.ifmt.options[i].selected=true;
  form.file2txt.checked=true;
  form.verbose.checked=false;
  form.smarts.value='';
  form.intxt.value='';
  form.hilight_smarts.checked=true;
  form.align_smarts.checked=false;
  form.showarom.checked=false;
  form.showh.checked=false;
  form.gen2d.checked=false;
  for (i=0;i<form.ncols.length;++i)
    if (form.ncols.options[i].value=='auto')
      form.ncols.options[i].selected=true;
  for (i=0;i<form.imgsize.length;++i)
    if (form.imgsize.options[i].value=='m')
      form.imgsize.options[i].selected=true;
  //form.transparent.checked=true;
  form.zoomable.checked=true;
}
function checkform(form)
{
  if (!form.intxt.value && !form.infile.value) {
    alert('ERROR: No input specified');
    return 0;
  }
  return true;
}
function go_depict(form)
{
  if (!checkform(form)) return;
  form.depict.value='TRUE'
  form.submit()
}
'''

#############################################################################
def PrintForm():
  ifmt_menu='<SELECT NAME="ifmt">\n'
  ifmt_menu+='<OPTION VALUE="auto">Automatic'
  ifmt_menu+='<OPTION VALUE="SMI">SMI'
  ifmt_menu+='<OPTION VALUE="MDL">MDL|SDF|MOL'
  ifmt_menu+='</SELECT>'
  ifmt_menu=re.sub('"%s">'%IFMT,'"%s" SELECTED>'%IFMT,ifmt_menu)

  imgsize_menu='<SELECT NAME="imgsize">'
  for size in ('xs','s','m','l','xl'):
    imgsize_menu+=('<OPTION VALUE="%s">imgsize:%s - '%(size,size.upper()))
    imgsize_menu+=('%dx%d'%(IMGSIZES[size][0],IMGSIZES[size][1]))
  imgsize_menu+='</SELECT>'
  imgsize_menu=re.sub('"%s">'%IMGSIZE,'"%s" SELECTED>'%IMGSIZE,imgsize_menu)

  ncols_menu='<SELECT NAME="ncols">'
  ncols_menu+='<OPTION VALUE="auto">columns: auto'
  for i in range(1,11):
    ncols_menu+=('<OPTION VALUE="%d">columns: %d'%(i,i))
  ncols_menu+='</SELECT>'
  ncols_menu=re.sub('"%s">'%NCOLS,'"%s" SELECTED>'%NCOLS,ncols_menu)

  htm = (('<FORM NAME="mainform" ACTION="%s" METHOD="POST" ENCTYPE="multipart/form-data">\n'%CGIURL)
    +('<INPUT TYPE="HIDDEN" NAME="depict">\n')
    +('<TABLE WIDTH="100%" CELLPADDING="0" CELLSPACING="0">\n')
    +('<TR><TD><H1>%s</H1></TD><TD>%s</TD>'%(APPNAME,'- depictions from RDKit\n'))
    +('<TD ALIGN=RIGHT>')
    +('<BUTTON TYPE=BUTTON onClick="window.location.replace(\'%s\')"><B>Reset</B></BUTTON>'%CGIURL )
    +('<BUTTON TYPE="BUTTON" onClick="void window.open(\'%s?help=TRUE\',\'helpwin\',\'width=600,height=400,scrollbars=1,resizable=1\')"><B>Help</B></BUTTON>'%CGIURL)
    +('</TD></TR></TABLE>\n')
    +('<TABLE WIDTH="100%" CELLPADDING="3" CELLSPACING="3">\n')
    +('<TR BGCOLOR="#CCCCCC">\n')
    +('<TD WIDTH="50%" VALIGN="top">')
    +('input format:%s<BR>\n'%ifmt_menu)
    +('upload:<INPUT TYPE="FILE" NAME="infile">&nbsp;')
    +('<INPUT TYPE="CHECKBOX" NAME="file2txt" VALUE="CHECKED" %s>file2txt<BR>'%FILE2TXT)
    +('or paste...<BR>')
    +('<TEXTAREA NAME="intxt" WRAP="OFF" ROWS="12" COLS="50">%s</TEXTAREA>'%INTXT)
    +('</TD>')
    +('<TD VALIGN="top">')
    +('<B>options:</B><BR/>')
    +('<TABLE WIDTH="100%" CELLPADDING="0" CELLSPACING="0"><TR><TD VALIGN="top">\n')
    +('<INPUT TYPE="CHECKBOX" NAME="showarom" VALUE="CHECKED" %s>showarom<BR>'%SHOWAROM)
    +('<INPUT TYPE="CHECKBOX" NAME="showh" VALUE="CHECKED" %s>show Hs<BR>'%SHOWH)
    +('</TD><TD VALIGN="top">')
    +('<INPUT TYPE=CHECKBOX NAME="gen2d" VALUE="CHECKED" %s>force gen2D<BR>'%GEN2D)
    +('</TD></TR></TABLE>\n')
    +('<HR>\n')
    +('smarts:<INPUT TYPE="TEXT" NAME="smarts" SIZE="25" VALUE="%s">'%SMARTS)
    +('<INPUT TYPE="CHECKBOX" NAME="hilight_smarts" VALUE="CHECKED" %s>highlight'%HILIGHT_SMARTS)
    +('<INPUT TYPE="CHECKBOX" NAME="align_smarts" VALUE="CHECKED" %s>align<BR>'%ALIGN_SMARTS)
    +('<HR>')
    +(ncols_menu)
    +(imgsize_menu)
    +('<HR>')
    ##  +('<INPUT TYPE="CHECKBOX" NAME="transparent" VALUE="CHECKED" %s>transparent<BR>'%TRANSPARENT)
    +('<INPUT TYPE="CHECKBOX" NAME="zoomable" VALUE="CHECKED" %s>zoomable<BR>'%('CHECKED' if ZOOMABLE else ''))
    +('<INPUT TYPE="CHECKBOX" NAME="verbose" VALUE="CHECKED" %s>verbose<BR>'%VERBOSE)
    +('</TD></TR>\n')
    +('<TR BGCOLOR="#CCCCCC">')
    +('<TD COLSPAN="2" ALIGN="CENTER">')
    +('<BUTTON TYPE="BUTTON" onClick="go_depict(this.form)"><B>Go %s</B></BUTTON>'%APPNAME)
    +('</TD></TR>\n')
    +('</TABLE>\n')
    +('</FORM>\n'))
  print(htm)

#############################################################################
def Depict():
  width,height=IMGSIZES[IMGSIZE]
  opts='mode=%s'%(MODE)
  opts += '&aromdash=true' if SHOWAROM else '&kekule=true'
  ##if TRANSPARENT: opts+='&transparent=true'
  if SHOWH: opts+='&showh=true'
  if SMARTS and HILIGHT_SMARTS:
    opts+='&smartscode=%s'%urllib.quote(SMARTS)
  if GEN2D: opts+='&gen2d=true'

  thtm='<CENTER><TABLE BORDER>\n<TR>\n'

  if NCOLS=='auto':
    ncols=int(900.0/(IMGSIZES[IMGSIZE][0]+20))
  else:
    ncols=int(NCOLS)
  n_mols=0

  for mol in MOLSUP:
    if not mol: continue
    n_mols+=1
    if not (n_mols-1)%ncols: thtm+='<TR>\n'
    molname=mol.GetProp('_Name')
    opts_this=opts
    if molname:
      opts_this+='&title=%s'%urllib.quote(molname)

    if IFMT in ('MDL','SDF','MOL'):
      mdl = rdkit.Chem.MolToMolBlock(mol)
      imghtm= util.http.mol2imghtm.Mol2ImgHtm(mdl,opts_this,height,width,MOL2IMG,ZOOMABLE,'go_zoom_mol2img',4,'','<CENTER>%s</CENTER>'%molname)
    else:
      smi = rdkit.Chem.MolToSmiles(mol,isomericSmiles=True)
      imghtm= util.http.mol2imghtm.Smi2ImgHtm(smi,opts_this,height,width,MOL2IMG,ZOOMABLE,'go_zoom_smi2img',4,'','<CENTER>%s</CENTER>'%molname)

    thtm+=('<TD ALIGN="center" BGCOLOR="white">%s<BR>%s</TD>\n'%(imghtm,molname))
    if not n_mols%ncols: thtm+='</TR>\n'
    if n_mols==NMAX:
      ERRORS.append('NOTE: max %d mols reached.'%NMAX)
      break

  if n_mols%ncols:
    if n_mols>ncols: thtm+=((ncols-n_mols%ncols)*'<TD ALIGN=CENTER>~</TD>\n')
    thtm+='</TR>\n'
  thtm+='</TABLE></CENTER>\n'
  OUTPUTS.append(thtm)
  ERRORS.append('RESULT: %d molecules displayed'%n_mols)
  LOG.write('%s\t%s\t%d\n' % (DATE,os.environ['REMOTE_ADDR'],n_mols))

#############################################################################
def Initialize():
  global FORM,VERBOSE,NMAX,ERRORS,OUTPUTS,CGIURL,PROG,LOG,DATE,APPNAME
  global TRANSPARENT,SHOWH,FILE2TXT,INTXT
  global GEN2D  #a.k.a. ignore2d
  global IMGSIZE,IMGSIZES,NCOLS,MODE,SMARTS,ZOOMABLE
  global SHOWAROM,KEKULE,NOTITLE
  ERRORS=[]; OUTPUTS=[];
  FORM=cgi.FieldStorage(keep_blank_values=1)
  global USER_AGENT
  #host=os.environ['HTTP_HOST'] 
  USER_AGENT=os.environ['HTTP_USER_AGENT']
  CGIURL=os.environ['SCRIPT_NAME']
  PROG=os.path.basename(os.environ['SCRIPT_FILENAME'])
  APPNAME='Depict'

  logohtm="<TABLE CELLSPACING=5 CELLPADDING=5>\n<TR><TD>"

  href="http://medicine.unm.edu/informatics/"
  imghtm=('<IMG SRC="/images/biocomp_logo_only.gif">')
  tiphtm=('%s web app from UNM Translational Informatics Divsion'%APPNAME)
  logohtm+=(util.http.Utils.HtmTipper(imghtm,'',tiphtm,href,200,"white"))

  href="http://www.rdkit.org"
  imghtm=('<IMG HEIGHT="60" BORDER="0" SRC="/images/rdkit_logo.png">')
  tiphtm=('RDKit')
  logohtm+=(util.http.Utils.HtmTipper(imghtm,'',tiphtm,href,200,"white"))

  logohtm+="</TD></TR></TABLE>"
  ERRORS.append(logohtm)

  DATE=time.strftime('%Y%m%d%H%M',time.localtime())
  VERBOSE=FORM.getvalue('verbose')

  ### These CGIs must be in same dir.
  global MOL2IMG
  urldir=os.path.dirname(re.sub(r'\?.*$','',os.environ['REQUEST_URI']))
  MOL2IMG=urldir+'/mol2img_rdk.cgi'

  global IFMT
  IFMT=FORM.getvalue('ifmt','SMI')

  KEKULE=FORM.getvalue('kekule')
  SHOWAROM=FORM.getvalue('showarom')
  ##TRANSPARENT=FORM.getvalue('transparent')
  SHOWH=FORM.getvalue('showh')
  GEN2D = 'CHECKED' if FORM.getvalue('gen2d') else ''
  FILE2TXT=FORM.getvalue('file2txt')
  INTXT=FORM.getvalue('intxt')
  IMGSIZE=FORM.getvalue('imgsize','m')
  MODE=FORM.getvalue('mode','cow')
  IMGSIZES={'xs':[96,96],'s':[160,160],'m':[260,180],'l':[380,280],'xl':[640,480]}
  NCOLS=FORM.getvalue('ncols','auto')
  global ALIGN_SMARTS,HILIGHT_SMARTS
  SMARTS=FORM.getvalue('smarts')
  ALIGN_SMARTS=FORM.getvalue('align_smarts')
  HILIGHT_SMARTS=FORM.getvalue('hilight_smarts')
  ZOOMABLE=FORM.getvalue('zoomable')
  NMAX=100

  LOGFILE='./logs/'+PROG+'.log'
  LOGFIELDS=['ip','nmols']
  if os.access(LOGFILE,os.R_OK):
    tmp=open(LOGFILE)
    lines=tmp.readlines()
    tmp.close()
    n=0
    for i,line in enumerate(lines):
      fields=line.split()
      if i>0 and len(fields)>2:
        n+=int(fields[2])
    if i>0:
      date=lines[1].split('\t')[0][:8]
      date='%s/%s/%s'%(date[4:6],date[6:8],date[:4])
      ERRORS.append('%s: since %s: used %d times; %d mols.'%(PROG,date,i,n))
  if os.access(LOGFILE,os.W_OK):
    LOG=open(LOGFILE,'a')
    if os.path.getsize(LOGFILE)==0:
      LOG.write('date\t' + '\t'.join(LOGFIELDS) + '\n')
  else:
    LOG=open(LOGFILE,'w+')
    if not LOG:
      ERRORS.append('ERROR: LOGFILE cannot be opened.')
      return False
    LOG.write('date\t' + '\t'.join(LOGFIELDS) + '\n')

  if VERBOSE:
    ERRORS.append('RDKit version: %s'%rdkit.rdBase.rdkitVersion)
    ERRORS.append('Boost version: %s'%rdkit.rdBase.boostVersion)
    ERRORS.append('Python version: %s'%platform.python_version())

  if not FORM.getvalue('depict'): return True

  ### STUFF FOR A RUN

  fdata=None; fext=None;
  if FORM.has_key('infile') and FORM['infile'].file and FORM['infile'].filename:
    fdata=FORM['infile'].value
    if FORM.getvalue('file2txt'):
      INTXT=fdata
    fext=os.path.splitext(FORM['infile'].filename)[1]
    #ERRORS.append('DEBUG: filename: "%s", ext="%s"'%(FORM['infile'].filename,fext))
  elif FORM.has_key('intxt') and FORM['intxt'].value:
    INTXT=FORM['intxt'].value
    fdata=INTXT
  else:
    ERRORS.append('ERROR: No input data.')
    return False

  global MOLSUP #Molecule supplier
  IFMT=SelectFormat(IFMT,fext)
  try:
    MOLSUP=SelectMolsupplier(IFMT)
    ### Fails to read 1st SMILES without "title line".
    MOLSUP.SetData('smiles name\n'+fdata)
  except:
    ERRORS.append('Format error: ifmt="%s", filename="%s", ext="%s"'%(IFMT,FORM['infile'].filename,fext))
    return False

  global SMARTSPAT
  if SMARTS:
    SMARTSPAT=rdkit.Chem.MolFromSmarts(SMARTS)
    if not SMARTSPAT:
      ERRORS.append('ERROR: Bad SMARTS: %s'%SMARTS)
      return False

  global TS,PREFIX
  TS=time.strftime('%Y%m%d%H%M%S',time.localtime())
  PREFIX=re.sub('\..*$','',PROG)+'.'+TS+'.'+str(random.randint(100,999))

  return True

#############################################################################
def SelectFormat(ifmt,fext):
  if ifmt == 'auto' and fext:
    if fext.upper() in ('.MDL','.SDF','.MOL'):
      ifmt='MDL'
    elif fext.upper() in ('.SMI','.SMILES','.ISM'):
      ifmt='SMI'
  return ifmt

#############################################################################
def SelectMolsupplier(ifmt):
  molsup=None
  if ifmt in ('MDL','SDF','MOL'):
    molsup=rdkit.Chem.AllChem.SDMolSupplier()
  elif ifmt == 'SMI':
    molsup=rdkit.Chem.AllChem.SmilesMolSupplier()
  else:
    ERRORS.append('Unknown format: %s'%(ifmt))
  return molsup

#############################################################################
def Help():
  htm='''\
<HR>
<H3>%(APPNAME)s help</H3>
</P><P>
This web app built with RDKit (Python API) and is fully open-source.
</P><P>
NMAX = %(NMAX)d
</P><P>
author/support: Jeremy Yang
</P>
'''%{'APPNAME':APPNAME,'NMAX':NMAX}
  print(htm)

#############################################################################
if __name__=='__main__':
  ok=Initialize()
  if not ok:
    util.http.Utils.PrintHeader(PROG,JavaScript())
    util.http.Utils.PrintFooter(ERRORS)
  elif FORM.getvalue('depict'):
    util.http.Utils.PrintHeader(PROG,JavaScript())
    Depict()
    PrintForm()
    util.http.Utils.PrintOutput(OUTPUTS)
    util.http.Utils.PrintFooter(ERRORS)
  elif FORM.getvalue('downloadtxt'):
    util.http.Utils.DownloadString(FORM.getvalue('downloadtxt'))
  elif FORM.has_key('help'):
    util.http.Utils.PrintHeader(PROG,JavaScript())
    Help()
    util.http.Utils.PrintFooter(ERRORS)
  else:
    util.http.Utils.PrintHeader(PROG,JavaScript())
    PrintForm()
    print('<SCRIPT>go_init(window.document.mainform)</SCRIPT>')
    util.http.Utils.PrintFooter(ERRORS)
