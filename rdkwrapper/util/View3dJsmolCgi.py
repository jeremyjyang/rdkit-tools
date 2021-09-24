#! /usr/bin/env python3
###
"""
Display JSMol w/ mol[s] from file[s].
http://wiki.jmol.org/index.php/Main_Page
http://wiki.jmol.org/index.php/JSmol
"""
###
import os,sys,cgi,urllib.parse,re

from .. import util

#############################################################################
def JavaScript():
  return """\
//These *_checkbox_* functions update JSmol based on changed checkbox input.
function update_checkbox_labels(form)
{
  if (form.labels.checked)
    Jmol.script(JMOL1,'select *;label %a;select *;');
    //Jmol.script(JMOL1,'select *;labels on');
  else
    Jmol.script(JMOL1,'select *;labels off');
}
function update_checkbox_antialias(form)
{
  if (form.antialias.checked)
    Jmol.script(JMOL1,'set antialiasDisplay true');
  else
    Jmol.script(JMOL1,'set antialiasDisplay false');
}
function update_checkbox_hydrogens(form)
{
  if (form.hydrogens.checked)
    Jmol.script(JMOL1,'set showHydrogens on');
  else
    Jmol.script(JMOL1,'set showHydrogens off');
}
function update_checkbox_dots(form)
{
  if (form.dots.checked)
    Jmol.script(JMOL1,'dots on');
  else
    Jmol.script(JMOL1,'dots off');
}
function update_checkbox_spin(form)
{
  if (form.spin.checked)
    Jmol.script(JMOL1,'spin on');
  else
    Jmol.script(JMOL1,'spin off');
}
//These *_radio_* functions update JSmol based on changed radio input.
function update_radio_color(form)
{
  var i;
  for (i=0;i<form.color.length;++i)     //radio
  {
    if (form.color[i].checked)
    {
      if (form.color[i].value=='cpk')
        Jmol.script(JMOL1,'color cpk');
      else if (form.color[i].value=='atomno')
        Jmol.script(JMOL1,'color property atomno');
      else if (form.color[i].value=='structure')
        Jmol.script(JMOL1,'color structure');
      break;
    }
  }
}
function update_radio_surface(form)
{
  var i;
  for (i=0;i<form.surface.length;++i)     //radio
  {
    if (form.surface[i].checked)
    {
      if (form.surface[i].value=='none')
        Jmol.script(JMOL1,'isosurface off');
      else if (form.surface[i].value=='vdw')
        Jmol.script(JMOL1,'select *;isosurface vdw');
      else if (form.surface[i].value=='mep')
        Jmol.script(JMOL1,'if ({atomno<10}.partialcharge==0){calculate partialcharge};isosurface vdw map mep');
      break;
    }
  }
}
function update_radio_psurface(form)
{
  var i;
  for (i=0;i<form.psurface.length;++i)     //radio
  {
    if (form.psurface[i].checked)
    {
      if (form.psurface[i].value=='none')
        Jmol.script(JMOL1,'isosurface ID psurf off;');
      else if (form.psurface[i].value=='vdw')
        Jmol.script(JMOL1,'select 1.1; isosurface ID psurf RESOLUTION 0 SOLVENT 0.0 COLOR yellow;');
      break;
    }
  }
}
function update_checkbox_translucent(form)
{
  if (form.translucent.checked)
    Jmol.script(JMOL1,'isosurface translucent');
  else
    Jmol.script(JMOL1,'isosurface opaque');
}
function update_radio_style(form)
{
  var i;
  for (i=0;i<form.style.length;++i)     //radio
  {
    if (form.style[i].checked)
    {
      if (form.style[i].value=='ballstick')
        Jmol.script(JMOL1,'select *;cartoons off;spacefill 23%;wireframe 0.15');
      else if (form.style[i].value=='wireframe')
        Jmol.script(JMOL1,'select *;cartoons off;wireframe -0.1');
      else if (form.style[i].value=='spacefill')
        Jmol.script(JMOL1,'select *;cartoons off;spacefill only');
      break;
    }
  }
}
function update_radio_cartoons(form)
{
  var i;
  for (i=0;i<form.cartoons.length;++i)     //radio
  {
    if (form.cartoons[i].checked)
    {
      if (form.cartoons[i].value=='none')
        Jmol.script(JMOL1,'select *;cartoons off;wireframe -0.1');
      else if (form.cartoons[i].value=='on')
        Jmol.script(JMOL1,'select protein or nucleic;cartoons only; set cartoonFancy false');
      else if (form.cartoons[i].value=='fancy')
        Jmol.script(JMOL1,'select protein or nucleic;cartoons only; set cartoonFancy true');
      break;
    }
  }
}
function update_checkbox_refmolgreen(form)
{
  if (form.refmolgreen.checked)
    Jmol.script(JMOL1,'select 1.1; color atoms green;');
  else
    Jmol.script(JMOL1,'select 1.1; color atoms CPK;');
}
function update_radio_molsvisible(form)
{
  var i;
  for (i=0;i<form.molsvisible.length;++i)     //radio
  {
    if (form.molsvisible[i].checked)
    {
      if (form.molsvisible[i].value=='refmol')
        Jmol.script(JMOL1,'model 1.1;');
      else if (form.molsvisible[i].value=='dbmol')
        Jmol.script(JMOL1,'model 1.2;');
      else if (form.molsvisible[i].value=='both')
        Jmol.script(JMOL1,'model 1.0;');
      break;
    }
  }
}
function update_radio_protein(form)
{
  var i;
  for (i=0;i<form.protein.length;++i)     //radio
  {
    if (form.protein[i].checked)
    {
      if (form.protein[i].value=='wireframe')
        Jmol.script(JMOL1,'select 1.1; wireframe on; spacefill off; ribbons off; cartoon off;');
      else if (form.protein[i].value=='ballstick')
        Jmol.script(JMOL1,'select 1.1; wireframe off; wireframe 0.15; spacefill 0.4; ribbons off; cartoon off;');
      else if (form.protein[i].value=='cartoon')
        Jmol.script(JMOL1,'select 1.1; wireframe off; spacefill off; ribbons off; cartoon on; color cartoon structure;');
      break;
    }
  }
}
function update_radio_ligand(form)
{
  var i;
  for (i=0;i<form.ligand.length;++i)     //radio
  {
    if (form.ligand[i].checked)
    {
      if (form.ligand[i].value=='ballstick')
        Jmol.script(JMOL1,'select 2.1;cartoons off;spacefill 23%;wireframe 0.15');
      else if (form.ligand[i].value=='wireframe')
        Jmol.script(JMOL1,'select 2.1;cartoons off;wireframe -0.1');
      else if (form.ligand[i].value=='spacefill')
        Jmol.script(JMOL1,'select 2.1;cartoons off;spacefill only');
      break;
    }
  }
}
function update_checkbox_animation(form)
{
  if (form.animation.checked)
    Jmol.script(JMOL1,'animation fps 2; animation mode loop; animation on;');
  else
    Jmol.script(JMOL1,'animation off;');
}
"""

#############################################################################
### https://chemapps.stolaf.edu/jmol/docs/#load
### "Loads the specified file or URL."
### MUST USE URL SINCE JS RUNS ON CLIENT!
### http://wiki.jmol.org/index.php/Jmol_JavaScript_Object/Functions#loadFile
### Jmol.loadFile = function(JmolObject, fileName, params)
### Use instead Jmol.script(myJmol, "load '" + fileName + "'....")
#############################################################################
def ViewMolecule(furl, molname=''):
  js_init = (f"""\
Jmol.loadFile(JMOL1, '{furl}');
var jmolscript='';
jmolscript+='set showHydrogens off; set measurements angstroms;';
jmolscript+='color echo yellow; ';
jmolscript+='set echo top left; ';
// jmolscript+='echo "{molname}";');
Jmol.script(JMOL1, jmolscript);
""")
  PrintHeader(title=TITLE, js=JavaScript(), initjs=js_init, css='')
  print("""
<DIV ID="JMOLDIV1" style="width:100%%;height:95%%;background-color:white"></DIV>
<FORM>
Hs:<INPUT TYPE="CHECKBOX" NAME="hydrogens" onChange="update_checkbox_hydrogens(this.form)">
dots:<INPUT TYPE="CHECKBOX" NAME="dots" onChange="update_checkbox_dots(this.form)">
labels:<INPUT TYPE="CHECKBOX" NAME="labels" onChange="update_checkbox_labels(this.form)">
spin:<INPUT TYPE="CHECKBOX" NAME="spin" onChange="update_checkbox_spin(this.form)">
&nbsp;&nbsp;&nbsp;
<b>style:</b>
<INPUT TYPE="RADIO" NAME="style" VALUE="ballstick" onChange="update_radio_style(this.form)" CHECKED>ball&amp;stick
<INPUT TYPE="RADIO" NAME="style" VALUE="wireframe" onChange="update_radio_style(this.form)">wireframe
<INPUT TYPE="RADIO" NAME="style" VALUE="spacefill" onChange="update_radio_style(this.form)">spacefill
</FORM>
</DIV>
""")

#############################################################################
def ViewMulticonformer(furl):
  js_extra = ("""
function JmolCBFunc(app, frameno, fileno, modelno, firstno, lastno) {{
  var n_conf=lastno-firstno+1;
  Jmol.script(JMOL1, 'echo "conformer: '+modelno+' / '+n_conf+'";');
}}
""")
  js_init = ("""\
Jmol.loadFile(JMOL1, '"""+furl+"""');
var jmolscript='';
jmolscript+='set showHydrogens off; set measurements angstroms;';
jmolscript+='color echo yellow; ';
jmolscript+='set echo top left; ';
jmolscript+='model 1; ';
//jmolscript+='echo "conformer: 1";');
jmolscript+='set AnimFrameCallback "JmolCBFunc"; ';
Jmol.script(JMOL1, jmolscript);
""")
  PrintHeader(title=TITLE, js=JavaScript()+js_extra, initjs=js_init, css='')
  print("""\
<DIV ID="JMOLDIV1" style="width:100%%;height:85%%;background-color:white"></DIV>
<br />
<FORM>
<BUTTON TYPE="BUTTON" onClick="Jmol.script(JMOL1, 'model NEXT')">next</BUTTON>
<BUTTON TYPE="BUTTON" onClick="Jmol.script(JMOL1, 'model PREVIOUS')">previous</BUTTON>
<BUTTON TYPE="BUTTON" onClick="Jmol.script(JMOL1, 'model 1')">1st</BUTTON>
<BUTTON TYPE="BUTTON" onClick="Jmol.script(JMOL1, 'model all')">all</BUTTON>

Hs:<INPUT TYPE="CHECKBOX" NAME="hydrogens" onChange="update_checkbox_hydrogens(this.form)">
dots:<INPUT TYPE="CHECKBOX" NAME="dots" onChange="update_checkbox_dots(this.form)">
labels:<INPUT TYPE="CHECKBOX" NAME="labels" onChange="update_checkbox_labels(this.form)">
spin:<INPUT TYPE="CHECKBOX" NAME="spin" onChange="update_checkbox_spin(this.form)">
showinfo:<INPUT TYPE="CHECKBOX" NAME="showinfo" onChange="Jmol.showInfo(JMOL1, this.form.showinfo.checked)">
slideshow:<INPUT TYPE="CHECKBOX" NAME="animation" onChange="update_checkbox_animation(this.form)">
<br />
<b>style:</b>
<INPUT TYPE="RADIO" NAME="style" VALUE="ballstick" onChange="update_radio_style(this.form)" CHECKED>ball&amp;stick
<INPUT TYPE="RADIO" NAME="style" VALUE="wireframe" onChange="update_radio_style(this.form)">wireframe
<INPUT TYPE="RADIO" NAME="style" VALUE="spacefill" onChange="update_radio_style(this.form)">spacefill
&nbsp;&nbsp;&nbsp;
<b>color:</b>
<INPUT TYPE="RADIO" NAME="color" VALUE="cpk" onChange="update_radio_color(this.form)" CHECKED>cpk
<INPUT TYPE="RADIO" NAME="color" VALUE="atomno" onChange="update_radio_color(this.form)">atomno
<INPUT TYPE="RADIO" NAME="color" VALUE="structure" onChange="update_radio_color(this.form)">structure
&nbsp;&nbsp;&nbsp;
<b>surface:</b>
<INPUT TYPE="RADIO" NAME="surface" VALUE="none" onChange="update_radio_surface(this.form)" CHECKED>none
<INPUT TYPE="RADIO" NAME="surface" VALUE="vdw" onChange="update_radio_surface(this.form)">vdw
<INPUT TYPE="RADIO" NAME="surface" VALUE="mep" onChange="update_radio_surface(this.form)">mep
<INPUT TYPE="CHECKBOX" NAME="translucent" onChange="update_checkbox_translucent(this.form)">translucent
<br />
</FORM>
""")

#############################################################################
def ViewComplex(furlP, furlL):
  js_init = ("""\
var jmolscript='';
jmolscript+='set showHydrogens off; set measurements angstroms;';
jmolscript+='load FILES "'''+furlP+'''" "'''+furlL+'''";';
jmolscript+='select 1.1; wireframe on; cpk off;';
//jmolscript+='isosurface ID psurf RESOLUTION 0 SOLVENT 0.0 COLOR yellow; ';
jmolscript+='model 2; select 2.1; cpk on; center 2.1; zoom 200; ';
jmolscript+='model ALL; display ALL; ';
jmolscript+='slab on; slab 80; depth 40; ';
jmolscript+='hover %%a %%D (%%n,%%c)|%%.2x,%%.2y,%%.2z;';
Jmol.script(JMOL1, jmolscript);
""")
  PrintHeader(title=TITLE, js=JavaScript(), initjs=js_init, css='')
  print("""\
<DIV ID="JMOLDIV1" style="width:100%%;height:85%%;background-color:white"></DIV>
<FORM>
<b>protein:</b>
<INPUT TYPE="RADIO" NAME="protein" VALUE="wireframe" onChange="update_radio_protein(this.form)" CHECKED>wireframe
<INPUT TYPE="RADIO" NAME="protein" VALUE="ballstick" onChange="update_radio_protein(this.form)">ballstick
<INPUT TYPE="RADIO" NAME="protein" VALUE="cartoon" onChange="update_radio_protein(this.form)">cartoon
&nbsp; &nbsp; &nbsp;
<b>surface:</b>
<INPUT TYPE="RADIO" NAME="psurface" VALUE="none" onChange="update_radio_psurface(this.form)" CHECKED>none
<INPUT TYPE="RADIO" NAME="psurface" VALUE="vdw" onChange="update_radio_psurface(this.form)">vdw
<INPUT TYPE="CHECKBOX" NAME="translucent" onChange="update_checkbox_translucent(this.form)">translucent
<br />
<b>ligand:</b>
<INPUT TYPE="RADIO" NAME="ligand" VALUE="ballstick" onChange="update_radio_ligand(this.form)">ballstick
<INPUT TYPE="RADIO" NAME="ligand" VALUE="spacefill" onChange="update_radio_ligand(this.form)" CHECKED>spacefill
<INPUT TYPE="RADIO" NAME="ligand" VALUE="wireframe" onChange="update_radio_ligand(this.form)">wireframe
<br />
showinfo:<INPUT TYPE="CHECKBOX" NAME="showinfo" onChange="Jmol.showInfo(JMOL1, this.form.showinfo.checked)">
</FORM>
</DIV>
""")

#############################################################################
def ViewOverlay(furl):
  js_init = ("""\
Jmol.loadFile(JMOL1, '"""+furl+"""');
var jmolscript='';
jmolscript+='set showHydrogens off; set measurements angstroms;';
jmolscript+='model 1.0;';
jmolscript+='dots ON; ';
Jmol.script(JMOL1, jmolscript);
""")
  PrintHeader(title=TITLE, js=JavaScript(), initjs=js_init, css='')
  print("""\
<DIV ID="JMOLDIV1" style="width:100%%;height:85%%;background-color:white"></DIV>
<FORM>
<b>show:</b>
<INPUT TYPE="RADIO" NAME="molsvisible" VALUE="refmol" onChange="update_radio_molsvisible(this.form)" >refmol
<INPUT TYPE="RADIO" NAME="molsvisible" VALUE="dbmol" onChange="update_radio_molsvisible(this.form)" >dbmol
<INPUT TYPE="RADIO" NAME="molsvisible" VALUE="both" onChange="update_radio_molsvisible(this.form)" CHECKED>both
<br />
refmol green:<INPUT TYPE="CHECKBOX" NAME="refmolgreen" onChange="update_checkbox_refmolgreen(this.form)">
<br />
Hs:<INPUT TYPE="CHECKBOX" NAME="hydrogens" onChange="update_checkbox_hydrogens(this.form)">
dots:<INPUT TYPE="CHECKBOX" NAME="dots" onChange="update_checkbox_dots(this.form)" CHECKED>
labels:<INPUT TYPE="CHECKBOX" NAME="labels" onChange="update_checkbox_labels(this.form)">
spin:<INPUT TYPE="CHECKBOX" NAME="spin" onChange="update_checkbox_spin(this.form)">
showinfo:<INPUT TYPE="CHECKBOX" NAME="showinfo" onChange="Jmol.showInfo(JMOL1, this.form.showinfo.checked)">
</FORM>
</DIV>
""")

#############################################################################
def ViewOverlay1xN(furlQ, furlD):
  js_init = ("""\
var jmolscript='';
jmolscript+='set showHydrogens off; set measurements angstroms;';
jmolscript+='load FILES "'''+furlQ+'''" "'''+furlD+'''";';
jmolscript+='center 1.1; zoom IN; ';
jmolscript+='dots ON; ';
jmolscript+='display 1.1,2.2; ';
jmolscript+='model 2; ';
jmolscript+='display ALL; ';
jmolscript+='set backgroundModel 1.1; ';
Jmol.script(JMOL1, jmolscript);
""")
  PrintHeader(title=TITLE, js=JavaScript(), initjs=js_init, css='')
  print("""\
<DIV ID="JMOLDIV1" style="width:100%%;height:85%%;background-color:white"></DIV>
<FORM>
<BUTTON TYPE="BUTTON" onClick="Jmol.script(JMOL1, 'model NEXT')">next</BUTTON>
<BUTTON TYPE="BUTTON" onClick="Jmol.script(JMOL1, 'model PREVIOUS')">previous</BUTTON>
<BUTTON TYPE="BUTTON" onClick="Jmol.script(JMOL1, 'model 2.2')">1st</BUTTON>
<BUTTON TYPE="BUTTON" onClick="Jmol.script(JMOL1, 'model 2; display ALL')">all</BUTTON>
<br />

<b>show:</b>
<INPUT TYPE="RADIO" NAME="molsvisible" VALUE="refmol" onChange="update_radio_molsvisible(this.form)" >refmol
<INPUT TYPE="RADIO" NAME="molsvisible" VALUE="dbmol" onChange="update_radio_molsvisible(this.form)" >dbmol
<INPUT TYPE="RADIO" NAME="molsvisible" VALUE="both" onChange="update_radio_molsvisible(this.form)" CHECKED>both
<br />
refmol green:<INPUT TYPE="CHECKBOX" NAME="refmolgreen" onChange="update_checkbox_refmolgreen(this.form)">
<br />
Hs:<INPUT TYPE="CHECKBOX" NAME="hydrogens" onChange="update_checkbox_hydrogens(this.form)">
dots:<INPUT TYPE="CHECKBOX" NAME="dots" onChange="update_checkbox_dots(this.form)" CHECKED>
labels:<INPUT TYPE="CHECKBOX" NAME="labels" onChange="update_checkbox_labels(this.form)">
spin:<INPUT TYPE="CHECKBOX" NAME="spin" onChange="update_checkbox_spin(this.form)">
showinfo:<INPUT TYPE="CHECKBOX" NAME="showinfo" onChange="Jmol.showInfo(JMOL1, this.form.showinfo.checked)">
</FORM>
</DIV>
""")

#############################################################################
def ViewEmpty():
  js_init = ("""\
var jmolscript='';
jmolscript+='set showHydrogens off; set measurements angstroms;';
Jmol.script(JMOL1, jmolscript);
""")
  PrintHeader(title=TITLE, js=JavaScript(), initjs=js_init, css='')
  print("""\
<DIV ID="JMOLDIV1" style="width:100%%;height:85%%;background-color:white"></DIV>
<FORM>
Hs:<INPUT TYPE="CHECKBOX" NAME="hydrogens" onChange="update_checkbox_hydrogens(this.form)">
dots:<INPUT TYPE="CHECKBOX" NAME="dots" onChange="update_checkbox_dots(this.form)" CHECKED>
labels:<INPUT TYPE="CHECKBOX" NAME="labels" onChange="update_checkbox_labels(this.form)">
spin:<INPUT TYPE="CHECKBOX" NAME="spin" onChange="update_checkbox_spin(this.form)">
showinfo:<INPUT TYPE="CHECKBOX" NAME="showinfo" onChange="Jmol.showInfo(JMOL1, this.form.showinfo.checked)">
</FORM>
</DIV>
""")

#############################################################################
def Initialize():
  global FORM,FILEA,FILEB,TITLE,ERRORS,MODE,APPNAME
  APPNAME = 'MolView3D'
  ERRORS=[];
  FORM = cgi.FieldStorage(keep_blank_values=1)
  TITLE = FORM.getvalue('title', 'JSMol')
  MODE = FORM.getvalue('mode')
  global WIDTH,HEIGHT
  if FORM.getvalue('width'):
    WIDTH = int(FORM.getvalue('width'))
  elif FORM.getvalue('size'):
    WIDTH = int(FORM.getvalue('size'))
  else:
    WIDTH=500
  if FORM.getvalue('height'):
    HEIGHT = int(FORM.getvalue('height'))
  elif FORM.getvalue('size'):
    HEIGHT = int(FORM.getvalue('size'))
  else:
    HEIGHT=500

  # FILEA and FILEB must be URLs, since JS runs on client!
  FILEA = FORM.getvalue('file')
  if not FILEA: FILEA = FORM.getvalue('filecode')
  if not FILEA: FILEA = FORM.getvalue('fileA')
  if not FILEA: FILEA = FORM.getvalue('filecodeA')
  if not FILEA:
    ERRORS.append('NOTE: no input specified.  Empty viewer mode.')
    return True

  FILEB = FORM.getvalue('fileB')
  if not FILEB: FILEB = FORM.getvalue('filecodeB')
  if FILEA: FILEA = urllib.parse.unquote(FILEA)
  if FILEB: FILEB = urllib.parse.unquote(FILEB)

  #Convert URLs to filepaths.
  fpath = util.http.env_cgi.scratchdir+re.sub(r'^.*/', '/', FILEA)
  if not os.access(fpath, os.R_OK):
    ERRORS.append('ERROR: file not found: '+fpath)
    return False
  if re.search('complex', MODE) or MODE=='overlay1xN':
    fpath = util.http.env_cgi.scratchdir+re.sub(r'^.*/', '/', FILEB)
    if not os.access(fpath, os.R_OK):
      ERRORS.append('ERROR: file not found: '+fpath)
      return False
  global MOLNAMEA, MOLNAMEB
  MOLNAMEA = FORM.getvalue('molnameA')
  if not MOLNAMEA:
    MOLNAMEA = os.path.basename(FILEA) if FILEA else 'MolA'
  MOLNAMEB = FORM.getvalue('molnameB')
  if not MOLNAMEB:
    MOLNAMEB = os.path.basename(FILEB) if FILEB else 'MolB'

  return True

#############################################################################
### The initjs arg is where the molecule file[s] must be loaded,
### since this must take place after the page and applet are loaded.
### This is needed to allow resizing, and relative (%) height/width.
#############################################################################
def PrintHeader(title='', js='', initjs='', css=''):
  print('Content-type: text/html\n\n<HTML STYLE="height:100%">')
  print('<HEAD><TITLE>'+title+'</TITLE>')
  print('<SCRIPT SRC="'+util.http.env_cgi.HTML_SUBDIR+'/js/biocomp.js"></SCRIPT>')
  print('<SCRIPT SRC="'+util.http.env_cgi.HTML_SUBDIR+'/jsmol/js/JSmoljQuery.js"></SCRIPT>')
  print('<SCRIPT SRC="'+util.http.env_cgi.HTML_SUBDIR+'/jsmol/js/JSmolCore.js"></SCRIPT>')
  print('<SCRIPT SRC="'+util.http.env_cgi.HTML_SUBDIR+'/jsmol/js/JSmolApplet.js"></SCRIPT>')
  print('<SCRIPT SRC="'+util.http.env_cgi.HTML_SUBDIR+'/jsmol/js/JSmolApi.js"></SCRIPT>')
  print('<SCRIPT SRC="'+util.http.env_cgi.HTML_SUBDIR+'/jsmol/js/j2sjmol.js"></SCRIPT>')
  print('<SCRIPT SRC="'+util.http.env_cgi.HTML_SUBDIR+'/jsmol/js/JSmol.js"></SCRIPT>')
  print('''\
<script type="text/javascript">

//Jmol.debugCode=true;
var JMOL1;

;(function() {

jmol_isReady = function(applet) {
  document.title = ("JSmol ("+applet._id+") is ready (UNM Translational Informatics)");
  Jmol._getElement(applet, "appletdiv").style.border="2px solid lightblue";
}
JSMolParams = {
        width: "100%%",
        height: "100%%",
        debug: false,
        color: "white",
        addSelectionOptions: false,
        serverURL: "http://chemapps.stolaf.edu/jmol/jsmol/jsmol.php",
        use: "HTML5",
        j2sPath: "'''+util.http.env_cgi.HTML_SUBDIR+'''/jsmol/j2s",
        readyFunction: jmol_isReady,
        script: "background black; set antialiasDisplay;"
};
})();

Jmol.setDocument(document);

//This renders the applet in the pre-sized div after the page loads.
$(document).ready(function(){
  Jmol.setDocument(0);
  JMOL1 = Jmol.getApplet("JMOL1", JSMolParams)    
  $("#JMOLDIV1").html(Jmol.getAppletHtml(JMOL1));
'''+initjs+'''
});
</script>
''')
  print('<SCRIPT SRC="'+util.http.env_cgi.HTML_SUBDIR+'/js/tip_main.js"></SCRIPT>')
  print('<SCRIPT SRC="'+util.http.env_cgi.HTML_SUBDIR+'/js/tip_style.js"></SCRIPT>')
  print('<SCRIPT SRC="'+util.http.env_cgi.HTML_SUBDIR+'/js/ddtip.js"></SCRIPT>')
  if js: print('<SCRIPT LANGUAGE="JavaScript">'+js+'</SCRIPT>')
  print('<LINK REL="stylesheet" type="text/css" HREF="'+util.http.env_cgi.HTML_SUBDIR+'/css/biocomp.css" />')
  if css: print('<STYLE TYPE="text/css">'+css+'</STYLE>')
  print('</HEAD>')
  print('<BODY HEIGHT="100%" BGCOLOR="#DDDDDD">')
  print('<DIV ID="TipLayer"></DIV>')

#############################################################################
if __name__=='__main__':
  ok=Initialize()
  if not ok:
    PrintHeader(title=TITLE)
    ERRORS.append('ERROR: no input; may be due to session timeout.')
  elif MODE=='multiconf':	##Used by omega.cgi
    ViewMulticonformer(FILEA)
  elif MODE=='complex':	##Used by fred.cgi
    ViewComplex(FILEA, FILEB)
  elif MODE=='overlay':	##Used by rocs.cgi, 1x1 mode.
    ViewOverlay(FILEA)
  elif MODE=='overlay1xN':	##Used by rocs.cgi, 1xN mode.
    ViewOverlay1xN(FILEA, FILEB)
  elif FILEA:
    ViewMolecule(FILEA, MOLNAMEA)	##Used by rdk_dgeom.cgi.
  else:
    ViewEmpty()
  util.http.Utils.PrintFooter(ERRORS)
