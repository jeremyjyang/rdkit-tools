#!/usr/bin/env python3
"""
Misc utils for use in web apps.
"""
import os,sys,urllib,base64,gzip,io

HTML_SUBDIR='/biocomp'	## Needed if subdir used.  Specify with leading '/'.

#############################################################################
def PrintHeader(title='', js='', css='', js_includes=[], jme=False, jmol=False):
  print('Content-type: text/html\n\n<HTML>')
  print('<HEAD><TITLE>%s</TITLE>'%title)
  print('<SCRIPT SRC="'+HTML_SUBDIR+'/js/biocomp.js"></SCRIPT>')
  if jme: print('<SCRIPT SRC="'+HTML_SUBDIR+'/js/JME_UNM.js"></SCRIPT>')
  if jmol: print('<SCRIPT SRC="'+HTML_SUBDIR+'/jmol/Jmol.js"></SCRIPT>')
  print('<SCRIPT SRC="'+HTML_SUBDIR+'/js/tip_main.js"></SCRIPT>')
  print('<SCRIPT SRC="'+HTML_SUBDIR+'/js/tip_style.js"></SCRIPT>')
  for js_include in js_includes:
    print('<SCRIPT TYPE="text/javascript" SRC="%s"></SCRIPT>'%js_include)
  if js: print('<SCRIPT LANGUAGE="JavaScript">%s</SCRIPT>'%js)
  print('<LINK REL="stylesheet" type="text/css" HREF="'+HTML_SUBDIR+'/css/biocomp.css" />')
  if css: print('<STYLE TYPE="text/css">%s</STYLE>'%css)
  print('</HEAD>')
  print('<BODY BGCOLOR="#DDDDDD">')
  print('<DIV ID="TipLayer"></DIV>')

#############################################################################
def PrintTestText(appname, t):
  print('Content-type: text/plain\n')
  print('APPNAME\t%s'%(appname))
  print('OK\t1')
  for key,val in t.items():
    print('%s\t%s'%(key, val))

#############################################################################
def PrintOutput(outputs):
  if len(outputs)>0:
    print('<HR>' + '<BR>\n'.join(outputs))

#############################################################################
def PrintFooter(errors):
  print('<HR>')
  if errors:
    print('<BR>\n'.join(errors))
  print('</BODY></HTML>')

#############################################################################
def Cleanup(delfiles):
  for file in delfiles:
    if file and os.access(file,os.F_OK):
      os.unlink(file)

#############################################################################
def HtmTipper(htm_in, tipnam, tiphtm, href=None, width=0, color="yellow"):
  htm=('<a onMouseOver="stm([\'{}\',\'{}\'],Style[1])"'.format(tipnam, tiphtm))
  htm+=(' onMouseOut="htm()"')
  if href: htm+=(' href="{}" target="_blank"'.format(href))
  htm+=('>{}</A>'.format(htm_in))
  return htm

#############################################################################
#Maybe better?
#def HtmTipper(htm_in,tipnam,tiphtm,width=150,color="yellow"):
#  htm=('<A onMouseOver="dd_tip(\'%s<BR>%s\',\'color\',width)"'%(tipnam,tiphtm))
#  htm+=(' onMouseOut="dd_hidetip()">%s</A>'%htm_in)
#  return htm

#############################################################################
def ChartImgHtm(url, values, xmaxes, labels, title, subtitle, imgfmt, hexcolor, gnuplotver, h, w, zoomable=True):
  valuesstr=','.join(map(lambda x:("%d"%x),values))
  xmaxsstr=','.join(map(lambda x:("%.2f"%x),xmaxes))
  labelsstr=','.join(labels)
  opts=("fmt=%s"%imgfmt)
  opts+=("&color=%s"%hexcolor)
  opts+=("&gpver=%s"%gnuplotver)
  opts+=("&title=%s"%urllib.quote(title))
  opts+=("&subtitle=%s"%urllib.quote(subtitle))
  opts+=("&values=%s"%urllib.quote(valuesstr))
  opts+=("&labels=%s"%urllib.quote(labelsstr))
  opts+=("&xmaxs=%s"%urllib.quote(xmaxsstr))
  htm='''\
<IMG SRC="%(CHART2IMG)s?%(OPTS)s&h=%(H)d&w=%(W)d">
'''%{'CHART2IMG':url,'OPTS':opts,'H':h,'W':w}
  if zoomable:
    htm='''
<A HREF="javascript:void(0)"
 onClick="javascript:go_zoom_chartimg('%(CHART2IMG)s','%(OPTS)s',480,640,'chartwin')"
 >%(IMGHTM)s</A>
'''%{'IMGHTM':htm,'CHART2IMG':url,'OPTS':opts}
  return htm

#############################################################################
def ChartLinkHtm(url, values, xmaxes, labels, title, subtitle, imgfmt, hexcolor, gnuplotver, txt):
  valuesstr=','.join(map(lambda x:("%d"%x),values))
  xmaxsstr=','.join(map(lambda x:("%.2f"%x),xmaxes))
  labelsstr=','.join(labels)
  opts=("fmt=%s"%imgfmt)
  opts+=("&color=%s"%hexcolor)
  opts+=("&gpver=%s"%gnuplotver)
  opts+=("&title=%s"%urllib.quote(title))
  opts+=("&subtitle=%s"%urllib.quote(subtitle))
  opts+=("&values=%s"%urllib.quote(valuesstr))
  opts+=("&labels=%s"%urllib.quote(labelsstr))
  opts+=("&xmaxs=%s"%urllib.quote(xmaxsstr))
  return '''\
<A HREF="javascript:void(0)"
 onClick="javascript:go_zoom_chartimg('%(CHART2IMG)s','%(OPTS)s',480,640,'chartwin')"
 >%(TXT)s</A>
'''%{'CHART2IMG':url,'OPTS':opts,'TXT':txt}

#############################################################################
def DownloadString(s, fname=None):
  #sys.stdout.write('Content-type: application/x-savefile\n')
  sys.stdout.write('Content-type: application/x-force-download\n')
  sys.stdout.write('Pragma: no-cache\n')
  sys.stdout.write('Cache-Control: no-cache\n')
  if fname:
    sys.stdout.write('Content-Disposition: attachment; filename={}\n'.format(fname))
  else:
    sys.stdout.write('Content-Disposition: attachment\n')
  s = base64.decodebytes(s.encode('utf-8'))
  sys.stdout.write('Content-length: %d\n\n'%len(s))
  sys.stdout.write(s.decode("utf-8"))
  sys.stdout.flush()

#############################################################################
def DownloadFile(fpath, fname=None):
  #sys.stdout.write('Content-type: application/x-savefile\n')
  sys.stdout.write('Content-type: application/x-force-download\n')
  sys.stdout.write('Pragma: no-cache\n')
  sys.stdout.write('Cache-Control: no-cache\n')
  if fname:
    sys.stdout.write('Content-Disposition: attachment; filename=%s\n'%fname)
  else:
    sys.stdout.write('Content-Disposition: attachment\n')
  sys.stdout.write('Content-length: %s\n\n'%str(os.path.getsize(fpath)))
  sys.stdout.flush()
  os.system('cat '+fpath)
  sys.stdout.flush()

#############################################################################
def NiceBytes(bytes):
  """ Express file size in human readable format. """
  if bytes<1e3:    return "%d bytes"%bytes
  elif bytes<1e6:  return "%.1f KB"%(bytes/1e3)
  elif bytes<1e9:  return "%.1f MB"%(bytes/1e6)
  elif bytes<1e12: return "%.1f GB"%(bytes/1e9)
  else:            return "%.1f TB"%(bytes/1e12)

#############################################################################
def NiceTime(secs):
  """ Express time in human readable format. """
  s=int(secs)
  if s<60: return '%ds'%s
  m,s = divmod(s,60)
  if m<60: return '%dm:%02ds'%(m,s)
  h,m = divmod(m,60)
  if h<24: return '%dh:%02dm:%02ds'%(h,m,s)
  d,h = divmod(h,24)
  return '%dd:%02dh:%02dm:%02ds'%(d,h,m,s)

#############################################################################
def GzipString(str):
  iob = io.BytesIO()
  fgz = gzip.GzipFile(fileobj=iob,mode='wb')
  fgz.write(str)
  fgz.close()
  str_gz = iob.getvalue()
  iob.close()
  return str_gz

#############################################################################
def GunzipBytes(str_gz):
  iob = io.BytesIO(str_gz)
  fgz = gzip.GzipFile(fileobj=iob, mode='rb')
  s = fgz.read()
  fgz.close()
  return s

