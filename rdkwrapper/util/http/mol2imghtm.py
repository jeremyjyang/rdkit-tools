#!/usr/bin/env python3
"""
These functions are supposed to be generic and work with multiple
depiction tools.  Not all parameters are implemented by all tools;
Unimplemented parameters are ignored.

Note that the zlib module compress() function does not produce a gzip
compatible file!   So we use gzip module.
"""
import re,urllib,base64

import htm_utils

#############################################################################
def Smi2ImgHtm(smi, opts, h, w, smi2img="/cgi-bin/smi2img.cgi", zoom=True, jsfunc='go_zoom_smi2img', zoomfactor=3, tipnam='', tiphtm=None):
  htm='<A HREF="javascript:void(0)"'
  if not opts: opts=''
  smicode=urllib.parse.quote(smi, safe='') #python3
  if tiphtm:
    htm+=''' onMouseOver="stm(['%(TIPNAM)s','%(TIPHTM)s'],Style[1])" onMouseOut="htm()"'''%{'TIPNAM':tipnam, 'TIPHTM':tiphtm}
  if zoom: htm+=''' onClick="%(jsfunc)s('%(SMI2IMG)s','%(SMICODE)s','%(OPTS)s',%(WZ)d,%(HZ)d)"'''%{'SMI2IMG':smi2img, 'SMICODE':smicode, 'OPTS':opts,
       'WZ':w*zoomfactor, 'HZ':h*zoomfactor, 'jsfunc':jsfunc}
  htm+='''\
><IMG BORDER="0" SRC="%(SMI2IMG)s?%(OPTS)s&h=%(H)d&w=%(W)d&smicode=%(SMICODE)s"></A>
'''%{'SMI2IMG':smi2img, 'SMICODE':smicode, 'OPTS':opts, 'W':w, 'H':h}
  return htm

#############################################################################
def Mdl2ImgHtm(mdl, opts, h, w, mdl2img="/cgi-bin/mdl2img.cgi", zoom=True, jsfunc='go_zoom_mdl2img', zoomfactor=3, tipnam='', tiphtm=None):
  htm='<A HREF="javascript:void(0)"'
  if not opts: opts=''
  opts+='&gz=true'
  mdl_gz=htm_utils.GzipString(mdl)
  mdl_gz_b64=base64.encodestring(mdl_gz)
  mdl_gz_b64=urllib.parse.quote(mdl_gz_b64.rstrip(), safe='')  ## handle newlines
  if tiphtm: htm+=''' onMouseOver="stm(['%(TIPNAM)s','%(TIPHTM)s'],Style[1])" onMouseOut="htm()"'''%{'TIPNAM':tipnam, 'TIPHTM':tiphtm}
  if zoom: htm+=''' onClick="%(jsfunc)s('%(MDL2IMG)s','%(MDLCODE)s','%(OPTS)s',%(WZ)d,%(HZ)d)"'''%{'MDL2IMG':mdl2img, 'MDLCODE':mdl_gz_b64, 'OPTS':opts,
       'WZ':w*zoomfactor, 'HZ':h*zoomfactor, 'jsfunc':jsfunc}
  htm+='''\
><IMG BORDER="0" SRC="%(MDL2IMG)s?%(OPTS)s&h=%(H)d&w=%(W)d&mdlcode=%(MDLCODE)s"></A>
'''%{'H':h, 'W':w, 'MDL2IMG':mdl2img, 'MDLCODE':mdl_gz_b64, 'OPTS':opts}
  return htm

#############################################################################
def Mol2ImgHtm(fdata, opts, h, w, mol2img="/cgi-bin/mol2img.cgi", zoom=True,
	jsfunc='go_zoom_mol2img', zoomfactor=3, tipnam='', tiphtm=None):
  htm='<A HREF="javascript:void(0)"'
  if not opts: opts='format=mdl'	## a guess; format expected
  opts+='&gz=true'
  fdata_gz=htm_utils.GzipString(fdata)
  fcode_gz_b64=base64.encodestring(fdata_gz)
  fcode_gz_b64=urllib.parse.quote(fcode_gz_b64.rstrip(), safe='')  ## handle newlines
  tip=''
  if tiphtm: htm+=''' onMouseOver="stm(['%(TIPNAM)s','%(TIPHTM)s'],Style[1])" onMouseOut="htm()"'''%{'TIPNAM':tipnam, 'TIPHTM':tiphtm}
  if zoom: htm+=''' onClick="%(jsfunc)s('%(MOL2IMG)s','%(MOLCODE)s','%(OPTS)s',%(WZ)d,%(HZ)d)"'''%{'MOL2IMG':mol2img, 'MOLCODE':fcode_gz_b64, 'OPTS':opts,
       'WZ':w*zoomfactor, 'HZ':h*zoomfactor, 'jsfunc':jsfunc}
  htm+='''\
><IMG BORDER="0" SRC="%(MOL2IMG)s?%(OPTS)s&h=%(H)d&w=%(W)d&fcode=%(MOLCODE)s"></A>
'''%{'H':h, 'W':w, 'MOL2IMG':mol2img, 'MOLCODE':fcode_gz_b64, 'OPTS':opts}
  return htm

#############################################################################
