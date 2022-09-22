#! /usr/bin/env python3
"""
Generates an inline image.
"""
### 
import os,sys,cgi,logging

from rdktools.depict import Mol2ImgCgi
 
#############################################################################
if __name__=='__main__':
  form = cgi.FieldStorage(keep_blank_values=1)
  logging.debug(f"sys.path:{sys.path}")
  Mol2ImgCgi.Mol2Img(form)
  sys.exit(0)
