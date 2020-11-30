#!/usr/bin/env python3
###
import os,sys,cgi

import mol2img_rdk

#############################################################################
if __name__=='__main__':
  form=cgi.FieldStorage(keep_blank_values=1)
  mol2img_rdk.Mol2Img(form)
