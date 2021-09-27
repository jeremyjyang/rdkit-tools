#!/usr/bin/env python3
#############################################################################
import sys,os,re,argparse,tempfile,time,logging
import pylatex

#############################################################################
class SAGPlus(pylatex.StandAloneGraphic):
  '''Simple subclass of StandAloneGraphic to handle name.'''
  def __init__(self, filename,image_options, image_name=None):
    pylatex.StandAloneGraphic.__init__(self,filename,image_options)
    self.name = image_name

#############################################################################
def PILImage2SAGPlus(img):
  ## Maybe use pylatex.utils.make_temp_dir() and .rm_temp_dir()
  fd,tmpfile=tempfile.mkstemp('.png',os.path.basename(__name__))
  img.save(tmpfile, "PNG")
  img_opts = 'width=%dpx,height=%dpx'%(img.width,img.height)
  sagp = SAGPlus(tmpfile, image_options=img_opts, image_name=img.info['name'])
  return sagp

#############################################################################
def WriteImageGrid(imgs, img_width, img_height, npr, npc, title, ofile):
  title = title if title else 'Depictions'

  geo_opts = {
	"includeheadfoot":False,
	"headheight":"10pt",
	"headsep":"8pt",
	"landscape":False,
	"tmargin": "2cm", 
	"lmargin": "2cm", 
	"rmargin":"2cm"}
  doc_opts = ["usletter", "12pt"]
  doc = pylatex.Document(geometry_options=geo_opts, document_options=doc_opts,
	documentclass="article",
	inputenc="utf8",
	indent=False, page_numbers=False)

  header = pylatex.PageStyle("header")
  with header.create(pylatex.Head("R")) as head:
    head.append('N = %d'%len(imgs))
  with header.create(pylatex.Foot("C")) as foot:
    foot.append('%s'%time.strftime('%Y-%m-%d ',time.localtime()))
    foot.append(pylatex.utils.italic('Powered by PyLaTeX (%s).'%os.path.basename(__name__)))
  doc.preamble.append(header)
  doc.change_document_style("header")

  with doc.create(pylatex.Section(title, numbering=False)) as section:
    i_img=0; i_page=0;
    while i_img<len(imgs): #page-loop
      if i_page>0: doc.append(pylatex.NewPage())
      i_page+=1
      page_row=0;
      with doc.create(pylatex.Tabular('|c'*npr +'|')) as table:
        table.add_hline()
        while i_img<len(imgs):
          page_row+=1
          imgs_this = imgs[i_img:i_img+npr]
          imgnames_this = [img.name for img in imgs_this]
          logging.debug('DEBUG: page %d ; row %d (%s)...'%(i_page,page_row,','.join(imgnames_this)))
          table.add_row(imgs_this+['' for j in range(npr-len(imgs_this))])
          #table.add_hline()
          table.add_row(imgnames_this+['' for j in range(npr-len(imgs_this))])
          table.add_hline()
          i_img+=len(imgs_this)
          if (page_row%npc)==0:
            break #page-break

  ofile=re.sub('\.pdf$','', ofile, re.I)
  for ext in ('tex', 'pdf'):
    logging.debug('Output file: %s.%s'%(ofile,ext))

  #Create .tex and .pdf files (adds extensions).
  doc.generate_pdf(ofile, clean_tex=False, silent=True)

  #Delete temp files:
  for img in imgs:
    tmpfile = str(img.arguments.__dict__['_positional_args'][0])
    #logging.debug('DEBUG: deleting %s...'%tmpfile)
    os.unlink(tmpfile)

  return doc

#############################################################################
def Test():
  FIMG_DIR = '/home/jjyang/rdkit/python/data'
  NPR = 4 # #-per-row
  NPC = 5 # #-per-col
  IMG_HEIGHT = 100
  IMG_WIDTH = 120
  OFILE = 'data/depictions' # implicit ".tex", ".pdf"

  fnames = os.listdir(FIMG_DIR)
  for i in range(len(fnames)-1,-1,-1):
    if not re.search('.png$', fnames[i], re.I):
      del fnames[i]
  fnames.sort()
  fimgs = [os.path.join(FIMG_DIR, fname) for fname in fnames]

  img_opts = 'width=%dpx,height=%dpx'%(IMG_WIDTH,IMG_HEIGHT)
  imgs = [SAGPlus(fimg, image_options=img_opts) for fimg in fimgs]

  for i in range(len(imgs)):
    imgs[i].name = fnames[i]

  for ext in ('tex', 'pdf'):
    logging.debug('DEBUG: output file: %s.%s'%(OFILE,ext))
  doc = WriteImageGrid(imgs, IMG_WIDTH, IMG_HEIGHT, NPR, NPC, "TEST", OFILE)


#############################################################################
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='PyLaTeX utility')
  parser.add_argument("op",choices=['test'],help='operation')
  args = parser.parse_args()

  if args.op=='test':
    Test()
  else:
    parser.print_help()

