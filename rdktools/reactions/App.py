#!/usr/bin/env python3
#
"""
https://www.rdkit.org/docs/source/rdkit.Chem.rdChemReactions.html
https://www.rdkit.org/docs/RDKit_Book.html
"""
import sys,os,re,math,argparse,logging

from PIL import Image, ImageDraw

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from .. import util
from .. import reactions

#############################################################################
# https://stackoverflow.com/questions/63671018/how-can-i-draw-an-arrow-using-pil
def ArrowedLine(img, ptA, ptB, width=1, color=(0, 0, 0)):
    """Draw line from ptA to ptB with arrowhead at ptB"""
    draw = ImageDraw.Draw(img)
    draw.line((ptA, ptB), width=width, fill=color)
    # Arrowhead - triangle with one vertex at ptB, extend 8*width either side of line
    x0, y0 = ptA
    x1, y1 = ptB
    # start at 95% of the length of the line
    xb = 0.95*(x1-x0)+x0
    yb = 0.95*(y1-y0)+y0
    # Other two vertices 
    if x0==x1: # if line is vertical
       vtx0 = (xb-5, yb)
       vtx1 = (xb+5, yb)
    elif y0==y1: # if line is horizontal
       vtx0 = (xb, yb+5)
       vtx1 = (xb, yb-5)
    else:
       alpha = math.atan2(y1-y0,x1-x0)-90*math.pi/180
       a = 8*width*math.cos(alpha)
       b = 8*width*math.sin(alpha)
       vtx0 = (xb+a, yb+b)
       vtx1 = (xb-a, yb-b)
    #draw.point((xb,yb), fill=(255,0,0))    # DEBUG: draw point of base in red - comment out draw.polygon() below if using this line
    # Now draw the arrowhead triangle
    draw.polygon([vtx0, vtx1, ptB], fill=color)
    return img

#############################################################################
if __name__=='__main__':
  EPILOG = """For 'react' operation, reactants are specified as disconnected
components of single input molecule record. For 'enumerateLibrary', reactants for
each role are specfied from separate input files, ordered as in the SMIRKS."""
  parser = argparse.ArgumentParser(description="RDKit chemical reactions utility", epilog=EPILOG)
  OPS = ["enumerateLibrary", "react", "demo", "demo2", "demo3", "demo4",]
  OUTPUT_MODES = ["products", "reactions",]
  parser.add_argument("op", choices=OPS, help="OPERATION")
  parser.add_argument("--i", dest="ifiles", help="input file[s] (SMILES/TSV or SDF)")
  parser.add_argument("--o", dest="ofile", help="output file (SMILES) [stdout]")
  parser.add_argument("--output_mode", choices=OUTPUT_MODES, default="reactions", help=f"{'|'.join(OUTPUT_MODES)} [products]")
  parser.add_argument("--o_depict", dest="ofile_depict", help="output depiction file (PNG) [display]")
  parser.add_argument("--smirks", help="SMIRKS reaction transform")
  parser.add_argument("--kekulize", action="store_true", help="Kekulize")
  parser.add_argument("--sanitize", action="store_true", help="Sanitize")
  parser.add_argument("--depict", action="store_true", help="Depict (1st reaction or product only)")
  parser.add_argument("--header", action="store_true", help="input SMILES/TSV file has header line")
  parser.add_argument("--delim", default="\t", help="delimiter for SMILES/TSV")
  parser.add_argument("--smilesColumn", type=int, default=0, help="input SMILES column")
  parser.add_argument("--nameColumn", type=int, default=1, help="input name column")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  molReaders=[];
  if args.ifiles is not None:
    for ifile in re.split(r'[\s*,\s*]', args.ifiles):
      molReader = util.File2Molreader(ifile, args.delim, args.smilesColumn, args.nameColumn, args.header)
      molReaders.append(molReader)
  fout = open(args.ofile, "w+") if args.ofile else sys.stdout

  if args.op == "react":
    if len(molReaders)!=1:
      parser.error(f"Exactly 1 input file required for operation: {args.op}")
    osmis = reactions.React(args.smirks, molReaders[0], args.output_mode, fout)
    if args.depict:
      deph = 300; depw = 400;
      #rxn = Chem.MolFromSmarts(osmis[0]) #RDKit can't depict?
      #rxn = Chem.rdChemReactions.ReactionFromSmarts(osmis[0]) #RDKit can't depict?
      rmixsmi,pmixsmi = re.split(r">>", osmis[0])
      rmol = Chem.MolFromSmiles(rmixsmi)
      Chem.AllChem.Compute2DCoords(rmol, clearConfs=True)
      pmol = Chem.MolFromSmiles(pmixsmi)
      Chem.AllChem.Compute2DCoords(pmol, clearConfs=True)
      rimg = Chem.Draw.MolToImage(rmol, size=(depw, deph), kekulize=True, fitImage=True)
      pimg = Chem.Draw.MolToImage(pmol, size=(depw, deph), kekulize=True, fitImage=True)
      arrow_img = Image.new("RGBA", (int(depw/4), deph), "white")
      arrow_img = ArrowedLine(arrow_img, (0, int(deph/2)), (int(depw/4), int(deph/2)), width=2, color=(0, 0, 0))
      rxnimg = Image.new("RGBA", (int(depw*2.25), deph), "white")
      rxnimg.paste(rimg, (0, 0))
      rxnimg.paste(arrow_img, (depw, 0))
      rxnimg.paste(pimg, (int(depw*1.25), 0))
      rxnimg.show()
      if args.ofile_depict:
        rxnimg.save(args.ofile_depict, 'PNG')

  elif args.op == "enumerateLibrary":
    if len(molReaders)<1:
      parser.error(f"1+ input file required for operation: {args.op}")
    reactions.EnumerateLibrary(args.smirks, molReaders, args.output_mode, fout)

  elif args.op == "demo": reactions.Demo()
  elif args.op == "demo2": reactions.Demo2()
  elif args.op == "demo3": reactions.Demo3()
  elif args.op == "demo4": reactions.Demo4()
  else:
    parser.error(f"Unsupported operation: {args.op}")

