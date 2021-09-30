#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/
https://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints
"""
###
#############################################################################
## RDKFingerprint # returns ExplicitBitVect
## args: minPath=1,maxPath=7,fpSize=2048,nBitsPerHash=4,useHs=True, tgtDensity=0.0,minSize=128
##
## GenMACCSKeys (no args) # returns SparseBitVect
##
## AllChem.GetMorganFingerprint(mol,2) # returns SparseBitVect
## AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
##
## When comparing the ECFP/FCFP fingerprints and Morgan fingerprints by
## RDKit, remember that 4 in ECFP4 corresponds to diameter of atom environments,
## while Morgan fingerprints take a radius parameter. So radius=2 roughly
## equivalent to ECFP4 and FCFP4.
#rdkit.DataStructs.BitVectToText
#rdkit.DataStructs.BitVectToBinaryText
#rdkit.DataStructs.BitVectToFPSText
#rdkit.DataStructs.CreateFromBinaryText
#rdkit.DataStructs.CreateFromBitString
#rdkit.DataStructs.CreateFromFPSText
#############################################################################
#https://www.rdkit.org/docs/source/rdkit.DataStructs.cDataStructs.html
#rdkit.DataStructs.cDataStructs.ExplicitBitVect
#############################################################################
import os,sys,re,json,time,inspect,argparse,logging,pickle,tempfile

import rdkit
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolToSmiles, MolFromSmiles
from rdkit import DataStructs
from rdkit.ML.Cluster import Murtagh

#from rdkit.Chem.Fingerprints import FingerprintMols #USING CUSTOM VERSION.
#from rdkit.Chem.Fingerprints import MolSimilarity #USING CUSTOM VERSION.
#from rdkit.Chem.Fingerprints import ClusterMols #USING CUSTOM VERSION.

from .. import fp as rdktools_fp
from .. import util as rdktools_util

#############################################################################
def ShowDetails(details):
  import inspect
  for key,val in inspect.getmembers(details):
    if not key.startswith('_'): # ignore private and protected functions
      if not inspect.ismethod(val): # ignore other methods 
        logging.debug(f"{key:>16}: {val}")

#############################################################################
def Demo():
  mols=[];
  for smi in rdktools_util.Utils.DEMOSMIS:
    mol = MolFromSmiles(re.sub(r'\s.*$', '', smi))
    mols.append(mol)
  molWriter = SmilesWriter("-", delimiter="\t", nameHeader='Name', includeHeader=True, isomericSmiles=True, kekuleSmiles=False)
  rdktools_fp.Utils.Mols2FPs_RDK(mols, molWriter)
 
#############################################################################
def DemoMorgan():
  mols=[];
  for smi in rdktools_util.Utils.DEMOSMIS:
    mol = MolFromSmiles(re.sub(r'\s.*$', '', smi))
    mols.append(mol)
  molWriter = SmilesWriter("-", delimiter="\t", nameHeader='Name', includeHeader=True, isomericSmiles=True, kekuleSmiles=False)
  rdktools_fp.Utils.Mols2FPs_Morgan(mols, 2, 1024, molWriter)
 
#############################################################################
def DemoMACCSKeys():
  mols=[];
  for smi in rdktools_util.Utils.DEMOSMIS:
    mol = MolFromSmiles(re.sub(r'\s.*$', '', smi))
    mols.append(mol)
  molWriter = SmilesWriter("-", delimiter="\t", nameHeader='Name', includeHeader=True, isomericSmiles=True, kekuleSmiles=False)
  rdktools_fp.Utils.Mols2FPs_MACCSKeys(mols, molWriter)
 
#############################################################################
def ParseArgs(args):
  """Based on Chem/Fingerprints/FingerprintMols.py"""
  details = rdktools_fp.FingerprintMols.FingerprinterDetails()
  if args.ifile: details.inFileName = args.ifile
  if args.ofile: details.outFileName = args.ofile
  if args.useHs: details.useHs = 1
  if args.useValence: details.useValence = 1
  if args.dbName: details.dbName = args.dbName
  if args.tableName: details.tableName = args.tableName
  if args.minSize: details.minSize = args.minSize
  if args.maxSize: details.fpSize = args.maxSize
  if args.density: details.tgtDensity = args.density
  if args.outTable: details.outTableName = args.outTable
  if args.outDbName: details.outDbName = args.outDbName
  if args.fpColName: details.fpColName = args.fpColName
  if args.minPath: details.minPath = args.minPath
  if args.maxPath: details.maxPath = args.maxPath
  if args.nBitsPerHash: details.bitsPerHash = args.nBitsPerHash
  if args.discrim: details.discrimHash = 1
  if args.smilesName: details.smilesName = args.smilesName
  if args.molPkl: details.molPklName = args.molPkl
  if args.useSD: details.useSmiles=False; details.useSD=True
  if args.idName: details.idName = args.idName
  if args.maxMols: details.maxMols = args.maxMols
  if args.fingerprinter=="MACCS":
    details.fingerprinter = rdkit.Chem.MACCSkeys.GenMACCSKeys
  else:  details.fingerprinter = rdkit.Chem.RDKFingerprint
  if args.keepTable: details.replaceTable = False
  if args.smilesTable: details.smilesTableName = args.smilesTable
  if args.topN: details.doThreshold = 0; details.topN = args.topN
  elif args.thresh: details.doThreshold = 1; details.screenThresh = args.thresh
  if args.smiles: details.probeSmiles = args.smiles
  if args.metric=="dice": details.metric = DataStructs.DiceSimilarity
  elif args.metric=="cosine": details.metric = DataStructs.CosineSimilarity
  if args.clusterAlgo=="SLINK": details.clusterAlgo = Murtagh.SLINK
  elif args.clusterAlgo=="CLINK": details.clusterAlgo = Murtagh.CLINK
  elif args.clusterAlgo=="UPGMA": details.clusterAlgo = Murtagh.UPGMA
  else: pass #(WARD)
  if args.actTable: details.actTableName = args.actTable
  if args.actName: details.actName = args.actName
  return details

#############################################################################
if __name__ == "__main__":
  MORGAN_NBITS=1024; MORGAN_RADIUS=2;
  parser = argparse.ArgumentParser(description="RDKit fingerprint generator", epilog="")
  OPS = [ "fpgen", "fpgen_morgan", "fpgen_maccs",
	"demo", "demo_maccs", "demo_morgan",
	"FingerprintMols", "MolSimilarity", "ClusterMols" ]
  
  parser.add_argument("op", choices=OPS, help="OPERATION")
  parser.add_argument("--scratchdir", default="/tmp")
  parser.add_argument("--smicol", type=int, default=1, help="SMILES column from TSV (counting from 1)")
  parser.add_argument("--namcol", type=int, default=2, help="name column from TSV (counting from 1)")
  parser.add_argument("--idelim", default="\t", help="delim for input TSV")
  parser.add_argument("--odelim", default="\t", help="delim for output TSV")
  parser.add_argument("--iheader", action="store_true", help="input TSV has header")
  parser.add_argument("--oheader", action="store_true", help="output TSV has header")
  parser.add_argument("--morgan_nbits", type=int, default=MORGAN_NBITS)
  parser.add_argument("--morgan_radius", type=int, default=MORGAN_RADIUS)
  parser.add_argument("--reportFreq", type=int, default=100)
  #FingerprintMols, MolSimilarity
  parser.add_argument("--i", dest="ifile", help="Input file; if provided and no tableName is specified, data will be read from the input file.  Text files delimited with either commas (extension .csv) or tabs (extension .txt) are supported.")
  parser.add_argument("--o", dest="ofile", help="Name of the output file (output will be a pickle file with one label,fingerprint entry for each molecule).")
  parser.add_argument("--useHs", action="store_true", help="Include Hs in the fingerprint Default is *false*.")
  parser.add_argument("--useValence", action="store_true", help="Include valence information in the fingerprints Default is *false*.")
  parser.add_argument("--dbName", help="Name of the database from which to pull input molecule information.  If output is going to a database, this will also be used for that unless the --outDbName option is used.")
  parser.add_argument("--tableName", help="Name of the database table from which to pull input molecule information")
  parser.add_argument("--minSize", type=int, default=64, help="Minimum size of the fingerprints to be generated (limits the amount of folding that happens).")
  parser.add_argument("--maxSize", type=int, default=2048, help="Base size of the fingerprints to be generated.")
  parser.add_argument("--density", type=float, default=0.3, help="Target bit density in the fingerprint.  The fingerprint will be folded until this density is reached.")
  parser.add_argument("--outTable", help="name of the output db table used to store fingerprints.  If this table already exists, it will be replaced.")
  parser.add_argument("--outDbName", help="name of output database, if it's being used.  Defaults to be the same as the input db.")
  parser.add_argument("--fpColName", default="AutoFragmentFP", help="name to use for the column which stores fingerprints (in pickled format) in the output db table.")
  parser.add_argument("--minPath", type=int, default=1, help="Minimum path length to be included in fragment-based fingerprints.")
  parser.add_argument("--maxPath", type=int, default=7, help="Maximum path length to be included in fragment-based fingerprints.")
  parser.add_argument("--nBitsPerHash", type=int, default=2, help="Number of bits to be set in the output fingerprint for each fragment.")
  parser.add_argument("--discrim", action="store_true", help="Use of path-based discriminators to hash bits.")
  parser.add_argument("--smilesName", default="#SMILES", help="Name of the SMILES column in the input database.")
  parser.add_argument("--molPkl", help="")
  parser.add_argument("--useSD", action="store_true", help="Assume that the input file is an SD file, not a SMILES table.")
  parser.add_argument("--idName", default="Name", help="Name of the id column in the input database.  Defaults to the first column for dbs.")
  parser.add_argument("--maxMols", type=int, help="Maximum number of molecules to be fingerprinted.")
  parser.add_argument("--fingerprinter", default="RDKIT", choices=["RDKIT", "MACCS"], help="RDKIT: daylight-type; MACCS: MACCS 166 keys")
  parser.add_argument("--keepTable", action="store_true", help="")
  # SCREENER:
  parser.add_argument("--smilesTable", help="")
  parser.add_argument("--topN", type=int, default=12, help="Top N similar; precedence over threshold.")
  parser.add_argument("--thresh", type=float, help="Similarity threshold.")
  parser.add_argument("--smiles", help="Query smiles for similarity screening.")
  parser.add_argument("--metric", choices=["dice", "cosine"], default="dice", help="Similarity algorithm")
  #CLUSTERS:
  parser.add_argument("--clusterAlgo", choices=["WARD", "SLINK", "CLINK", "UPGMA"],
default="WARD", help="Clustering algorithm: WARD = Ward's minimum variance; SLINK = single-linkage clustering algorithm; CLINK = complete-linkage clustering algorithm; UPGMA = group-average clustering algorithm")
  parser.add_argument("--actTable", help="name of table containing activity values (used to color points in the cluster tree).")
  parser.add_argument("--actName", help="name of column with activities in the activity table. The values in this column should either be integers or convertible into integers.")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f"RDKit version: {rdkit.rdBase.rdkitVersion}")
  #logging.info(f"MatplotLib version: {mpl.__version__}")

  t0=time.time()

  if args.op=="FingerprintMols": 
    details = ParseArgs(args)
    rdktools_fp.FingerprintMols.FingerprintsFromDetails(details, reportFreq=args.reportFreq)
    sys.exit()

  elif args.op=="MolSimilarity":
    details = ParseArgs(args)
    ftmp = tempfile.NamedTemporaryFile(suffix='.pkl', delete=True)
    details.outFileName = ftmp.name
    logging.debug(f"Temporary file: {ftmp.name}")
    ftmp.close()
    ShowDetails(details)
    queryMol = MolFromSmiles(re.sub(r'\s.*$', '', args.smiles))
    queryName = re.sub(r'^[\S]*\s', '', args.smiles)
    logging.info("MolSimilarity ({0}): Query: {1}".format((f"TopN:{args.topN}" if args.topN else f"threshold: {args.thresh}"), (f"{queryName}" if queryName else f"{querySmiles}")))
    rdktools_fp.FingerprintMols.FingerprintsFromDetails(details, reportFreq=args.reportFreq)
    n_fp=0; bvs=[];
    with open(ftmp.name, "rb") as fin:
      while True:
        try:
          ID,bv = pickle.load(fin)
          bvs.append((ID, bv))
        except Exception as e:
          break
        n_fp+=1
        #logging.debug(f"{n_fp}. {ID} ({len(bv)}): {bv.ToBitString()}")
    os.remove(ftmp.name)
    details.outFileName = args.ofile
    results = rdktools_fp.MolSimilarity.ScreenFingerprints(details, bvs, mol=queryMol, probeFp=None)
    n_hit=0; 
    results.reverse()
    for ID,score in results:
      n_hit+=1
      logging.debug(f"{n_hit:2d}. {ID:>16}: {score:.3f}")
    sys.exit()

  elif args.op=="ClusterMols": 
    details = ParseArgs(args)
    ftmp = tempfile.NamedTemporaryFile(suffix='.pkl', delete=True)
    details.outFileName = ftmp.name
    logging.debug(f"Temporary file: {ftmp.name}")
    ftmp.close()
    rdktools_fp.FingerprintMols.FingerprintsFromDetails(details, reportFreq=args.reportFreq)
    n_fp=0; bvs=[];
    with open(ftmp.name, "rb") as fin:
      while True:
        try:
          ID,bv = pickle.load(fin)
          bvs.append((ID, bv))
        except Exception as e:
          break
        #logging.debug(f"{n_fp}. {ID}; type(bv): {type(bv)}")
        n_fp+=1
    os.remove(ftmp.name)
    details.outFileName = args.ofile

    logging.info(f"ClusterMols ({args.fingerprinter}, {args.metric}, {args.clusterAlgo})")
    clustTree = rdktools_fp.ClusterMols.ClusterPoints(bvs, details.metric, details.clusterAlgo, haveLabels=0, haveActs=1)
    if args.ofile:
      pickle.dump(clustTree, args.ofile)
    sys.exit()

  elif args.op=="demo":
    Demo()
    sys.exit()

  elif args.op=="demo_morgan":
    DemoMorgan()
    sys.exit()

  elif args.op=="demo_maccs":
    DemoMACCSKeys()
    sys.exit()

  #if not (args.ifile): parser.error('--i required.')

  molReader = rdktools_util.Utils.File2Molreader(args.ifile, args.idelim, args.smicol-1, args.namcol-1, args.iheader)
  molWriter = rdktools_util.Utils.File2Molwriter(args.ofile, args.odelim, args.oheader)
  mols = rdktools_util.Utils.ReadMols(molReader)

  if args.op=="fpgen":
    rdktools_fp.Utils.Mols2FPs_RDK(mols, molWriter)

  elif args.op=="fpgen_morgan":
    rdktools_fp.Utils.Mols2FPs_Morgan(mols, args.morgan_radius, args.morgan_nbits, molWriter)

  elif args.op=="fpgen_maccs":
    rdktools_fp.Utils.Mols2FPs_MACCS(mols, molWriter)

  else:
    parser.error(f"Unsupported operation: {args.op}")

  logging.info('Elapsed: %s'%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))))

