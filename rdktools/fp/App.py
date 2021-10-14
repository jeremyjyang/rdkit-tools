#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/

When comparing the ECFP/FCFP fingerprints and Morgan fingerprints by
RDKit, remember that 4 in ECFP4 corresponds to diameter of atom environments,
while Morgan fingerprints take a radius parameter. So radius=2 roughly
equivalent to ECFP4 and FCFP4.

https://www.rdkit.org/docs/source/rdkit.DataStructs.cDataStructs.html
rdkit.DataStructs.cDataStructs.ExplicitBitVect
"""
#############################################################################
import os,sys,re,json,time,argparse,logging,pickle,tempfile

import rdkit
import rdkit.Chem.AllChem
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolToSmiles, MolFromSmiles
from rdkit import DataStructs
from rdkit.ML.Cluster import Murtagh
from rdkit.ML.Cluster import Butina
from rdkit.ML.Cluster import ClusterVis
from rdkit.ML.Cluster import ClusterUtils

#from rdkit.Chem.Fingerprints import FingerprintMols #USING CUSTOM VERSION.
#from rdkit.Chem.Fingerprints import MolSimilarity #USING CUSTOM VERSION.
#from rdkit.Chem.Fingerprints import ClusterMols #USING CUSTOM VERSION.

from .. import fp as rdktools_fp
from .. import util as rdktools_util

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
  if args.smilesColumn: details.smilesName = args.smilesColumn
  if args.molPkl: details.molPklName = args.molPkl
  details.useSmiles = bool(args.input_format=="SMILES")
  details.useSD = bool(args.input_format=="SD")
  if args.idColumn: details.idName = args.idColumn
  if args.maxMols: details.maxMols = args.maxMols
  details.fingerprinter = args.fpAlgo
  details.morgan_radius = args.morgan_radius
  details.morgan_nbits = args.morgan_nbits
  details.replaceTable = args.replaceTable
  if args.smilesTable: details.smilesTableName = args.smilesTable
  if args.topN: details.doThreshold = 0; details.topN = args.topN
  elif args.thresh: details.doThreshold = 1; details.screenThresh = args.thresh
  if args.querySmiles: details.probeSmiles = args.querySmiles
  if args.metric=="ALLBIT": details.metric = DataStructs.AllBitSimilarity
  elif args.metric=="ASYMMETRIC": details.metric = DataStructs.AsymmetricSimilarity
  elif args.metric=="DICE": details.metric = DataStructs.DiceSimilarity
  elif args.metric=="COSINE": details.metric = DataStructs.CosineSimilarity
  elif args.metric=="KULCZYNSKI": details.metric = DataStructs.KulczynskiSimilarity
  elif args.metric=="MCCONNAUGHEY": details.metric = DataStructs.McConnaugheySimilarity
  elif args.metric=="ONBIT": details.metric = DataStructs.OnBitSimilarity
  elif args.metric=="RUSSEL": details.metric = DataStructs.RusselSimilarity
  elif args.metric=="SOKAL": details.metric = DataStructs.SokalSimilarity
  elif args.metric=="TANIMOTO": details.metric = DataStructs.TanimotoSimilarity
  elif args.metric=="TVERSKY": details.metric = DataStructs.TverskySimilarity
  details.tversky_alpha = args.tversky_alpha
  details.tversky_beta = args.tversky_beta
  if args.clusterAlgo=="SLINK": details.clusterAlgo = Murtagh.SLINK
  elif args.clusterAlgo=="CLINK": details.clusterAlgo = Murtagh.CLINK
  elif args.clusterAlgo=="UPGMA": details.clusterAlgo = Murtagh.UPGMA
  elif args.clusterAlgo=="BUTINA": details.clusterAlgo = Butina.ClusterData
  else: pass #(WARD)
  if args.actTable: details.actTableName = args.actTable
  if args.actName: details.actName = args.actName
  return details

#############################################################################
if __name__ == "__main__":
  EPILOG="""This app employs custom, updated versions of RDKit
FingerprintMols.py, MolSimilarity.py, ClusterMols.py,
with enhanced command-line functionality for molecular
fingerprint-based analytics."""
  FPALGOS=["RDKIT", "MACCS", "MORGAN"]
  MORGAN_NBITS=1024; MORGAN_RADIUS=2;
  METRICS = ["ALLBIT", "ASYMMETRIC", "DICE", "COSINE", "KULCZYNSKI", "MCCONNAUGHEY", "ONBIT", "RUSSEL", "SOKAL", "TANIMOTO", "TVERSKY"]
  parser = argparse.ArgumentParser(description="RDKit fingerprint-based analytics", epilog=EPILOG)
  OPS=[ "FingerprintMols", "MolSimilarity", "ClusterMols" ]
  parser.add_argument("op", choices=OPS, help="OPERATION")
  #FingerprintMols, MolSimilarity
  parser.add_argument("--i", dest="ifile", help="input file; if provided and no tableName is specified, data will be read from the input file.  Text files delimited with either commas (extension .csv) or tabs (extension .txt) are supported.")
  parser.add_argument("--o", dest="ofile", help="name of the output file (output will be a pickle file with one label,fingerprint entry for each molecule).")
  parser.add_argument("--useHs", action="store_true", help="include Hs in the fingerprint Default is *false*.")
  parser.add_argument("--useValence", action="store_true", help="include valence information in the fingerprints Default is *false*.")
  parser.add_argument("--dbName", help="name of the database from which to pull input molecule information.  If output is going to a database, this will also be used for that unless the --outDbName option is used.")
  parser.add_argument("--tableName", help="name of the database table from which to pull input molecule information")
  parser.add_argument("--minSize", type=int, default=64, help="minimum size of the fingerprints to be generated (limits the amount of folding that happens).")
  parser.add_argument("--maxSize", type=int, default=2048, help="base size of the fingerprints to be generated.")
  parser.add_argument("--density", type=float, default=0.3, help="target bit density in the fingerprint.  The fingerprint will be folded until this density is reached.")
  parser.add_argument("--outTable", help="name of the output db table used to store fingerprints.  If this table already exists, it will be replaced.")
  parser.add_argument("--outDbName", help="name of output database, if it's being used.  Defaults to be the same as the input db.")
  parser.add_argument("--fpColName", default="AutoFragmentFP", help="name to use for the column which stores fingerprints (in pickled format) in the output db table.")
  parser.add_argument("--minPath", type=int, default=1, help="minimum path length to be included in fragment-based fingerprints.")
  parser.add_argument("--maxPath", type=int, default=7, help="maximum path length to be included in fragment-based fingerprints.")
  parser.add_argument("--nBitsPerHash", type=int, default=2, help="number of bits to be set in the output fingerprint for each fragment.")
  parser.add_argument("--discrim", action="store_true", help="use of path-based discriminators to hash bits.")
  parser.add_argument("--smilesColumn", default="#SMILES", help="name of the SMILES column in the input database.")
  parser.add_argument("--molPkl", help="")
  parser.add_argument("--input_format", choices=["SMILES", "SD"], default="SMILES", help="SMILES table or SDF file.")
  parser.add_argument("--idColumn", "--nameColumn", default="Name", help="name of the id column in the input database.  Defaults to the first column for dbs.")
  parser.add_argument("--maxMols", type=int, help="maximum number of molecules to be fingerprinted.")
  parser.add_argument("--fpAlgo", default="RDKIT", choices=FPALGOS, help="RDKIT = Daylight path-based; MACCS = MDL MACCS 166 keys")
  parser.add_argument("--morgan_nbits", type=int, default=MORGAN_NBITS)
  parser.add_argument("--morgan_radius", type=int, default=MORGAN_RADIUS)
  parser.add_argument("--replaceTable", action="store_true", help="")
  parser.add_argument("--smilesTable", help="name of database table which contains SMILES for the input fingerprints.  If provided with --smilesName, output will contain SMILES data.")
  parser.add_argument("--topN", type=int, default=12, help="top N similar; precedence over threshold.")
  parser.add_argument("--thresh", type=float, help="similarity threshold.")
  parser.add_argument("--querySmiles", help="query smiles for similarity screening.")
  parser.add_argument("--metric", choices=METRICS, default="TANIMOTO", help="similarity algorithm")
  parser.add_argument("--tversky_alpha", type=float, default=.8, help="Tversky alpha parameter, weights query molecule features")
  parser.add_argument("--tversky_beta", type=float, default=.2, help="Tversky beta parameter, weights target molecule features")
  parser.add_argument("--clusterAlgo", choices=["WARD", "SLINK", "CLINK", "UPGMA", "BUTINA"], default="WARD", help="clustering algorithm: WARD = Ward's minimum variance; SLINK = single-linkage clustering algorithm; CLINK = complete-linkage clustering algorithm; UPGMA = group-average clustering algorithm; BUTINA = Butina JCICS 39 747-750 (1999)")
  parser.add_argument("--actTable", help="name of table containing activity values (used to color points in the cluster tree).")
  parser.add_argument("--actName", help="name of column with activities in the activity table. The values in this column should either be integers or convertible into integers.")
  parser.add_argument("--reportFreq", type=int, default=100)
  parser.add_argument("--showVis", action="store_true", help="show visualization if available.")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f"RDKit version: {rdkit.rdBase.rdkitVersion}")

  t0=time.time()

  if args.op=="FingerprintMols": 
    logging.info("FingerprintMols ({0})".format(f"{args.fpAlgo}({args.morgan_radius},{args.morgan_nbits})" if args.fpAlgo=="MORGAN" else args.fpAlgo))
    details = ParseArgs(args)
    rdktools_fp.ShowDetails(details)
    rdktools_fp.FingerprintMols.FingerprintsFromDetails(details, reportFreq=args.reportFreq)
    if args.ofile is not None: logging.info(f"FPs written to: {args.ofile}")

  elif args.op=="MolSimilarity":
    logging.info("MolSimilarity ({0}, {1})".format(args.fpAlgo, f"{args.metric}({args.tversky_alpha},{args.tversky_beta})" if args.metric=="TVERSKY" else args.metric))
    details = ParseArgs(args)
    ftmp = tempfile.NamedTemporaryFile(suffix='.pkl', delete=True)
    details.outFileName = ftmp.name
    logging.debug(f"Temporary file: {ftmp.name}")
    ftmp.close()
    rdktools_fp.ShowDetails(details)
    queryMol = MolFromSmiles(re.sub(r'\s.*$', '', args.querySmiles))
    if queryMol is None:
      logging.error(f"Failed to parse SMILES: {args.querySmiles}")
    queryName = re.sub(r'^[\S]*\s', '', args.querySmiles)
    logging.info("{0}; Query: {1}".format((f"TopN:{args.topN}" if args.topN else f"threshold: {args.thresh}"), (f"{queryName}" if queryName else f"{querySmiles}")))
    rdktools_fp.FingerprintMols.FingerprintsFromDetails(details, reportFreq=args.reportFreq)
    n_fp=0; fps=[];
    with open(ftmp.name, "rb") as fin:
      while True:
        try:
          id_this,fp = pickle.load(fin)
          fps.append((id_this, fp))
        except Exception as e:
          break
        n_fp+=1
    os.remove(ftmp.name)
    details.outFileName = args.ofile
    queryFp = rdktools_fp.FingerprintMols.FingerprintMol(queryMol, **details.__dict__)
    results = rdktools_fp.MolSimilarity.ScreenFingerprints(details, fps, probeFp=queryFp)
    n_hit=0; 
    results.reverse()
    fout = open(args.ofile, "w+") if args.ofile else sys.stdout
    for id_this,score in results:
      n_hit+=1
      fout.write(f"{n_hit}\t{queryName}\t{id_this}\t{score:.3f}\n")

  elif args.op=="ClusterMols": 
    details = ParseArgs(args)
    ftmp = tempfile.NamedTemporaryFile(suffix='.pkl', delete=True)
    details.outFileName = ftmp.name
    logging.debug(f"Temporary file: {ftmp.name}")
    ftmp.close()
    rdktools_fp.FingerprintMols.FingerprintsFromDetails(details, reportFreq=args.reportFreq)
    n_fp=0; fps=[];
    with open(ftmp.name, "rb") as fin:
      while True:
        try:
          id_this,fp = pickle.load(fin)
          fps.append((id_this, fp))
        except Exception as e:
          break
        n_fp+=1
    os.remove(ftmp.name)
    details.outFileName = args.ofile

    logging.info(f"ClusterMols ({args.fpAlgo}, {args.metric}, {args.clusterAlgo})")
    clustTree, dMat = rdktools_fp.ClusterMols.ClusterPoints(fps, details.metric, details.clusterAlgo, haveLabels=0, haveActs=1, returnDistances=True)

    logging.info(f"dMat.ndim:{dMat.ndim}; dMat.shape:{dMat.shape}; dMat.size:{dMat.size}")

    fout = open(args.ofile, "w+") if args.ofile else sys.stdout

    nodes = ClusterUtils.GetNodeList(clustTree)
    for node in nodes:
      node.Print(level=0, showData=1, offset="  ")
    for node in nodes:
      logging.debug(f"{node.GetIndex()}; IsTerminal:{node.IsTerminal()}; len(GetChildren()):{len(node.GetChildren())}")
      #for child in node.GetChildren():
      #  fout.write(f"{node.GetIndex()}\t{node.getName()}\t{child.GetIndex()}\t{child.getName()}\n")

    if args.showVis:
      from PIL import Image
      ftmp = tempfile.NamedTemporaryFile(suffix='.png', delete=True)
      ftmp.close()
      ClusterVis.ClusterToImg(clustTree, ftmp.name, (400,600), ptColors=[], lineWidth=None, showIndices=0, stopAtCentroids=0, logScale=0)
      img = Image.open(ftmp.name)
      logging.debug(f"img.size: {img.size}")
      img.show()

    if args.ofile:
      pickle.dump(clustTree, args.ofile)

  else:
    parser.error(f"Unsupported operation: {args.op}")

  logging.info('Elapsed: %s'%(time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))))

