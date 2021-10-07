###
#  Copyright (c) 2006-2015, Rational Discovery LLC, Greg Landrum, and Julie Penzotti and others
#
#  All Rights Reserved.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" utility functionality for molecular similarity
 includes a command line app for screening databases


Sample Usage:

  python MolSimilarity.py  -d data.gdb -t daylight_sig --idName="Mol_ID" \
      --topN=100 --smiles='c1(C=O)ccc(Oc2ccccc2)cc1' --smilesTable=raw_dop_data \
      --smilesName="structure" -o results.csv

"""
import pickle,logging

from rdkit import Chem
from rdkit import DataStructs
from rdkit.DataStructs.TopNContainer import TopNContainer
from rdkit.Dbase import DbModule
from rdkit.Dbase.DbConnection import DbConnect

from rdkit.Chem.Fingerprints import DbFpSupplier
#from rdkit.Chem.Fingerprints import FingerprintMols #USE CUSTOM VERSION.
from .. import fp as rdktools_fp

try:
  from rdkit.VLib.NodeLib.DbPickleSupplier import _lazyDataSeq as _dataSeq
except ImportError:
  _dataSeq = None


def _ConstructSQL(details, extraFields=''):
  fields = f"{details.tableName}.{details.idName}"
  join = ""
  if details.smilesTableName:
    if details.smilesName:
      fields = fields + f",{details.smilesName}"
    join = f"JOIN {details.smilesTableName} smi ON smi.{details.idName}={details.tableName}.{details.idName}"
                                            
  if details.actTableName:
    if details.actName:
      fields = fields + f",{details.actName}"
    join = join + f"JOIN {details.actTableName} act ON act.{details.idName}={details.tableName}.{details.idName}"

  # data = conn.GetData(fields=fields,join=join)
  if extraFields:
    fields += "," + extraFields
  cmd = f"SELECT {fields} FROM {details.tableName} {join}"
  return cmd


def ScreenInDb(details, mol):
  try:
    probeFp = rdktools_fp.FingerprintMols.FingerprintMol(mol, **details.__dict__)
  except Exception:
    import traceback
    logging.error("Problems fingerprinting molecule.")
    traceback.print_exc()
    return []
  if details.dbName and details.tableName:
    try:
      conn = DbConnect(details.dbName, details.tableName)
      if hasattr(details, "dbUser"):
        conn.user = details.dbUser
      if hasattr(details, "dbPassword"):
        conn.password = details.dbPassword
    except Exception:
      import traceback
      logging.error(f"Problems establishing connection to database: {details.dbName}|{details.tableName}\n")
      traceback.print_exc()

  if details.metric not in (DataStructs.TanimotoSimilarity, DataStructs.DiceSimilarity, DataStructs.CosineSimilarity):
    data = GetFingerprints(details)
    res = ScreenFingerprints(details, data, mol)
  else:
    res = []
    if details.metric == DataStructs.TanimotoSimilarity:
      func = "rd_tanimoto"
      pkl = probeFp.ToBitString()
    elif details.metric == DataStructs.DiceSimilarity:
      func = "rd_dice"
      pkl = probeFp.ToBitString()
    elif details.metric == DataStructs.CosineSimilarity:
      func = "rd_cosine"
      pkl = probeFp.ToBitString()
    extraFields = f"{func}({DbModule.placeHolder},{details.fpColName}) as tani"
    cmd = _ConstructSQL(details, extraFields=extraFields)

    if details.doThreshold:
      # we need to do a subquery here:
      cmd = f"SELECT * FROM ({cmd}) tmp WHERE tani>{details.screenThresh}"
    cmd += " ORDER BY tani DESC"
    if not details.doThreshold and details.topN > 0:
      cmd += f" LIMIT {details.topN}"
    curs = conn.GetCursor()
    curs.execute(cmd, (pkl, ))
    res = curs.fetchall()
  return res


def GetFingerprints(details):
  """ returns an iterable sequence of fingerprints
  each fingerprint will have a _fieldsFromDb member whose first entry is
  the id.
  """
  if details.dbName and details.tableName:
    try:
      conn = DbConnect(details.dbName, details.tableName)
      if hasattr(details, 'dbUser'):
        conn.user = details.dbUser
      if hasattr(details, 'dbPassword'):
        conn.password = details.dbPassword
    except Exception:
      import traceback
      logging.error(f"Problems establishing connection to database: {details.dbName}|{details.tableName}")
      traceback.print_exc()
    cmd = _ConstructSQL(details, extraFields=details.fpColName)
    curs = conn.GetCursor()
    # curs.execute(cmd)
    # print "CURSOR:",curs,curs.closed
    if _dataSeq:
      suppl = _dataSeq(curs, cmd, depickle=not details.noPickle, klass=DataStructs.ExplicitBitVect)
      _dataSeq._conn = conn
    else:
      suppl = DbFpSupplier.ForwardDbFpSupplier(data, fpColName=details.fpColName)
  elif details.inFileName:
    conn = None
    try:
      inF = open(details.inFileName, "r")
    except IOError:
      import traceback
      logging.error(f"Problems reading from file {details.inFileName}")
      traceback.print_exc()

    suppl = []
    done = 0
    while not done:
      try:
        ID, fp = pickle.load(inF)
      except Exception:
        done = 1
      else:
        fp._fieldsFromDb = [ID]
        suppl.append(fp)
  else:
    suppl = None

  return suppl


def ScreenFingerprints(details, data, mol=None, probeFp=None):
  """ Returns a list of results
  """
  if probeFp is None:
    try:
      probeFp = rdktools_fp.FingerprintMols.FingerprintMol(mol, **details.__dict__)
    except Exception:
      import traceback
      logging.error('Problems fingerprinting molecule.')
      traceback.print_exc()
      return []
  if not probeFp:
    return []
  res = [];
  if not details.doThreshold and details.topN > 0:
    topN = TopNContainer(details.topN)
  else:
    topN = [];
  res=[]; count=0;
  for pt in data:
    fp1 = probeFp
    if not details.noPickle:
      if type(pt) in (tuple, list):
        ID, fp = pt
      else:
        fp = pt
        ID = pt._fieldsFromDb[0]
      if details.metric is DataStructs.TverskySimilarity:
        score = DataStructs.TverskySimilarity(fp1, fp, details.tversky_alpha, details.tversky_beta)
      else:
        score = DataStructs.FingerprintSimilarity(fp1, fp, details.metric)
    else:
      ID, pkl = pt
      if details.metric is DataStructs.TverskySimilarity:
        score = DataStructs.TverskySimilarity(fp1, str(pkl), details.tversky_alpha, details.tversky_beta)
      else:
        score = details.metric(fp1, str(pkl))
    if topN:
      topN.Insert(score, ID)
    elif not details.doThreshold or (details.doThreshold and score >= details.screenThresh):
      res.append((ID, score))
    count += 1
    if hasattr(details, 'stopAfter') and count>=details.stopAfter:
      break
  for score, ID in topN:
    res.append((ID, score))
  return res


def ScreenFromDetails(details, mol=None):
  """ Returns a list of results
  """
  if not mol:
    if not details.probeMol:
      smi = details.probeSmiles
      try:
        mol = Chem.MolFromSmiles(smi)
      except Exception:
        import traceback
        logging.error(f"Problems generating molecule for smiles: {smi}")
        traceback.print_exc()
        return
    else:
      mol = details.probeMol
  if not mol:
    return

  if details.outFileName:
    try:
      outF = open(details.outFileName, 'w+')
    except IOError:
      logging.error(f"Could not open output file {details.outFileName} for writing")
      return None
  else:
    outF = None

  if not hasattr(details, 'useDbSimilarity') or not details.useDbSimilarity:
    data = GetFingerprints(details)
    res = ScreenFingerprints(details, data, mol)
  else:
    res = ScreenInDb(details, mol)
  if outF:
    for pt in res:
      outF.write(','.join([str(x) for x in pt]))
      outF.write('\n')
  return res


_usageDoc = """
Usage: MolSimilarity.py [args] <fName>

  If <fName> is provided and no tableName is specified (see below),
  data will be read from the pickled file <fName>.  This file should
  contain a series of pickled (ID,fingerprint) tuples.

  NOTE: at the moment the user is responsible for ensuring that the
  fingerprint parameters given at run time (used to fingerprint the
  probe molecule) match those used to generate the input fingerprints.

  Command line arguments are:
    - --smiles=val: sets the SMILES for the input molecule.  This is
      a required argument.

    - -d _dbName_: set the name of the database from which
      to pull input fingerprint information.

    - -t _tableName_: set the name of the database table
      from which to pull input fingerprint information

    - --smilesTable=val: sets the name of the database table
      which contains SMILES for the input fingerprints.  If this
      information is provided along with smilesName (see below),
      the output file will contain SMILES data

    - --smilesName=val: sets the name of the SMILES column
      in the input database.  Default is *SMILES*.

    - --topN=val: sets the number of results to return.
      Default is *10*.

    - --thresh=val: sets the similarity threshold.

    - --idName=val: sets the name of the id column in the input
      database.  Default is *ID*.

    - -o _outFileName_:  name of the output file (output will
      be a CSV file with one line for each of the output molecules

    - --dice: use the DICE similarity metric instead of Tanimoto

    - --cosine: use the cosine similarity metric instead of Tanimoto

    - --fpColName=val: name to use for the column which stores
      fingerprints (in pickled format) in the output db table.
      Default is *AutoFragmentFP*

    - --minPath=val:  minimum path length to be included in
      fragment-based fingerprints. Default is *1*.

    - --maxPath=val:  maximum path length to be included in
      fragment-based fingerprints. Default is *7*.

    - --nBitsPerHash: number of bits to be set in the output
      fingerprint for each fragment. Default is *4*.

    - --discrim: use of path-based discriminators to hash bits.
      Default is *false*.

    - -V: include valence information in the fingerprints
      Default is *false*.

    - -H: include Hs in the fingerprint
      Default is *false*.

    - --useMACCS: use the public MACCS keys to do the fingerprinting
      (instead of a daylight-type fingerprint)


"""
if __name__ == '__main__':
  logging.info("This is MolSimilarity")
  rdktools_fp.FingerprintMols._usageDoc = _usageDoc
  details = rdktools_fp.FingerprintMols.ParseArgs()
  ScreenFromDetails(details)
