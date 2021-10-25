# `RDKIT-TOOLS`

Tools for use with RDKit.

* <https://rdkit.org>
* <https://www.rdkit.org/docs/Install.html>

## Dependencies

* RDKit Python package (possibly via conda).

```
$ conda create -n rdkit -c conda-forge rdkit ipykernel
$ conda activate rdkit
(rdkit) $ conda install -c conda-forge pyvis 
(rdkit) $ conda install -c conda-forge networkx=2.5 
```

## Scaffold analysis

```
$ conda activate rdkit
(rdkit) $ python3 -m rdktools.scaffold.App -h
usage: App.py [-h] [--i IFILE] [--o OFILE] [--o_html OFILE_HTML]
              [--scratchdir SCRATCHDIR] [--smicol SMICOL] [--namcol NAMCOL]
              [--idelim IDELIM] [--odelim ODELIM] [--iheader] [--oheader]
              [--brics] [-v]
              {bmscaf,scafnet,demobm,demonet,demonetvis}

RDKit scaffold analysis

positional arguments:
  {bmscaf,scafnet,demobm,demonet,demonetvis}
                        OPERATION

optional arguments:
  -h, --help            show this help message and exit
  --i IFILE             input file, TSV or SDF
  --o OFILE             output file, TSV|SDF
  --o_html OFILE_HTML   output file, HTML
  --scratchdir SCRATCHDIR
  --smicol SMICOL       SMILES column from TSV (counting from 0)
  --namcol NAMCOL       name column from TSV (counting from 0)
  --idelim IDELIM       delim for input TSV
  --odelim ODELIM       delim for output TSV
  --iheader             input TSV has header
  --oheader             output TSV has header
  --brics               BRICS fragmentation rules (Degen, 2008)
  -v, --verbose
```

## Standardization

```
(rdkit) $ python3 -m rdktools.standard.App
usage: App.py [-h] [--i IFILE] [--o OFILE] [--norms {default,unm}]
              [--i_norms IFILE_NORMS] [--remove_isomerism] [-v]
              {standardize,list_norms,show_params,demo}
App.py: error: the following arguments are required: op
lengua$ python3 -m rdktools.standard.App -h
usage: App.py [-h] [--i IFILE] [--o OFILE] [--norms {default,unm}]
              [--i_norms IFILE_NORMS] [--remove_isomerism] [-v]
              {standardize,list_norms,show_params,demo}

RDKit chemical standardizer

positional arguments:
  {standardize,list_norms,show_params,demo}
                        operation

optional arguments:
  -h, --help            show this help message and exit
  --i IFILE             input file, SMI or SDF
  --o OFILE             output file, SMI or SDF
  --norms {default,unm}
                        normalizations
  --i_norms IFILE_NORMS
                        input normalizations file, format: SMIRKS<space>NAME
  --remove_isomerism    if true, output SMILES isomerism removed
  -v, --verbose
```

## Conformation generation:

```
(rdkit) $ python3 -m rdktools.conform.App -h
usage: App.py [-h] [--i IFILE] [--o OFILE] [--ff {UFF,MMFF}] [--optiters OPTITERS]
              [--nconf NCONF] [--etol ETOL] [--title_in_header] [-v]

RDKit Conformer Generation

optional arguments:
  -h, --help           show this help message and exit
  --i IFILE            input file, SMI or SDF
  --o OFILE            output SDF with 3D
  --ff {UFF,MMFF}      force-field
  --optiters OPTITERS  optimizer iterations per conf
  --nconf NCONF        # confs per mol
  --etol ETOL          energy tolerance
  --title_in_header    title line in header
  -v, --verbose

Based on distance geometry method by Blaney et al.
```

## Fingerprint-based analytics

```
(rdkit) $ python3 -m rdktools.fp.App MolSimilarity -h
usage: App.py [-h] [--i IFILE] [--o OFILE] [--useHs] [--useValence]
              [--dbName DBNAME] [--tableName TABLENAME] [--minSize MINSIZE]
              [--maxSize MAXSIZE] [--density DENSITY] [--outTable OUTTABLE]
              [--outDbName OUTDBNAME] [--fpColName FPCOLNAME]
              [--minPath MINPATH] [--maxPath MAXPATH]
              [--nBitsPerHash NBITSPERHASH] [--discrim]
              [--smilesColumn SMILESCOLUMN] [--molPkl MOLPKL]
              [--input_format {SMILES,SD}] [--idColumn IDCOLUMN]
              [--maxMols MAXMOLS] [--fpAlgo {RDKIT,MACCS,MORGAN}]
              [--morgan_nbits MORGAN_NBITS] [--morgan_radius MORGAN_RADIUS]
              [--keepTable] [--smilesTable SMILESTABLE] [--topN TOPN]
              [--thresh THRESH] [--querySmiles QUERYSMILES]
              [--metric {ALLBIT,ASYMMETRIC,DICE,COSINE,KULCZYNSKI,MCCONNAUGHEY,ONBIT,RUSSEL,SOKAL,TANIMOTO,TVERSKY}]
              [--tversky_alpha TVERSKY_ALPHA] [--tversky_beta TVERSKY_BETA]
              [--clusterAlgo {WARD,SLINK,CLINK,UPGMA,BUTINA}]
              [--actTable ACTTABLE] [--actName ACTNAME]
              [--reportFreq REPORTFREQ] [-v]
              {FingerprintMols,MolSimilarity,ClusterMols}

RDKit fingerprint-based analytics

positional arguments:
  {FingerprintMols,MolSimilarity,ClusterMols}
                        OPERATION

optional arguments:
  -h, --help            show this help message and exit
  --i IFILE             Input file; if provided and no tableName is specified,
                        data will be read from the input file. Text files
                        delimited with either commas (extension .csv) or tabs
                        (extension .txt) are supported.
  --o OFILE             Name of the output file (output will be a pickle file
                        with one label,fingerprint entry for each molecule).
  --useHs               Include Hs in the fingerprint Default is *false*.
  --useValence          Include valence information in the fingerprints
                        Default is *false*.
  --dbName DBNAME       Name of the database from which to pull input molecule
                        information. If output is going to a database, this
                        will also be used for that unless the --outDbName
                        option is used.
  --tableName TABLENAME
                        Name of the database table from which to pull input
                        molecule information
  --minSize MINSIZE     Minimum size of the fingerprints to be generated
                        (limits the amount of folding that happens).
  --maxSize MAXSIZE     Base size of the fingerprints to be generated.
  --density DENSITY     Target bit density in the fingerprint. The fingerprint
                        will be folded until this density is reached.
  --outTable OUTTABLE   name of the output db table used to store
                        fingerprints. If this table already exists, it will be
                        replaced.
  --outDbName OUTDBNAME
                        name of output database, if it's being used. Defaults
                        to be the same as the input db.
  --fpColName FPCOLNAME
                        name to use for the column which stores fingerprints
                        (in pickled format) in the output db table.
  --minPath MINPATH     Minimum path length to be included in fragment-based
                        fingerprints.
  --maxPath MAXPATH     Maximum path length to be included in fragment-based
                        fingerprints.
  --nBitsPerHash NBITSPERHASH
                        Number of bits to be set in the output fingerprint for
                        each fragment.
  --discrim             Use of path-based discriminators to hash bits.
  --smilesColumn SMILESCOLUMN
                        Name of the SMILES column in the input database.
  --molPkl MOLPKL
  --input_format {SMILES,SD}
                        SMILES table or SDF file.
  --idColumn IDCOLUMN   Name of the id column in the input database. Defaults
                        to the first column for dbs.
  --maxMols MAXMOLS     Maximum number of molecules to be fingerprinted.
  --fpAlgo {RDKIT,MACCS,MORGAN}
                        RDKIT = Daylight path-based; MACCS = MDL MACCS 166
                        keys
  --morgan_nbits MORGAN_NBITS
  --morgan_radius MORGAN_RADIUS
  --keepTable
  --smilesTable SMILESTABLE
  --topN TOPN           Top N similar; precedence over threshold.
  --thresh THRESH       Similarity threshold.
  --querySmiles QUERYSMILES
                        Query smiles for similarity screening.
  --metric {ALLBIT,ASYMMETRIC,DICE,COSINE,KULCZYNSKI,MCCONNAUGHEY,ONBIT,RUSSEL,SOKAL,TANIMOTO,TVERSKY}
                        Similarity algorithm
  --tversky_alpha TVERSKY_ALPHA
                        Tversky alpha parameter, weights query molecule
                        features
  --tversky_beta TVERSKY_BETA
                        Tversky beta parameter, weights target molecule
                        features
  --clusterAlgo {WARD,SLINK,CLINK,UPGMA,BUTINA}
                        Clustering algorithm: WARD = Ward's minimum variance;
                        SLINK = single-linkage clustering algorithm; CLINK =
                        complete-linkage clustering algorithm; UPGMA = group-
                        average clustering algorithm; BUTINA = Butina JCICS 39
                        747-750 (1999)
  --actTable ACTTABLE   name of table containing activity values (used to
                        color points in the cluster tree).
  --actName ACTNAME     name of column with activities in the activity table.
                        The values in this column should either be integers or
                        convertible into integers.
  --reportFreq REPORTFREQ
  -v, --verbose

This app employs custom, updated versions of RDKit FingerprintMols.py,
MolSimilarity.py, ClusterMols.py, with enhanced command-line functionality for
molecular fingerprint-based analytics.
```

Examples:

```
(rdkit) $ python3 -m rdktools.fp.App FingerprintMols --i drugcentral.smi --smilesColumn "smiles" --idColumn "name" --fpAlgo MORGAN --morgan_nbits 2048
```

```
(rdkit) $ python3 -m rdktools.fp.App MolSimilarity --i drugcentral.smi --smilesColumn "smiles" --idColumn "name" --querySmiles "NCCc1ccc(O)c(O)c1 dopamine" --fpAlgo MORGAN --morgan_nbits 512 --metric TVERSKY --tversky_alpha 0.8 --tversky_beta 0.2
```

```
(rdkit) $ python3 -m rdktools.fp.App ClusterMols --i drugcentral.smi --smilesColumn "smiles" --idColumn "name" --fpAlgo MORGAN --morgan_nbits 512 --clusterAlgo BUTINA --metric TANIMOTO
```

## Molecular properties

```
$ python3 -m rdktools.properties.App -h
usage: App.py [-h] --i IFILE [--o OFILE] [--iheader] [--oheader] [--kekulize]
              [--sanitize] [--delim DELIM] [--smilesColumn SMILESCOLUMN]
              [--nameColumn NAMECOLUMN] [-v]
              {descriptors,descriptors3d,lipinski,logp,estate,freesasa,demo}

RDKit molecular properties utility

positional arguments:
  {descriptors,descriptors3d,lipinski,logp,estate,freesasa,demo}
                        OPERATION

optional arguments:
  -h, --help            show this help message and exit
  --i IFILE             input molecule file
  --o OFILE             output file with data (TSV)
  --iheader             input file has header line
  --oheader             include TSV header line with smiles output
  --kekulize            Kekulize
  --sanitize            Sanitize
  --delim DELIM         SMILES/TSV delimiter
  --smilesColumn SMILESCOLUMN
                        input SMILES column
  --nameColumn NAMECOLUMN
                        input name column
  -v, --verbose
```


