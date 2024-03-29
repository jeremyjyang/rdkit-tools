# `RDKIT-TOOLS`

Tools for use with RDKit. Motivated and intended for use with
[CFDE](https://nih-cfde.org/) and CFChemDb, developed by the IDG-CFDE team.

See also:

* [CFChemDb](https://github.com/unmtransinfo/CFChemDb) (repository)
* [CFChemDb_UI](https://github.com/jeremyjyang/CFChemDb_UI) (repository)
* [rdktools](https://pypi.org/project/rdktools/) (Pypi package)
* [CFDE: Common Fund Data Ecosystem](https://nih-cfde.org/)

RDKit:

* <https://rdkit.org>
* <https://www.rdkit.org/docs/Install.html>

## Dependencies

* RDKit Python package (via conda recommended).

```
$ conda create -n rdktools -c conda-forge rdkit ipykernel
$ conda activate rdktools
(rdktools) $ conda install -c conda-forge pyvis 
(rdktools) $ conda install -c conda-forge networkx=2.5 
```

See also: [conda/environment.yml](conda/environment.yml)

## Contents

* [Formats](#Formats) - chemical file format conversion
* [Depictions](#Depictions) - 2D molecular depictions
* [Standardization](#Standardization) - molecular standardization 
* [Fingerprints](#Fingerprints) - molecular path and pattern based binary feature vectors, similarity, and clustering tools
* [Conformations](#Conformations) - distance geometry based 3D conformation generation
* [Properties](#Properties) - molecular property calculation: Lipinsky, Wildman-Crippen LogP, Kier-Hall electrotopological descriptors, solvent accessible surface area (SASA), and more.
* [Scaffolds](#Scaffolds) - Bemis-Murcko and BRICS scaffold analysis, rdScaffoldNetworks.
* [SMARTS](#SMARTS) - molecular pattern matching (subgraph isomorphism)
* [Reactions](#Reactions) - Reaction SMILES, SMARTS, and SMIRKS based reaction analytics
* [util.sklearn](#util.sklearn) - Scikit-learn utilities for processing molecular fingerprints and other feature vectors.


## Formats

```
(rdktools) $ python3 -m rdktools.formats.App -h
usage: App.py [-h] [--i IFILE] [--o OFILE] [--kekulize] [--sanitize] [--header]
              [--delim DELIM] [--smilesColumn SMILESCOLUMN] [--nameColumn NAMECOLUMN]
              [-v]
              {mdl2smi,mdl2tsv,smi2mdl,smiclean,mdlclean,mol2inchi,mol2inchikey,demo}

RDKit chemical format utility

positional arguments:
  {mdl2smi,mdl2tsv,smi2mdl,smiclean,mdlclean,mol2inchi,mol2inchikey,demo}
                        operation

optional arguments:
  -h, --help            show this help message and exit
  --i IFILE             input file (SMILES/TSV or SDF)
  --o OFILE             output file (specify '-' for stdout)
  --kekulize            Kekulize
  --sanitize            Sanitize
  --header              input SMILES/TSV file has header line
  --delim DELIM         delimiter for SMILES/TSV
  --smilesColumn SMILESCOLUMN
                        input SMILES column
  --nameColumn NAMECOLUMN
                        input name column
  -v, --verbose
```

## Depictions

```
(rdktools) $ python3 -m rdktools.depict.App -h
usage: App.py [-h] [--i IFILE] [--ifmt {AUTO,SMI,MDL}] [--ofmt {PNG,JPEG,PDF}]
              [--smilesColumn SMILESCOLUMN] [--nameColumn NAMECOLUMN] [--header]
              [--delim DELIM] [--height HEIGHT] [--width WIDTH] [--kekulize]
              [--wedgebonds] [--pdf_title PDF_TITLE] [--batch_dir BATCH_DIR]
              [--batch_prefix BATCH_PREFIX] [--o OFILE] [-v]
              {single,batch,pdf,demo,demo2}

RDKit molecule depiction utility

positional arguments:
  {single,batch,pdf,demo,demo2}
                        OPERATION

optional arguments:
  -h, --help            show this help message and exit
  --i IFILE             input molecule file
  --ifmt {AUTO,SMI,MDL}
                        input file format
  --ofmt {PNG,JPEG,PDF}
                        output file format
  --smilesColumn SMILESCOLUMN
  --nameColumn NAMECOLUMN
  --header              SMILES/TSV file has header
  --delim DELIM         SMILES/TSV field delimiter
  --height HEIGHT       height of image
  --width WIDTH         width of image
  --kekulize            display Kekule form
  --wedgebonds          stereo wedge bonds
  --pdf_title PDF_TITLE
                        PDF doc title
  --batch_dir BATCH_DIR
                        destination for batch files
  --batch_prefix BATCH_PREFIX
                        prefix for batch files
  --o OFILE             output file
  -v, --verbose

Modes: single = one image; batch = multiple images; pdf = multi-page
```

```
python3 -m rdktools.depict.App single -height 500 --width 600 --i valium.smiles --o valium.png 
```

<img src="data/valium.png" height="400">

## Scaffolds
Tools for processing SMILES inputs and performing scaffold analysis. 
```
(rdktools) $ python3 -m rdktools.scaffold.App -h
usage: App.py [-h] {bmscaf,scafnet,scafnet_rings,demobm,demonet_img,demonet_html} ...

RDKit scaffold analysis

positional arguments:
  {bmscaf,scafnet,scafnet_rings,demobm,demonet_img,demonet_html}
                        operation
    bmscaf              Generate scaffolds using Bemis-Murcko clustering
    scafnet             Generate a scaffold network using the given SMILES
    scafnet_rings       Generate a scaffold network using the given SMILES, with output
                        containing unique ringsystems only
    demobm              Demo scaffold generated using Bemis-Murcko clustering
    demonet_img         Demo generating scaffold network image
    demonet_html        Demo generating interactive scaffold network using pyvis

options:
  -h, --help            show this help message and exit
```
Additional information for a specific operation can be found by using the `-h` flag after providing the operation. For example:
```
(rdktools) $ python3 -m rdktools.scaffold.App bmscaf -h
usage: App.py bmscaf [-h] [--log_fname LOG_FNAME] [-v] [--o_png OFILE_PNG]
                     [--mols_per_row MOLS_PER_ROW] --i IFILE [--o OFILE]
                     [--smiles_column SMILES_COLUMN] [--name_column NAME_COLUMN]
                     [--idelim IDELIM] [--odelim ODELIM] [--iheader] [--oheader]

options:
  -h, --help            show this help message and exit
  --log_fname LOG_FNAME
                        File to save logs to. If not given will log to stdout. (default:
                        None)
  -v, --verbose         verbosity of logging (default: 0)
  --o_png OFILE_PNG     visualization output file, PNG (default: None)
  --mols_per_row MOLS_PER_ROW
                        Mols per row in output PNG (default: 8)
  --i IFILE             input file, SMI or SDF (default: None)
  --o OFILE             output file, SMI or SDF. Will use stdout if not specified (default:
                        None)
  --smiles_column SMILES_COLUMN
                        (integer) column where SMILES are located (for SMI file) (default: 0)
  --name_column NAME_COLUMN
                        (integer) column where molecule names are located (for SMI file)
                        (default: 1)
  --idelim IDELIM       delim for input SMI/TSV file (default is tab) (default: )
  --odelim ODELIM       delim for output SMI/TSV file (default is tab) (default: )
  --iheader             input SMILES/TSV has header line (default: False)
  --oheader             output TSV has header (default: False)
```

### Examples
1) Viewing scaffolds generated from SMILES file using Bemis-Murcko clustering. Labels for each scaffold correspond to the input SMILES.
```
python3 -m rdktools.scaffold.App bmscaf --i data/drugs.smiles --o_png data/drugs_bmscaf.png
```

<img src="data/drugs_bmscaf.png" height="400"> <br />

2) Viewing scaffold network generated from SMILES file. In addition to the HTML file, the networks nodes and edges are stored in `mcs_net.tsv`.

```
python3 -m rdktools.scaffold.App scafnet --i data/mcs_example.smiles --o mcs_net.tsv --o_html mcs_net.html
```

Screenshot of output network (will be interactive if ran locally)
<img src="data/scafnet_out.png" height="400"> <br />  

## Standardization

```
(rdktools) $ python3 -m rdktools.standard.App -h
usage: App.py [-h] [--i IFILE] [--o OFILE] [--delim DELIM] [--smilesColumn SMILESCOLUMN]
              [--nameColumn NAMECOLUMN] [--nameSDField NAMESDFIELD] [--header]
              [--sanitize] [--kekuleSmiles] [--normset {DEFAULT,UNM}]
              [--i_normset IFILE_NORMSET] [--isomericSmiles] [--metalRemove]
              [--largestFragment] [--neutralize] [-v]
              {standardize,canonicalize,saltremove,list_norms,show_params,demo}

RDKit chemical standardizer

positional arguments:
  {standardize,canonicalize,saltremove,list_norms,show_params,demo}
                        OPERATION

options:
  -h, --help            show this help message and exit
  --i IFILE             input file, SMI or SDF
  --o OFILE             output file, SMI or SDF
  --delim DELIM         SMILES/TSV delimiter
  --smilesColumn SMILESCOLUMN
  --nameColumn NAMECOLUMN
  --nameSDField NAMESDFIELD
                        SD field to use as name
  --header              SMILES/TSV has header line
  --sanitize            Sanitize molecules as read.
  --kekuleSmiles        Kekule SMILES output.
  --normset {DEFAULT,UNM}
                        normalization sets
  --i_normset IFILE_NORMSET
                        input normalizations file, format: SMIRKS<space>NAME
  --isomericSmiles      If false, output SMILES isomerism removed
  --metalRemove         Remove disconnected metals like salts charges (use with
                        saltremove).
  --largestFragment     Remove non-largest fragments (use with saltremove).
  --neutralize          Neutralize charges (use with saltremove).
  -v, --verbose

```

For documentation on RDKit Molecular Sanitization,
see [The RDKit Book](https://www.rdkit.org/docs/RDKit_Book.html). Briefly:

  The idea is to generate useful computed properties (like hybridization, ring
  membership, etc.) for the rest of the code and to ensure that the molecules are
  "reasonable": that they can be represented with octet-complete Lewis dot
  structures.

## Conformations

```
(rdktools) $ python3 -m rdktools.conform.App -h
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

## Fingerprints

By default, RDKit and Morgan fingerprints are generated length 2048 bits,
by the following methods:

**RDKit path-based, Daylight-like**:
```
Chem.RDKFingerprint(mol, minPath=1, maxPath=7, fpSize=2048, nBitsPerHash=2, useHs=False, minSize=2048)
```

**Morgan ECFP-like**:
```
AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
```

```
(rdktools) $ python3 -m rdktools.fp.App -h

usage: App.py [-h] [--i IFILE] [--iheader] [--o OFILE] [--output_as_dataframe]
              [--output_as_tsv] [--useHs] [--useValence] [--dbName DBNAME]
              [--tableName TABLENAME] [--minSize MINSIZE] [--maxSize MAXSIZE]
              [--density DENSITY] [--outTable OUTTABLE] [--outDbName OUTDBNAME]
              [--fpColName FPCOLNAME] [--minPath MINPATH] [--maxPath MAXPATH]
              [--nBitsPerHash NBITSPERHASH] [--discrim] [--smilesColumn SMILESCOLUMN]
              [--molPkl MOLPKL] [--input_format {SMILES,SD}] [--idColumn IDCOLUMN]
              [--maxMols MAXMOLS] [--fpAlgo {RDKIT,MACCS,MORGAN}]
              [--morgan_nbits MORGAN_NBITS] [--morgan_radius MORGAN_RADIUS]
              [--replaceTable] [--smilesTable SMILESTABLE] [--topN TOPN]
              [--thresh THRESH] [--querySmiles QUERYSMILES]
              [--metric {ALLBIT,ASYMMETRIC,DICE,COSINE,KULCZYNSKI,MCCONNAUGHEY,ONBIT,RUSSEL,SOKAL,TANIMOTO,TVERSKY}]
              [--tversky_alpha TVERSKY_ALPHA] [--tversky_beta TVERSKY_BETA]
              [--clusterAlgo {WARD,SLINK,CLINK,UPGMA,BUTINA}] [--actTable ACTTABLE]
              [--actName ACTNAME] [--reportFreq REPORTFREQ] [--showVis] [-v]
              {FingerprintMols,MolSimilarity,ClusterMols}

RDKit fingerprint-based analytics

positional arguments:
  {FingerprintMols,MolSimilarity,ClusterMols}
                        OPERATION

optional arguments:
  -h, --help            show this help message and exit
  --i IFILE             input file; if provided and no tableName is specified, data will
                        be read from the input file. Text files delimited with either
                        commas (extension .csv) or tabs (extension .txt) are supported.
  --iheader             input file has header line
  --o OFILE             output file (pickle file with one label,fingerprint entry for
                        each molecule).
  --output_as_dataframe
                        Output FPs as Pandas dataframe (pickled) with names as index,
                        columns as feature names, if available.
  --output_as_tsv       Output FPs as TSV with names as index, columns as feature names,
                        if available.
  --useHs               include Hs in the fingerprint Default is *false*.
  --useValence          include valence information in the fingerprints Default is
                        *false*.
  --dbName DBNAME       name of the database from which to pull input molecule
                        information. If output is going to a database, this will also be
                        used for that unless the --outDbName option is used.
  --tableName TABLENAME
                        name of the database table from which to pull input molecule
                        information
  --minSize MINSIZE     minimum size of the fingerprints to be generated (limits the
                        amount of folding that happens) [64].
  --maxSize MAXSIZE     base size of the fingerprints to be generated [2048].
  --density DENSITY     target bit density in the fingerprint. The fingerprint will be
                        folded until this density is reached [0.3].
  --outTable OUTTABLE   name of the output db table used to store fingerprints. If this
                        table already exists, it will be replaced.
  --outDbName OUTDBNAME
                        name of output database, if it's being used. Defaults to be the
                        same as the input db.
  --fpColName FPCOLNAME
                        name to use for the column which stores fingerprints (in pickled
                        format) in the output db table [AutoFragmentFP].
  --minPath MINPATH     minimum path length to be included in fragment-based
                        fingerprints [1].
  --maxPath MAXPATH     maximum path length to be included in fragment-based
                        fingerprints [7].
  --nBitsPerHash NBITSPERHASH
                        number of bits to be set in the output fingerprint for each
                        fragment [2].
  --discrim             use of path-based discriminators to hash bits.
  --smilesColumn SMILESCOLUMN
                        name of the SMILES column in the input database [#SMILES].
  --molPkl MOLPKL
  --input_format {SMILES,SD}
                        SMILES table or SDF file [{DEFAULTS['input_format']}].
  --idColumn IDCOLUMN, --nameColumn IDCOLUMN
                        name of the id column in the input database. Defaults to the
                        first column for dbs [Name].
  --maxMols MAXMOLS     maximum number of molecules to be fingerprinted.
  --fpAlgo {RDKIT,MACCS,MORGAN}
                        RDKIT = Daylight path-based; MACCS = MDL MACCS 166 keys [RDKIT]
  --morgan_nbits MORGAN_NBITS
                        [1024]
  --morgan_radius MORGAN_RADIUS
                        [2]
  --replaceTable
  --smilesTable SMILESTABLE
                        name of database table which contains SMILES for the input
                        fingerprints. If provided with --smilesName, output will contain
                        SMILES data.
  --topN TOPN           top N similar; precedence over threshold [12].
  --thresh THRESH       similarity threshold.
  --querySmiles QUERYSMILES
                        query smiles for similarity screening.
  --metric {ALLBIT,ASYMMETRIC,DICE,COSINE,KULCZYNSKI,MCCONNAUGHEY,ONBIT,RUSSEL,SOKAL,TANIMOTO,TVERSKY}
                        similarity algorithm [TANIMOTO]
  --tversky_alpha TVERSKY_ALPHA
                        Tversky alpha parameter, weights query molecule features [0.8]
  --tversky_beta TVERSKY_BETA
                        Tversky beta parameter, weights target molecule features [0.2]
  --clusterAlgo {WARD,SLINK,CLINK,UPGMA,BUTINA}
                        clustering algorithm: WARD = Ward's minimum variance; SLINK =
                        single-linkage clustering algorithm; CLINK = complete-linkage
                        clustering algorithm; UPGMA = group-average clustering
                        algorithm; BUTINA = Butina JCICS 39 747-750 (1999) [WARD]
  --actTable ACTTABLE   name of table containing activity values (used to color points
                        in the cluster tree).
  --actName ACTNAME     name of column with activities in the activity table. The values
                        in this column should either be integers or convertible into
                        integers.
  --reportFreq REPORTFREQ
                        [100]
  --showVis             show visualization if available.
  -v, --verbose

This app employs custom, updated versions of RDKit FingerprintMols.py, MolSimilarity.py,
ClusterMols.py, with enhanced command-line functionality for molecular fingerprint-based
analytics.
```

Examples:

```
(rdktools) $ python3 -m rdktools.fp.App FingerprintMols --i drugcentral.smiles --smilesColumn "smiles" --idColumn "name" --fpAlgo MORGAN --morgan_nbits 2048 --output_as_tsv --o drugcentral_morganfp.tsv
```

```
(rdktools) $ python3 -m rdktools.fp.App MolSimilarity --i drugcentral.smiles --smilesColumn "smiles" --idColumn "name" --querySmiles "NCCc1ccc(O)c(O)c1 dopamine" --fpAlgo MORGAN --morgan_nbits 512 --metric TVERSKY --tversky_alpha 0.8 --tversky_beta 0.2
```

```
(rdktools) $ python3 -m rdktools.fp.App ClusterMols --i drugcentral.smiles --smilesColumn "smiles" --idColumn "name" --fpAlgo MORGAN --morgan_nbits 512 --clusterAlgo BUTINA --metric TANIMOTO
```

## SMARTS
Tools for processing SMILES inputs based on SMARTS queries. Includes options for filtering and counting matches. Outputs are in TSV format.
```
(rdktools) $ python3 -m rdktools.smarts.App -h
usage: App.py [-h]
              {matchCounts,matchFilter,matchCountsMulti,matchFilterMulti,filterPAINS,demo}
              ...

RDKit SMARTS utility

positional arguments:
  {matchCounts,matchFilter,matchCountsMulti,matchFilterMulti,filterPAINS,demo}
                        operation
    matchCounts         Count matches of a single SMARTS in each molecule
    matchFilter         Filter molecules that match a single SMARTS
    matchCountsMulti    Count matches of multiple SMARTS (from file) in each molecule
    matchFilterMulti    Filter molecules that match multiple SMARTS (from file)
    filterPAINS         Filter molecules that match PAINS
    demo                Demo of matchCounts

options:
  -h, --help            show this help message and exit
```
Additional information for a specific operation can be found by using the `-h` flag after providing the operation. For example:
```
(rdktools) $ python3 -m rdktools.smarts.App matchCountsMulti -h

usage: App.py matchCountsMulti [-h] [--log_fname LOG_FNAME] [-v] --smartsfile SMARTSFILE
                               [--strict] [--usa] [--nonzero_rows] --i IFILE [--o OFILE]
                               [--delim DELIM] [--smiles_column SMILES_COLUMN]
                               [--name_column NAME_COLUMN] [--iheader] [--exclude_mol_props]

options:
  -h, --help            show this help message and exit
  --log_fname LOG_FNAME
                        File to save logs to. If not given will log to stdout. (default:
                        None)
  -v, --verbose         verbosity of logging (default: 0)
  --smartsfile SMARTSFILE
                        input SMARTS file (for multi-ops)
  --strict              raise error if any SMARTS cannot be parsed. If not set, will ignore
                        invalid SMARTS. (default: False)
  --usa                 unique set-of-atoms match counts (default: False)
  --nonzero_rows        only include rows (ie molecules) with at least one match in output
                        (default: False)
  --i IFILE             input file, SMI or SDF
  --o OFILE             output file, TSV. Will use stdout if not specified. (default: None)
  --delim DELIM         delimiter for SMILES/TSV (default is tab) (default: )
  --smiles_column SMILES_COLUMN
                        (integer) column where SMILES are located (for SMI file) (default: 0)
  --name_column NAME_COLUMN
                        (integer) column where molecule names are located (for SMI file)
                        (default: 1)
  --iheader             input SMILES/TSV has header line (default: False)
  --exclude_mol_props   exclude molecular properties present in input SMILES/SDF in output
                        (i.e., only include SMILES & Name properties) (default: False)
```
### Examples
1) Counting matches from `ursu_pains.sma` in SMILES from `hscaf_testset.smi`. Note that `hscaf_testset.smi` uses a seperator of space, hence the choice of `--delim`. The `-vvv` option gives verbose logs which will be saved to `out.log`. Output is saved to `mcs_out.tsv`.

```
python3 -m rdktools.smarts.App matchCountsMulti --i rdktools/smarts/data/hscaf_testset.smi --smartsfile rdktools/smarts/data/ursu_pains.sma --o mcs_out.tsv -vvv --log_fname out.log --delim " "
```

2) Similar to above, but using a csv SMILES file where the ordering of columns is different than the default and the file includes a header. Note that in this case the `--nonzero_rows` flag is also given, which means that only molecules which had at least one match will be included in the output.

```
python3 -m rdktools.smarts.App matchCountsMulti --i rdktools/smarts/data/mcs_demo2.csv --smartsfile rdktools/smarts/data/ursu_pains.sma --o mcs_out2.tsv --nonzero_rows -vvv --log_fname out2.log --delim "," --smiles_column 1 --name_column 0 --iheader
```

3) Testing if the molecule from `example.sdf` passes the PAINS filtering. The `--exclude_mol_props` flag is passed to exclude molecular properties that are present in `example.sdf` (e.g., IR.FREQUENCIES) from the output. 

```
python3 -m rdktools.smarts.App filterPAINS --i rdktools/smarts/data/example.sdf --o sdf_out.tsv -vvv --log_fname out3.log --exclude_mol_props
```

4) Filtering molecules from `mcs_demo2.csv` that match a given SMARTS query. Note that the SMARTS string is captured with single quotes (''). This is important because "!" can cause conflicts with the command line due to history expansion (see [here](https://stackoverflow.com/questions/33685239/in-a-string-makes-it-unusable-in-a-command-line-error-message-event-not-fou)).

```
python3 -m rdktools.smarts.App matchFilter --i rdktools/smarts/data/mcs_demo2.csv --smarts '[$([N;!H0]-[#6]);!$(N-C=[O,N,S])]' --delim "," --smiles_column 1 --name_column 0 --iheader -vv
```

## Reactions

```
$ python3 -m rdktools.reactions.App -h
usage: App.py [-h] [--i IFILES] [--o OFILE] [--output_mode {products,reactions}]
              [--o_depict OFILE_DEPICT] [--smirks SMIRKS] [--kekulize] [--sanitize]
              [--depict] [--header] [--delim DELIM] [--smilesColumn SMILESCOLUMN]
              [--nameColumn NAMECOLUMN] [-v]
              {enumerateLibrary,react,demo,demo2,demo3,demo4}

RDKit chemical reactions utility

positional arguments:
  {enumerateLibrary,react,demo,demo2,demo3,demo4}
                        OPERATION

optional arguments:
  -h, --help            show this help message and exit
  --i IFILES            input file[s] (SMILES/TSV or SDF)
  --o OFILE             output file (SMILES) [stdout]
  --output_mode {products,reactions}
                        products|reactions [products]
  --o_depict OFILE_DEPICT
                        output depiction file (PNG) [display]
  --smirks SMIRKS       SMIRKS reaction transform
  --kekulize            Kekulize
  --sanitize            Sanitize
  --depict              Depict (1st reaction or product only)
  --header              input SMILES/TSV file has header line
  --delim DELIM         delimiter for SMILES/TSV
  --smilesColumn SMILESCOLUMN
                        input SMILES column
  --nameColumn NAMECOLUMN
                        input name column
  -v, --verbose

For 'react' operation, reactants are specified as disconnected components of single
input molecule record. For 'enumerateLibrary', reactants for each role are specfied from
separate input files, ordered as in the SMIRKS.
```

```
python3 -m rdktools.reactions.App react --smirks '[O:2]=[C:1][OH].[N:3]>>[O:2][C:1][N:3]' --i reactants.smiles --nameColumn 0 --depict --o_depict reaction.png
```

<img src="data/reaction.png" height="300">

## Properties

```
(rdktools) $ python3 -m rdktools.properties.App -h
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

### util.sklearn

Scikit-learn utilities for processing molecular fingerprints and other feature vectors.

```
(rdktools) lengua$ python3 -m rdktools.util.sklearn.ClusterFingerprints -h
usage: ClusterFingerprints.py [-h] [--i IFILE] [--o OFILE] [--o_vis OFILE_VIS]
                              [--scratchdir SCRATCHDIR] [--idelim IDELIM]
                              [--odelim ODELIM]
                              [--affinity {euclidean,l1,l2,manhattan,cosine,precomputed}]
                              [--linkage {ward,complete,average,single}]
                              [--truncate_level TRUNCATE_LEVEL] [--iheader] [--oheader]
                              [--dendrogram_orientation {left,top,right,bottom}]
                              [--display] [-v]
                              {cluster,demo}

Hierarchical, agglomerative clustering by Scikit-learn

positional arguments:
  {cluster,demo}        OPERATION

optional arguments:
  -h, --help            show this help message and exit
  --i IFILE             input file, TSV
  --o OFILE             output file, TSV
  --o_vis OFILE_VIS     output file, PNG or HTML
  --scratchdir SCRATCHDIR
  --idelim IDELIM       delim for input TSV
  --odelim ODELIM       delim for output TSV
  --affinity {euclidean,l1,l2,manhattan,cosine,precomputed}
  --linkage {ward,complete,average,single}
  --truncate_level TRUNCATE_LEVEL
                        Level from root of hierarchy for clusters and dendrogram.
  --iheader             input TSV has header
  --oheader             output TSV has header
  --dendrogram_orientation {left,top,right,bottom}
  --display             display dendrogram
  -v, --verbose
```

```
(rdktools) $ python3 -m rdktools.util.sklearn.ClusterFingerprints cluster --i drugcentral_morganfp.tsv --truncate_level 5 --o_vis drugcentral_morganfp_ward-clusters_dendrogram.png
```

<img src="data/drugcentral_morganfp_ward-clusters_dendrogram.png" height="500">
