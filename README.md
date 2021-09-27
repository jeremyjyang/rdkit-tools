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

## Execution

Scaffold analysis:

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

Standardization:

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

Conformation generation:

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

