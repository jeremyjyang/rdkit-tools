#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.MurckoScaffold.html
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.rdScaffoldNetwork.html
rdScaffoldNetwork available RDKit 2020.03.1+.
"""
import argparse
import logging
import os
import sys
import time

import matplotlib as mpl
import pyvis
import rdkit

from .. import scaffold, util

# allowable OPS
BMSCAF = "bmscaf"
HIERARCHICAL_SCAFFOLDS = "hier_scaf"
SCAFNET = "scafnet"
SCAFNET_RINGS = "scafnet_rings"
DEMOBM = "demobm"
DEMONET_IMG = "demonet_img"
DEMONET_HTML = "demonet_html"


def parse_args(parser: argparse.ArgumentParser):
    # mapping of operation to help string
    OPS = {
        BMSCAF: "Generate scaffolds using Bemis-Murcko clustering",
        SCAFNET: "Generate a scaffold network using the given SMILES",
        SCAFNET_RINGS: "Generate a scaffold network using the given SMILES, with output containing unique ringsystems only",
        HIERARCHICAL_SCAFFOLDS: "Derive hierarchical scaffolds for the given SMILES",
        DEMOBM: "Demo scaffold generated using Bemis-Murcko clustering",
        DEMONET_IMG: "Demo generating scaffold network image",
        DEMONET_HTML: "Demo generating interactive scaffold network using pyvis",
    }
    subparsers = parser.add_subparsers(dest="op", help="operation")
    # create a subparser for each operation
    parsers = []
    for op, op_help in OPS.items():
        parsers.append(
            subparsers.add_parser(
                op, help=op_help, formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
        )
    for sub_parser in parsers:
        # sub_parser.prog == "App.py <PROG_NAME>"
        prog_name = sub_parser.prog.split()[1]
        sub_parser.add_argument(
            "--log_fname",
            help="File to save logs to. If not given will log to stdout.",
            default=None,
        )
        sub_parser.add_argument(
            "-v", "--verbose", action="count", default=0, help="verbosity of logging"
        )
        if prog_name in [DEMONET_IMG, DEMONET_HTML, SCAFNET, HIERARCHICAL_SCAFFOLDS]:
            sub_parser.add_argument(
                "--scratchdir",
                default=os.path.join(os.environ["HOME"], "tmp", "rdktools"),
                help="Directory where demo and other temporary files will be stored. Note that the dir and its contents will be created but not destroyed",
            )
        if prog_name in [BMSCAF, SCAFNET]:
            sub_parser.add_argument(
                "--o_png",
                dest="ofile_png",
                help="visualization output file, PNG",
            )
            default_per_row = 8 if prog_name == BMSCAF else 4
            sub_parser.add_argument(
                "--mols_per_row",
                type=int,
                default=default_per_row,
                help="Mols per row in output PNG",
            )
        if prog_name in [BMSCAF, SCAFNET, SCAFNET_RINGS, HIERARCHICAL_SCAFFOLDS]:
            sub_parser.add_argument(
                "--i", dest="ifile", required=True, help="input file, SMI or SDF"
            )
            if prog_name == HIERARCHICAL_SCAFFOLDS:
                # two output files
                sub_parser.add_argument(
                    "--o_scaf",
                    dest="o_scaf",
                    required=True,
                    default=argparse.SUPPRESS,
                    help="output file with detected scaffolds and their ids",
                )
                sub_parser.add_argument(
                    "--o_mol",
                    dest="o_mol",
                    required=True,
                    default=argparse.SUPPRESS,
                    help="output file with molecules and their ids",
                )
                sub_parser.add_argument(
                    "--o_mol2scaf",
                    dest="o_mol2scaf",
                    required=True,
                    default=argparse.SUPPRESS,
                    help="output file mapping mol ids to associated scaffold ids and their relative depth",
                )
            else:
                sub_parser.add_argument(
                    "--o",
                    dest="ofile",
                    help="output file, SMI or SDF. Will use stdout if not specified",
                )
            sub_parser.add_argument(
                "--smiles_column",
                type=int,
                default=0,
                help="(integer) column where SMILES are located (for input SMI file)",
            )
            sub_parser.add_argument(
                "--name_column",
                type=int,
                default=1,
                help="(integer) column where molecule names are located (for input SMI file)",
            )
            sub_parser.add_argument(
                "--idelim",
                default="\t",
                help="delim for input SMI/TSV file (default is tab)",
            )
            sub_parser.add_argument(
                "--odelim",
                default="\t",
                help="delim for output SMI/TSV file (default is tab)",
            )
            sub_parser.add_argument(
                "--iheader",
                action="store_true",
                help="input SMILES/TSV has header line",
            )
            sub_parser.add_argument(
                "--oheader", action="store_true", help="output TSV has header"
            )
        if prog_name in [SCAFNET, SCAFNET_RINGS, HIERARCHICAL_SCAFFOLDS]:
            sub_parser.add_argument(
                "--brics",
                action="store_true",
                help="use BRICS fragmentation rules (Degen, 2008)",
            )
        if prog_name == SCAFNET:
            sub_parser.add_argument(
                "--o_html",
                dest="ofile_html",
                help="location to store HTML visualization output file. Will show in browser once execution is complete",
            )
            sub_parser.add_argument(
                "--scafname", default="RDKit Scaffold Analysis", help="title for output"
            )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RDKit scaffold analysis", epilog="")
    args = parse_args(parser)

    level = logging.WARNING  # default
    if args.verbose == 1:
        level = logging.INFO
    elif args.verbose >= 2:
        level = logging.DEBUG

    logging.basicConfig(
        filename=args.log_fname,
        filemode="a",
        format="%(levelname)s:%(message)s",
        level=level,
    )

    logging.info(f"RDKit version: {rdkit.rdBase.rdkitVersion}")
    logging.info(f"MatplotLib version: {mpl.__version__}")
    logging.info(f"Pyvis version: {pyvis.__version__}")

    t0 = time.time()

    if hasattr(args, "scratchdir"):
        os.makedirs(args.scratchdir, exist_ok=True)

    if args.op == DEMOBM:
        scaffold.DemoBM()
        sys.exit()

    elif args.op == DEMONET_IMG:
        scaffold.DemoNetImg(args.scratchdir)
        sys.exit()

    elif args.op == DEMONET_HTML:
        logging.debug(f"scratchdir: {args.scratchdir}")
        scaffold.DemoNetHtml(args.scratchdir)
        sys.exit()

    molReader = util.File2Molreader(
        args.ifile, args.idelim, args.smiles_column, args.name_column, args.iheader
    )
    if args.op == BMSCAF:
        molWriter = util.File2Molwriter(args.ofile, args.odelim, args.oheader)
        mols = util.ReadMols(molReader)
        scafmols, legends = scaffold.Mols2BMScaffolds(mols, molWriter)
        if args.ofile_png:
            img = rdkit.Chem.Draw.MolsToGridImage(
                scafmols,
                molsPerRow=args.mols_per_row,
                legends=legends,
            )
            img.save(args.ofile_png, format="PNG")
    elif args.op == SCAFNET:
        ofile = args.ofile if args.ofile else sys.stdout
        scafnet, _ = scaffold.Mols2ScafNet(
            molReader, args.brics, ofile, args.odelim, args.oheader
        )
        if args.ofile_png:
            scaffold.Scafnet2Img(scafnet, args.ofile_png, args.mols_per_row)
        if args.ofile_html:
            scaffold.Scafnet2Html(
                scafnet, args.scafname, args.scratchdir, args.ofile_html
            )
    elif args.op == SCAFNET_RINGS:
        molWriter = util.File2Molwriter(args.ofile, args.odelim, args.oheader)
        n_mol = 0
        for mol in molReader:
            n_mol += 1
            molname = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Mol_{n_mol:d}"
            scafnet = scaffold.Mols2ScafNet([mol], args.brics, None)
            rings = scaffold.ScafNet2Rings(scafnet, molname, molWriter)
        logging.info(f"Mols processed: {n_mol}")
    elif args.op == HIERARCHICAL_SCAFFOLDS:
        scaffold.HierarchicalScaffolds(
            molReader,
            args.o_mol,
            args.o_scaf,
            args.o_mol2scaf,
            args.odelim,
            args.oheader,
        )
    else:
        parser.error(f"Unsupported operation: {args.op}")

    logging.info(
        f"""Elapsed: {time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))}"""
    )
