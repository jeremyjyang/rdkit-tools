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


def parse_args(parser: argparse.ArgumentParser):
    # mapping of operation to help string
    OPS = {
        "bmscaf": "",
        "scafnet": "",
        "scafnet_rings": "",
        "demobm": "",
        "demonet_img": "",
        "demonet_html": "",
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
        if prog_name in ["demonet_img", "demonet_html", "scafnet"]:
            sub_parser.add_argument(
                "--scratchdir",
                default=os.path.join(os.environ["HOME"], "tmp", "rdktools"),
            )
        if prog_name in ["bmscaf", "scafnet", "scafnet_rings"]:
            sub_parser.add_argument("--i", dest="ifile", help="input file, TSV or SDF")
            sub_parser.add_argument("--o", dest="ofile", help="output file, TSV or SDF")
        if prog_name in ["bmscaf", "scafnet"]:
            sub_parser.add_argument(
                "--o_png",
                dest="ofile_png",
                default="--",
                help="visualization output file, PNG",
            )
        if prog_name == "scafnet":
            sub_parser.add_argument(
                "--o_html",
                dest="ofile_html",
                default="--",
                help="visualization output file, HTML",
            )
            sub_parser.add_argument(
                "--scafname", default="RDKit Scaffold Analysis", help="title for output"
            )
        if prog_name in ["bmscaf", "scafnet", "scafnet_rings"]:
            sub_parser.add_argument(
                "--smilesColumn",
                type=int,
                default=0,
                help="SMILES column from TSV (counting from 0)",
            )
            sub_parser.add_argument(
                "--nameColumn",
                type=int,
                default=1,
                help="name column from TSV (counting from 0)",
            )
            sub_parser.add_argument(
                "--idelim", default="\t", help="delim for input TSV"
            )
            sub_parser.add_argument(
                "--odelim", default="\t", help="delim for output TSV"
            )
            sub_parser.add_argument(
                "--iheader", action="store_true", help="input TSV has header"
            )
            sub_parser.add_argument(
                "--oheader", action="store_true", help="output TSV has header"
            )
        if sub_parser.prog in ["scafnet", "scafnet_rings"]:
            sub_parser.add_argument(
                "--brics",
                action="store_true",
                help="BRICS fragmentation rules (Degen, 2008)",
            )
        sub_parser.add_argument(
            "--display", action="store_true", help="Display scafnet interactively."
        )
    args = parser.parse_args()
    # use given scratchdir as default for o_png and o_html
    if hasattr(args, "ofile_png") and not args.ofile_png:
        args.ofile_png = os.path.join(args.scratchdir, "rdktools_scafnet.png")
    if hasattr(args, "ofile_html") and not args.ofile_html:
        args.ofile_html = os.path.join(args.scratchdir, "rdktools_scafnet.html")
    return args


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RDKit scaffold analysis", epilog="")
    args = parse_args(parser)

    logging.basicConfig(
        format="%(levelname)s:%(message)s",
        level=(logging.DEBUG if args.verbose > 1 else logging.INFO),
    )

    logging.info(f"RDKit version: {rdkit.rdBase.rdkitVersion}")
    logging.info(f"MatplotLib version: {mpl.__version__}")
    logging.info(f"Pyvis version: {pyvis.__version__}")

    t0 = time.time()

    if hasattr(args, "scratchdir"):
        os.makedirs(args.scratchdir, exist_ok=True)

    if args.op == "demobm":
        scaffold.DemoBM()
        sys.exit()

    elif args.op == "demonet_img":
        scaffold.DemoNetImg(args.scratchdir)
        sys.exit()

    elif args.op == "demonet_html":
        logging.debug(f"scratchdir: {args.scratchdir}")
        scaffold.DemoNetHtml(args.scratchdir)
        sys.exit()

    if not (args.ifile):
        parser.error("--i required.")

    if args.op == "bmscaf":
        molReader = util.File2Molreader(
            args.ifile, args.idelim, args.smilesColumn, args.nameColumn, args.iheader
        )
        molWriter = util.File2Molwriter(args.ofile, args.odelim, args.oheader)
        mols = util.ReadMols(molReader)
        scafmols = scaffold.Mols2BMScaffolds(mols, molWriter)
        if args.ofile_png:
            img = rdkit.Chem.Draw.MolsToGridImage(scafmols, molsPerRow=8)
            img.save(args.ofile_png, format="PNG")

    elif args.op == "scafnet":
        molReader = util.File2Molreader(
            args.ifile, args.idelim, args.smilesColumn, args.nameColumn, args.iheader
        )
        fout = open(args.ofile, "w") if args.ofile else sys.stdout
        mols = util.ReadMols(molReader)
        scafnet = scaffold.Mols2ScafNet(mols, args.brics, fout)
        if args.ofile_png:
            scaffold.Scafnet2Img(scafnet, args.ofile_png)
        if args.ofile_html:
            scaffold.Scafnet2Html(
                scafnet, args.scafname, args.scratchdir, args.ofile_html
            )

    elif args.op == "scafnet_rings":
        molReader = util.File2Molreader(
            args.ifile, args.idelim, args.smilesColumn, args.nameColumn, args.iheader
        )
        molWriter = util.File2Molwriter(args.ofile, args.odelim, args.oheader)
        n_mol = 0
        for mol in molReader:
            n_mol += 1
            molname = (
                mol.GetProp("_Name") if mol.HasProp("_Name") else f"Mol_{n_mol:03d}"
            )
            scafnet = scaffold.Mols2ScafNet([mol], args.brics, None)
            rings = scaffold.ScafNet2Rings(scafnet, molname, molWriter)
        logging.info(f"Mols processed: {n_mol}")

    else:
        parser.error(f"Unsupported operation: {args.op}")

    logging.info(
        f"""Elapsed: {time.strftime('%Hh:%Mm:%Ss', time.gmtime(time.time()-t0))}"""
    )
