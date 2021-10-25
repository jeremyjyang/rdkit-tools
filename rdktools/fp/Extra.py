#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.MurckoScaffold.html
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.rdScaffoldNetwork.html
rdScaffoldNetwork available RDKit 2020.03.1+.
"""
#############################################################################
import os,sys,re,logging,argparse

import rdkit

from .. import fp as rdktools_fp

#############################################################################
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="RDKit fp utils")
  OPS=[ "demopath", "demomorgan", "demomaccs", "list_maccskeys" ]
  parser.add_argument("op", choices=OPS, help="OPERATION")
  parser.add_argument("-v", "--verbose", action="count", default=0)
  args = parser.parse_args()

  logging.basicConfig(format='%(levelname)s:%(message)s', level=(logging.DEBUG if args.verbose>1 else logging.INFO))

  logging.info(f"RDKit version: {rdkit.rdBase.rdkitVersion}")

  if args.op == "list_maccskeys":
    rdktools_fp.ListMACCSkeys()

  elif args.op == "demomaccs":
    rdktools_fp.DemoMACCSKeys()

  elif args.op == "demopath":
    rdktools_fp.DemoPath()

  elif args.op == "demomorgan":
    rdktools_fp.DemoMorgan()

  else:
    parser.error(f"Unsupported operation: {args.op}")
