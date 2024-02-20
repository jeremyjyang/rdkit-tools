#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.MurckoScaffold.html
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.rdScaffoldNetwork.html
rdScaffoldNetwork available RDKit 2020.03.1+.

https://matplotlib.org/
https://pyvis.readthedocs.io/en/latest/
"""
import csv
import inspect
import json
import logging

#############################################################################
import os
import re
import sys
import tempfile
from typing import Optional

import pyvis
import rdkit
import rdkit.Chem
import rdkit.Chem.AllChem

# fingerprints:
from rdkit.Chem import MolFromSmiles, MolToSmiles
from rdkit.Chem.AllChem import Compute2DCoords

# scaffolds:
from rdkit.Chem.Scaffolds import MurckoScaffold, rdScaffoldNetwork

from .. import util


# from matplotlib import pyplot as plt
def ensure_path_separator(dir: str):
    """
    Ensure that given path is a directory by appending
    '/' or whatever path separator is appropriate.

    :param str dir: path to directory
    :return _type_: dir with path separator appended (will be same as dir if separator already present)
    """
    return os.path.join(dir, "")


def fix_pyvis_header_html(html_file_path: str):
    """
    Fix issue with double-heading in pyvis-generated HTML file.
    :param str html_file_path: path to HTML file
    """
    with open(html_file_path, "r+") as f:
        html_txt = f.read()
        html_rev = re.sub(r"<center>.+?<\/h1>\s+<\/center>", "", html_txt, 1, re.DOTALL)
        f.seek(0)
        f.write(html_rev)
        f.truncate()


def center_align_pyvis_html(html_file_path: str):
    """
    Center pyvis-generated graph in HTML file (rather than left-align).
    :param str html_file_path: path to HTML file
    """
    # Open the HTML file, add a CSS rule to center the graph, and save the changes.
    with open(html_file_path, "r+") as f:
        html_txt = f.read()
        html_txt = html_txt.replace("float: left;", "margin: auto;")
        f.seek(0)
        f.write(html_txt)
        f.truncate()


def write_scaffold_net(
    scafnet, ofile: str, odelimeter: str = ",", oheader: bool = False
):
    if ofile is sys.stdout:
        f = ofile
    else:
        f = open(ofile, "w")
    net_writer = csv.writer(f, delimiter=odelimeter)
    if oheader:
        net_writer.writerow(["element_type", "index", "info"])
    for i in range(len(scafnet.nodes)):
        info = {
            "SMILES": scafnet.nodes[i],
            "Counts": scafnet.counts[i],
        }
        net_writer.writerow(["node", i, json.dumps(info)])
    for i in range(len(scafnet.edges)):
        info = {
            "beginIdx": scafnet.edges[i].beginIdx,
            "endIdx": scafnet.edges[i].endIdx,
            "edgeType": str(scafnet.edges[i].type),
        }
        net_writer.writerow(
            [
                "edge",
                i,
                json.dumps(info),
            ]
        )
    if ofile is not sys.stdout:
        f.close()


#############################################################################
def Mols2BMScaffolds(mols, molWriter):
    scafmols = []
    legends = []
    for i, mol in enumerate(mols):
        molname = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        logging.debug(f"{i+1}. {molname}")
        scafmol = MurckoScaffold.GetScaffoldForMol(mol)
        Compute2DCoords(scafmol, clearConfs=True)
        scafmols.append(scafmol)
        legends.append(molname)
        molWriter.write(scafmol)
    logging.info(f"{len(mols)} mols written to {molWriter}")
    return scafmols, legends


#############################################################################
def Mols2ScafNet(
    mols,
    brics=False,
    ofile: Optional[str] = None,
    odelimeter: str = ",",
    oheader: bool = False,
):
    if brics:
        params = rdScaffoldNetwork.BRICSScaffoldParams()
    else:
        params = rdScaffoldNetwork.ScaffoldNetworkParams()
        params.flattenChirality = True
        params.flattenIsotopes = True
        params.flattenKeepLargest = True
        params.includeGenericBondScaffolds = False
        params.includeGenericScaffolds = False
        params.includeScaffoldsWithAttachments = False
        params.includeScaffoldsWithoutAttachments = True
        params.keepOnlyFirstFragment = False
        params.pruneBeforeFragmenting = True

    attrs = [a for a in inspect.getmembers(params) if not (a[0].startswith("__"))]
    for a in attrs:
        logging.debug(f"{a[0]}: {a[1]}")

    for i, mol in enumerate(mols):
        molname = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        logging.debug(f"{i+1}. {molname}:")

    scafnet = rdScaffoldNetwork.CreateScaffoldNetwork(mols, params)
    if ofile is not None:
        write_scaffold_net(scafnet, ofile, odelimeter, oheader)
    logging.info(f"nodes: {len(scafnet.nodes)}; edges:{len(scafnet.edges)}")
    return scafnet


#############################################################################
def ScafNet2Rings(scafnet, name, molWriter):
    """Output unique ringsystems only."""
    ringsmis = set()
    rings = []
    pat = rdkit.Chem.MolFromSmarts("*!@-*")  # Non-ring single bond.
    for smi in scafnet.nodes:
        logging.debug(f"node: {smi}")
        if smi not in ringsmis:
            ring = MolFromSmiles(smi)  # Ring-maybe.
            if len(smi) == 0 or len(ring.GetSubstructMatches(pat)) > 0:
                continue
            else:
                ringsmis.add(smi)
                rings.append(ring)
    rings_mol = MolFromSmiles(".".join(sorted(list(ringsmis))))
    rings_mol.SetProp("_Name", f"{name}_RINGS")
    molWriter.write(rings_mol)
    logging.info(f"{name}: rings: {len(ringsmis)}")
    return rings


#############################################################################
def DemoBM():
    scafmols = []
    for smi in util.DEMOSMIS:
        mol = MolFromSmiles(re.sub(r"\s.*$", "", smi))
        scafmol = MurckoScaffold.GetScaffoldForMol(mol) if mol else None
        scafmols.append(scafmol)
        smi_std = MolToSmiles(scafmol, isomericSmiles=False) if scafmol else None
        logging.info(f"{smi:>28s} >> {smi_std}")
    img = rdkit.Chem.Draw.MolsToGridImage(scafmols, molsPerRow=4)
    img.show()


#############################################################################
def DemoNetImg(scratchdir):
    fout = tempfile.NamedTemporaryFile(
        prefix=ensure_path_separator(scratchdir),
        suffix=".png",
        mode="w+b",
        delete=False,
    )
    ofile = fout.name
    fout.close()
    brics = True
    logging.debug(f"DemoNetImg({brics}, {fout.name})")
    smi = "Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12 flucloxacillin"
    mols = [MolFromSmiles(re.sub(r"\s.*$", "", smi))]
    scafnet = Mols2ScafNet(mols, False)
    logging.info(f"Scafnet nodes: {len(scafnet.nodes)}; edges: {len(scafnet.edges)}")
    # scafmols = [MolFromSmiles(m) for m in scafnet.nodes]
    scafmols = []
    for i, m in enumerate(scafnet.nodes[:]):
        logging.debug(f"{i+1}. MolFromSmiles({m})...")
        scafmols.append(MolFromSmiles(m))
    logging.info(f"Scafmols: {len(scafmols)}")
    img = Scafnet2Img(scafnet, ofile)
    img.show()


#############################################################################
def Scafnet2Img(scafnet, ofile, molsPerRow: int = 4):
    # title="RDKit_ScafNet:"+re.sub(r'^[^\s]*\s+(.*)$', r'\1', smi)) #How to add title?
    scafmols = [MolFromSmiles(m) for m in scafnet.nodes]
    img = rdkit.Chem.Draw.MolsToGridImage(
        scafmols,
        legends=[f"Idx: {i} , Counts: {c}" for i, c in enumerate(scafnet.counts)],
        molsPerRow=molsPerRow,
    )
    logging.debug(f"Writing scafnet PNG to: {ofile}")
    img.save(ofile, format="PNG")
    return img


#############################################################################
def DemoNetHtml(scratchdir):
    logging.debug(f"scratchdir: {scratchdir}")
    fout = tempfile.NamedTemporaryFile(
        prefix=ensure_path_separator(scratchdir),
        suffix=".html",
        mode="w+",
        delete=False,
    )
    ofile = fout.name
    logging.debug(f"ofile: {ofile}")
    logging.debug(f"DemoNetHtml({scratchdir}, {ofile})")
    demosmi = "Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12 flucloxacillin"
    mols = [MolFromSmiles(re.sub(r"\s.*$", "", demosmi))]
    scafnet = Mols2ScafNet(mols, False)
    logging.info(f"Scafnet nodes: {len(scafnet.nodes)}; edges: {len(scafnet.edges)}")
    g = Scafnet2Html(
        scafnet,
        "RDKit_ScafNet: " + re.sub(r"^[^\s]*\s+(.*)$", r"\1", demosmi),
        scratchdir,
        ofile,
    )
    fout.close()


#############################################################################
def Scafnet2Html(scafnet, scafname, scratchdir, ofile):
    logging.debug(f"pyvis.network.Network()...")
    g = pyvis.network.Network(
        notebook=False, height="800px", width="1000px", heading=scafname
    )
    logging.debug(f"pyvis.network.Network()... Done.")

    for i, smiles in enumerate(scafnet.nodes):
        logging.debug(f"{i+1}. util.moltosvg(rdkit.Chem.MolFromSmiles({smiles})...")
        svg = util.moltosvg(rdkit.Chem.MolFromSmiles(smiles))
        image_path = os.path.join(scratchdir, f"{i}.svg")
        with open(image_path, "w") as outf:
            outf.write(svg)
        logging.debug(f"g.add_node()...")
        g.add_node(
            i,
            shape="image",
            label=" ",
            image=image_path,
            title=f"SMILES: {smiles}\nIndex: {i}\nCounts: {scafnet.counts[i]}",
            size=60,
        )
        # Segmentation fault after last node.
    logging.debug(f"util.moltosvg()... Done.")
    for i, e in enumerate(scafnet.edges):
        logging.debug(f"{i+1}. g.add_edge()...")
        g.add_edge(e.beginIdx, e.endIdx, label=str(e.type))
    logging.debug(f"g.add_edge()... Done.")

    VisJS_options = {
        "edges": {"font": {"size": 20}},
        "nodes": {"font": {"color": "rgba(214,47,66,1)", "size": 16, "face": "tahoma"}},
        "physics": {
            "forceAtlas2Based": {
                "gravitationalConstant": -120,
                "springLength": 200,
                "avoidOverlap": 0.42,
            },
            "minVelocity": 0.75,
            "solver": "forceAtlas2Based",
        },
    }
    logging.debug(f"g.set_options()...")
    g.set_options(options=json.dumps(VisJS_options))
    if ofile:
        logging.info(f"Writing SCAFNET HTML to: {ofile}")
        # have to change working directory g.save_graph() only handles local files (for some reason)
        current_dir = os.getcwd()
        ofile_dir, ofile_name = os.path.dirname(ofile), os.path.basename(ofile)
        if ofile_dir != "":
            os.chdir(ofile_dir)
        g.save_graph(ofile_name)
        fix_pyvis_header_html(ofile)
        center_align_pyvis_html(ofile)
        # go back to original dir
        os.chdir(current_dir)
    return g


#############################################################################
