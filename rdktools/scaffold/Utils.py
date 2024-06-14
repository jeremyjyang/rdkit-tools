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
import os
import re
import sys
import tempfile
from typing import Optional, Union

import networkx as nx
import pyvis
import rdkit
import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.Chem.rdmolops

# fingerprints:
from rdkit.Chem import MolFromSmiles, MolToSmiles, rdchem, rdChemReactions
from rdkit.Chem.AllChem import Compute2DCoords

# scaffolds:
from rdkit.Chem.Scaffolds import MurckoScaffold, rdScaffoldNetwork
from rdkit.Chem.Scaffolds.rdScaffoldNetwork import EdgeType, NetworkEdge

from .. import util


def get_csv_writer(file_path: str, delimiter: str):
    if file_path is sys.stdout:
        f = file_path
    else:
        f = open(file_path, "w")
    csv_writer = csv.writer(f, delimiter=delimiter)
    return csv_writer, f


def close_file(f):
    if f is not sys.stdout:
        f.close()


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
    scafnet: rdScaffoldNetwork, ofile: str, odelimeter: str = ",", oheader: bool = False
):
    net_writer, f = get_csv_writer(ofile, odelimeter)
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
    close_file(f)


#############################################################################
def Mols2BMScaffolds(mols: list[rdchem.Mol], molWriter):
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
# TODO: add option to specify scafnet params via file
def _get_default_scafnet_params(brics: bool):
    if brics:
        params = rdScaffoldNetwork.BRICSScaffoldParams()
    else:
        params = rdScaffoldNetwork.ScaffoldNetworkParams()
        params.flattenChirality = True
        params.flattenIsotopes = True
        params.flattenKeepLargest = True
        params.includeGenericBondScaffolds = False
        params.includeGenericScaffolds = False
        params.includeScaffoldsWithAttachments = True
        params.includeScaffoldsWithoutAttachments = False
        params.keepOnlyFirstFragment = False
        params.pruneBeforeFragmenting = True
    return params


def is_valid_scaf(can_smiles: str):
    if len(can_smiles) == 0:
        # empty string given
        return False
    elif can_smiles == "c1ccccc1" or can_smiles == "C1=CC=CC=C1":
        # benzene excluded from scaffolds
        return False
    mol = MolFromSmiles(can_smiles)
    if mol is None:
        # invalid mol
        return False
    ring_system = rdkit.Chem.MolFromSmarts("[R]")
    if not (mol.HasSubstructMatch(ring_system)):
        # scaffolds should contain at least one ring-system
        return False
    single_terminal_bond = rdkit.Chem.MolFromSmarts("*-[D1]")
    if mol.HasSubstructMatch(single_terminal_bond):
        # scaffolds should not include single terminal bonds
        return False
    return True


def update_mol_indices(scafnet, mol_indices: list[int], self_scafs: list[str]):
    node_map = {smi: i for i, smi in enumerate(scafnet.nodes)}
    n = len(mol_indices)
    for mol_smiles in self_scafs:
        if mol_smiles in node_map:
            idx = node_map[mol_smiles]
            mol_indices.append(idx)
        else:
            # not expected to occur
            raise ValueError(f"Mol not found for {mol_smiles}")
    is_self_scaf = [False] * len(mol_indices)
    is_self_scaf[n:] = [True] * (len(mol_indices) - n)
    return mol_indices, is_self_scaf


def Mols2ScafNet(
    mol_supplier: Union[list[rdchem.Mol], rdkit.Chem.rdmolfiles.SmilesMolSupplier],
    brics: bool = False,
    ofile: Optional[str] = None,
    odelimeter: str = ",",
    oheader: bool = False,
    params: Optional[rdScaffoldNetwork.ScaffoldNetworkParams] = None,
):
    if params is None:
        params = _get_default_scafnet_params(brics)
    attrs = [a for a in inspect.getmembers(params) if not (a[0].startswith("__"))]
    for a in attrs:
        logging.debug(f"{a[0]}: {a[1]}")

    scafnet = None
    mol_idx = 0
    mol_indices = []  # node indices for input molecules
    seen_mol_smiles = {}
    # molecules which are scaffolds of themselves OR are scaffolds of another mol
    # we need to track these for downstream tasks
    self_scafs = []
    prev_edge_i = 0
    for i, mol in enumerate(mol_supplier):
        mol_smiles = MolToSmiles(mol)
        if mol_smiles in seen_mol_smiles:
            # molecule is a duplicate
            logging.warn(f"Mol {i, mol_smiles} is a duplicate; skipping over")
            continue
        seen_mol_smiles[mol_smiles] = True
        if rdkit.Chem.rdMolDescriptors.CalcNumRings(mol) > 10:
            logging.warn(f"Not including {i, mol_smiles} in graph as it has > 10 rings")
            continue
        if scafnet is None:
            scafnet = rdScaffoldNetwork.CreateScaffoldNetwork([mol], params)
        else:
            try:
                rdScaffoldNetwork.UpdateScaffoldNetwork([mol], scafnet, params)
            except rdkit.Chem.rdchem.KekulizeException:
                logging.warn(f"Could not kekulize: {i, mol_smiles}")
        if (
            mol_idx == len(scafnet.nodes)
            or prev_edge_i == len(scafnet.edges)
            or scafnet.edges[prev_edge_i].type != EdgeType.Initialize
        ):
            # molecule is a self-scaf or is the scaffold of another mol
            self_scafs.append(mol_smiles)
        else:
            mol_indices.append(mol_idx)
        mol_idx = len(scafnet.nodes)
        prev_edge_i = len(scafnet.edges)
    # include indices of mols that are self-scafs
    mol_indices, is_self_scaf = update_mol_indices(scafnet, mol_indices, self_scafs)
    if ofile is not None:
        write_scaffold_net(scafnet, ofile, odelimeter, oheader)
    logging.info(f"nodes: {len(scafnet.nodes)}; edges:{len(scafnet.edges)}")
    return scafnet, mol_indices, is_self_scaf


#############################################################################
# convert a scafnet to networkx representation
def ScafNet2NetworkX(scafnet):
    edges = list(
        map(
            lambda e: (e.beginIdx, e.endIdx),
            scafnet.edges,
        )
    )
    nodes_with_smiles = [(i, {"smiles": smi}) for i, smi in enumerate(scafnet.nodes)]
    nx_graph = nx.DiGraph()
    nx_graph.add_nodes_from(nodes_with_smiles)
    nx_graph.add_edges_from(edges)
    return nx_graph


def write_hier_scafs(
    o_mol: str,
    o_scaf: str,
    o_mol2scaf: str,
    fragment_maps: list[dict],
    mol_indices: list[int],
    is_self_scaf: list[bool],
    nodes: list,
    odelimeter: str,
    oheader: bool,
) -> None:
    # idx == ids
    mol_writer, f_mol = get_csv_writer(o_mol, odelimeter)
    scaf_writer, f_scaf = get_csv_writer(o_scaf, odelimeter)
    mol2scaf_writer, f_mol2scaf = get_csv_writer(o_mol2scaf, odelimeter)
    if oheader:
        mol_writer.writerow(["mol_id", "smiles"])
        scaf_writer.writerow(["scaffold_id", "smiles"])
        mol2scaf_writer.writerow(["mol_id", "scaffold_id", "scaffold_depth"])
    seen_scafs = {}
    seen_invalid_scafs = {}
    for fragment_map, mol_idx, iss in zip(fragment_maps, mol_indices, is_self_scaf):
        mol_smile = nodes[mol_idx]
        mol_writer.writerow([mol_idx, mol_smile])
        for fragment_idx in fragment_map:
            fragment_smile = nodes[fragment_idx]
            if fragment_idx == mol_idx and not (iss):
                # don't need to include the molecule itself (redundant)
                # (unless molecule is a self-scaffold)
                pass
            elif fragment_idx in seen_invalid_scafs or not (
                is_valid_scaf(fragment_smile)
            ):
                seen_invalid_scafs[fragment_idx] = True
            else:
                scaf_depth = fragment_map[fragment_idx]
                mol2scaf_writer.writerow([mol_idx, fragment_idx, scaf_depth])
                if fragment_idx not in seen_scafs:
                    scaf_writer.writerow([fragment_idx, fragment_smile])
                    seen_scafs[fragment_idx] = True
    close_file(f_mol)
    close_file(f_scaf)
    close_file(f_mol2scaf)


def _get_hier_scafnet_params(fragmentation_rules: Optional[list[str]] = None):
    if fragmentation_rules:
        params = rdScaffoldNetwork.ScaffoldNetworkParams(fragmentation_rules)
    else:
        params = rdScaffoldNetwork.ScaffoldNetworkParams()
    params.collectMolCounts = False
    params.flattenChirality = True
    params.flattenIsotopes = True
    params.flattenKeepLargest = True
    params.includeGenericBondScaffolds = False
    params.includeGenericScaffolds = False
    params.includeScaffoldsWithAttachments = False
    params.includeScaffoldsWithoutAttachments = True
    params.keepOnlyFirstFragment = False
    params.pruneBeforeFragmenting = True
    return params


def trim_graph(g: nx.DiGraph, mol_indices: list[int]):
    # remove invalid scaffolds from the graph
    to_remove = []
    # convert to dict for faster lookup time
    mol_indices_dict = {i: True for i in mol_indices}
    for node_idx, node_data in g.nodes(data=True):
        if node_idx in mol_indices_dict:
            # node is a molecule
            pass
        elif not (is_valid_scaf(node_data["smiles"])):
            to_remove.append(node_idx)
    g.remove_nodes_from(to_remove)
    # reconnecting graph where necessary
    for node_idx in to_remove:
        for in_src, _ in g.in_edges(node_idx):
            if in_src in to_remove:
                continue
            for _, out_dst in g.out_edges(node_idx):
                if out_dst in to_remove:
                    continue
                g.add_edge(in_src, out_dst)
    return g


def HierarchicalScaffolds(
    molReader,
    o_mol: str,
    o_scaf: str,
    o_mol2scaf: str,
    odelim: str = None,
    oheader: bool = False,
):
    params = _get_hier_scafnet_params()
    scafnet, mol_indices, is_self_scaf = Mols2ScafNet(molReader, params=params)
    nx_graph = ScafNet2NetworkX(scafnet)
    # nx_graph = trim_graph(nx_graph, mol_indices)  # remove invalid scaffolds
    fragment_maps = []
    for mol_idx in mol_indices:
        bfs_tree = nx.bfs_tree(nx_graph, source=mol_idx)
        fragment_map = nx.shortest_path_length(bfs_tree, source=mol_idx)
        fragment_maps.append(fragment_map)
    if o_mol is not None:
        nodes = nx_graph.nodes.data("smiles")
        write_hier_scafs(
            o_mol,
            o_scaf,
            o_mol2scaf,
            fragment_maps,
            mol_indices,
            is_self_scaf,
            nodes,
            odelim,
            oheader,
        )
    return fragment_maps


############################################################################
def ScafNet2Rings(scafnet: rdScaffoldNetwork, name: str, molWriter):
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
def DemoNetImg(scratchdir: str):
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
    scafnet, _, _ = Mols2ScafNet(mols, False)
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
def Scafnet2Img(scafnet: rdScaffoldNetwork, ofile: str, molsPerRow: int = 4):
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
def DemoNetHtml(scratchdir: str):
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
    scafnet, _, _ = Mols2ScafNet(mols, False)
    logging.info(f"Scafnet nodes: {len(scafnet.nodes)}; edges: {len(scafnet.edges)}")
    g = Scafnet2Html(
        scafnet,
        "RDKit_ScafNet: " + re.sub(r"^[^\s]*\s+(.*)$", r"\1", demosmi),
        scratchdir,
        ofile,
    )
    fout.close()


#############################################################################
def Scafnet2Html(
    scafnet: rdScaffoldNetwork,
    scafname: str,
    scratchdir: str,
    ofile: Optional[str] = None,
):
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
