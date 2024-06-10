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
import queue
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
from rdkit.Chem import MolFromSmiles, MolToSmiles, rdchem
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
    elif mol.HasSubstructMatch("[!R]-[!R;D1]"):
        # scaffolds should not include single bonds that are not part of linker between two rings
        return False
    return True


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
    # [!R]-[!R;D1] -> single, non-ring bond
    SMARTS_pat = ["[!#0;R:1]-!@[!#0:2]>>[*:1]-[#0].[#0]-[*:2]"]
    params = _get_hier_scafnet_params(SMARTS_pat)  # TODO: remove

    attrs = [a for a in inspect.getmembers(params) if not (a[0].startswith("__"))]
    for a in attrs:
        logging.debug(f"{a[0]}: {a[1]}")

    # for i, mol in enumerate(mols):
    #    molname = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
    #    logging.debug(f"{i+1}. {molname}:")
    scafnet = None
    mol_idx = 0
    mol_indices = []  # node indices for input molecules
    nx_graph = nx.empty_graph(default=nx.DiGraph)
    prev_edge_i = 0
    for mol in mol_supplier:
        if rdkit.Chem.rdMolDescriptors.CalcNumRings(mol) > 10:
            logging.warn(
                f"Not including {MolToSmiles(mol)} in graph as it has > 10 rings"
            )
            continue
        if scafnet is None:
            scafnet = rdScaffoldNetwork.CreateScaffoldNetwork([mol], params)
        else:
            try:
                rdScaffoldNetwork.UpdateScaffoldNetwork([mol], scafnet, params)
            except rdkit.Chem.rdchem.KekulizeException:
                logging.debug(f"Could not kekulize: {MolToSmiles(mol)}")
        mol_indices.append(mol_idx)
        # have to convert edges to tuple for networkx compatibility
        new_edges = list(
            map(
                lambda e: (e.beginIdx, e.endIdx, {"type": e.type}),
                scafnet.edges[prev_edge_i:],
            )
        )
        if len(new_edges) == 0:
            # molecule is its own scaffold
            new_edges = [(mol_idx, mol_idx, {"type": EdgeType.Initialize})]
        nx_graph.add_edges_from(new_edges)
        mol_idx = len(scafnet.nodes)
        prev_edge_i = len(scafnet.edges)
    if ofile is not None:
        write_scaffold_net(scafnet, ofile, odelimeter, oheader)
    logging.info(f"nodes: {len(scafnet.nodes)}; edges:{len(scafnet.edges)}")
    return scafnet, mol_indices, nx_graph


#############################################################################
class CustomNetworkEdge:
    def __init__(self, edge_type: EdgeType, beginIdx: int, endIdx: int) -> None:
        self.type = edge_type
        self.beginIdx = beginIdx
        self.endIdx = endIdx


def _get_node_edges(node_idx: int, edges: list[NetworkEdge]):
    """
    Get the edges which originate from a given node.
    """
    node_edges = []
    for e in edges:
        if e.beginIdx == node_idx:
            node_edges.append(e)
    return node_edges


def _construct_adjacency_list(scafnet: rdScaffoldNetwork) -> list[list]:
    n_nodes = len(scafnet.nodes)
    adjacency_list = list(
        map(lambda x: _get_node_edges(x, scafnet.edges), range(n_nodes))
    )
    return adjacency_list


def _insert_init_edge(adj_list: list[list], idx: int) -> int:
    """
    Insert initial edge from dummy node to molecule if molecule is a
    scaffold of itself. This is necessary because when we write out scaffolds
    we're looking for edges that have edge type "EdgeType.Initialize".
    :param list[list] adj_list: adjacency list constructed by _construct_adjacency_list.
        adj_list will be updated in-place
    :return int: index of dummy node in adjacency list
    """
    dummy_idx = len(adj_list)  # counts start from 0, so this is max_idx+1
    adj_list.append([CustomNetworkEdge(EdgeType.Initialize, dummy_idx, idx)])
    return dummy_idx


def _insert_init_edge_nx(g: nx.Graph, idx: int) -> int:
    """
    Insert initial edge from dummy node to molecule if molecule is a
    scaffold of itself. This is necessary because when we write out scaffolds
    we're looking for edges that have edge type "EdgeType.Initialize".
    :param list[list] adj_list: adjacency list constructed by _construct_adjacency_list.
        adj_list will be updated in-place
    :return int: index of dummy node in adjacency list
    """
    dummy_idx = len(g.nodes())  # counts start from 0, so this is max_idx+1
    g.add_edge(dummy_idx, idx, type=EdgeType.Initialize)
    return dummy_idx


# TODO: may be better to use floyd-warshall or related algo here
def _get_fragment_map(init_edge: NetworkEdge, adj_list: list[list]) -> dict:
    """
    Perform BFS to generate a mapping from all nodes which are connected
    to init_edge.beginIdx (ie the original molecule) to their depth
    relative to the first molecule. first molecule has depth of 0.
    :param NetworkEdge init_edge: initial edge from input molecule.
    :param list[list] adj_list: adjacency list representation of a ScaffoldNet
    :return dict: mapping from node indices to their depth relative to
        init_edge.beginIdx
    """
    visited = {}
    fragment_nodes = {}  # map nodes to depth
    visited[init_edge] = True
    fragment_nodes[init_edge.beginIdx] = 0
    q = queue.Queue()
    q.put(init_edge)
    while not (q.empty()):
        edge = q.get()
        next_idx = edge.endIdx
        fragment_nodes[edge.endIdx] = fragment_nodes[edge.beginIdx] + 1
        next_edges = adj_list[next_idx]
        for e in next_edges:
            if e.type == EdgeType.Fragment and e not in visited:
                visited[e] = True
                q.put(e)
    return fragment_nodes


def write_hier_scafs(
    fragment_maps: list[dict],
    mol_indices: list[int],
    nodes: list,
    o_mol: str,
    o_scaf: str,
    odelimeter: str,
    oheader: bool,
) -> None:
    # idx == ids
    mol_writer, f_mol = get_csv_writer(o_mol, odelimeter)
    scaf_writer, f_scaf = get_csv_writer(o_scaf, odelimeter)
    if oheader:
        mol_writer.writerow(["mol_id", "mol_smiles", "scaffold_id", "scaffold_depth"])
        scaf_writer.writerow(["scaffold_id", "scaffold_smiles"])
    seen_scafs = [False] * len(nodes)
    for fragment_map, mol_idx in zip(fragment_maps, mol_indices):
        mol_smile = nodes[mol_idx]
        for fragment_idx in fragment_map:
            fragment_smile = nodes[fragment_idx]
            if fragment_idx == mol_idx:
                # don't need to include the molecule itself (redundant)
                continue
            scaf_depth = fragment_map[fragment_idx]
            mol_writer.writerow([mol_idx, mol_smile, fragment_idx, scaf_depth])
            if not (seen_scafs[fragment_idx]):
                scaf_writer.writerow([fragment_idx, fragment_smile])
                seen_scafs[fragment_idx] = True
    close_file(f_mol)
    close_file(f_scaf)


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
    params.keepOnlyFirstFragment = True
    params.pruneBeforeFragmenting = True
    return params


def get_init_edge_idx(edges: list):
    for i, e in enumerate(edges):
        if e.type == EdgeType.Initialize:
            return i
    return -1  # not found


def HierarchicalScaffolds(
    molReader,
    brics: bool = False,
    o_mol: str = None,
    o_scaf: str = None,
    odelim: str = None,
    oheader: bool = False,
):
    params = _get_hier_scafnet_params()
    scafnet, mol_indices, nx_graph = Mols2ScafNet(molReader, params=params)
    nodes = list(scafnet.nodes)
    fragment_maps = []
    for i, mol_idx in enumerate(mol_indices):
        bfs_tree = nx.bfs_tree(nx_graph, source=mol_idx)
        fragment_map = nx.shortest_path_length(bfs_tree, source=mol_idx)
        fragment_maps.append(fragment_map)
    """
    adjacency_list = _construct_adjacency_list(scafnet)
    for i, mol_idx in enumerate(mol_indices):
        edges = adjacency_list[mol_idx]
        init_edge_idx = get_init_edge_idx(edges)
        if init_edge_idx != -1:
            e = edges[init_edge_idx]
        else:
            # molecule is its own scaffold, need to insert dummy edge
            new_mol_idx = _insert_init_edge(adjacency_list, mol_idx)
            mol_indices[i] = new_mol_idx
            mol_idx = new_mol_idx
            nodes.append(new_mol_idx)
            e = adjacency_list[new_mol_idx][0]
        fragment_map = _get_fragment_map(e, adjacency_list)
        fragment_maps.append(fragment_map)
    """
    if o_mol is not None:
        write_hier_scafs(
            fragment_maps, mol_indices, nodes, o_mol, o_scaf, odelim, oheader
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
