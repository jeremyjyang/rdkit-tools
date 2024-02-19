#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.MurckoScaffold.html
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.rdScaffoldNetwork.html
rdScaffoldNetwork available RDKit 2020.03.1+.

https://matplotlib.org/
https://pyvis.readthedocs.io/en/latest/
"""
import inspect
import json
import logging

#############################################################################
import os
import re
import stat
import tempfile

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


#############################################################################
def Mols2BMScaffolds(mols, molWriter):
    scafmols = []
    for i, mol in enumerate(mols):
        molname = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        logging.debug(f"{i+1}. {molname}")
        scafmol = MurckoScaffold.GetScaffoldForMol(mol)
        Compute2DCoords(scafmol, clearConfs=True)
        scafmols.append(scafmol)
        molWriter.write(scafmol)
    logging.info(f"{len(mols)} mols written to {molWriter}")
    return scafmols


#############################################################################
def Mols2ScafNet(mols, brics=False, fout=None):
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
    if fout is not None:
        for i in range(len(scafnet.nodes)):
            fout.write(f"node\t{i}\t{scafnet.nodes[i]}\t{scafnet.counts[i]}\n")
        for i in range(len(scafnet.edges)):
            fout.write(
                f"edge\t{i}\t{scafnet.edges[i].beginIdx}\t{scafnet.edges[i].endIdx}\t{scafnet.edges[i].type}\n"
            )
        fout.flush()
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
            if len(ring.GetSubstructMatches(pat)) > 0:
                continue
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
        prefix=scratchdir + "/", suffix=".png", mode="w+b", delete=False
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
def Scafnet2Img(scafnet, ofile):
    # title="RDKit_ScafNet:"+re.sub(r'^[^\s]*\s+(.*)$', r'\1', smi)) #How to add title?
    scafmols = [MolFromSmiles(m) for m in scafnet.nodes]
    img = rdkit.Chem.Draw.MolsToGridImage(
        scafmols,
        legends=[f"{i}, counts: {c}" for i, c in enumerate(scafnet.counts)],
        molsPerRow=4,
    )
    logging.debug(f"Writing scafnet PNG to: {ofile}")
    img.save(ofile, format="PNG")
    return img


#############################################################################
def DemoNetHtml(scratchdir):
    logging.debug(f"scratchdir: {scratchdir}")
    fout = tempfile.NamedTemporaryFile(
        prefix=scratchdir + "/", suffix=".html", mode="w+", delete=False
    )
    ofile = fout.name
    logging.debug(f"ofile: {ofile}")
    fout.close()
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
    os.chmod(
        ofile,
        stat.S_IRUSR
        | stat.S_IWUSR
        | stat.S_IRGRP
        | stat.S_IWGRP
        | stat.S_IROTH
        | stat.S_IWOTH,
    )
    g.show(ofile)


#############################################################################
def Scafnet2Html(scafnet, scafname, scratchdir, ofile):
    logging.debug(f"pyvis.network.Network()...")
    g = pyvis.network.Network(
        notebook=False, height="800px", width="1000px", heading=scafname
    )
    logging.debug(f"pyvis.network.Network()... Done.")

    for i, n in enumerate(scafnet.nodes):
        logging.debug(f"{i+1}. util.moltosvg(rdkit.Chem.MolFromSmiles({n})...")
        svg = util.moltosvg(rdkit.Chem.MolFromSmiles(n))
        with open(f"{scratchdir}/{i}.svg", "w") as outf:
            outf.write(svg)
        logging.debug(f"g.add_node()...")
        g.add_node(
            i,
            shape="image",
            label=" ",
            image=f"{scratchdir}/{i}.svg",
            title=svg,
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
    # g.show_buttons()
    if ofile:
        logging.info(f"Writing SCAFNET HTML to: {ofile}")
        g.save_graph(ofile)
    return g


#############################################################################
