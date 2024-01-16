#!/usr/bin/env python3
#############################################################################
import logging

import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.rdBase
from rdkit.Chem import FilterCatalog

from rdktools.smarts.SmartsFile import SmartsFile


#############################################################################
### DeduplicateMatches() - remove matches not USA (Unique Sets of Atoms)
### uumatches - unique set of matches, sorted atom order
### umatches - unique set of matches, original atom order
#############################################################################
def DeduplicateMatches(matches):
    uumatches = []
    umatches = []
    for match in matches:
        uumatch = list(match)
        uumatch.sort()
        if uumatch not in uumatches:
            uumatches.append(uumatch)
            umatches.append(match)
    return tuple(umatches)


#############################################################################
def clearMolProps(mol):
    for prop in mol.GetPropNames():
        if prop != "_Name":
            mol.ClearProp(prop)


#############################################################################
def MatchFilter(smarts, molReader, molWriter, exclude_mol_props: bool):
    pat = rdkit.Chem.MolFromSmarts(smarts)
    n_mol = 0
    n_mol_matched = 0
    for mol in molReader:
        if exclude_mol_props:
            clearMolProps(mol)
        if not mol.HasProp("_Name") or mol.GetProp("_Name") == "":
            mol.SetProp("_Name", f"mol_{n_mol+1}")
        matches = mol.GetSubstructMatches(pat, uniquify=True, useChirality=False)
        if len(matches) > 0:
            n_mol_matched += 1
            molWriter.write(mol)
        n_mol += 1
    logging.info(f"n_mol: {n_mol}; n_mol_matched: {n_mol_matched}")


#############################################################################
def MatchFilterMulti(smarts_file_path, molReader, molWriter, exclude_mol_props: bool):
    """All SMARTS must match, or mol is filtered."""
    sf = SmartsFile(smarts_file_path)
    logging.info(f"SMARTS read from file {smarts_file_path}: {len(sf.smarts_strs)}")
    querys = [
        {"smarts": smarts, "pat": mol, "n_mol_matched": 0}
        for smarts, mol in zip(sf.smarts_strs, sf.smarts_mols)
    ]
    n_mol = 0
    n_mol_matched = 0
    for mol in molReader:
        if exclude_mol_props:
            clearMolProps(mol)
        if not mol.HasProp("_Name") or mol.GetProp("_Name") == "":
            mol.SetProp("_Name", f"mol_{n_mol+1}")
        for j, query in enumerate(querys):
            matches = mol.GetSubstructMatches(
                query["pat"], uniquify=True, useChirality=False
            )
            if len(matches) > 0:
                query["n_mol_matched"] += 1
        matchcounts = [q["n_mol_matched"] for q in querys]
        if min(matchcounts) > 0:
            molWriter.write(mol)
            n_mol_matched += 1
        n_mol += 1
    logging.info(f"n_mol: {n_mol}; n_mol_matched: {n_mol_matched}")


#############################################################################
def MatchCounts(smarts: str, usa: bool, molReader, molWriter, exclude_mol_props: bool):
    pat = rdkit.Chem.MolFromSmarts(smarts)
    n_mol = 0
    n_mol_matched = 0
    for mol in molReader:
        if exclude_mol_props:
            clearMolProps(mol)
        if not mol.HasProp("_Name") or mol.GetProp("_Name") == "":
            mol.SetProp("_Name", f"mol_{n_mol+1}")
        # name = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
        # smi = MolToSmiles(mol, isomericSmiles=True)
        matches = mol.GetSubstructMatches(pat, uniquify=True, useChirality=False)
        if usa:
            matches = DeduplicateMatches(matches)
        n_matches = len(matches)
        n_mol_matched += 1 if n_matches > 0 else 0
        mol.SetProp("n_matches", f"{n_matches}")
        if n_mol == 0:
            molWriter.SetProps(mol.GetPropNames())
        molWriter.write(mol)
        n_mol += 1
    logging.info(f"n_mol: {n_mol}; n_mol_matched: {n_mol_matched}")


#############################################################################
def MatchCountsMulti(
    smarts_file_path, usa, molReader, molWriter, exclude_mol_props: bool
):
    sf = SmartsFile(smarts_file_path)
    logging.info(f"SMARTS read from file {smarts_file_path}: {len(sf.smarts_strs)}")
    querys = [
        {"smarts": smarts, "pat": mol, "n_mol_matched": 0}
        for smarts, mol in zip(sf.smarts_strs, sf.smarts_mols)
    ]
    n_mol = 0
    for mol in molReader:
        if exclude_mol_props:
            clearMolProps(mol)
        if not mol.HasProp("_Name") or mol.GetProp("_Name") == "":
            mol.SetProp("_Name", f"mol_{n_mol+1}")
        print("NAME:", mol.GetProp("_Name"))
        for j, query in enumerate(querys):
            matches = mol.GetSubstructMatches(
                query["pat"], uniquify=True, useChirality=False
            )
            if usa:
                matches = DeduplicateMatches(matches)
            n_matches = len(matches)
            if n_matches > 0:
                query["n_mol_matched"] += 1
            mol.SetProp(
                f"""n_matches(query_{j+1:02d} = "{query['smarts']}")""", f"{n_matches}"
            )
        if n_mol == 0:
            molWriter.SetProps(mol.GetPropNames())
        molWriter.write(mol)
        n_mol += 1
    logging.info(
        f"n_mol: {n_mol}; mols matched: "
        + (
            ",".join(
                [f"query_{j+1:02d}:{q['n_mol_matched']}" for j, q in enumerate(querys)]
            )
        )
    )


#############################################################################
# https://github.com/rdkit/rdkit/tree/master/Code/GraphMol/FilterCatalog
def FilterPAINS(molReader, molWriter, exclude_mol_props: bool):
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
    catalog = FilterCatalog.FilterCatalog(params)
    n_mol = 0
    n_filtered = 0
    n_out = 0
    for mol in molReader:
        if exclude_mol_props:
            clearMolProps(mol)
        if not mol.HasProp("_Name") or mol.GetProp("_Name") == "":
            mol.SetProp("_Name", f"mol_{n_mol+1}")
        name = mol.GetProp("_Name")
        if catalog.HasMatch(mol):
            entry = catalog.GetFirstMatch(mol)
            # logging.debug(f"{n_mol}. matched filter; first match: {entry.GetDescription() if entry else ''}")
            for filterMatch in entry.GetFilterMatches(mol):
                logging.debug(
                    f'{n_mol}. "{name}": PAINS filter matched: {filterMatch.filterMatch}'
                )
                mol.SetProp(f"{filterMatch.filterMatch}", "TRUE")
            mol.SetProp(f"PAINS", "FAIL")
            n_filtered += 1
        else:
            mol.SetProp(f"PAINS", "PASS")
            if n_mol == 0:
                molWriter.SetProps(mol.GetPropNames())
            molWriter.write(mol)
            n_out += 1
        n_mol += 1
    logging.info(f"n_mol: {n_mol}; n_out: {n_out}; n_filtered: {n_filtered}")


#############################################################################
