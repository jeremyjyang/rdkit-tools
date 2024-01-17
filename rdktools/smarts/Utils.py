#!/usr/bin/env python3
#############################################################################
import logging

import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.rdBase
from rdkit.Chem import FilterCatalog

from rdktools.smarts.SmartsFile import SmartsFile


def get_smarts_queries(smarts_file_path: str, strict_smarts: bool):
    # get SMARTS queries from file
    sf = SmartsFile(smarts_file_path, strict_smarts)
    logging.info(f"SMARTS read from file {smarts_file_path}: {len(sf.smarts_strs)}")
    if len(sf.failed_smarts) > 0:
        failed_smarts = [
            (fs, fsr) for fs, fsr in zip(sf.failed_smarts, sf.failed_smarts_raw)
        ]
        logging.error(
            "The following SMARTS could not be interpreted and will be ignored:"
        )
        logging.error(f"(parsed string, raw string): {failed_smarts}")
    queries = [
        {"smarts": smarts, "pat": mol, "n_mol_matched": 0}
        for smarts, mol in zip(sf.smarts_strs, sf.smarts_mols)
    ]
    return queries


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
def MatchFilter(smarts: str, molReader, molWriter, exclude_mol_props: bool):
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
def MatchFilterMulti(
    smarts_file_path: str,
    strict_smarts: bool,
    molReader,
    molWriter,
    exclude_mol_props: bool,
):
    """All SMARTS must match, or mol is filtered."""
    queries = get_smarts_queries(smarts_file_path, strict_smarts)
    n_mol = 0
    n_mol_matched = 0
    for mol in molReader:
        if exclude_mol_props:
            clearMolProps(mol)
        if not mol.HasProp("_Name") or mol.GetProp("_Name") == "":
            mol.SetProp("_Name", f"mol_{n_mol+1}")
        for j, query in enumerate(queries):
            matches = mol.GetSubstructMatches(
                query["pat"], uniquify=True, useChirality=False
            )
            if len(matches) > 0:
                query["n_mol_matched"] += 1
        matchcounts = [q["n_mol_matched"] for q in queries]
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
    smarts_file_path: str,
    strict_smarts: bool,
    usa: bool,
    molReader,
    molWriter,
    exclude_mol_props: bool,
):
    queries = get_smarts_queries(smarts_file_path, strict_smarts)
    n_mol = 0
    for mol in molReader:
        if exclude_mol_props:
            clearMolProps(mol)
        if not mol.HasProp("_Name") or mol.GetProp("_Name") == "":
            mol.SetProp("_Name", f"mol_{n_mol+1}")
        for j, query in enumerate(queries):
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
                [f"query_{j+1:02d}:{q['n_mol_matched']}" for j, q in enumerate(queries)]
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
    n_fail = 0
    n_empty = 0
    n_pass = 0
    n_err = 0
    for mol in molReader:
        if mol is None:
            n_empty += 1
            n_err += 1
            continue
        if exclude_mol_props:
            clearMolProps(mol)
        if not mol.HasProp("_Name") or mol.GetProp("_Name") == "":
            mol.SetProp("_Name", f"mol_{n_mol+1}")
        name = mol.GetProp("_Name")
        if catalog.HasMatch(mol):
            entry = catalog.GetFirstMatch(mol)
            for filterMatch in entry.GetFilterMatches(mol):
                logging.debug(
                    f'{n_mol}. "{name}": PAINS filter matched: {filterMatch.filterMatch}'
                )
                mol.SetProp(f"{filterMatch.filterMatch}", "TRUE")
            mol.SetProp(f"PAINS", "FAIL")
            n_fail += 1
        else:
            mol.SetProp(f"PAINS", "PASS")
            if n_mol == 0:
                molWriter.SetProps(mol.GetPropNames())
            molWriter.write(mol)
            n_pass += 1
        n_mol += 1
    logging.info(f"n_mol: {n_mol}; n_pass: {n_pass}; n_fail: {n_fail}; n_empty: {n_empty}; n_err: {n_err}")

#############################################################################
