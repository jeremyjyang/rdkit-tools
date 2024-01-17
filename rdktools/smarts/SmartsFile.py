"""
@author Jack Ringer
Date: 1/11/2024
Description:
Class definition for a smarts (".sma" or ".smarts") file.

Acknowledgements:
This code is based on:
https://github.com/unmtransinfo/unmatch_biocomp_smarts/blob/master/src/main/java/edu/unm/health/biocomp/smarts/SmartsFile.java
"""
import re

from rdkit.Chem import MolFromSmarts


class SmartsFile:
    def __init__(self, file_path: str, strict: bool = False):
        """
        Initialize SmartsFile object. Will parse the given file.
        :param str file_path: path to smarts file
        :param bool strict: strictly enforce file correctness
            (throw errors if there are issues with SMARTS), defaults to False
        """
        self.file_path = file_path
        self.strict = strict
        self.smarts_strs = []
        self.smarts_mols = []
        self.smarts_names = []
        self.defines = {}
        self.failed_smarts = []
        self.failed_smarts_raw = []  # for debugging
        self._parse()

    def _resolve_smarts(self, smarts_str: str) -> str:
        """
        Resolve smarts string using defines.
        :param str smarts_str: raw SMART string from file (possibly using definitions)
        :return str: SMART string with definitions resolved
        :raises ValueError: if self.strict and a definition is not found
        """
        pat_tag = re.compile("\\$[a-zA-Z0-9_]+")
        match_tag = pat_tag.search(smarts_str)
        resolved_smarts_str = ""
        tag = ""
        i, j, k = 0, 0, 0
        while match_tag:
            j = match_tag.start()
            k = match_tag.end()
            resolved_smarts_str += smarts_str[i:j]
            tag = smarts_str[j + 1 : k]
            smarts_def = self.defines[tag]
            if smarts_def is not None:
                if pat_tag.search(smarts_def):  # nested define
                    smarts_def = self._resolve_smarts(smarts_def)
                resolved_smarts_str += "$(" + smarts_def + ")"
            elif self.strict:
                raise ValueError(f"Undefined definition: {tag}")
            else:
                return None
            match_tag = pat_tag.search(smarts_str, k)
            i = k
        resolved_smarts_str += smarts_str[i:]
        return resolved_smarts_str

    def _parse(self):
        """
        Parses the smarts file.
        :raises ValueError: if self.strict and a SMARTS cannot be parsed
        """
        # define patterns
        pat_comment = re.compile("^\\s*#.*$", re.IGNORECASE)
        pat_blank = re.compile("^\\s*$", re.IGNORECASE)
        pat_define = re.compile("^define\\s+\\$(\\S+)\\s+(\\S+)\\s*$", re.IGNORECASE)
        pat_smarts_name = re.compile("^(\\S+)\\s+(.+)$", re.IGNORECASE)

        # parse file using patterns above
        with open(self.file_path, "r") as f:
            for line in f:
                raw_smarts_str = ""
                smarts_name = ""
                match_comment = pat_comment.match(line)
                match_blank = pat_blank.match(line)
                match_define = pat_define.match(line)
                match_smarts_name = pat_smarts_name.match(line)
                if match_comment or match_blank:
                    continue
                elif match_define:
                    tag = match_define.group(1)
                    val = match_define.group(2)
                    self.defines[tag] = val
                    continue
                elif match_smarts_name:
                    raw_smarts_str = match_smarts_name.group(1)
                    smarts_name = match_smarts_name.group(2).strip()
                else:
                    raw_smarts_str = line.strip()
                resolved_smarts_str = self._resolve_smarts(raw_smarts_str)
                if resolved_smarts_str is None:
                    self.failed_smarts.append("")
                    self.failed_smarts_raw.append(raw_smarts_str)
                    continue
                smarts_mol = MolFromSmarts(resolved_smarts_str)
                if smarts_mol is None:
                    if self.strict:
                        raise ValueError(f"Invalid SMARTS: {resolved_smarts_str}")
                    else:
                        self.failed_smarts.append(resolved_smarts_str)
                        self.failed_smarts_raw.append(raw_smarts_str)
                        continue
                smarts_mol.SetProp("Name", smarts_name)  # will just be "" if not given
                self.smarts_strs.append(resolved_smarts_str)
                self.smarts_mols.append(smarts_mol)
                self.smarts_names.append(smarts_name)


if __name__ == "__main__":
    fpath = "data/pains_unm.sma"
    sf = SmartsFile(fpath, False)
    print("*" * 80)
    for i in range(len(sf.failed_smarts)):
        print(f"{sf.failed_smarts_raw[i]}")
        print(f"{sf.failed_smarts[i]}")
        print()
