#!/usr/bin/env python
#
from rdkit import Chem
from rdkit.Chem import MACCSkeys

for i,val in MACCSkeys.smartsPatts.items():
  smarts,n = val
  print i, smarts, n
