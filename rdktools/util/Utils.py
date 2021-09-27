#!/usr/bin/env python3
"""
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.MurckoScaffold.html
https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.rdScaffoldNetwork.html
rdScaffoldNetwork available RDKit 2020.03.1+.
"""
#############################################################################
import os,sys,re,logging,json,time,inspect

import matplotlib as mpl
#from matplotlib import pyplot as plt

import rdkit
import rdkit.Chem
import rdkit.Chem.AllChem
from rdkit.Chem import SmilesMolSupplier, SDMolSupplier, SDWriter, SmilesWriter, MolStandardize, MolToSmiles, MolFromSmiles, Draw, rdDepictor

# fingerprints:
from rdkit.Chem import RDKFingerprint, PatternFingerprint, LayeredFingerprint, LayeredFingerprint_substructLayers
from rdkit.Chem.MACCSkeys import GenMACCSKeys
from rdkit.Chem.AllChem import GetMorganFingerprint,GetMorganFingerprintAsBitVect

from rdkit.Chem.Draw import rdMolDraw2D

import pyvis
from pyvis.network import Network

from rdkit.Chem.Draw import rdMolDraw2D

DEMOSMIS = [
'C[C@H]1CN(CCCN1[S+]([O-])(=O)C2=CC=CC3=C2C(=CN=C3)C)C(=O)CN Rho Kinase Inhibitor IV',
'N[S+]([O-])(=O)C1=C(Cl)C=C2NC(N[S+]([O-])(=O)C2=C1)C(Cl)Cl trichlormethiazide',
'C[S+]([O-])C1=CC=C(C=C1)\C=C2\C(=C(\CC(O)=O)C3=C2C=CC(=C3)F)C Sulindac',
'O=[N+]([O-])c1cc(S(=O)(=O)N2CCCC2)ccc1NN=Cc1ccccn1 PUBCHEM_CID:4036736',
'Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12 flucloxacillin',
'CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O ampicillin',
'CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O.[Na] penicillin',
'Cc1onc(-c2ccccc2)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12 oxacillin'
	]

#############################################################################
def moltosvg(mol, molSize=(450,250), kekulize=True):
    mc = rdMolDraw2D.PrepareMolForDrawing(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
    opts = drawer.drawOptions()
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg.replace('svg:', '')

#############################################################################
def ReadMols(molReader):
  mols=[]; 
  for mol in molReader:
    if mol is not None: mols.append(mol)
  logging.info(f'{len(mols)} mols read from {molReader}')
  return mols

#############################################################################
def File2Molreader(ifile, idelim, smicol, namcol, iheader):
  if re.sub(r'.*\.', '', ifile).lower()in ('smi', 'smiles', 'csv', 'tsv'):
    molReader = SmilesMolSupplier(ifile, delimiter=idelim, smilesColumn=smicol, nameColumn=namcol, titleLine=iheader, sanitize=True)
  elif re.sub(r'.*\.', '', ifile).lower() in ('sdf','sd','mdl','mol'):
    molReader = SDMolSupplier(ifile, sanitize=True, removeHs=True)
  else:
    molReader = None
    logging.error(f'Invalid file extension: {ifile}')
  return molReader

#############################################################################
def File2Molwriter(ofile, odelim, oheader):
  if not ofile:
    molWriter = SmilesWriter("-", delimiter=odelim, nameHeader='Name', includeHeader=oheader, isomericSmiles=True, kekuleSmiles=False)
  elif re.sub(r'.*\.', '', ofile).lower() in ('sdf','sd','mdl','mol'):
    molWriter = SDWriter(ofile)
  elif re.sub(r'.*\.', '', ofile).lower() in ('smi', 'smiles', 'csv', 'tsv'):
    molWriter = SmilesWriter(ofile, delimiter=odelim, nameHeader='Name', includeHeader=oheader, isomericSmiles=True, kekuleSmiles=False)
  else:
    logging.error(f'Invalid file extension: {ofile}')
    molWriter = None
  return molWriter

#############################################################################
def Mols2FPs_RDK(mols, molWriter=None):
  for molA in mols:
    fpA = RDKFingerprint(molA)
    fpstxt = rdkit.DataStructs.BitVectToFPSText(fpA)
    logging.debug(fpstxt)
    molA.SetProp("_Name", (molA.GetProp("_Name") if molA.HasProp("_Name") else "")+"\t"+fpstxt)
    count=0;
    for i in range(fpA.GetNumBits()):
      if fpA.GetBit(i): count+=1
    logging.debug(f"count={count}; size={fpA.GetNumBits()}")
    if molWriter is not None:
      molWriter.write(molA)

#############################################################################
def Mols2FPs_Morgan(mols, radius=2, nbits=1024, molWriter=None):
  for molA in mols:
    fpA = GetMorganFingerprintAsBitVect(molA, radius, nbits)
    fpstxt = rdkit.DataStructs.BitVectToFPSText(fpA)
    logging.debug(fpstxt)
    molA.SetProp("_Name", (molA.GetProp("_Name") if molA.HasProp("_Name") else "")+"\t"+fpstxt)
    count=0;
    for i in range(fpA.GetNumBits()):
      if fpA.GetBit(i): count+=1
    logging.debug(f"count={count}; size={fpA.GetNumBits()}")
    if molWriter is not None:
      molWriter.write(molA)

#############################################################################
def Mols2FPs_MACCSKeys(mols, molWriter=None):
  for molA in mols:
    fpA = GenMACCSKeys(molA)
    fpstxt = rdkit.DataStructs.BitVectToFPSText(fpA)
    logging.debug(fpstxt)
    molA.SetProp("_Name", (molA.GetProp("_Name") if molA.HasProp("_Name") else "")+"\t"+fpstxt)
    count=0;
    for i in range(fpA.GetNumBits()):
      if fpA.GetBit(i): count+=1
    logging.debug(f"count={count}; size={fpA.GetNumBits()}")
    onbits = tuple(fpA.GetOnBits())
    smarts = {i:rdkit.Chem.MACCSkeys.smartsPatts[i][0] for i in onbits}
    logging.debug(str(smarts))
    if molWriter is not None:
      molWriter.write(molA)

#############################################################################
def GenerateConformations(mol, nconf, ffname, optiters, etol):
  '''verbose=2 means MMFF generates much data, unfortunately to stdout.'''
  verbose = 2 if logging.getLogger().getEffectiveLevel()==logging.DEBUG else 1 if logging.getLogger().getEffectiveLevel()==logging.INFO else 0
  molh = rdkit.Chem.AddHs(mol)
  confIds=rdkit.Chem.AllChem.EmbedMultipleConfs(molh, numConfs=nconf)
  ok_opts = {confId:False for confId in confIds} #optimization status
  molh_0 = rdkit.Chem.Mol(molh)
  i_conf=0;
  for confId in confIds:
    i_conf+=1
    if ffname.upper()=='MMFF':
      ff = rdkit.Chem.AllChem.MMFFGetMoleculeForceField(molh, rdkit.Chem.AllChem.MMFFGetMoleculeProperties(molh, "MMFF94", verbose), confId=confId)
    else:
      ff = rdkit.Chem.AllChem.UFFGetMoleculeForceField(molh, confId=confId)
    try:
      e_i = ff.CalcEnergy()
      rval = ff.Minimize(maxIts=optiters, energyTol=etol)
      ok_opts[confId] = bool(rval==0)
      e_f = ff.CalcEnergy()
    except Exception as e:
      logging.info('%s'%str(e))
      e_i,e_f,ok_opts[confId] = None,None,False
    try:
      rmsd = rdkit.Chem.AllChem.AlignMol(molh, molh_0, prbCid=confId, refCid=confId)
    except Exception as e:
      rmsd = 0.0
      logging.info('%s'%str(e))
    logging.debug('\t%3d. RMSD = %.2f ; Ei = %s, Ef = %s, converged = %s'%(i_conf, rmsd, ('%.2f'%e_i if e_i else 'None'), ('%.2f'%e_f if e_f else 'None'), str(ok_opts[confId]))) 
  if not confIds:
    logging.info('ERROR: EmbedMultipleConfs() failed.') 

  return molh, list(confIds)

#############################################################################
