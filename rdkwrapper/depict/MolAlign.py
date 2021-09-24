#!/usr/bin/env python3
#
import sys,os,logging

from rdkit import Chem
from rdkit.Chem import AllChem
#from rdkit.Chem.rdMolAlign import AlignMol  ##3D
from rdkit.Chem.ChemUtils.AlignDepict import AlignDepict 

from rdkit.Chem import rdDepictor 
from rdkit import Geometry 

#############################################################################
def MyAlignDepict(mol,core,corepat=None,acceptFailure=False):
  if core and corepat: 
    coreMatch = core.GetSubstructMatch(corepat) 
    if not coreMatch and not acceptFailure: 
      raise ValueError("Core smarts does not map to core.")
  else: 
    coreMatch = range(core.GetNumAtoms(onlyHeavy=True)) 
  if corepat: 
    match = mol.GetSubstructMatch(corepat) 
  else: 
    match = mol.GetSubstructMatch(core) 
   
  if not match: 
    if not acceptFailure: 
      raise ValueError("Substructure match with core not found.")
    else: 
      coordMap={} 
  else: 
    conf = core.GetConformer() 
    coordMap={} 
    for i,idx in enumerate(match): 
      pt3 = conf.GetAtomPosition(coreMatch[i]) 
      pt2 = Geometry.Point2D(pt3.x, pt3.y) 
      coordMap[idx]=pt2 
  rdDepictor.Compute2DCoords(mol, clearConfs=True, coordMap=coordMap) 


#############################################################################
if __name__=='__main__':
  fpath_in=sys.argv[1]
  fpath_out=sys.argv[2]
  smarts=sys.argv[3]

  corepat=Chem.MolFromSmarts(smarts)
  if not corepat:
    logging.info('Bad smarts: %s'%(smarts))
    sys.exit()

  sdreader = Chem.SDMolSupplier(fpath_in, sanitize=True, removeHs=True)

  sdwriter=Chem.SDWriter(fpath_out)

  i_mol=0
  i_fail=0
  core = sdreader.next()
  sdreader.reset()
  for mol in sdreader:
    i_mol+=1
    logging.info('%d. %s' %(i_mol, mol.GetProp('_Name')))

    try:
      MyAlignDepict(mol, core, corepat, acceptFailure=False)
    except ValueError:
      i_fail+=1
      logging.info('%d. alignment failed.' %(i_mol))

    sdwriter.write(mol)

  logging.info('%d mols written to %s' %(i_mol, fpath_out))
