#!/usr/bin/env python3
#############################################################################
import os,sys,re,logging
import rdkit.Chem
import rdkit.Chem.AllChem

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
