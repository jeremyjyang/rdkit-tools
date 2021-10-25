#!/usr/bin/env python3
#
import sys,os,logging,tempfile

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.EState.EState import EStateIndices
from rdkit.Chem.Lipinski import *
from rdkit.Chem.Descriptors import *
from rdkit.Chem.Descriptors3D import *
from rdkit.Chem.rdFreeSASA import *

#############################################################################
def CalcCrippenLogP(molReader, molWriter):
  i_mol=0
  for mol in molReader:
    i_mol+=1
    name = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    logging.info(f"{i_mol}. {name}")
    logp = MolLogP(mol)
    mol.SetProp('WildmanCrippenLogP', f"{logp:.3f}")
    if i_mol==1:
      molWriter.SetProps(mol.GetPropNames())
    molWriter.write(mol)
  logging.info(f"n_out: {i_mol}")

#############################################################################
def CalcEStateIndices(molReader, molWriter):
  """
Calculate Kier-Hall electrotopological descriptors
  [x] smiles or molfile input and output with data
  [ ] stdout and stdin
  [ ] sumDeltaI -- can it be done with RDKit?
 
References:
  1. L.H. Hall, B. Mohney and L.B. Kier., "The electrotopological state:
     structure information at the atomic level for molecular graphs",
     JCICS 31 (1) 76-82 (1991)
  2. L.B. Kier and L.H. Hall _Molecular Structure Description:
     The Electrotopological State"_  Academic Press (1999)
  """

  i_mol=0
  for mol in molReader:
    i_mol+=1
    esvals = EStateIndices(mol, force=1)
    mol.SetProp('RDKit_EStates', (','.join(map(lambda x:(f"{x:.3f}"), esvals))))
    mol.SetProp('RDKit_MaxEState', (f"{max(esvals):.3f}"))
    if i_mol==1:
      molWriter.SetProps(mol.GetPropNames())
    molWriter.write(mol)
  logging.info(f"n_out: {i_mol}")

#############################################################################
def CalcLipinski(molReader, molWriter):
  i_mol=0
  for mol in molReader:
    i_mol+=1
    name = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    logging.info(f"{i_mol}. {name}")

    mol.SetProp("FractionCSP3", f"{FractionCSP3(mol):.3f}") # fraction of Cs SP3-hybridized
    mol.SetProp("HeavyAtomCount", f"{HeavyAtomCount(mol)}") # heavy atoms 
    mol.SetProp("NHOHCount", f"{NHOHCount(mol)}") # NHs or OHs
    mol.SetProp("NOCount", f"{NOCount(mol)}") # Ns and Os
    mol.SetProp("NumAliphaticCarbocycles", f"{NumAliphaticCarbocycles(mol)}") # aliphatic (containing at least one non-aromatic bond) carbocycles
    mol.SetProp("NumAliphaticHeterocycles", f"{NumAliphaticHeterocycles(mol)}") # aliphatic (containing at least one non-aromatic bond) heterocycles
    mol.SetProp("NumAliphaticRings", f"{NumAliphaticRings(mol)}") # aliphatic (containing at least one non-aromatic bond) rings
    mol.SetProp("NumAromaticCarbocycles", f"{NumAromaticCarbocycles(mol)}") # aromatic carbocycles
    mol.SetProp("NumAromaticHeterocycles", f"{NumAromaticHeterocycles(mol)}") # aromatic heterocycles
    mol.SetProp("NumAromaticRings", f"{NumAromaticRings(mol)}") # aromatic rings
    mol.SetProp("NumHAcceptors", f"{NumHAcceptors(mol)}") # H Bond Acceptors
    mol.SetProp("NumHDonors", f"{NumHDonors(mol)}") # H Bond Donors
    mol.SetProp("NumHeteroatoms", f"{NumHeteroatoms(mol)}") # Heteroatoms
    mol.SetProp("NumRotatableBonds", f"{NumRotatableBonds(mol)}") # Rotatable Bonds
    mol.SetProp("NumSaturatedCarbocycles", f"{NumSaturatedCarbocycles(mol)}") # saturated carbocycles
    mol.SetProp("NumSaturatedHeterocycles", f"{NumSaturatedHeterocycles(mol)}") # saturated heterocycles
    mol.SetProp("NumSaturatedRings", f"{NumSaturatedRings(mol)}") # saturated rings
    mol.SetProp("RingCount", f"{RingCount(mol)}") #
  
    if i_mol==1:
      molWriter.SetProps(mol.GetPropNames())
    molWriter.write(mol)
  logging.info(f"n_out: {i_mol}")
  
#############################################################################
def CalcDescriptors(molReader, molWriter):
  i_mol=0
  for mol in molReader:
    i_mol+=1
    name = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    logging.info(f"{i_mol}. {name}")

    mol.SetProp("MolWt", f"{MolWt(mol):.3f}") # The average molecular weight of the molecule
    mol.SetProp("ExactMolWt", f"{ExactMolWt(mol):.3f}") # The exact molecular weight of the molecule
    mol.SetProp("HeavyAtomMolWt", f"{HeavyAtomMolWt(mol):.3f}")
    mol.SetProp("MaxAbsPartialCharge", f"{MaxAbsPartialCharge(mol):.3f}")
    mol.SetProp("MaxPartialCharge", f"{MaxPartialCharge(mol):.3f}")
    mol.SetProp("MinAbsPartialCharge", f"{MinAbsPartialCharge(mol):.3f}")
    mol.SetProp("MinPartialCharge", f"{MinPartialCharge(mol):.3f}")
    mol.SetProp("NumRadicalElectrons", f"{NumRadicalElectrons(mol)}") # The number of radical electrons the molecule has (says nothing about spin state)
    mol.SetProp("NumValenceElectrons", f"{NumValenceElectrons(mol)}") # The number of valence electrons the molecule has
  
    if i_mol==1:
      molWriter.SetProps(mol.GetPropNames())
    molWriter.write(mol)
  logging.info(f"n_out: {i_mol}")

#############################################################################
def CalcDescriptors3D(molReader, molWriter):
  """Requires 3D conformations."""
  i_mol=0
  for mol in molReader:
    i_mol+=1
    name = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    logging.info(f"{i_mol}. {name}")


  # molecular asphericity
  # from Todeschini and Consoni “Descriptors from Molecular Geometry” Handbook of Chemoinformatics https://doi.org/10.1002/9783527618279.ch37
  # Definition: 0.5 * ((pm3-pm2)**2 + (pm3-pm1)**2 + (pm2-pm1)**2)/(pm1**2+pm2**2+pm3**2)
    mol.SetProp("Asphericity", f"{Asphericity(mol):.3f}")

  # molecular eccentricity
  # from Todeschini and Consoni “Descriptors from Molecular Geometry” Handbook of Chemoinformatics https://doi.org/10.1002/9783527618279.ch37
  # Definition: sqrt(pm3**2 -pm1**2) / pm3**2
    mol.SetProp("Eccentricity", f"{Eccentricity(mol):.3f}")

  # Inertial shape factor
  # from Todeschini and Consoni “Descriptors from Molecular Geometry” Handbook of Chemoinformatics https://doi.org/10.1002/9783527618279.ch37
  # Definition: pm2 / (pm1*pm3)
    mol.SetProp("InertialShapeFactor", f"{InertialShapeFactor(mol):.3f}")

  # Normalized principal moments ratio 1 (=I1/I3)
  # from Sauer and Schwarz JCIM 43:987-1003 (2003) https://doi.org/10.1021/ci025599w
    mol.SetProp("NPR1", f"{NPR1(mol):.3f}")

  # Normalized principal moments ratio 2 (=I2/I3)
  # from Sauer and Schwarz JCIM 43:987-1003 (2003) https://doi.org/10.1021/ci025599w
    mol.SetProp("NPR2", f"{NPR2(mol):.3f}")

  # First (smallest) principal moment of inertia
    mol.SetProp("PMI1", f"{PMI1(mol):.3f}")

  # Second principal moment of inertia
    mol.SetProp("PMI2", f"{PMI2(mol):.3f}")

  # Third (largest) principal moment of inertia
    mol.SetProp("PMI3", f"{PMI3(mol):.3f}")

  # Radius of gyration
  # from Todeschini and Consoni “Descriptors from Molecular Geometry” Handbook of Chemoinformatics https://doi.org/10.1002/9783527618279.ch37
  # Definition: for planar molecules: sqrt( sqrt(pm3*pm2)/MW ) for nonplanar molecules: sqrt( 2*pi*pow(pm3*pm2*pm1,1/3)/MW )
    mol.SetProp("RadiusOfGyration", f"{RadiusOfGyration(mol):.3f}")

  # Molecular spherocityIndex
  # from Todeschini and Consoni “Descriptors from Molecular Geometry” Handbook of Chemoinformatics https://doi.org/10.1002/9783527618279.ch37
  # Definition: 3 * pm1 / (pm1+pm2+pm3) where the moments are calculated without weights
    mol.SetProp("SpherocityIndex", f"{SpherocityIndex(mol):.3f}")

    if i_mol==1:
      molWriter.SetProps(mol.GetPropNames())
    molWriter.write(mol)
  logging.info(f"n_out: {i_mol}")

#############################################################################
def CalcFreeSASA(molReader, molWriter):
  """Requires 3D conformations."""
  i_mol=0
  sasa_opts = SASAOpts()
  atomQuery = MakeFreeSasaPolarAtomQuery()
  #atomQuery = MakeFreeSasaAPolarAtomQuery()
  sasa_opts.algorithm = SASAAlgorithm.LeeRichards
  #sasa_opts.algorithm = SASAAlgorithm.ShrakeRupley
  #sasa_opts.classifier = SASAClassifier.NACCESS
  sasa_opts.classifier = SASAClassifier.OONS
  #sasa_opts.classifier = SASAClassifier.Protor
  sasa_opts.probeRadius = 1.4 #Angstroms

  for mol in molReader:
    i_mol+=1
    name = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
    radii = classifyAtoms(mol, sasa_opts)
    val = CalcSASA(mol,  radii, -1, atomQuery, sasa_opts)
    mol.SetProp("SASA", f"{val:.3f}")
    logging.info(f"{i_mol}. {name}; SASA: {val:f}")
    if i_mol==1:
      molWriter.SetProps(mol.GetPropNames())
    molWriter.write(mol)
  logging.info(f"n_out: {i_mol}")

#############################################################################
def Demo():
  molWriter = Chem.SDWriter("-")
  f = tempfile.NamedTemporaryFile(delete=False)
  f.write(b"""\
dopamine
  -OEChem-03200608123D

 22 22  0     0  0  0  0  0  0999 V2000
    0.1459   -1.4217    0.2165 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2188   -1.5581   -0.0373 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1068    0.9801    0.2844 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7020   -0.1525    0.3774 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0276   -0.4254   -0.1304 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4715    0.8437    0.0305 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1611   -0.0069    0.6483 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0113    0.1181   -0.6221 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4148    0.2369   -0.3309 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3588   -0.5713   -0.3787 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2515    1.9570   -0.0577 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7685   -2.3096    0.2875 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6457   -2.5496   -0.1618 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.3244    1.9703    0.4094 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.5156   -0.8799    1.2219 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.3505    0.8630    1.3007 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.8935   -0.7696   -1.2516 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.7206    0.9947   -1.2103 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.7033    0.7051    0.5132 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.0852    0.0389   -1.0565 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5069   -1.3096   -0.9904 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1215    1.7289   -0.4214 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  1  4  1  0  0  0  0
  2  5  1  0  0  0  0
  3  4  2  0  0  0  0
  3  6  1  0  0  0  0
  4  7  1  0  0  0  0
  5  6  2  0  0  0  0
  5 10  1  0  0  0  0
  6 11  1  0  0  0  0
  7  8  1  0  0  0  0
  8  9  1  0  0  0  0
  1 12  1  0  0  0  0
  2 13  1  0  0  0  0
  3 14  1  0  0  0  0
  7 15  1  0  0  0  0
  7 16  1  0  0  0  0
  8 17  1  0  0  0  0
  8 18  1  0  0  0  0
  9 19  1  0  0  0  0
  9 20  1  0  0  0  0
 10 21  1  0  0  0  0
 11 22  1  0  0  0  0
M  END
$$$$
""")
  f.close()

  molReader = Chem.SDMolSupplier(f.name)
  CalcCrippenLogP(molReader, molWriter)

  molReader = Chem.SDMolSupplier(f.name)
  CalcEStateIndices(molReader, molWriter)

  molReader = Chem.SDMolSupplier(f.name)
  CalcLipinski(molReader, molWriter)

  molReader = Chem.SDMolSupplier(f.name)
  CalcDescriptors(molReader, molWriter)

  molReader = Chem.SDMolSupplier(f.name)
  CalcDescriptors3D(molReader, molWriter)

  molReader = Chem.SDMolSupplier(f.name)
  CalcFreeSASA(molReader, molWriter)

  os.unlink(f.name)

#############################################################################
