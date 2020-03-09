See: 
https://rdkit.readthedocs.org/en/latest/Cartridge.html#reference-guide
https://github.com/rdkit/rdkit/issues/1762

Types:
	mol : rdkit molecule
	qmol : rdkit query molecule
	sfp  : sparse fp
	bfp : bitvector fp

Parameters:

	rdkit.tanimoto_threshold : threshold value for the Tanimoto similarity
	rdkit.dice_threshold : threshold value for the Dice similiarty
	rdkit.do_chiral_sss : whether or not stereo used in substructure matching

Operators:

  Similarity search

    % : operator used for similarity searches using Tanimoto similarity. Returns whether or not the Tanimoto
    % similarity between two fps (either two sfp or two bfp values) exceeds
    % rdkit.tanimoto_threshold.
    # : operator used for similarity searches using Dice similarity. Returns whether or not the Dice
    # similarity between two fps (either two sfp or two bfp values) exceeds rdkit.dice_threshold.
    <%> : used for Tanimoto KNN searches (to return ordered lists of neighbors).
    <#> : used for Dice KNN searches (to return ordered lists of neighbors).

  Substructure and exact structure search

    @> : substructure operator. mol on the right is substructure query
    <@ : substructure operator. mol on the left is substructure query
    @= : whether or not two molecules are the same.

  Molecule comparison

    < : returns whether or not the left mol is less than the right mol
    > : returns whether or not the left mol is greater than the right mol
    = : returns whether or not the left mol is equal to the right mol
    <= : returns whether or not the left mol is less than or equal to the right mol
    >= : returns whether or not the left mol is greater than or equal to the right mol

Functions
  Fingerprint Related
    Generating fingerprints

    morgan_fp(mol,int default 2): sfp, count-based; second argument radius; ECFP-like.
    morganbv_fp(mol,int default 2): bfp, bitvector; second argument radius; ECFP-like.
    featmorgan_fp(mol,int default 2): sfp, count-based, chemical-feature-based, second argument radius; FCFP-like.
    featmorganbv_fp(mol,int default 2): bfp, bitvector; chemical-feature-based, second argument radius; FCFP-like.
    rdkit_fp(mol): bfp, RDKit fp; daylight-like using hashed molecular subgraphs.
    atompair_fp(mol): sfp, count-based atom-pair.
    atompairbv_fp(mol): bfp, bitvector atom-pair.
    torsion_fp(mol): sfp, count-based topological-torsion .
    torsionbv_fp(mol): bfp, bitvector topological-torsion .
    layered_fp(mol): bfp, layered, experimental substructure using hashed molecular subgraphs.
    maccs_fp(mol): bfp, MACCS fp.


    Working with fingerprints

    tanimoto_sml(fp,fp): Tanimoto similarity between two fps, same type (sfp or bfp).
    dice_sml(fp,fp): Dice similarity between two fps, same type (sfp or bfp).
    size(bfp): length of (number of bits) a bfp.
    add(sfp,sfp): returns an sfp formed by the element-wise addition of the two sfps.
    subtract(sfp,sfp): returns an sfp formed by the element-wise subtraction of the two sfps.
    all_values_lt(sfp,int): whether or not all elements of the sfp argument are less than the int argument.
    all_values_gt(sfp,int): whether or not all elements of the sfp argument are greater than the int argument.

Fingerprint I/O

    bfp_to_binary_text(bfp): returns a bytea with the binary string of the fp that can be converted back into an RDKit fp in other software. 
    bfp_from_binary_text(bytea): constructs a bfp from a binary string of the fp.


Molecule Related
Molecule I/O and Validation

    is_valid_smiles(smiles) : returns whether or not a SMILES string produces a valid RDKit molecule.
    is_valid_ctab(ctab) : returns whether or not a CTAB (mol block) string produces a valid RDKit molecule.
    is_valid_smarts(smarts) : returns whether or not a SMARTS string produces a valid RDKit molecule.
    is_valid_mol_pkl(bytea) : returns whether or not a binary string (bytea) can be converted into an RDKit molecule. 
    mol_from_smiles(smiles) : returns a molecule for a SMILES string, NULL if the molecule construction fails.
    mol_from_smarts(smarts) : returns a molecule for a SMARTS string, NULL if the molecule construction fails.
    mol_from_ctab(ctab, bool default false) : returns a molecule for a CTAB (mol block) string, NULL if the molecule construction fails. The optional second argument controls whether or not the molecule’s coordinates are saved.
    mol_from_pkl(bytea) : returns a molecule for a binary string (bytea), NULL if the molecule construction fails. 
    mol_to_smiles(mol) : returns the canonical SMILES for a molecule.
    mol_to_smarts(mol) : returns SMARTS string for a molecule.
    mol_to_pkl(mol) : returns binary string (bytea) for a molecule.
    mol_to_ctab(mol,bool default true) : returns a CTAB (mol block) string for a molecule. The optional second argument controls whether or not 2D coordinates will be generated for molecules that don’t have coordinates.

Substructure operations

    substruct(mol,mol) : returns whether or not the second mol is a substructure of the first.
    substruct_count(mol,mol,bool default true) : returns the number of substructure matches between the second molecule and the first. The third argument toggles whether or not the matches are uniquified. 

Descriptors

    mol_amw(mol) : returns the AMW for a molecule.
    mol_logp(mol) : returns the MolLogP for a molecule.
    mol_tpsa(mol) : returns the topological polar surface area for a molecule 
    mol_fractioncsp3(mol) : returns the fraction of carbons that are sp3 hybridized 
    mol_hba(mol) : returns the number of Lipinski H-bond acceptors (i.e. number of Os and Ns) for a molecule.
    mol_hbd(mol) : returns the number of Lipinski H-bond donors (i.e. number of Os and Ns that have at least one H) for a molecule.
    mol_numatoms(mol) : returns the total number of atoms in a molecule.
    mol_numheavyatoms(mol) : returns the number of heavy atoms in a molecule.
    mol_numrotatablebonds(mol) : returns the number of rotatable bonds in a molecule 
    mol_numheteroatoms(mol) : returns the number of heteroatoms in a molecule 
    mol_numrings(mol) : returns the number of rings in a molecule 
    mol_numaromaticrings(mol) : returns the number of aromatic rings in a molecule
    mol_numaliphaticrings(mol) : returns the number of aliphatic (at least one non-aromatic bond) rings in a molecule 
    mol_numsaturatedrings(mol) : returns the number of saturated rings in a molecule 
    mol_numaromaticheterocycles(mol) : returns the number of aromatic heterocycles in a molecule 
    mol_numaliphaticheterocycles(mol) : returns the number of aliphatic (at least one non-aromatic bond) heterocycles in a molecule 
    mol_numsaturatedheterocycles(mol) : returns the number of saturated heterocycles in a molecule
    mol_numaromaticcarbocycles(mol) : returns the number of aromatic carbocycles in a molecule 
    mol_numaliphaticcarbocycles(mol) : returns the number of aliphatic (at least one non-aromatic bond) carbocycles in a molecule 
    mol_numsaturatedcarbocycles(mol) : returns the number of saturated carbocycles in a molecule 
    mol_inchi(mol) : returns an InChI for the molecule. (available from the 2011_06 release, requires that the RDKit be built with InChI support).
    mol_inchikey(mol) : returns an InChI key for the molecule. (available from the 2011_06 release, requires that the RDKit be built with InChI support).
    mol_formula(mol,bool default false, bool default true) : returns a string with the molecular formula.  The second argument controls whether isotope information is included in the formula; the third argument controls whether “D” and “T” are used instead of [2H] and [3H].

Connectivity Descriptors

    mol_chi0v(mol) - mol_chi4v(mol) : returns the ChiXv value for a molecule for X=0-4 
    mol_chi0n(mol) - mol_chi4n(mol) : returns the ChiXn value for a molecule for X=0-4 
    mol_kappa1(mol) - mol_kappa3(mol) : returns the kappaX value for a molecule for X=1-3 

Other

    rdkit_version() : returns a string with the cartridge version number.

