#############################################################################
# "New Substructure Filters for Removal of Pan Assay Interference Compounds (PAINS) from Screening
# Libraries and for their Exclusion in Bioassays", J. B. Baell and G. A. Holloway,
# J. Med. Chem, 2010, Vol. 53, No 7, 2719-2740.
#############################################################################
# Translated to SMARTS by: Jeremy Yang ; 30 Nov 2010
#############################################################################
# TABLE S1. Functional Group Filters (#) Used For the WEHI 93K HTS
# Library to Remove Inappropriate Compounds. Subsequent Groups Used for
# the CTX 136K HTS Library are double hashed (##). These are currently
# blocked out but it is recommended that they be activated (unblocked)
# prior to the purchase of new libraries as a more expansive set
# of filters
#############################################################################
##Many definitions are dependent on others (e.g. if allow sulfonates, others defs need modification) 
#
define $Hal [F,Cl,Br,I]
define $Het [!#6&!#1]
define $Hev [!#1]
#
C(=O)[$Hal]	acid halide and related
C(=S)[$Hal] 	acid halide and related
#
[$Het][$Hal]	any het halide (includes sulfonyl halides) 
#
O=[C,S]C=[CH2]	acrylamide and related, including vinylpyridines etc
S=CC=[CH2] 	acrylamide and related, including vinylpyridines etc
N#CC=[CH2] 	acrylamide and related, including vinylpyridines etc
C(=[CH2])c:1:n:[$Hev]:[$Hev]:[$Hev]:[$Hev]1 
#
##acyl carbamate and related hot linear aryl carbonates/carbamates
C(=O)NC(=O)*[O,S]	acyl carbamate and related hot linear aryl carbonates/carbamates
[$Het]C(=O)O[!R]c:1:[$Hev](*[!O;!SX2;!CX4;!NX3]):[$Hev]:[$Hev](*[$Hal,$(C#N),$(C(F)(F)F),$(S(=O)=O),$(C=O);!$(C(=O)[OH1])]):[$Hev]:[$Hev](*[!O;!SX2;!CX4;!NX3])1	acyl carbamate and related hot linear aryl carbonates/carbamates (1)
[$Het]C(=O)O[!R]c:1:[$Hev](*[$Hal,$(C#N),$(C(F)(F)F),$(S(=O)=O),$(C(=O));!$(C(=O)[OH1])]):[$Hev]:[$Hev](*[!O;!*X2;!CX4;!NX3]):[$Hev]:[$Hev](*[!O;!SX2;!CX4;!NX3])1	acyl carbamate and related hot linear aryl carbonates/carbamates (2)
##acyl hydrazone (and oxime) not in ring: 
##leave in (as well as carbazides, semicarbazides, 
##acylsemicarbazides, hydrazones, acyl hydrazides) 
##despite some reservation of Metal chelating for some of these 
O=C[$Het]N=[!R]C(C)[$(*C),$(*H)]	acyl hydrazone (and oxime) not in ring
O=C[NH1][!R][NH1]C(=O)[$Het]	acyl hydrazide not in ring leave in
O=CN(C=O)C=O ##acyl imide taken care of by imide filter 
O=C[!R]N(C)S(=O)(=O)C ##Nalkyl acyl sulfonamide 
[C,S](=O)C#N #acyl and sulfo cyanide 
#
#C(=S)C#N taken care of by thiocarbonyl filter 
*[!O&!N]C(*[!O&!N])=[!R]NC=O	acylimines 
[CH1]=O	aldehyde 
##C(=S)H taken care of by thiocarbonyl filter 
C[!RH2][CH2][CH2]C[!RH2][CH2][CH3]	aliphatic long chain 
[Br,Cl][CX4H1]	alkyl halide 
N1(C=NN=N1)C=O	amidotetrazole
[$Het][CX4][!R][$Het]	aminals etc not in ring unless basic amine
[$Het][CX4]N([$(*H),$(*[CX4])])[$(*H),$(*[CX4])]	aminals etc not in ring unless basic amine
C(=*[O,S])*[O,S]C=*[O,S]	anhydride/thioanhydrides
##anilines not removed; medchem issue: specifics dealt with in frequent hitter filter
##azepanes not removed
N~[NX1]	azido/diazonium/diazo specifics also in frequent hitters:
N[!R]=N	azo (and carbazones etc)
#N1CC1 ... made redundant
[$Hev]1[$Het]C1	aziridines, epoxides etc
##benzidine not removed; medchem issue
[$Hev]~1C(=O)[$Het]C1	betalactams/lactones, not thietanes, oxetanes, azetidines etc
##beta-keto in diketos
##beta-amino ketones and related
CC[!R](=O)[CH2][CH2][$Het]
CS[!R](=O)[CH2][CH2][$Het]
S1[CH2][CH1]2[CH1]([CH1]1[CH2][CH2][CH2][CH2]C=O)[NH1]C([NH1]2)=O	biotin not now
B	boron
NC(=O)[OH1]	carbamic acids
[NH1]C(=O)[!R][NH1][NH1]	carbazide leave in
[C+1,S+1,O+1,OX3,$Hal&+1,$(S(*)(*)(*)(*)*)]	carbocation/anion and other wrongly charged atoms
c:1:c(:c:c:c(:c1)[OH1])[OH1]	catechol
##crown ethers main ones not now
[$Het]C~C[$Het]C~C[$Het]C~C[$Het]C~C[$Het]C~C[$Het]C~C	chromones reactive; in Michael acceptors
##coumarins reactive; in Michael Acceptors; other specific ones in frequent hitters
[CH2]NC#N	cyanamides
[O,S]C#N	cyanate/thiocyanate
N#CC[OH1]	cyanohydrin
##cycloalkanes not removed
[CH1]1C=CC=C[CH1]1	cyclohexadienes
[CH1]1C=C[CH2]C=C1	cyclohexadienes
C#CC#C	dialkynes
##diaminobenzenes specific ones in frequent hitters
##diketo (includes oquinones) and betadiketo and related hot keto
CC(=[$Het][!R])C(=[$Het][!R])C
CC(=O)C(=O)*[O,N]
CC(=[$Het][!R])S(=O)(=O)
#CC(=O)C(=S)C covered by thiocarbonyl
CC(=O)[CX4]C(=O)C
CC(=O)[CX4]S(=O)C
CC(=O)[CX4]C#N
CC(=O)[CX4]C([$Hal])([$Hal])[$Hal]
[SX2&!R][SX2&!R]	disulfide and related
SSS	disulfide and related
##enamines in frequent hitters
C(=O)O[!R][$(n:[$Hev]),$(NC=O)]	ester of hobt and su (and carbamates/carbonates) A
CC(=S)O[!R]C	esters thiotype; some accounted for by thiocarbonyl removal
CC(=O)S[!R]C	esters thiotype; some accounted for by thiocarbonyl removal
C(=S)S	esters thiotype; some accounted for by thiocarbonyl removal
CC(=O)O[!R][CX3]	esters aryl B
[CX4]C(=O)O[!R][CH2]c:1:[$Hev](*[!O;!SX2;!CX4;!NX3]):[$Hev]:[$Hev](*[$Hal,$(C#N),$(C(F)(F)F),$(S(=O)=O),$(C(=O)[$Hev][!OH1])]):[$Hev]:[$Hev](*[!O;!SX2;!CX4;!NX3])1	ester hot benzyl, non aryl C
[CX4]C(=O)O[!R][CH2]c:1:[$Hev](*[$Hal,$(C#N),$(C(F)(F)F),$(S(=O)=O),$(C(=O)[$Hev][!OH1])]):[$Hev]:[$Hev](*[!O;!SX2;!CX4;!NX3]):[$Hev]:[$Hev](*[!O;!SX2;!CX4;!NX3])1	ester hot benzyl, non aryl C
CC(=O)O[CH1]C([$Hal])([$Hal])[$Hal]	esters other labile D
[CX4]C(=O)O[!R][CH1]C(=O)*[!OH1]	esters other labile D
[CX4]C(=O)O[!R][CH1]C#N	esters other labile D
COC[!R](=O)C(=O)	esters other labile D
COC(=O)S(=O)(=O)	esters other labile D
COC[!R](=O)[CX4]S(=O)=O	esters other labile D
COC[!R](=O)[CX4]C#N	esters other labile D
COC[!R](=O)[CX4]C([$Hal])([$Hal])[$Hal]	esters other labile D
COC[!R](=O)[CX4]N(C(=O))C(=O)	esters other labile D

c:1:c(:c(:[$Hev]:[$Hev]:[$Hev]1)Br)*[Cl,Br]	over halogenated rings; 1,2
c:1:c(:[$Hev]:c(:[$Hev]:[$Hev]1)Br)*[Cl,Br]	over halogenated rings; 1,3
c:1:c(:[$Hev]:[$Hev]:c(:[$Hev]1)Br)*[Cl,Br]	over halogenated rings; 1,4
[$Hev]:1:c(:c(:[$Hev]:c(:[$Hev]1)[$Hal])[$Hal])[$Hal]	over halogenated rings; 1,2,4
c:1(:c(:c(:[$Hev]:[$Hev]:[$Hev]1)[$Hal])[$Hal])[$Hal]	over halogenated rings; 1,2,3
c:1(:c:c(:c:c(:c1)Cl)Cl)Cl	over halogenated rings; 1,3,5
*[O,N][NH2]	hydrazine and other nucleophilic NH2
[OH1][NX3][CX4]	hydroxylamine
c:1:c:c:c(:c(:c1)[OH1])[OH1]	hydroquinone more related ones in frequent hitters
C(=O)N[OH1]	hydroxamic acid
CC(=O)N[!R]C(=O)C	imide
N=C[$Hal]	imidoyl chlorides
*[!N;!O]N=C[!R](*[!N])*[!N]	imines reactive ones
*[!O;!N]N=[!R]C([$(*H),$(*C)])[$(*H),$(*C)]	imines 
*[!O;!N]N=[!R]C([CX4])C	imines 
I	iodine
#isocyanate, isothiocyanate dealt with in ketenes
*=C=*	ketene includes allenes, carbodimides, isocyanates etc
C=1C(NC([CH1]1)=O)=O	maleimides and surrogates
C1([CH2]C(NC1=O)=O)[$Het][CX3]	maleimides and surrogates
C[CH1]=[CH1][!R]C(=O)C	michael acceptors more in frequent hitter file
C[CH1]=C[!R](C#N)C#N	michael acceptors more in frequent hitter file
C[CH1]=C[!R](C(=O)[$Hev]~[$Hev])C#N	michael acceptors more in frequent hitter file
C[CH1]=[!R]C([$(*H),$(*C)])*[$(C#N),$(C(=O));!$(C(=O)*[N,O])]	first, alphabeta unsubstituted, unsaturated nitriles and ketones
*[C,O][CH1]=[!R]C(*[$(S=C(=O)),$(S(=O)),$(C#N),$Hal,$(C([$Hal])([$Hal])[$Hal]);!$(C(=O)[OH1])])*[$(C(=O)),$(S(=O)),$(C#N);!$(C(=O)[OH1])]	second, linear diactivated A
[$Hev][!OH1]C(=O)[CH1]=[!R][CH1]C(=O)[$Hev][!OH1]	linear diactivated B
*[$(C#N),$(C(=O)[$Hev][!OH1]),$(S=O)]CH=[!R]C([$(*H),$(*C)])c:1:[$Hev](*[!O;!SX2;!CX4;!NX3]):[$Hev]:[$Hev](*[$Hal,$(C#N),$(C(F)(F)F),$(S(=O)=O),$(C(=O)[$Hev][!OH1])]):[$Hev]:[$Hev](*[!O;!SX2;!CX4;!NX3])1	linear diactivated C
*[$(C#N),$(C(=O)[$Hev][!OH1]),$(S=O)][CH1]=[!R]C([$(*H),$(*C)])c:1:[$Hev](*[$Hal,$(C#N),$(C(F)(F)F),$(S(=O)=O),$(C(=O)[$Hev][!OH1])]):[$Hev]:[$Hev](*[!O;!SX2;!CX4;!NX3]):[$Hev]:[$Hev](*[!O;!SX2;!CX4;!NX3])1	linear diactivated C
C=1(C(=O)Oc:2:c:c:c:c:c2[CH1]1)*[$(C=O),$(C=S),$(S=O),$(C#N),$Hal,$(C([$Hal])([$Hal])[$Hal]);!$(C(=O)[OH1])]	thirdly, cyclic: activated coumarins; could revisit some of the less activated amides
C1(=[CH1]Oc:2:c:c:c:c:c2C1=O)*[$(C=O),$(C=S),$(S=O),$(C#N),$Hal,$(C([$Hal])([$Hal])[$Hal]);!$(C(=O)[OH1])]	activated chromones
[SX2]c:1:s:c:*:n1	mercaptothiazoles and related
#C[n+1]:[$Hev] redundant by definition below
##Broaden to include Nnitros, N oxides etc: can review
[n+1]:[$Hev]	Nacyl and quaternary alkyl pyridines etc
C[!R](=O)N(:[$Hev]):[$Hev]	Nacyl pyrroles, imidazoles etc not in ring, but not ureas, carbamates
##naphthylamines alpha and beta; leave in; medchem issue
#c:1:c:c:3:c(:c:c1):c:c:c:c3[!R]N
#c:1:c:c:3:c(:c:c1):c:c:c(:c3)[!R]N
#[cH1]:1:[cH1]:c:5:c(:[cH1]:[cH1]1):[cH1]:[cH1]:[cH1]:c5[NH1]
#[cH1]:1:[cH1]:c:5:c(:[cH1]:[cH1]1):[cH1]:[cH1]:[cH1]:c5N([$(*H),$(*C(=O)C)])[$(*H),$(*C(=O)C)]
#[cH1]:1:[cH1]:c:5:c(:[cH1]:[cH1]1):[cH1]:[cH1]:c(:[cH1]5)[NH1]
[$Hev]N[!R](~[OX1])~[OX1]	nitro, including NNO2 etc:
CN[!R](~[OX1])~[OX1]	nitro, including NNO2 etc:
[NX2]=O	nitroso
[$Hev]:[n+1]~[!R]O	N oxide and Nhydroxypyridine
N[OX1]	Noxide: broader
##oxetanes, thietanes (see beta lactams)
###oximes: not removed
##ON single bond: not removed
OO	peroxide
##phenyl carbonate and carbamate hot ones removed above
P	phosphorous
##phthalimides: not removed (except hot ones above in triaryls)
##polycyclic aromatic: not removed
##polyene: not removed
C[N+1](C)(C)C	quaternary nitrogen
##saponins not removed
[Si]	silicon
##stilbenes not removed
S[O!X1]	sulfur oxygen single bond (includes sulfinic/sulfonic acids/esters)
N[SX2]	sulfurnitrogen single bond, exluding sulfonamides
c:1(C(=O)N):c(:c:c:c:c1)SCCC(=O)*[!OH1]	sulfurnitrogen single bond surrogate
C=S	thiocarbonyl
[S&!H0]	thiols
C(Cl)(Cl)Cl	trichloromethyl
C(=O)C[!R]([$Hal])[$Hal]	tri (and di) haloketones and amides (+esters) and related
OC(=[NH1])C([$Hal])([$Hal])[$Hal]	trihaloimidate
[NX2]~[CX1]	isonitrile
##sulfoneactivated rings and related below
*[$([$Hal]),$(S(=O)(=O))]c:1:n:*(*[!N]):*:*(*[!N]):n1	2halo and 2sulfonylpyrimidines, but deactivated triazines OK
c:1(*[!N]):n:*(*[!N]):*:*(*[$([$Hal]),$(S(=O)(=O))]):n1	4halo and 4sulfonylpyrimidines not too deactivated
c:1:n:*(*[$Hal,$(S(=O)(=O)C)]):*:*(*[$Hal,$(C(C)=NO),$(C#N),$(C=O),$(C(F)(F)F),$(S(=O)=O);!$(C(=O)[OH1])]):*1	4activated 2halo and 2sulfonyl pyridines
#
c:1:n:*(*[$Hal,$(C(C)=NO),$(C#N),$(C=O),$(C(F)(F)F),$(S(=O)=O);!$(C(=O)[OH1])]):*:*(*[$Hal,$(S(=O)(=O)C)]):*1	2activated 4halo and 4sulfonyl pyridines
#
c:1:n:*(*[$Hal,$(S(=O)(=O)C)]):*(*[$Hal,$(C(C)=NO),$(C#N),$(C=O),$(C(F)(F)F),$(S(=O)=O);!$(C(=O)[OH1])]):*:*1	3activated 2halo and 2sulfonyl pyridines
#
c:1:n:*(*[$Hal,$(S(=O)(=O)C)]):*:*:*(*[$Hal,$(C(C)=NO),$(C#N),$(C=O),$(C(F)(F)F),$(S(=O)=O);!$(C(=O)[OH1])])1	5activated 2halo and 2sulfonyl pyridines
#
c:1(*[$Hal,$(C(C)=NO),$(C#N),$(C(=O)),$(C(F)(F)F),$(S(=O)=O);!$(C(=O)[OH1])]):n:*(*[$Hal,$(S(=O)(=O)C)]):*(*[!N]):*(*[!N]):*1	6activated 2halo and 2sulfonyl pyridines and not too deactivated
#
c:1:n:*:*(*[$Hal,$(C(C)=NO),$(C#N),$(C=O),$(C(F)(F)F),$(S(=O)=O);!$(C(=O)[OH1])]):*(*[$Hal,$(S(=O)(=O)C)]):*1	3activated 4halo and 4sulfonyl pyridines
#
c:1(*[!N&!O]):c(*[$Hal,$(C#N),$(C(=O)),$(C(F)(F)F),$(S(=O)=O);!$(C(=O)[OH1])]):c(*[$Hal,$(S(=O)(=O)C),$(C[R](=O)NC)]):c(*[!N&!O]):c(*[!N&!O]):c(*[$Hal,$(C#N),$(C(=O)),$(C(F)(F)F),$(S(=O)=O);!$(C(=O)[OH1])])1	1,2,4activated aryls (including activated phthalimides
#
c:1(*[$Hal,$(C#N),$(C=O),$(C(F)(F)F),$(S(=O)=O);!$(C(=O)[OH1])]):c(*[$Hal,$(S(=O)(=O)C),$(C[R](=O)NC)]):c(*[$Hal,$(C#N),$(C=O),$(C(F)(F)F),$(S(=O)=O);!$(C(=O)[OH1])]):c(*[!N&!O]):c(*[!N&!O]):c(*[!N&!O])1	1,2,3activated aryls
#
c:1(*[!OX2;!NX3]):n:c(:[$Hev]:[$Hev]:n1)[!R]*[O,S]c:8:[$Hev](*[!O;!SX2;!CX4;!NX3]):[$Hev]:[$Hev](*[$Hal,$(C#N),$(C(F)(F)F),$(S(=O)=O)]):[$Hev]:[$Hev](*[!O;!SX2;!CX4;!NX3])8	activated pyrimidine phenol ethers (1)
#
c:1(*[!OX2;!NX3]):n:c(:[$Hev]:[$Hev]:n1)[!R]*[O,S]c:8:[$Hev](*[$Hal,$(C#N),$(C(F)(F)F),$(S(=O)=O)]):[$Hev]:[$Hev](*[!O;!SX2;!CX4;!NX3]):[$Hev]:[$Hev](*[!O;!SX2;!CX4;!NX3])8	activated pyrimidine phenol ethers (2)
#
c:1(*[!OX2;!NX3]):n:c(:n:[$Hev]:[$Hev][!n]1)[!R]*[O,S]c:8:[$Hev](*[!O;!SX2;!CX4;!NX3]):[$Hev]:[$Hev](*[$Hal,$(C#N),$(C(F)(F)F),$(S(=O)=O)]):[$Hev]:[$Hev](*[!O;!SX2;!CX4;!NX3])8	activated pyrimidine phenol ethers (3)
#
c:1(*[!OX2;!NX3]):n:c(:n:[$Hev]:[$Hev][!n]1)[!R]*[O,S]c:8:[$Hev](*[$Hal,$(C#N),$(C(F)(F)F),$(S(=O)=O)]):[$Hev]:[$Hev](*[!O;!SX2;!CX4;!NX3]):[$Hev]:[$Hev](*[!O;!SX2;!CX4;!NX3])8	activated pyrimidine phenol ethers (4)
c:1(S(=O)(=O)Nc:n):[cH1]:[cH1]:c(F):[cH1]:[cH1]1	activated 4fluorosulphonamides
[$Het]:1:c([$Hal]):n:c:c1	thiazolinium halides and related
CC(C)=C([$Hal])N([$(*H),$(*[CX4])])[$(*H),$(*[CX4])]	reactive haloenamines
C#[CH1]	terminal acetylenes
##################################
##Limiting number of certain slns#
##################################
s:1:c:[cH1]:[cH1]:[cH1]1.s:1:c:[cH1]:[cH1]:[cH1]1	thiophenes <max=1>
o:1:c:[cH1]:[cH1]:[cH1]1.o:1:c:[cH1]:[cH1]:[cH1]1	furans <max=1>
#
C#N.C#N.C#N	nitriles <max=2>
C(=O)[OH1].C(=O)[OH1].C(=O)[OH1]	acids <max=2>
#
C(=O)OC.C(=O)OC.C(=O)OC	esters <max=2>
#
[CX4][SX2][CX4].[CX4][SX2][CX4].[CX4][SX2][CX4]	thioethers <max=2>
#
Br.Br bromine <max=1>
Cl.Cl.Cl.Cl chlorine <max=3>
F.F.F.F.F.F fluorine <max=5>
#############################################################################
