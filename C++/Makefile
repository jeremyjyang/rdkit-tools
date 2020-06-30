#RDCXXFLAGS = -I$(RDBASE)/Code
#RDLFLAGS = -L$(RDBASE)/lib

###
# Ubuntu 20.04 (apt install rdkit python3-rdkit rdkit-data)
RDCXXFLAGS = -I/usr/include/rdkit
RDLFLAGS = -L/usr/lib

CC       = /usr/bin/gcc
CXX      = /usr/bin/g++
#CFLAGS = -m64 -W -Wall -pthread -O3 -fomit-frame-pointer -ffast-math
CFLAGS = -m64 -pthread -O3 -fomit-frame-pointer -ffast-math
#CXXFLAGS = -W -Wall -O3 -fomit-frame-pointer -ffast-math
CXXFLAGS = -O3 -fomit-frame-pointer -ffast-math
LINKER   = $(CXX)
LFLAGS   = -m64 -s
LSOFLAGS = -shared
#LIBS = -lboost_python -lpython2.6 -lboost_thread -lz -lm -lpthread
#LIBS = -lboost_python -lpython2.7 -lboost_thread -lz -lm -lpthread

# Ubuntu 20.04 (apt install rdkit python3-rdkit rdkit-data)
LIBS = -lboost_python38 -lpython3.8 -lboost_thread -lz -lm -lpthread

#LIBS = -lboost_python -lpython -lboost_thread -lz -lm -lpthread
#LIBS = /usr/lib64/libboost_python.so /usr/lib64/libboost_thread.so /usr/lib64/libpython2.7.so -lz -lm -lpthread
RM = /bin/rm -f


RDLIBS = \
	-lRDKitAlignment \
	-lRDKitCatalogs \
	-lRDKitChemicalFeatures \
	-lRDKitChemReactions \
	-lRDKitChemTransforms \
	-lRDKitDataStructs \
	-lRDKitDepictor \
	-lRDKitDescriptors \
	-lRDKitDistGeometry \
	-lRDKitDistGeomHelpers \
	-lRDKitEigenSolvers \
	-lRDKitFileParsers \
	-lRDKitFingerprints \
	-lRDKitForceFieldHelpers \
	-lRDKitForceField \
	-lRDKitFragCatalog \
	-lRDKitGraphMol \
	-lRDKitMolAlign \
	-lRDKitMolCatalog \
	-lRDKitMolChemicalFeatures \
	-lRDKitMolTransforms \
	-lRDKitOptimizer \
	-lRDKitPartialCharges \
	-lRDKitRDBoost \
	-lRDKitRDGeneral \
	-lRDKitRDGeometryLib \
	-lRDKitShapeHelpers \
	-lRDKitSimDivPickers \
	-lRDKitSLNParse \
	-lRDKitSmilesParse \
	-lRDKitSubgraphs \
	-lRDKitSubstructMatch

PROGS = substruct_test1 sample smartstest_usa

default: $(PROGS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(RDCXXFLAGS) -g -c $<
.o :
	$(LINKER) $(LFLAGS) $(RDLFLAGS) -g -o $@ $@.o \
	$(RDLIBS) $(LIBS)
.cpp:
	@make $@.o
	@make $@

clean:
	$(RM) *.o *.a $(PROGS) core
