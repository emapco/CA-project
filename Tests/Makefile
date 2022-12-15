# Chem 274B - Cellular Automata Final Project
# Trevor Oldham
# Contributor: Emmanuel Cortes
# Created: 12/03/2022
# This creates object files for datatypes in CAdatatypes.h


# GNU C++ Compiler
CPP         = g++      

# compiler flags -g debug, -O3 optimized version -c create a library object
CPPFLAGS    =-O3
OMPFLAGS    =-fopenmp -DENABLE_OMP
LDFLAGS     =

# Add additional flags for Mac OS X
ifeq ($(detected_OS),Darwin)
OMPFLAGS    =-Xclang -DENABLE_OMP
LDFLAGS     +=-L/usr/local/opt/llvm/lib -lomp
endif

# The directory where the include files needed to create the library objects are
INC_DIR = ../Include

# The directory where all object files are going to be moved after creation
LIB_DIR     = ../Libdir
BIN_DIR     = ../Bindir

# The next line contains the list of object files created by this Makefile.
EXECS = test_CA test_CA_omp unit_test_CA_utils

test_CA:
	$(CPP) $(CPPFLAGS) -I$(INC_DIR) test_CA.cpp -o test_CA \
	 $(LIB_DIR)/cellularautomata.a
	mv test_CA $(BIN_DIR)

test_CA_omp:
	$(CPP) $(CPPFLAGS) $(OMPFLAGS) $(LDFLAGS) -I$(INC_DIR) test_CA.cpp \
	-o test_CA_omp $(LIB_DIR)/cellularautomata_omp.a
	mv test_CA_omp $(BIN_DIR)

unit_test_CA_utils:
	$(CPP) $(CPPFLAGS) -I$(INC_DIR) unit_test_CA_utils.cpp $(LIB_DIR)/cellularautomata.a -o unit_test_CA_utils
	mv unit_test_CA_utils $(BIN_DIR)

sequential: test_CA unit_test_CA_utils

parallel: test_CA_omp

all: $(EXECS)

cleanall:
	cd $(BIN_DIR); rm -f $(EXECS)
