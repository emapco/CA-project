# Chem 274B - Cellular Automata Final Project
# Trevor Oldham
# Contributor: Emmanuel Cortes
# Created: 12/03/2022
# This creates object files for datatypes in CAdatatypes.h


# GNU C++ Compiler
CPP         = g++      

# compiler flags -g debug, -O3 optimized version -c create a library object
CPPFLAGS    =-O0 -g 
OMPFLAGS    =-fopenmp -DENABLE_OMP
LDFLAGS     =

# Add additional flags for Mac OS X
ifeq ($(detected_OS),Darwin)
OMPFLAGS    +=-Xclang 
LDFLAGS     +=-L/usr/local/opt/llvm/lib -lomp
endif

# The directory where the include files needed to create the library objects are
INC_DIR = ../Include

# The directory where all object files are going to be moved after creation
LIB_DIR     = ../Libdir
BIN_DIR     = ../Bindir

# The next line contains the list of object files created by this Makefile.
EXECS = test_CA test_CA_omp

test_CA:
	$(CPP) $(CPPFLAGS) -I$(INC_DIR) test_CA.cpp -o test_CA \
	 $(LIB_DIR)/cellularautomata.a
	mv test_CA $(BIN_DIR)

test_CA_omp:
	$(CPP) $(CPPFLAGS) $(OMPFLAGS) $(LDFLAGS) -I$(INC_DIR) test_CA.cpp \
	-o test_CA_omp $(LIB_DIR)/cellularautomata_omp.a
	mv test_CA_omp $(BIN_DIR)

all: $(EXECS)

cleanall:
	cd $(BIN_DIR); rm $(EXECS)
