# Chem 274B - Cellular Automata Final Project
# Trevor Oldham
# Contributor: Emmanuel Cortes
# Created: 12/03/2022
# This creates object files for datatypes in CAdatatypes.h


# GNU C++ Compiler
CPP         = g++      

# compiler flags -g debug, -O3 optimized version -c create a library object
CPPFLAGS    =-O0 -g 

# The directory where the include files needed to create the library objects are
INC_DIR = ../Include

# The directory where all object files are going to be moved after creation
LIB_DIR     = ../Libdir
BIN_DIR     = ../Bindir

# The next line contains the list of object files created by this Makefile.
EXECS = test_CA

test_CA:
	$(CPP) $(CPPFLAGS) -I$(INC_DIR) test_CA.cpp -o test_CA \
	 $(LIB_DIR)/cellularautomata.o $(LIB_DIR)/CAutils.o
	mv test_CA $(BIN_DIR)

all: $(EXECS)

cleanall:
	cd $(BIN_DIR); rm $(EXECS)
