# Chem 274B - Cellular Automata Final Project
# Trevor Oldham
# Created: 12/03/2022
# This creates object files for datatypes in CAdatatypes.h


# GNU C++ Compiler
CPP         = g++      

# compiler flags -g debug, -O3 optimized version -c create a library object
CPPFLAGS    =-O3 -o

# The directory where the include files needed to create the library objects are
INC_DIR = ../Include

# The directory where all object files are going to be moved after creation
LIB_DIR     = ../Libdir
BIN_DIR     = ../Bindir

# The next line contains the list of object files created by this Makefile.
EXECS = test_CA

test_CA: $(INC_DIR)/CAdatatypes.h
	$(CPP) $(CPPFLAGS) test_CA test_CA.cpp \
	-I$(INC_DIR) $(LIB_DIR)/cellularautomata.o
	mv test_CA $(BIN_DIR)

all: $(EXECS)

cleanall:
	cd $(BIN_DIR); rm $(EXECS)
