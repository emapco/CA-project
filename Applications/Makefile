# Chem 274B - Cellular Automata Final Project
# Emmanuel Cortes
# Contributor(s): Trevor Oldham
# Created: 12/06/2022
# This creates an executable for our galaxy CA model


# GNU C++ Compiler
CPP         = g++      

# compiler flags -g debug, -O3 optimized version
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
SRC_DIR     = ../Source/Datatypes
UTIL_DIR    = ../Utils

# The next line contains the list of object files created by this Makefile.
EXECS = galaxy_model galaxy_model_omp

galaxy_model:
	$(CPP) $(CPPFLAGS) -I$(INC_DIR) \
	$(LIB_DIR)/galaxy.o $(LIB_DIR)/cellularautomata.a \
	galaxy_app.cpp -o galaxy_model
	mv galaxy_model $(BIN_DIR)

galaxy_model_omp:
	$(CPP) $(CPPFLAGS) $(OMPFLAGS) $(LDFLAGS) -I$(INC_DIR) \
	$(LIB_DIR)/galaxy.o $(LIB_DIR)/cellularautomata_omp.a \
	galaxy_app.cpp -o galaxy_model_omp
	mv galaxy_model_omp $(BIN_DIR)

sequential: galaxy_model

parallel: galaxy_model_omp

all: $(EXECS)

cleanall:
	cd $(BIN_DIR); rm -f $(EXECS)
