# Chem 274B - Cellular Automata Final Project
# Emmanuel Cortes
# Contributor(s): Trevor Oldham
# Created: 12/06/2022
# This creates object files for utility functions in CAutils.h


# GNU C++ Compiler
CPP         = g++

# compiler flags -g debug, -O3 optimized version -c create a library object
CPPFLAGS    =-Wall -Wextra -O3 -c
OMPFLAGS    =-fopenmp -DENABLE_OMP

# Add additional flags for Mac OS X
ifeq ($(detected_OS),Darwin)
OMPFLAGS    =-Xclang -DENABLE_OMP
endif

# The directory where the include files needed to create the library objects are
INC_DIR = ../Include

# The directory where all object files are going to be moved after creation
LIB_DIR     = ../Libdir

# The next line contains the list of object files created by this Makefile.
OBJS = CA_utils.o CA_utils_omp.o

CA_utils.o:
	$(CPP) $(CPPFLAGS) -I$(INC_DIR) CA_utils.cpp 
	mv CA_utils.o $(LIB_DIR)

CA_utils_omp.o:
	$(CPP) $(CPPFLAGS) $(OMPFLAGS) $(LDFLAGS) -I$(INC_DIR) CA_utils.cpp \
	-o CA_utils_omp.o
	mv CA_utils_omp.o $(LIB_DIR)

sequential: CA_utils.o

parallel: CA_utils_omp.o

all: $(OBJS)

cleanall:
	cd $(LIB_DIR); rm -f $(OBJS)
