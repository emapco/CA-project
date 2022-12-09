# Chem 274B - Cellular Automata Final Project
# Emmanuel Cortes
# Created: 12/06/2022
# This creates object files for utility functions in CAutils.h


# GNU C++ Compiler
CPP         = g++      

# compiler flags -g debug, -O3 optimized version -c create a library object
CPPFLAGS    =-Wall -Wextra -O0 -c -g

# The directory where the include files needed to create the library objects are
INC_DIR = ../Include

# The directory where all object files are going to be moved after creation
LIB_DIR     = ../Libdir

# The next line contains the list of object files created by this Makefile.
OBJS = CAutils.o

CAutils.o:  $(INC_DIR)/CAutils.h
	$(CPP) $(CPPFLAGS) CAutils.cpp -I$(INC_DIR)
	mv CAutils.o $(LIB_DIR)

all: $(OBJS)

cleanall:
	cd $(LIB_DIR); rm $(OBJS)
