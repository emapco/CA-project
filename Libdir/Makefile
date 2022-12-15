# Chem 274B - Cellular Automata Final Project
# Emmanuel Cortes
# Contributor(s): Trevor Oldham
# Created: 12/10/2022
# This makefile creates library object files that containing all
# cellular automata functionality. 


# cellular automata object files (sequential and parallelized)
CA_OBJS = cellularautomata.o CA_utils.o
CA_OMP_OBJS = cellularautomata_omp.o CA_utils_omp.o
# shared library files
CA_LIB = cellularautomata.a
CA_OMP_LIB = cellularautomata_omp.a

cellularautomata.a: cleanall
	ar rU $(CA_LIB) $(CA_OBJS)
	ranlib $(CA_LIB) 
	rm -f $(CA_OBJS)

cellularautomata_omp.a: cleanall
	ar rU $(CA_OMP_LIB) $(CA_OMP_OBJS)
	ranlib $(CA_OMP_LIB) 
	rm -f $(CA_OMP_OBJS)

sequential: $(CA_LIB)

parallel: $(CA_OMP_LIB)

all: $(CA_LIB) $(CA_OMP_LIB)

cleanall:
	rm -f $(CA_LIB)
	rm -f $(CA_OMP_LIB)
