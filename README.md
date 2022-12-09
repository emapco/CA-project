# Chem 274B - Cellular Automata Final Project
// Created: 12/03/2022

## updates:
    12/3/2022:  Trevor: added basic directory structure, added CAdatatypes.h, added cellularautomata.cpp class functionality, added test_CA.cpp in Tests, added makefiles in /Tests and /Source/Datatypes

    12/4/2022:  Emmanuel: Added root directory Makefile to call all makefiles in the project. Added README.md files to /Bindir, /Include, and /Libdir subdirectories.

    12/5/2022:  Emmanuel: Added Doxyfile for generating documentation.

    12/5/2022:  Group Meeting/Trevor: Updated mydatatypes.h to CAdatatypes.h. Added doxygen comment briefs to cellularautomata.cpp methods. Updated test_CA.cpp to test 3d case. And Updated CellularAutomata API.

    12/5/2022:  Emmanuel: Updated Doxyfile and renamed Docs/ to docs/ for GitHub Pages support. 
    Resolved bug encountered in test_CA.cpp when initializing a vector or tensor. 
    Miscellaneous changes to existing docstrings.

    12/8/2022:  Emmanuel: Added step function and related functions for computing the next cellular automata state. The step function supports the ability to pass a custom rule function to allow the user to implement their own rules. Added guard clauses to various methods so that the class properly handles failure modes. Also defined a method for printing an error message depending on the error code given.
