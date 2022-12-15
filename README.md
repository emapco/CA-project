# Chem 274B - Cellular Automata Final Project

The repository contains a general-purpose library for create cellular automata models. The cellular automata cell states can be `ints` or `class`/`struct` types. See CAdatatypes.h for information on the `class`/`struct` requirements. The general-purpose library also supports parallelization with the aid of OpenMP.

Provided in this repository is an example model that simulates the formation of galaxies using GalaxyCell `class` instances as the cell type. Once the model is compiled (see below for information), you will find an executable called `galaxy_app` (sequential) or optionally `galaxy_app_omp` (parallel). To run the model from root project directory run: `Bindir/galaxy_model`.


Requirements for compilation:
- `g++`
- `make`
- `OpenMP` (parallelized code)

To compile the sequential implementation, run:
- `make sequential`.

To compile the parallel implementation, run:
- `make parallel`.

To compile both, run:
- `make all`.

To clean all object files and executables, run:
- `make cleanall`.


## Directories:
- Applications: contains galaxy model program source code.
- Bindir: contains executable binary files
- docs: contains doxygen generated html site
- Include: contains Header files
- Libdir: contains compiled object code
- Source: contains implementation source code
- Tests: contains unit tests and test programs
- Utils: contains utility source code for `*util.h` files in the `Include` directory

## Files:
- Makefile: contains targets for compile the project source code
- Doxyfile: doxygen documentation generation settings 


## Updates:
    12/3/2022:  Trevor: added basic directory structure, added CAdatatypes.h, added cellularautomata.cpp class functionality, added test_CA.cpp in Tests, added makefiles in /Tests and /Source/Datatypes

    12/4/2022:  Emmanuel: Added root directory Makefile to call all makefiles in the project. Added README.md files to /Bindir, /Include, and /Libdir subdirectories.

    12/5/2022:  Emmanuel: Added Doxyfile for generating documentation.

    12/5/2022:  Trevor: Updated mydatatypes.h to CAdatatypes.h. Added doxygen comment briefs to cellularautomata.cpp methods. 
    Updated test_CA.cpp to test 3d case. And Updated CellularAutomata API.

    12/5/2022:  Emmanuel: Updated Doxyfile and renamed Docs/ to docs/ for GitHub Pages support. 
    Resolved bug encountered in test_CA.cpp when initializing a vector or tensor. 
    Miscellaneous changes to existing docstrings.

    12/8/2022:  Emmanuel: Added step function and related functions for computing the next cellular automata state. 
    The step function supports the ability to pass a custom rule function to allow the user to implement their own rules.
    Added guard clauses to various methods so that the class properly handles failure modes. 
    Also defined a method for printing an error message depending on the error code given.

    12/10/2022: Emmanuel: Added parallelization with the aid of OpenMP directives. 
    Updated makefiles to create new parallelized targets. The makefiles support Linux and Mac OSX OpenMP systems.
    Class implementation and utility object files are now combined into a single library object file.

    12/11/2022: Emmanuel: Added templates to CellularAutomata class. Also created a specialized template class for `int` data type.
    Added utility functions to CAutils.h that are general and will useful for our galaxy model. Add unit tests for CA_utils.cpp.

    12/14/2022: Emmanuel: Added `galaxydatatypes.h`, `galaxy.cpp`, and an application that runs our model `galaxy_app.cpp` using user input. 
    Modified makefiles based on Trevor's input so targets compile for MasOS systems and can handle other edge cases. Updated doc.
