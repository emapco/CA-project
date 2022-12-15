# Chem 274B - Cellular Automata Final Project
# Creator: Emmanuel Cortes
# Contributor(s): 
# Created: 12/04/2022
# This makes all targets in the project directory

SOURCE_DIR = Source/Datatypes
TEST_DIR = Tests
UTILS_DIR = Utils
LIB_DIR = Libdir
APP_DIR = Applications

sequential:
	cd $(SOURCE_DIR); make sequential
	cd $(UTILS_DIR); make sequential
	cd $(LIB_DIR); make sequential
	cd $(TEST_DIR); make sequential
	cd $(APP_DIR); make sequential

parallel:
	cd $(SOURCE_DIR); make parallel
	cd $(UTILS_DIR); make parallel
	cd $(LIB_DIR); make parallel
	cd $(TEST_DIR); make parallel
	cd $(APP_DIR); make parallel

all:                       
	cd $(SOURCE_DIR); make all
	cd $(UTILS_DIR); make all
	cd $(LIB_DIR); make all
	cd $(TEST_DIR); make all
	cd $(APP_DIR); make all

cleanall:
	cd $(SOURCE_DIR); make cleanall
	cd $(UTILS_DIR); make cleanall
	cd $(LIB_DIR); make cleanall
	cd $(TEST_DIR); make cleanall
	cd $(APP_DIR); make cleanall

