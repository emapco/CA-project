/**
 * @file test_CA.cpp
 * @author Trevor Oldham (trevoldham@berkeley.edu)
 *
 * <b>Contributor(s)</b> <br> &emsp;&emsp; Emmanuel Cortes
 * @brief this file creates a new CA object and
 * tests the basic functionality to prove correctness
 * @date 2022-12-03
 */
#include "CAdatatypes.h"
#include <iostream>
#include <array>

/**
 * @brief Test custom rule function for setting a new cell state for the cell at cell_index
 *
 * @param cell_index array of cell indices that we are going to update it state for
 * @param index_size number of indices need to address the cell
 * @param neighborhood_cells array of neighboring cells
 * @param neighborhood_size size of neighborhood_cells array
 * @param new_cell_state reference to the new cell state
 */
void custom_rule(int *cell_index, const int index_size,
                 int *neighborhood_cells, const int neighborhood_size,
                 int &new_cell_state, CellularAutomata<int>&)
{
    new_cell_state = 1;
}

/**
 * @brief Tests the initialization of CellularAutomata instances.
 *
 * @param CA the CellularAutomata instance
 */
void test_CA_init(CellularAutomata<int> &CA)
{
    // initialize the first grid state using a value for x_state and a probability
    std::cout << "Initializing " << CA.axis1_dim
              << " x " << CA.axis2_dim
              << " x " << CA.axis3_dim << " grid of cells" << std::endl;
    CA.init_condition(1, 0.5);

    // print the current grid of cell states
    std::cout << "Printing the current grid of cell states" << std::endl;
    CA.print_grid();

    // print other relevant class attributes of the CA object
    std::cout << "Printing other relevant info about CellularAutomata object" << std::endl;
    std::cout << "Number of States: " << CA.num_states << std::endl;
    std::cout << "Boundary Radius: " << CA.boundary_radius << std::endl;
}

/**
 * @brief Tests the CellularAutomata instances's step function
 *
 * @param CA the CellularAutomata instance
 */
void test_CA_step(CellularAutomata<int> &CA)
{
    for (int i = 0; i < 1; i++)
    {
        CA.step(custom_rule);
        CA.print_grid();
    }
}

/**
 * @brief Main testing function that calls all other types of tests
 *
 * @param CA the CellularAutomata instance
 */
void test_CA(CellularAutomata<int> &CA)
{
    test_CA_init(CA);
    test_CA_step(CA);
}

/**
 * @brief Tests the implementation of the CA with vector of cells
 *
 */
void test_CA_vector()
{
    {
        std::cout << "**** Testing Vector CellularAutomata: PERIODIC PARITY ****\n";
        // initialize a CA object with default constructor
        CellularAutomata<int> CA = CellularAutomata<int>();
        CA.setup_dimensions_1d(30);
        CA.setup_rule(CAEnums::Rule::Parity);
        test_CA(CA);
    }
    {
        std::cout << "**** Testing Vector CellularAutomata: PERIODIC MAJORITY ****\n";
        CellularAutomata<int> CA = CellularAutomata<int>();
        CA.setup_dimensions_1d(30);
        test_CA(CA);
    }
    {
        std::cout << "**** Testing Vector CellularAutomata: PERIODIC CUSTOM ****\n";
        CellularAutomata<int> CA = CellularAutomata<int>();
        CA.setup_dimensions_1d(30);
        CA.setup_rule(CAEnums::Rule::Custom);
        test_CA(CA);
    }
}

/**
 * @brief Tests the implementation of the CA with a matrix of cells
 *
 */
void test_CA_matrix()
{
    {
        std::cout << "\n**** Testing Matrix CellularAutomata: PERIODIC PARITY ****\n";
        // initialize a CA object with default constructor
        CellularAutomata<int> CA = CellularAutomata<int>();
        CA.setup_dimensions_2d(5, 10);
        CA.setup_rule(CAEnums::Rule::Parity);
        test_CA(CA);
    }
    {
        std::cout << "\n**** Testing Matrix CellularAutomata: PERIODIC MAJORITY ****\n";
        CellularAutomata<int> CA = CellularAutomata<int>();
        CA.setup_dimensions_2d(5, 10);
        test_CA(CA);
    }
    {
        std::cout << "**** Testing Vector CellularAutomata: PERIODIC CUSTOM ****\n";
        CellularAutomata<int> CA = CellularAutomata<int>();
        CA.setup_dimensions_2d(5, 10);
        CA.setup_rule(CAEnums::Rule::Custom);
        test_CA(CA);
    }
}

/**
 * @brief Tests the implementation of the CA with a tensor of cells
 *
 */
void test_CA_tensor()
{
    {
        std::cout << "\n**** Testing Tensor CellularAutomata: PERIODIC PARITY ****\n";
        // initialize a CA object with default constructor
        CellularAutomata<int> CA = CellularAutomata<int>();
        CA.setup_dimensions_3d(5, 5, 10);
        CA.setup_rule(CAEnums::Rule::Parity);
        test_CA(CA);
    }
    {
        std::cout << "\n**** Testing Tensor CellularAutomata: PERIODIC MAJORITY ****\n";
        CellularAutomata<int> CA = CellularAutomata<int>();
        // CA.setup_dimensions(50, 50, 50);
        CA.setup_dimensions_3d(5, 5, 10);
        test_CA(CA);
    }
    {
        std::cout << "**** Testing Vector CellularAutomata: PERIODIC CUSTOM ****\n";
        CellularAutomata<int> CA = CellularAutomata<int>();
        CA.setup_dimensions_3d(5, 5, 10);
        CA.setup_rule(CAEnums::Rule::Custom);
        test_CA(CA);
    }
}

int main()
{
    test_CA_vector();
    test_CA_matrix();
    test_CA_tensor();
    return 0;
}