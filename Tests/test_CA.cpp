/**
 * @file test_CA.cpp
 * @author Trevor Oldham (trevoldham@berkeley.edu)
 *
 * <b>Contributor(s)</b> <br> &emsp;&emsp; Emmanuel Cortes
 * @brief this file creates a new CA object and
 * tests the basic functionality to prove correctness
 * @date 2022-12-03
 */
#include <iostream>
#include <CAdatatypes.h>
#include <array>

/**
 * @brief Tests the initialization of CellularAutomata instances.
 *
 * @param CA the CellularAutomata instance
 */
void test_CA_init(CellularAutomata &CA)
{
    // initialize the first grid state using a value for x_state and a probability
    std::cout << "Initializing " << CA.axis1_dim
              << " x " << CA.axis2_dim 
              << " x " << CA.axis3_dim << " grid of cells" << std::endl;
    CA.init_condition(1, 0.4);

    // print the current grid of cell states
    std::cout << "Printing the current grid of cell states" << std::endl;
    CA.print_grid();

    // print other relevant class attributes of the CA object
    std::cout << "Printing other relevant info about CellularAutomata object" << std::endl;
    std::cout << "Number of States: " << CA.num_states << std::endl;
    std::cout << "Boundary Radius: " << CA.boundary_radius << std::endl;
//     std::cout << "Short Radius: " << CA.shortr_weight << std::endl;
//     std::cout << "Long Radius: " << CA.longr_weight << std::endl;
}

/**
 * @brief Tests the CellularAutomata instances's step function
 *
 * @param CA the CellularAutomata instance
 */
void test_CA_step(CellularAutomata &CA)
{
    for (int i = 0; i < 2; i++)
    {
        CA.step();
        CA.print_grid();
    }
}

/**
 * @brief Main testing function that calls all other types of tests
 *
 * @param CA the CellularAutomata instance
 */
void test_CA(CellularAutomata &CA)
{
    test_CA_init(CA);
    test_CA_step(CA);
}

/**
 * @brief Tests the implementation of the vector CA
 *
 */
void test_CA_vector()
{
    {
        std::cout << "**** Testing Vector CellularAutomata: PERIODIC PARITY ****\n";
        // initialize a CA object with default constructor
        CellularAutomata CA = CellularAutomata();
        CA.setup_dimensions(10);
        CA.setup_rule(Parity);
        test_CA(CA);
    }
    {
        std::cout << "**** Testing Vector CellularAutomata: PERIODIC MAJORITY ****\n";
        // initialize a CA object with default constructor
        CellularAutomata CA = CellularAutomata();
        CA.setup_dimensions(10);
        test_CA(CA);
    }
}

/**
 * @brief Tests the implementation of the matrix CA
 *
 */
void test_CA_matrix()
{
    {
        std::cout << "\n**** Testing Matrix CellularAutomata: PERIODIC PARITY ****\n";
        // initialize a CA object with default constructor
        CellularAutomata CA = CellularAutomata();
        CA.setup_dimensions(5, 10);
        CA.setup_rule(Parity);
        test_CA(CA);
    }
    {
        std::cout << "\n**** Testing Matrix CellularAutomata: PERIODIC MAJORITY ****\n";
        // initialize a CA object with default constructor
        CellularAutomata CA = CellularAutomata();
        CA.setup_dimensions(5, 10);
        test_CA(CA);
    }
}

/**
 * @brief Tests the implementation of the tensor CA
 *
 */
void test_CA_tensor()
{
    {
        std::cout << "\n**** Testing Tensor CellularAutomata: PERIODIC PARITY ****\n";
        // initialize a CA object with default constructor
        CellularAutomata CA = CellularAutomata();
        CA.setup_dimensions(5, 5, 10);
        test_CA(CA);
    }
    {
        std::cout << "\n**** Testing Tensor CellularAutomata: PERIODIC MAJORITY ****\n";
        // initialize a CA object with default constructor
        CellularAutomata CA = CellularAutomata();
        CA.setup_dimensions(5, 5, 10);
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