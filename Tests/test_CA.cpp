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
 * @param CA 
 */
void test_CA_init(CellularAutomata &CA)
{
    // initialize the first grid state using a value for x_state and a probability
    std::cout << "Initializing " << CA.axis1_dim << " x " << CA.axis2_dim << " grid of cells" << std::endl;
    CA.init_condition(1, 0.1);

    // print the current grid of cell states
    std::cout << "Printing the current grid of cell states" << std::endl;
    CA.print_grid();

    // print other relevant class attributes of the CA object
    std::cout << "Printing other relevant info about CellularAutomata object" << std::endl;
    std::cout << "Number of States: " << CA.num_states << std::endl;
    std::cout << "Boundary Radius: " << CA.boundary_radius << std::endl;
    std::cout << "Short Radius: " << CA.shortr_weight << std::endl;
    std::cout << "Long Radius: " << CA.longr_weight << std::endl;
}

int main()
{
    int error;
    {
        std::cout << "**** Testing Vector CellularAutomata ****\axis1_dim";
        // initialize a CA object with default constructor, pass dimensions as arguments
        CellularAutomata CA = CellularAutomata();
        error = CA.setup_dimensions(10);
        test_CA_init(CA);
    }
    {
        std::cout << "\axis1_dim\axis1_dim**** Testing Matrix CellularAutomata ****\axis1_dim";
        // initialize a CA object with default constructor, pass dimensions as arguments
        CellularAutomata CA = CellularAutomata();
        error = CA.setup_dimensions(10, 20);
        test_CA_init(CA);
    }
    {
        std::cout << "\axis1_dim\axis1_dim**** Testing Tensor CellularAutomata ****\axis1_dim";
        // initialize a CA object with default constructor, pass dimensions as arguments
        CellularAutomata CA = CellularAutomata();
        error = CA.setup_dimensions(5, 10, 20);
        test_CA_init(CA);
    }

    return 0;
}