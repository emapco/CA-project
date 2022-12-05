// Chem 274B - Cellular Automata Final Project
// Trevor Oldham
// Created: 12/03/2022
    
//this file creates a new CA object and tests the basic functionality to prove correctness

#include<iostream>
#include<mydatatypes.h>
#include<array>

int main()
{
    //initialize a CA object with default constructor, pass dimensions as arguments
    CellularAutomata CA = CellularAutomata(20, 20);
    
    //initialize the first grid state using a value for x_state and a probability
    std::cout << "Initializing " << CA.n << " x " << CA.m << " grid of cells" << std::endl;
    CA.init_condition(1, 0.1);

    //print the current grid of cell states
    std::cout << "Printing the current grid of cell states" << std::endl;
    CA.print_grid();

    //print other relevant class attributes of the CA object
    std::cout << "Printing other relevant info about CellularAutomata object" << std::endl;
    std::cout << "Number of States: " << CA.num_states << std::endl;
    std::cout << "Boundary Radius: " << CA.boundary_radius << std::endl;
    std::cout << "Short Radius: " << CA.shortr_weight << std::endl;
    std::cout << "Long Radius: " << CA.longr_weight << std::endl;

    return 0;

}