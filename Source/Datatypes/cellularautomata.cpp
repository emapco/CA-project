// Chem 274B - Cellular Automata Final Project
// Trevor Oldham
// Created: 12/03/2022
    
// This file contains the class CellularAutomata functions declared in mydatatypes.h

#include "mydatatypes.h"
#include <iostream>
#include<array>
#include<random>
#include<ctime>

//default constructor for object - sets the default value to all class attributes
//inputs:
//      n - number of rows in first dimension
//      m - number of columns in second dimension
CellularAutomata::CellularAutomata(int n, int m)
{
    this->n = n;
    this->m = m;
    num_states = 2;
    boundary b_type = None;
    b_radius = 1;
    neighborhood n_type = VonNeumann;
    rule r_type = Majority;
    this->shortr_weight = 1;
    this->longr_weight = 2;

    //matrix = new (nothrow) int[n][m]
    matrix = new (nothrow) int*[n];
    for (int i = 0; i < n; i++)
    {
        matrix[i] = new (nothrow) int[m];
    }

    //initialize matrix filled with zeros
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < m; k++)
        {
            matrix[j][k] = 0;
        }
    }
}

//default destructor: deallocates memory reserved for matrix
CellularAutomata::~CellularAutomata()
{
    //Free each sub-array
    for(int i = 0; i < n; ++i) {
        delete[] matrix[i];   
    }
    //Free the array of pointers
    delete[] matrix;
    
}

//function: setup neighborhood
//inputs:
//      neighboorhood_type:     type of neighborhood, described as enum with choices (VonNeumann or Moore)
int CellularAutomata::setup_neighborhood(neighborhood neighborhood_type)
{
    this->n_type = neighborhood_type;
}

//function: setup boundary
//inputs:
//      bound_type:     type of boundary described as enum with choices (None, Periodic, Walled, CutOff)
int CellularAutomata::setup_boundary(boundary bound_type, int radius)
{
    this->b_type = bound_type;
    this->b_radius = radius;
}

//function:  setup_cell_states
//describes the range of cell states to be used in the CA object
//inputs:
//      n_states:   int describing the range of numbers to use for cell states
int CellularAutomata::setup_cell_states(int n_states)
{
    this->num_states = n_states;
}

//function: init_condition
//initializes the first state of the grid using random numbers
//inputs:
//      x_state:    int: choose the cell state to initialize the grid with. 
//                      run this function multiple times with different numbers for x_state to create a grid with more than two states
//      prob:       double: the probability of a cell to turn to state given from x_state
int CellularAutomata::init_condition(int x_state, double prob)
{
    //initialize values in matrix
    //initialize matrix filled with random cell values
    srand(time(NULL));
    double r;
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < m; k++)
        {
            r = (double) rand()/RAND_MAX;
            if (r < prob)
            {
                matrix[j][k] = x_state;
            }
        }
    }
}

int CellularAutomata::init_cond_rewrite(int x_state, double prob)
{


}

//function:     setup_rule_short_long
//function to set the short and long weights used in the majority rule
//inputs:
//      shortr_weight:  double: the short weight
//      longr_weight:   double: the long weight
int CellularAutomata::setup_rule_short_long(double shortr_weight, double longr_weight)
{
    this->shortr_weight = shortr_weight;
    this->longr_weight = longr_weight;
}

//function:     setup_rule
//used to specifiy the rule type to be used in CA object
//inputs:
//      rule_type:      enum rule representing the rule type. choose from (Majority, Parity)
int CellularAutomata::setup_rule(rule rule_type)
{
    this->r_type = rule_type;
}

//function:     print_grid
//prints the current values in each cell in a formatted grid
int CellularAutomata::print_grid()
{
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < m; k++)
        {
            std::cout << matrix[j][k] << " ";
        }
        std::cout << std::endl;
    }
}

