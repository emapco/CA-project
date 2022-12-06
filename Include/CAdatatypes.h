/**
 * @file CAdatatypes.h
 * @author Trevor Oldham (trevoldham@berkeley.edu)
 *
 * <b>Contributor(s)</b> <br> &emsp;&emsp; Emmanuel Cortes
 * @brief  This file contains the custom datatypes used
 * in the Cellular Automata project
 * @date 2022-12-03
 */

#pragma once // Ensures that this file is only included once
             // during compilation

#include <array>
using namespace std;

enum neighborhood
{
    VonNeumann,
    Moore
};

enum boundary
{
    None,
    Periodic,
    Walled,
    CutOff
};

enum rule
{
    Majority,
    Parity
};

enum error_codes
{
    CellsAlreadyInitialized = -1,
    CellsAreNull = -2,
};

class CellularAutomata // The CA data structure
{
public:
    boundary boundary_type;         //!< enum code to hold boundary
    int boundary_radius;            //!< declare a radius for the boundary
    neighborhood neighborhood_type; //!< enum code to hold neighborhood type
    int num_states;                 //!< integer code for number of states

    rule rule_type;       //!< enum code for rule type
    double shortr_weight; //!< double for short radius
    double longr_weight;  //!< double for long radius
    int ndims;            //!< int for number of dimensions
    int n;                //!< count of cells in first dimension
    int m;                //!< count of cells in second dimension
    int p;                //!< count of cells in third dimension
    int *vector;          //!< vector for one dimensional grid of cells
    int **matrix;         //!< pointer to 2d array for grid cells holding a state
    int ***tensor;        //!< tensor for three dimensional grid of cells

    ///////////////////////////////////////////
    CellularAutomata();
    ~CellularAutomata();

    // provide a setup function to choose from two neighborhood types
    int setup_neighborhood(neighborhood neighborhood_type);

    // provide a setup function for setting the vector dimension
    int setup_dimensions(int n);

    // provide a setup function for setting the matrix's dimensions
    int setup_dimensions(int n, int m);

    // provide a setup function for setting the tensor's dimensions
    int setup_dimensions(int n, int m, int p);

    // provide a setup function to choose from four boundary types
    int setup_boundary(boundary bound_type, int radius);

    // provide a setup function to choose the number of states
    int setup_cell_states(int n_states);

    // provide a function to initialize the beginning state of the grid
    int init_condition(int x_state, double prob);

    // provide a function for setting the Majority rule radii
    int setup_rule_short_long(double shortr_weight, double longr_weight);

    // provide a function for setting the rule type
    int setup_rule(rule rule_type);

    // provide a function to print the current grid
    int print_grid();
};
