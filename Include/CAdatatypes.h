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
#include <unordered_map>
using namespace std;

/**
 * @brief enum containing the different type of neighborhoods
 * the CellularAutomata class supports
 *
 */
enum neighborhood
{
    VonNeumann,
    Moore
};

/**
 * @brief enum containing the different type of boundaries
 * the CellularAutomata class supports
 *
 */
enum boundary
{
    Periodic,
    Walled,
    CutOff
};

/**
 * @brief enum containing the different type of rules
 * the CellularAutomata class supports for it's transitions
 */
enum rule
{
    Majority,
    Parity
};

/**
 * @brief enum containing the various error codes
 * the CellularAutomata class can return
 *
 */
enum error_codes
{
    CellsAlreadyInitialized = -1,
    CellsAreNull = -2,
    InvalidNeighborhood = -3,
};

class CellularAutomata // The CA data structure
{
private:
    int *vector;        //!< vector for one dimensional grid of cells
    int *next_vector;   //!< vector for one dimensional grid of cells holding the next state
    int **matrix;       //!< pointer to 2d array for grid cells holding a state
    int **next_matrix;  //!< pointer to 2d array for grid cells holding the next state
    int ***tensor;      //!< tensor for three dimensional grid of cells
    int ***next_tensor; //!< tensor for three dimensional grid of cells holding the next state
    int steps_taken;    //!< the number of steps the CA has taken
public:
    boundary boundary_type;         //!< enum code to hold boundary
    int boundary_radius;      //!< declare a radius for the boundary
    int long_boundary_radius;       //!< declare a radius for a long range boundary
    neighborhood neighborhood_type; //!< enum code to hold neighborhood type
    int num_states;                 //!< integer code for number of states

    rule rule_type;       //!< enum code for rule type
    double shortr_weight; //!< double for short radius
    double longr_weight;  //!< double for long radius

    int axis1_dim; //!< count of cells in first dimension
    int axis2_dim; //!< count of cells in second dimension
    int axis3_dim; //!< count of cells in third dimension

    ///////////////////////////////////////////
    CellularAutomata();
    ~CellularAutomata();

    // provide a setup function to choose from two neighborhood types
    int setup_neighborhood(neighborhood neighborhood_type);

    // provide a setup function for setting the vector dimension
    int setup_dimensions(int axis1_dim);

    // provide a setup function for setting the matrix's dimensions
    int setup_dimensions(int axis1_dim, int axis2_dim);

    // provide a setup function for setting the tensor's dimensions
    int setup_dimensions(int axis1_dim, int axis2_dim, int axis3_dim);

    // provide a setup function to choose from three boundary types
    int setup_boundary(boundary bound_type, int radius);

    // provide a setup function to choose from three boundary types
    // with a long and short radius
    int setup_boundary(boundary bound_type, int long_radius, int short_radius);

    // provide a setup function to choose the number of states
    int setup_cell_states(int n_states);

    // provide a function to initialize the beginning state of the grid
    int init_condition(int x_state, double prob);

    // provide a function for setting the Majority rule radii
    int setup_rule_short_long(double shortr_weight, double longr_weight);

    // provide a function for setting the rule type
    int setup_rule(rule rule_type);

    int get_cell_state(int sum, const unordered_map<int, int> &votes_counter);

    int get_state_from_neighborhood(int i, int &new_cell_state);

    int get_state_from_neighborhood(int i, int j, int &new_cell_state);

    int get_state_from_neighborhood(int i, int j, int k, int &new_cell_state);

    int step();

    int step(int(func)(int, int));

    // provide a function to print the current grid
    int print_grid();
};
