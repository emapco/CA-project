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

/**
 * @brief counter type for determining which cell state is the majority
 *
 */
using MajorityCounter = std::unordered_map<int, int>;

/*! A CellularAutomata class for simulating cellular automata models */
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

    // function for setting the new cell state based on the specified rule
    int set_new_cell_state(int *cell_index, int index_size,
                           int *neighborhood_cells, int neighborhood_size,
                           int &new_cell_state, void(custom_rule)(int *, int, int *, int, int &));

    // function for setting the cell state with generated neighborhood array; vector case
    int get_state_from_neighborhood_1d(int *cell_index, int index_size, int &new_cell_state,
                                       void(custom_rule)(int *, int, int *, int, int &));

    // function for setting the cell state with generated neighborhood array; matrix case
    int get_state_from_neighborhood_2d(int *cell_index, int index_size, int &new_cell_state,
                                       void(custom_rule)(int *, int, int *, int, int &));

    // function for setting the cell state with generated neighborhood array; tensor case
    int get_state_from_neighborhood_3d(int *cell_index, int index_size, int &new_cell_state,
                                       void(custom_rule)(int *, int, int *, int, int &));

    // method for generating periodic neighborhood around cell_index (matrix case)
    void generate_periodic_neighborhood_1d(int *cell_index, int *neighborhood_cells, int &neighborhood_index);

    // method for generating cutoff neighborhood around cell_index (matrix case)
    void generate_cutoff_neighborhood_1d(int *cell_index, int *neighborhood_cells, int &neighborhood_index);

    // method for generating periodic neighborhood around cell_index (matrix case)
    void generate_periodic_neighborhood_2d(int *cell_index, int *neighborhood_cells, int &neighborhood_index);

    // method for generating cutoff neighborhood around cell_index (matrix case)
    void generate_cutoff_neighborhood_2d(int *cell_index, int *neighborhood_cells, int &neighborhood_index);

    // method for generating periodic neighborhood around cell_index (tensor case)
    void generate_periodic_neighborhood_3d(int *cell_index, int *neighborhood_cells, int &neighborhood_index);

    // method for generating cutoff neighborhood around cell_index (tensor case)
    void generate_cutoff_neighborhood_3d(int *cell_index, int *neighborhood_cells, int &neighborhood_index);

    // returns a pointer to dynamic neighborhood array
    int *malloc_neighborhood_array(int rank);

public:
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
        Parity,
        Custom
    };

    /**
     * @brief enum containing the various error codes
     * the CellularAutomata class can return
     *
     */
    enum error_code
    {
        CellsAlreadyInitialized = -1,
        CellsAreNull = -2,
        CellsMalloc = -3,
        InvalidCellState = -4,
        InvalidCellStateCondition = -5,
        InvalidRadius = -6,
        InvalidNumStates = -7,
        NeighborhoodCellsMalloc = -8,
        CustomRuleIsNull = -9
    };

    boundary boundary_type;         //!< enum code to hold boundary
    int boundary_radius;            //!< declare a radius for the boundary
    neighborhood neighborhood_type; //!< enum code to hold neighborhood type
    int num_states;                 //!< integer code for number of states
    rule rule_type;                 //!< enum code for rule type
    int axis1_dim;                  //!< count of cells in first dimension
    int axis2_dim;                  //!< count of cells in second dimension
    int axis3_dim;                  //!< count of cells in third dimension

    CellularAutomata();  // default constructor
    ~CellularAutomata(); // destructor

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

    // provide a setup function to choose the number of states
    int setup_cell_states(int n_states);

    // provide a function to initialize the beginning state of the grid
    int init_condition(int x_state, double prob);

    // provide a function for setting the rule type
    int setup_rule(rule rule_type);

    // simulates a CA step
    int step();

    // simulates a CA step with a custom_rule
    int step(void(custom_rule)(int *, int, int *, int, int &));

    // provide a function to print the current grid
    int print_grid();

    // provide a printing function that displays a message for the given error code
    void print_error_status(error_code error);
};
