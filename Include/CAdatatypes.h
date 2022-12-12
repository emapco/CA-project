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

/**
 * @brief Enum namespace that contains various enum definitions
 * used by CellularAutomata.
 */
namespace CAEnums
{
    /**
     * @brief enum containing the different type of neighborhoods
     * the CellularAutomata class supports
     *
     */
    enum Neighborhood
    {
        VonNeumann,
        Moore
    };

    /**
     * @brief enum containing the different type of boundaries
     * the CellularAutomata class supports
     *
     */
    enum Boundary
    {
        Periodic,
        Walled,
        CutOff
    };

    /**
     * @brief enum containing the different type of rules
     * the CellularAutomata class supports for it's transitions
     */
    enum Rule
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
    enum ErrorCode
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
}

/**
 * @brief A base CellularAutomata class that contains untemplated member variables and method definitions
 * from which templated classes can inherit.
 */
class BaseCellularAutomata
{
public:
    CAEnums::Boundary boundary_type;         //!< enum code to hold boundary
    int boundary_radius;                     //!< declare a radius for the boundary
    CAEnums::Neighborhood neighborhood_type; //!< enum code to hold neighborhood type
    int num_states;                          //!< integer code for number of states
    CAEnums::Rule rule_type;                 //!< enum code for rule type
    int axis1_dim;                           //!< count of cells in first dimension
    int axis2_dim;                           //!< count of cells in second dimension
    int axis3_dim;                           //!< count of cells in third dimension

    /**
     * @brief Construct a new Cellular Automata:: Cellular Automata object.
     * Sets the default value to all class attributes.
     *
     */
    BaseCellularAutomata();

    /**
     * @brief Destroy the Cellular Automata:: Cellular Automata object.
     * Deallocates memory reserved for vector/matrix/tensor.
     *
     */
    virtual ~BaseCellularAutomata(){};

    /**
     * @brief Setup neighborhood with values from enum neighborhood.
     *
     * @param neighborhood_type enum value for neighborhood type (VonNeumann or Moore)
     * @return int error code
     */
    int setup_neighborhood(CAEnums::Neighborhood neighborhood_type);

    /**
     * @brief Setup boundary with enum values from boundary and set the boundary radius.
     *
     * @param bound_type enum value for boundary (None, Periodic, Walled, CutOff)
     * @param radius radius for the boundary
     * @return int - error code\n
     * InvalidRadius: radius can't be less than equal to 0\n
     * 0: no error
     */
    int setup_boundary(CAEnums::Boundary bound_type, int radius);

    /**
     * @brief Defines the range of cell states to be used in the CA object.
     *
     * @param num_states describes the range of numbers to use for cell states
     * @return int - error code\n
     * InvalidNumStates: num_states can't be less than to 2\n
     * 0: no error
     */
    int setup_cell_states(int n_states);

    /**
     * @brief Setup the rule choice to specify the rule type to be used in CA object.
     *
     * @param rule_type enum rule representing the rule type
     * @return int error code
     */
    int setup_rule(CAEnums::Rule rule_type);

    /**
     * @brief Prints an error message for the given error code
     *
     * @param error ErrorCode enum type
     */
    void print_error_status(CAEnums::ErrorCode error);
};

/**
 * @brief A CellularAutomata class for simulating cellular automata models.
 * This templated class supports structs/class objects.
 * The object must contain a .state property that the CellularAutomata class updates.
 * The object must also define a move operator to ensure well defined swapping of states (current state and next state).
 * The object must also define an assignment operator to ensure proper copying happen when assigning one object to another.
 * If these requirements are not satisfied then the class will produce undefined behavior.
 *
 * @tparam T : struct/class with a .state property, move operator and assignment operator.
 */
template <typename T>
class CellularAutomata : public BaseCellularAutomata
{
private:
    T *vector;        //!< vector for one dimensional grid of cells
    T *next_vector;   //!< vector for one dimensional grid of cells holding the next state
    T **matrix;       //!< pointer to 2d array for grid cells holding a state
    T **next_matrix;  //!< pointer to 2d array for grid cells holding the next state
    T ***tensor;      //!< tensor for three dimensional grid of cells
    T ***next_tensor; //!< tensor for three dimensional grid of cells holding the next state
    int steps_taken;  //!< the number of steps the CA has taken

    /**
     * @brief Sets the new_cell_state variable based on the specified rule.
     *
     * @param cell_index cell of interest's index; can be modified for dynamic models
     * @param index_size size of cell_index array
     * @param neighborhood_cells flatten array of all neighboring cells
     * @param neighborhood_size size of neighborhood_cells array
     * @param new_cell_state reference variable for setting the new state
     * @param custom_rule function that is called when a Custom rule type is specified
     *@return int - error code\n
     * CellsAreNull: neither vector, matrix, nor tensor are initialized\n
     * CustomRuleIsNull: given custom rule function is null\n
     * 0: no error
     */
    int set_new_cell_state(int *cell_index, int index_size,
                           T *neighborhood_cells, int neighborhood_size,
                           T &new_cell_state, void(custom_rule)(int *, int, T *, int, T &));

    /**
     * @brief Generates an array of neighboring cells and then calls set_new_cell_state to set the state.
     * This method supports a vector of cell states.
     *
     * @param cell_index cell of interest's index; can be modified for dynamic models
     * @param index_size size of cell_index array
     * @param new_cell_state reference variable for setting the new state
     * @param custom_rule function that is called when a Custom rule type is specified
     *@return int - error code\n
     * Error codes returned by set_new_cell_state\n
     * 0: no error
     */
    int get_state_from_neighborhood_1d(int *cell_index, int index_size, T &new_cell_state,
                                       void(custom_rule)(int *, int, T *, int, T &));

    /**
     * @brief Generates an array of neighboring cells and then calls set_new_cell_state to set the state.
     * This method supports a matrix of cell states.
     *
     * @param cell_index cell of interest's index; can be modified for dynamic models
     * @param index_size size of cell_index array
     * @param new_cell_state reference variable for setting the new state
     * @param custom_rule function that is called when a Custom rule type is specified
     *@return int - error code\n
     * Error codes returned by set_new_cell_state\n
     * 0: no error
     */
    int get_state_from_neighborhood_2d(int *cell_index, int index_size, T &new_cell_state,
                                       void(custom_rule)(int *, int, T *, int, T &));

    /**
     * @brief Generates an array of neighboring cells and then calls set_new_cell_state to set the state.
     * This method supports a tensor of cell states.
     *
     * @param cell_index cell of interest's index; can be modified for dynamic models
     * @param index_size size of cell_index array
     * @param new_cell_state reference variable for setting the new state
     * @param custom_rule function that is called when a Custom rule type is specified
     *@return int - error code\n
     * Error codes returned by set_new_cell_state\n
     * 0: no error
     */
    int get_state_from_neighborhood_3d(int *cell_index, int index_size, T &new_cell_state,
                                       void(custom_rule)(int *, int, T *, int, T &));

    /**
     * @brief Adds periodic neighboring cell state to neighborhood_cells array.
     * Handles vector (1d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_periodic_neighborhood_1d(int *cell_index, T *neighborhood_cells, int &neighborhood_index);

    /**
     * @brief Adds cutoff neighboring cell state to neighborhood_cells array.
     * Handles vector (1d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_cutoff_neighborhood_1d(int *cell_index, T *neighborhood_cells, int &neighborhood_index);

    /**
     * @brief Adds periodic neighboring cell state to neighborhood_cells array.
     * Handles matrix (2d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_periodic_neighborhood_2d(int *cell_index, T *neighborhood_cells, int &neighborhood_index);

    /**
     * @brief Adds cutoff neighboring cell state to neighborhood_cells array.
     * Handles matrix (2d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_cutoff_neighborhood_2d(int *cell_index, T *neighborhood_cells, int &neighborhood_index);

    /**
     * @brief Adds periodic neighboring cell state to neighborhood_cells array.
     * Handles tensor (3d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_periodic_neighborhood_3d(int *cell_index, T *neighborhood_cells, int &neighborhood_index);

    /**
     * @brief Adds cutoff neighboring cell state to neighborhood_cells array.
     * Handles tensor (3d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_cutoff_neighborhood_3d(int *cell_index, T *neighborhood_cells, int &neighborhood_index);

    /**
     * @brief Dynamically allocates an array that will contain the neighborhood cell states.
     *
     * @param rank rank of the cells data
     * @return T* pointer to the dynamic array
     */
    T *malloc_neighborhood_array(int rank);

public:
    /**
     * @brief Construct a new Cellular Automata object.
     * Sets the default value to all class attributes.
     *
     */
    CellularAutomata();

    /**
     * @brief Destroy the Cellular Automata object.
     * Deallocates memory reserved for vector/matrix/tensor.
     *
     */
    ~CellularAutomata();

    /**
     * @brief Set up 1d array of cell states.
     *
     * @param axis1_dim the size of a one dimensional vector
     * @param fill_value the value to set every cell state to
     * @return int - error code\n
     * CellsAlreadyInitialized: vector was already allocated\n
     * CellsMalloc: couldn't allocate memory for the specified vector size\n
     * 0: no error
     */
    int setup_dimensions_1d(int axis1_dim, int fill_value = 0);

    /**
     * @brief Set up 2d matrix of cell states.
     *
     * @param axis1_dim  The size of the first dimension
     * @param axis2_dim  The size of the second dimension
     * @param fill_value the value to set every cell state to
     * @return int - error code\n
     * CellsAlreadyInitialized: matrix was already allocated\n
     * CellsMalloc: couldn't allocate memory for the specified matrix size\n
     * 0: no error
     */
    int setup_dimensions_2d(int axis1_dim, int axis2_dim, int fill_value = 0);

    /**
     * @brief Set up 3d tensor of cell states.
     *
     * @param axis1_dim The size of the first dimension
     * @param axis2_dim The size of the second dimension
     * @param axis3_dim The size of the third dimension
     * @param fill_value the value to set every cell state to
     * @return int - error code\n
     * CellsAlreadyInitialized: tensor was already allocated\n
     * CellsMalloc: couldn't allocate memory for the specified tensor size\n
     * 0: no error
     */
    int setup_dimensions_3d(int axis1_dim, int axis2_dim, int axis3_dim, int fill_value = 0);

    /**
     * @brief Initializes the first state of the grid using random numbers.
     *
     * @param x_state choose the cell state to initialize the grid with.
     * @param prob the probability of a cell to turn to state given from x_state
     *@return int - error code\n
     * CellsAreNull: tensor not initialized\n
     * InvalidCellStateCondition: x_state must be less than num_states\n
     * 0: no error
     */
    int init_condition(int x_state, double prob);

    /**
     * @brief Simulates a cellular automata step.
     * A new state is generated and stored stored as the new state for subsequent calls to step method.
     *
     * This method supports the use of a custom rule type.
     *
     * @param custom_rule function that is called when a Custom rule type is specified
     *@return int - error code\n
     * Error codes returned by get_state_from_neighborhood_Xd methods\n
     * 0: no error
     */
    int step(void(custom_rule)(int *, int, T *, int, T &));

    /**
     * @brief Simulates a cellular automata step.
     * A new state is generated and stored as the new state for subsequent calls to step method.
     *
     *@return int - error code\n
     * Error codes returned by get_state_from_neighborhood_Xd methods\n
     * 0: no error
     */
    int step();

    /**
     * @brief Print the current state of the grid.
     *
     * @return int - error code\n
     * CellsAreNull: neither vector, matrix, nor tensor are initialized\n
     * 0: no error
     */
    int print_grid();
};

/**
 * @brief A CellularAutomata class for simulating cellular automata models.
 * This specialized template handles integer cell states.
 *
 * @tparam int
 */
template <>
class CellularAutomata<int> : public BaseCellularAutomata
{
private:
    int *vector;        //!< vector for one dimensional grid of cells
    int *next_vector;   //!< vector for one dimensional grid of cells holding the next state
    int **matrix;       //!< pointer to 2d array for grid cells holding a state
    int **next_matrix;  //!< pointer to 2d array for grid cells holding the next state
    int ***tensor;      //!< tensor for three dimensional grid of cells
    int ***next_tensor; //!< tensor for three dimensional grid of cells holding the next state
    int steps_taken;    //!< the number of steps the CA has taken

    /**
     * @brief Sets the new_cell_state variable based on the specified rule.
     *
     * @param cell_index cell of interest's index; can be modified for dynamic models
     * @param index_size size of cell_index array
     * @param neighborhood_cells flatten array of all neighboring cells
     * @param neighborhood_size size of neighborhood_cells array
     * @param new_cell_state reference variable for setting the new state
     * @param custom_rule function that is called when a Custom rule type is specified
     *@return int - error code\n
     * CellsAreNull: neither vector, matrix, nor tensor are initialized\n
     * CustomRuleIsNull: given custom rule function is null\n
     * 0: no error
     */
    int set_new_cell_state(int *cell_index, int index_size,
                           int *neighborhood_cells, int neighborhood_size,
                           int &new_cell_state, void(custom_rule)(int *, int, int *, int, int &));

    /**
     * @brief Generates an array of neighboring cells and then calls set_new_cell_state to set the state.
     * This method supports a vector of cell states.
     *
     * @param cell_index cell of interest's index; can be modified for dynamic models
     * @param index_size size of cell_index array
     * @param new_cell_state reference variable for setting the new state
     * @param custom_rule function that is called when a Custom rule type is specified
     *@return int - error code\n
     * Error codes returned by set_new_cell_state\n
     * 0: no error
     */
    int get_state_from_neighborhood_1d(int *cell_index, int index_size, int &new_cell_state,
                                       void(custom_rule)(int *, int, int *, int, int &));

    /**
     * @brief Generates an array of neighboring cells and then calls set_new_cell_state to set the state.
     * This method supports a matrix of cell states.
     *
     * @param cell_index cell of interest's index; can be modified for dynamic models
     * @param index_size size of cell_index array
     * @param new_cell_state reference variable for setting the new state
     * @param custom_rule function that is called when a Custom rule type is specified
     *@return int - error code\n
     * Error codes returned by set_new_cell_state\n
     * 0: no error
     */
    int get_state_from_neighborhood_2d(int *cell_index, int index_size, int &new_cell_state,
                                       void(custom_rule)(int *, int, int *, int, int &));

    /**
     * @brief Generates an array of neighboring cells and then calls set_new_cell_state to set the state.
     * This method supports a tensor of cell states.
     *
     * @param cell_index cell of interest's index; can be modified for dynamic models
     * @param index_size size of cell_index array
     * @param new_cell_state reference variable for setting the new state
     * @param custom_rule function that is called when a Custom rule type is specified
     *@return int - error code\n
     * Error codes returned by set_new_cell_state\n
     * 0: no error
     */
    int get_state_from_neighborhood_3d(int *cell_index, int index_size, int &new_cell_state,
                                       void(custom_rule)(int *, int, int *, int, int &));

    /**
     * @brief Adds periodic neighboring cell state to neighborhood_cells array.
     * Handles vector (1d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_periodic_neighborhood_1d(int *cell_index, int *neighborhood_cells, int &neighborhood_index);

    /**
     * @brief Adds cutoff neighboring cell state to neighborhood_cells array.
     * Handles vector (1d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_cutoff_neighborhood_1d(int *cell_index, int *neighborhood_cells, int &neighborhood_index);

    /**
     * @brief Adds periodic neighboring cell state to neighborhood_cells array.
     * Handles matrix (2d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_periodic_neighborhood_2d(int *cell_index, int *neighborhood_cells, int &neighborhood_index);

    /**
     * @brief Adds cutoff neighboring cell state to neighborhood_cells array.
     * Handles matrix (2d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_cutoff_neighborhood_2d(int *cell_index, int *neighborhood_cells, int &neighborhood_index);

    /**
     * @brief Adds periodic neighboring cell state to neighborhood_cells array.
     * Handles tensor (3d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_periodic_neighborhood_3d(int *cell_index, int *neighborhood_cells, int &neighborhood_index);

    /**
     * @brief Adds cutoff neighboring cell state to neighborhood_cells array.
     * Handles tensor (3d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_cutoff_neighborhood_3d(int *cell_index, int *neighborhood_cells, int &neighborhood_index);

    /**
     * @brief Dynamically allocates an array that will contain the neighborhood cell states.
     *
     * @param rank rank of the cells data
     * @return T* pointer to the dynamic array
     */
    int *malloc_neighborhood_array(int rank);

public:
    /**
     * @brief Construct a new Cellular Automata object.
     * Sets the default value to all class attributes.
     *
     */
    CellularAutomata();

    /**
     * @brief Destroy the Cellular Automata object.
     * Deallocates memory reserved for vector/matrix/tensor.
     *
     */
    ~CellularAutomata();

    /**
     * @brief Set up 1d array of cell states.
     *
     * @param axis1_dim the size of a one dimensional vector
     * @param fill_value the value to set every cell state to
     * @return int - error code\n
     * CellsAlreadyInitialized: vector was already allocated\n
     * CellsMalloc: couldn't allocate memory for the specified vector size\n
     * 0: no error
     */
    int setup_dimensions_1d(int axis1_dim, int fill_value = 0);

    /**
     * @brief Set up 2d matrix of cell states.
     *
     * @param axis1_dim  The size of the first dimension
     * @param axis2_dim  The size of the second dimension
     * @param fill_value the value to set every cell state to
     * @return int - error code\n
     * CellsAlreadyInitialized: matrix was already allocated\n
     * CellsMalloc: couldn't allocate memory for the specified matrix size\n
     * 0: no error
     */
    int setup_dimensions_2d(int axis1_dim, int axis2_dim, int fill_value = 0);

    /**
     * @brief Set up 3d tensor of cell states.
     *
     * @param axis1_dim The size of the first dimension
     * @param axis2_dim The size of the second dimension
     * @param axis3_dim The size of the third dimension
     * @param fill_value the value to set every cell state to
     * @return int - error code\n
     * CellsAlreadyInitialized: tensor was already allocated\n
     * CellsMalloc: couldn't allocate memory for the specified tensor size\n
     * 0: no error
     */
    int setup_dimensions_3d(int axis1_dim, int axis2_dim, int axis3_dim, int fill_value = 0);

    /**
     * @brief Initializes the first state of the grid using random numbers.
     *
     * @param x_state choose the cell state to initialize the grid with.
     * @param prob the probability of a cell to turn to state given from x_state
     *@return int - error code\n
     * CellsAreNull: tensor not initialized\n
     * InvalidCellStateCondition: x_state must be less than num_states\n
     * 0: no error
     */
    int init_condition(int x_state, double prob);

    /**
     * @brief Simulates a cellular automata step.
     * A new state is generated and stored stored as the new state for subsequent calls to step method.
     *
     * This method supports the use of a custom rule type.
     *
     * @param custom_rule function that is called when a Custom rule type is specified
     *@return int - error code\n
     * Error codes returned by get_state_from_neighborhood_Xd methods\n
     * 0: no error
     */
    int step(void(custom_rule)(int *, int, int *, int, int &));

    /**
     * @brief Simulates a cellular automata step.
     * A new state is generated and stored as the new state for subsequent calls to step method.
     *
     *@return int - error code\n
     * Error codes returned by get_state_from_neighborhood_Xd methods\n
     * 0: no error
     */
    int step();

    /**
     * @brief Print the current state of the grid.
     *
     * @return int - error code\n
     * CellsAreNull: neither vector, matrix, nor tensor are initialized\n
     * 0: no error
     */
    int print_grid();
};
