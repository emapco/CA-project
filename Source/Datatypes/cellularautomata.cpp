/**
 * @file cellularautomata.cpp
 * @author Trevor Oldham (trevoldham@berkeley.edu)
 *
 * <b>Contributor(s)</b> <br> &emsp;&emsp; Emmanuel Cortes
 * @brief This file contains the class CellularAutomata functions declared in CAdatatypes.h
 * @date 2022-12-03
 */
#include "CAdatatypes.h"
#include "CAutils.h"
#include <iostream>
#include <array>
#include <random>        // srand, rand
#include <ctime>         // time
#include <unordered_map> // unordered_map
#include <utility>       // make_pair
#include <algorithm>     // max_element

/**
 * @brief Construct a new Cellular Automata:: Cellular Automata object.
 * Sets the default value to all class attributes.
 *
 */
CellularAutomata::CellularAutomata()
{
    axis1_dim = 0;
    axis2_dim = 0;
    axis3_dim = 0;
    num_states = 2;
    boundary_type = Periodic;
    boundary_radius = 1;
    neighborhood_type = Moore;
    rule_type = Majority;
    steps_taken = 0;
    vector = nullptr;
    next_vector = nullptr;
    matrix = nullptr;
    next_matrix = nullptr;
    tensor = nullptr;
    next_tensor = nullptr;
}

/**
 * @brief Destroy the Cellular Automata:: Cellular Automata object.
 * Deallocates memory reserved for vector/matrix/tensor.
 *
 */
CellularAutomata::~CellularAutomata()
{
    // Free each sub-array
    if (vector != nullptr)
    {
        // Free the array of pointers
        delete[] vector;
        delete[] next_vector;
    }

    if (matrix != nullptr)
    {
        // Free each sub-array
        for (int i = 0; i < axis1_dim; ++i)
        {
            delete[] matrix[i];
            delete[] next_matrix[i];
        }
        // Free the array of pointers
        delete[] matrix;
        delete[] next_matrix;
    }

    if (tensor != nullptr)
    {
        // Free each sub-array
        for (int i = 0; i < axis1_dim; ++i)
        {
            for (int j = 0; j < axis2_dim; j++)
            {
                delete[] tensor[i][j];
                delete[] next_tensor[i][j];
            }

            delete[] tensor[i];
            delete[] next_tensor[i];
        }

        // Free the array of pointers
        delete[] tensor;
        delete[] next_tensor;
    }
}

/**
 * @brief Construct a new CellularAutomata::setup_dimensions object.
 *
 * @param axis1_dim the size of a one dimensional vector
 * @return int - error code\n
 * CellsAlreadyInitialized: vector was already allocated\n
 * CellsMalloc: couldn't allocate memory for the specified vector size\n
 * 0: no error
 */
int CellularAutomata::setup_dimensions(int axis1_dim)
{
    if (vector != nullptr)
    {
        return CellsAlreadyInitialized;
    }

    this->axis1_dim = axis1_dim;
    vector = new (std::nothrow) int[axis1_dim];
    next_vector = new (std::nothrow) int[axis1_dim];

    if (vector == nullptr || next_vector == nullptr)
    {
        return CellsMalloc;
    }

    // initialize vector filled with zeros
    for (int j = 0; j < axis1_dim; j++)
    {
        vector[j] = 0;
        next_vector[j] = 0;
    }

    return 0;
}

/**
 * @brief Construct a new CellularAutomata::setup_dimensions object.
 *
 * @param axis1_dim  The size of the first dimension
 * @param axis2_dim  The size of the second dimension
 * @return int - error code\n
 * CellsAlreadyInitialized: matrix was already allocated\n
 * CellsMalloc: couldn't allocate memory for the specified matrix size\n
 * 0: no error
 */
int CellularAutomata::setup_dimensions(int axis1_dim, int axis2_dim)
{
    if (matrix != nullptr)
    {
        return CellsAlreadyInitialized;
    }

    this->axis1_dim = axis1_dim;
    this->axis2_dim = axis2_dim;
    matrix = new (std::nothrow) int *[axis1_dim];
    next_matrix = new (std::nothrow) int *[axis1_dim];

    if (matrix == nullptr || next_matrix == nullptr)
    {
        return CellsMalloc;
    }

    for (int i = 0; i < axis1_dim; i++)
    {
        matrix[i] = new (std::nothrow) int[axis2_dim];
        next_matrix[i] = new (std::nothrow) int[axis2_dim];
    }

    // initialize matrix filled with zeros
    for (int j = 0; j < axis1_dim; j++)
    {
        for (int k = 0; k < axis2_dim; k++)
        {
            matrix[j][k] = 0;
            next_matrix[j][k] = 0;
        }
    }

    return 0;
}

/**
 * @brief Construct a new CellularAutomata::setup_dimensions object.
 *
 * @param axis1_dim The size of the first dimension
 * @param axis2_dim The size of the second dimension
 * @param axis3_dim The size of the third dimension
 * @return int - error code\n
 * CellsAlreadyInitialized: tensor was already allocated\n
 * CellsMalloc: couldn't allocate memory for the specified tensor size\n
 * 0: no error
 */
int CellularAutomata::setup_dimensions(int axis1_dim, int axis2_dim, int axis3_dim)
{
    if (tensor != nullptr)
    {
        return CellsAlreadyInitialized;
    }

    this->axis1_dim = axis1_dim;
    this->axis2_dim = axis2_dim;
    this->axis3_dim = axis3_dim;
    tensor = new (std::nothrow) int **[axis1_dim];
    next_tensor = new (std::nothrow) int **[axis1_dim];

    if (tensor == nullptr || next_tensor == nullptr)
    {
        return CellsMalloc;
    }

    for (int i = 0; i < axis1_dim; i++)
    {
        tensor[i] = new (std::nothrow) int *[axis2_dim];
        next_tensor[i] = new (std::nothrow) int *[axis2_dim];

        for (int j = 0; j < axis2_dim; j++)
        {
            tensor[i][j] = new (std::nothrow) int[axis3_dim];
            next_tensor[i][j] = new (std::nothrow) int[axis3_dim];
        }
    }

    // initialize matrix filled with zeros
    for (int i = 0; i < axis1_dim; i++)
    {
        for (int j = 0; j < axis2_dim; j++)
        {
            for (int k = 0; k < axis3_dim; k++)
            {
                tensor[i][j][k] = 0;
                next_tensor[i][j][k] = 0;
            }
        }
    }

    return 0;
}

/**
 * @brief Setup neighborhood with values from enum neighborhood.
 *
 * @param neighborhood_type enum value for neighborhood type (VonNeumann or Moore)
 * @return int error code
 */
int CellularAutomata::setup_neighborhood(neighborhood neighborhood_type)
{
    this->neighborhood_type = neighborhood_type;
    return 0;
}

/**
 * @brief Setup boundary with enum values from boundary and set the boundary radius.
 *
 * @param bound_type enum value for boundary (None, Periodic, Walled, CutOff)
 * @param radius radius for the boundary
 * @return int - error code\n
 * InvalidRadius: radius can't be less than equal to 0\n
 * 0: no error
 */
int CellularAutomata::setup_boundary(boundary bound_type, int radius)
{
    if (radius <= 0)
    {
        return InvalidRadius;
    }
    this->boundary_type = bound_type;
    this->boundary_radius = radius;
    return 0;
}

/**
 * @brief Defines the range of cell states to be used in the CA object.
 *
 * @param num_states describes the range of numbers to use for cell states
 * @return int - error code\n
 * InvalidNumStates: num_states can't be less than to 2\n
 * 0: no error
 */
int CellularAutomata::setup_cell_states(int num_states)
{
    if (num_states < 2)
    {
        return InvalidNumStates;
    }
    this->num_states = num_states;
    return 0;
}

/**
 * @brief Setup the rule choice sed to specify the rule type to be used in CA object.
 *
 * @param rule_type enum rule representing the rule type
 * @return int error code
 */
int CellularAutomata::setup_rule(rule rule_type)
{
    this->rule_type = rule_type;
    return 0;
}

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
int CellularAutomata::init_condition(int x_state, double prob)
{
    if (!(x_state < num_states))
    {
        return InvalidCellStateCondition;
    }

    srand(time(NULL));
    double random_cell_state;

    if (vector != nullptr)
    {
        for (int i = 0; i < axis1_dim; i++)
        {
            random_cell_state = (double)rand() / RAND_MAX;
            if (random_cell_state < prob)
            {
                vector[i] = x_state;
            }
        }
    }
    else if (matrix != nullptr)
    {
        for (int j = 0; j < axis1_dim; j++)
        {
            for (int k = 0; k < axis2_dim; k++)
            {
                random_cell_state = (double)rand() / RAND_MAX;
                if (random_cell_state < prob)
                {
                    matrix[j][k] = x_state;
                }
            }
        }
    }
    else if (tensor != nullptr)
    {
        for (int i = 0; i < axis1_dim; i++)
        {
            for (int j = 0; j < axis2_dim; j++)
            {
                for (int k = 0; k < axis3_dim; k++)
                {
                    random_cell_state = (double)rand() / RAND_MAX;
                    if (random_cell_state < prob)
                    {
                        tensor[i][j][k] = x_state;
                    }
                }
            }
        }
    }
    else
    {
        return CellsAreNull;
    }
    return 0;
}

/**
 * @brief Sets the new_cell_state variable based on the specified rule.
 *
 * @param cell_index cell of interest's index
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
int CellularAutomata::set_new_cell_state(int *cell_index, int index_size,
                                         int *neighborhood_cells, int neighborhood_size,
                                         int &new_cell_state, void(custom_rule)(int *, int, int *, int, int &))
{
    int sum = 0;                         // sum of cells within boundary_radius for Parity rule
    MajorityCounter state_votes_counter; // counter to keep track of votes for Majority rule
    initialize_majority_rule_counter(state_votes_counter, num_states);

    switch (rule_type)
    {
    case Custom:
        if (vector == nullptr && matrix == nullptr && tensor == nullptr)
        {
            new_cell_state = -1;
            return CellsAreNull;
        }
        if (custom_rule == nullptr)
        {
            new_cell_state = -1;
            return CustomRuleIsNull;
        }
        else
        {
            // custom_rule should set the new_cell_state
            custom_rule(cell_index, index_size, neighborhood_cells, neighborhood_size, new_cell_state);
        }
        break;
    case Parity:
        for (int i = 0; i < neighborhood_size; i++)
        {
            // update sum with current cell value
            sum += neighborhood_cells[i];
        }
        new_cell_state = sum % num_states; // store the parity state as the new state
        break;
    case Majority:
        for (int i = 0; i < neighborhood_size; i++)
        {
            // increment the cell state's number of votes
            auto it = state_votes_counter.find(neighborhood_cells[i]);
            if (it != state_votes_counter.end())
            {
                it->second++;
            }
        }
        // find the max_element in the counter based on the pairs' second variable (number of "votes")
        auto max_elem = max_element(state_votes_counter.begin(), state_votes_counter.end(),
                                    [](const std::pair<int, int> &a, const std::pair<int, int> &b)
                                    { return a.second < b.second; });
        new_cell_state = max_elem->first; // set the majority state as the new state
        break;
    }
    return 0;
}

/**
 * @brief Generates an array of neighboring cells and then calls set_new_cell_state to set the state.
 * This method supports a vector of cell states.
 *
 * @param cell_index cell of interest's index
 * @param index_size size of cell_index array
 * @param new_cell_state reference variable for setting the new state
 * @param custom_rule function that is called when a Custom rule type is specified
 *@return int - error code\n
 * Error codes returned by set_new_cell_state\n
 * 0: no error
 */
int CellularAutomata::get_state_from_neighborhood_1d(int *cell_index, int index_size, int &new_cell_state,
                                                     void(custom_rule)(int *, int, int *, int, int &))
{
    int error_code = 0;         // store error code return by other methods
    int i = cell_index[0];      // get i-th index from array
    int periodic_index;         // used by Periodic boundary type
    int neighborhood_size;      // number of neighbors in neighborhood
    int neighborhood_index = 0; // keep track of neighborhood array as we iterate through CA cells

    /*
     * Generate a flatten array of the cell's neighborhood.
     * The neighborhood array can then be utilized for Majority, Parity, or Custom rule
     */
    neighborhood_size = get_neighborhood_size(index_size, boundary_radius, neighborhood_type);
    // allocate memory and check if operation was successful
    int *neighborhood_cells = new (std::nothrow) int[neighborhood_size];
    if (neighborhood_cells == nullptr)
    {
        return NeighborhoodCellsMalloc;
    }

    // VonNeumann and Moore do not differ for 1d (vector) case
    switch (boundary_type)
    {
    case Periodic:
        for (int di = -boundary_radius; di <= boundary_radius; di++)
        {
            periodic_index = get_periodic_index_axis(i, di, axis1_dim);
            // add state to flattened array
            neighborhood_cells[neighborhood_index] = vector[periodic_index];
            neighborhood_index++;
        }
        error_code = set_new_cell_state(cell_index, index_size, neighborhood_cells, neighborhood_size, new_cell_state, custom_rule);
        break;
    case Walled: // with walled boundaries the edge cells never change
        // check if i is a boundary cell
        if (i == 0 || i == axis1_dim - 1)
        {
            new_cell_state = vector[i]; // keep boundary cell state
            break;
        }
        // fall-through to CutOff case when cell isn't a boundary cell
    case CutOff:
        for (int di = -boundary_radius; di <= boundary_radius; di++)
        {
            // exclude cells that are out of bounds
            if (di <= 0 || di >= axis1_dim)
            {
                // outside bounds; don't include cell state in the sum/counter
                continue;
            }
            // add state to flattened array
            neighborhood_cells[neighborhood_index] = vector[di];
            neighborhood_index++;
        }
        error_code = set_new_cell_state(cell_index, index_size, neighborhood_cells, neighborhood_size, new_cell_state, custom_rule);
        break;
    }
    delete[] neighborhood_cells;
    return error_code;
}

/**
 * @brief Generates an array of neighboring cells and then calls set_new_cell_state to set the state.
 * This method supports a matrix of cell states.
 *
 * @param cell_index cell of interest's index
 * @param index_size size of cell_index array
 * @param new_cell_state reference variable for setting the new state
 * @param custom_rule function that is called when a Custom rule type is specified
 *@return int - error code\n
 * Error codes returned by set_new_cell_state\n
 * 0: no error
 */
int CellularAutomata::get_state_from_neighborhood_2d(int *cell_index, int index_size, int &new_cell_state,
                                                     void(custom_rule)(int *, int, int *, int, int &))
{
    int error_code = 0;         // store error code return by other methods
    int i = cell_index[0];      // get i-th index from array
    int j = cell_index[1];      // get j-th index from array
    int periodic_index1;        // axis1_dim index used by Periodic boundary type
    int periodic_index2;        // axis2_dim index used by Periodic boundary type
    int neighborhood_size;      // number of neighbors in neighborhood
    int neighborhood_index = 0; // keep track of neighborhood array as we iterate through CA cells

    /*
     * Generate a flatten array of the cell's neighborhood.
     * The neighborhood array can then be utilized for Majority, Parity, or Custom rule
     */
    neighborhood_size = get_neighborhood_size(index_size, boundary_radius, neighborhood_type);
    // allocate memory and check if operation was successful
    int *neighborhood_cells = new (std::nothrow) int[neighborhood_size];
    if (neighborhood_cells == nullptr)
    {
        return NeighborhoodCellsMalloc;
    }

    switch (boundary_type)
    {
    case Periodic:
        for (int di = -boundary_radius; di <= boundary_radius; di++)
        {
            periodic_index1 = get_periodic_index_axis(i, di, axis1_dim);
            for (int dj = -boundary_radius; dj <= boundary_radius; dj++)
            {
                // exclude diagonal cells from Moore neighborhood when VonNeumann is selected
                if ((neighborhood_type == VonNeumann) && is_diagonal_neighboring_cell_2d(di, dj))
                {
                    // current di,dj cell is a diagonal neighbor so exclude it
                    continue;
                }
                periodic_index2 = get_periodic_index_axis(j, dj, axis2_dim);
                // add state to flattened array
                neighborhood_cells[neighborhood_index] = matrix[periodic_index1][periodic_index2];
                neighborhood_index++;
            }
        }
        error_code = set_new_cell_state(cell_index, index_size, neighborhood_cells, neighborhood_size, new_cell_state, custom_rule);
        break;
    case Walled: // with walled boundaries the edge cells never change
        // check if i,j is a boundary cell
        if ((i == 0 || i == axis1_dim - 1) || (j == 0 || j == axis2_dim - 1))
        {
            new_cell_state = matrix[i][j]; // keep boundary cell state
            break;
        }
        // fall-through to CutOff case when cell isn't a boundary cell
    case CutOff:
        for (int di = -boundary_radius; di <= boundary_radius; di++)
        {
            for (int dj = -boundary_radius; dj <= boundary_radius; dj++)
            {
                // exclude diagonal cells from Moore neighborhood when VonNeumann is selected
                if ((neighborhood_type == VonNeumann) && is_diagonal_neighboring_cell_2d(di, dj))
                {
                    // current di,dj cell is a diagonal neighbor so exclude it
                    continue;
                }
                // exclude cells that are out of bounds
                if ((di <= 0 || di >= axis1_dim) || (dj <= 0 || dj >= axis2_dim))
                {
                    // outside bounds; don't include cell state in the sum/counter
                    continue;
                }
                // add state to flattened array
                neighborhood_cells[neighborhood_index] = matrix[di][dj];
                neighborhood_index++;
            }
        }
        error_code = set_new_cell_state(cell_index, index_size, neighborhood_cells, neighborhood_size, new_cell_state, custom_rule);
        break;
    }
    delete[] neighborhood_cells;
    return error_code;
}

/**
 * @brief Generates an array of neighboring cells and then calls set_new_cell_state to set the state.
 * This method supports a tensor of cell states.
 *
 * @param cell_index cell of interest's index
 * @param index_size size of cell_index array
 * @param new_cell_state reference variable for setting the new state
 * @param custom_rule function that is called when a Custom rule type is specified
 *@return int - error code\n
 * Error codes returned by set_new_cell_state\n
 * 0: no error
 */
int CellularAutomata::get_state_from_neighborhood_3d(int *cell_index, int index_size, int &new_cell_state,
                                                     void(custom_rule)(int *, int, int *, int, int &))
{
    int error_code = 0;         // store error code return by other methods
    int i = cell_index[0];      // get i-th index from array
    int j = cell_index[1];      // get j-th index from array
    int k = cell_index[2];      // get k-th index from array
    int periodic_index1;        // axis1_dim index used by Periodic boundary type
    int periodic_index2;        // axis2_dim index used by Periodic boundary type
    int periodic_index3;        // axis3_dim index used by Periodic boundary type
    int neighborhood_size;      // number of neighbors in neighborhood
    int neighborhood_index = 0; // keep track of neighborhood array as we iterate through CA cells

    /*
     * Generate a flatten array of the cell's neighborhood.
     * The neighborhood array can then be utilized for Majority, Parity, or Custom rule
     */
    neighborhood_size = get_neighborhood_size(index_size, boundary_radius, neighborhood_type);
    // allocate memory and check if operation was successful
    int *neighborhood_cells = new (std::nothrow) int[neighborhood_size];
    if (neighborhood_cells == nullptr)
    {
        return NeighborhoodCellsMalloc;
    }

    switch (boundary_type)
    {
    case Periodic:
        for (int di = -boundary_radius; di <= boundary_radius; di++)
        {
            periodic_index1 = get_periodic_index_axis(i, di, axis1_dim);
            for (int dj = -boundary_radius; dj <= boundary_radius; dj++)
            {
                periodic_index2 = get_periodic_index_axis(j, dj, axis2_dim);
                for (int dk = -boundary_radius; dk <= boundary_radius; dk++)
                {
                    // exclude diagonal cells from Moore neighborhood when VonNeumann is selected
                    if ((neighborhood_type == VonNeumann) && is_diagonal_neighboring_cell_3d(di, dj, dk))
                    {
                        // current di,dj,dk cell is a diagonal neighbor so exclude it
                        continue;
                    }
                    periodic_index3 = get_periodic_index_axis(k, dk, axis3_dim);
                    // add state to flattened array
                    neighborhood_cells[neighborhood_index] = tensor[periodic_index1][periodic_index2][periodic_index3];
                    neighborhood_index++;
                }
            }
        }
        error_code = set_new_cell_state(cell_index, index_size, neighborhood_cells, neighborhood_size, new_cell_state, custom_rule);
        break;
    case Walled:
        // check if i,j,k is a boundary cell
        if ((i == 0 || i == axis1_dim - 1) || (j == 0 || j == axis2_dim - 1) || (k == 0 || j == axis3_dim - 1))
        {
            // with walled boundaries the edge cells never change
            new_cell_state = tensor[i][j][k]; // keep boundary cell state
            break;
        }
        // fall-through to CutOff case when cell isn't a boundary cell
    case CutOff:
        for (int di = -boundary_radius; di <= boundary_radius; di++)
        {
            for (int dj = -boundary_radius; dj <= boundary_radius; dj++)
            {
                for (int dk = -boundary_radius; dk <= boundary_radius; dk++)
                {
                    // exclude diagonal cells from Moore neighborhood when VonNeumann is selected
                    if ((neighborhood_type == VonNeumann) && is_diagonal_neighboring_cell_3d(di, dj, dk))
                    {
                        // current di,dj,dk cell is a diagonal neighbor so exclude it
                        continue;
                    }
                    // exclude cells that are out of bounds
                    if ((di <= 0 || di >= axis1_dim) || (dj <= 0 || dj >= axis2_dim) || (dk <= 0 || dk >= axis3_dim))
                    {
                        // outside bounds; don't include cell state in the sum/counter
                        continue;
                    }
                    // add state to flattened array
                    neighborhood_cells[neighborhood_index] = tensor[di][dj][dk];
                    neighborhood_index++;
                }
            }
        }
        error_code = set_new_cell_state(cell_index, index_size, neighborhood_cells, neighborhood_size, new_cell_state, custom_rule);
        break;
    }
    delete[] neighborhood_cells;
    return error_code;
}

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
int CellularAutomata::step(void(custom_rule)(int *, int, int *, int, int &))
{
    int error_code = 0; // store error code return by other methods
    int new_cell_state; // stores the cell's new state
    int index_size;     // number of indices required to address the cell

    if (vector != nullptr)
    {
        // initialize index_size for a vector and declare the cell_index variable
        index_size = 1;
        int cell_index[index_size];
        for (int i = 0; i < axis1_dim; i++)
        {
            cell_index[0] = i; // store the i-th index
            // store the main cell's index (for custom rule type)
            error_code = get_state_from_neighborhood_1d(cell_index, index_size, new_cell_state, custom_rule);
            if (error_code < 0)
            {
                return error_code;
            }
            next_vector[i] = new_cell_state;
        }
        // store next cell state to the current cell state for the next time step
        swap_states(vector, next_vector, axis1_dim);
    }
    else if (matrix != nullptr)
    {
        // initialize index_size for a matrix and declare the cell_index variable
        index_size = 2;
        int cell_index[index_size];
        for (int i = 0; i < axis1_dim; i++)
        {
            cell_index[0] = i; // store the i-th index
            for (int j = 0; j < axis2_dim; j++)
            {
                cell_index[1] = j; // store the j-th index
                error_code = get_state_from_neighborhood_2d(cell_index, index_size, new_cell_state, custom_rule);
                if (error_code < 0)
                {
                    return error_code;
                }
                next_matrix[i][j] = new_cell_state;
            }
        }
        // store next cell state to the current cell state for the next time step
        swap_states(matrix, next_matrix, axis1_dim, axis2_dim);
    }
    else if (tensor != nullptr)
    {
        // initialize index_size for a tensor and declare the cell_index variable
        index_size = 3;
        int cell_index[index_size];
        for (int i = 0; i < axis1_dim; i++)
        {
            cell_index[0] = i; // store the i-th index
            for (int j = 0; j < axis2_dim; j++)
            {
                cell_index[1] = j; // store the j-th index
                for (int k = 0; k < axis3_dim; k++)
                {
                    cell_index[2] = k; // store the k-th index
                    error_code = get_state_from_neighborhood_3d(cell_index, index_size, new_cell_state, custom_rule);
                    if (error_code < 0)
                    {
                        return error_code;
                    }
                    next_tensor[i][j][k] = new_cell_state;
                }
            }
        }
        // store next cell state to the current cell state for the next time step
        swap_states(tensor, next_tensor, axis1_dim, axis2_dim, axis3_dim);
    }
    else
    {
        return CellsAreNull;
    }

    steps_taken++;
    return error_code;
}

/**
 * @brief Simulates a cellular automata step.
 * A new state is generated and stored stored as the new state for subsequent calls to step method.
 *
 *@return int - error code\n
 * Error codes returned by get_state_from_neighborhood_Xd methods\n
 * 0: no error
 */
int CellularAutomata::step()
{
    return step(nullptr); // return step(func) error code
}

/**
 * @brief Print the current state of the grid.
 *
 * @return int - error code\n
 * CellsAreNull: neither vector, matrix, nor tensor are initialized\n
 * 0: no error
 */
int CellularAutomata::print_grid()
{
    if (vector != nullptr)
    {
        for (int i = 0; i < axis1_dim; i++)
        {
            std::cout << vector[i] << " ";
        }
        std::cout << std::endl;
    }
    else if (matrix != nullptr)
    {
        for (int j = 0; j < axis1_dim; j++)
        {
            for (int k = 0; k < axis2_dim; k++)
            {
                std::cout << matrix[j][k] << " ";
            }
            std::cout << std::endl;
        }
    }
    else if (tensor != nullptr)
    {
        for (int i = 0; i < axis1_dim; i++)
        {
            std::cout << "Printing " << i << "'th slice of Tensor" << std::endl;
            for (int j = 0; j < axis2_dim; j++)
            {
                for (int k = 0; k < axis3_dim; k++)
                {
                    std::cout << tensor[i][j][k] << " ";
                }
                std::cout << std::endl;
            }
        }
    }
    else
    {
        return CellsAreNull;
    }
    return 0;
}

/**
 * @brief Prints an error message for the given error code
 *
 * @param error error_code enum type
 */
void CellularAutomata::print_error_status(error_code error)
{
    switch (error)
    {
    case CellsAlreadyInitialized:
        std::cout << "ERROR [" << error << "]: Can't reinitialize vector, matrix, nor tensor. \n";
        break;
    case CellsAreNull:
        std::cout << "ERROR [" << error << "]: The vector, matrix, and tensor are null. \n";
        break;
    case CellsMalloc:
        std::cout << "ERROR [" << error << "]: Could not allocate memory for either vector, matrix, or tensor. \n";
        break;
    case InvalidCellState:
        std::cout << "ERROR [" << error << "]: Invalid cell state given. Must be greater than equal to 2. \n";
        break;
    case InvalidCellStateCondition:
        std::cout << "ERROR [" << error << "]: Invalid cell state condition given. Must be less than the set number of states. \n";
        break;
    case InvalidRadius:
        std::cout << "ERROR [" << error << "]: Invalid boundary radius given. \n";
        break;
    case InvalidNumStates:
        std::cout << "ERROR [" << error << "]: Invalid number of states given. \n";
        break;
    case NeighborhoodCellsMalloc:
        std::cout << "ERROR [" << error << "]: Could not allocate memory for neighborhood array. \n";
        break;
    case CustomRuleIsNull:
        std::cout << "ERROR [" << error << "]: Custom rule function is null (none given). \n";
        break;
    }
}
