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

#include "CAutils.h"
#include <array>
#include <unordered_map>
#include <iostream>
#include <algorithm> // max_element
#include <utility>   // pair
#include <cmath>     // pow

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
        CustomRuleIsNull = -9,
        RadiusLargerThanDimensions = -10
    };
}

/**
 * @brief A base CellularAutomata class that contains non-templated member variables and method definitions
 * from which templated and specialized template classes can inherit.
 * Default values:<br>
 * &emsp;&emsp; boundary_type = CAEnums::Periodic <br>
 * &emsp;&emsp; neighborhood_type = CAEnums::Moore <br>
 * &emsp;&emsp; rule_type = CAEnums::Majority <br>
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
     * @brief Defines the range of cell states to be used in the CA object.
     *
     * @param num_states describes the range of numbers to use for cell states
     * @return int - error code\n
     * InvalidNumStates: num_states can't be less than to 2\n
     * 0: no error
     */
    int setup_cell_states(int num_states);

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
 * This templated class supports structs/class objects or integer states.
 *
 * The object must contain a .state property that the CellularAutomata class updates.<br>
 * The object must have a default constructor or default initialzed values.<br>
 * The object must meet the requirements of CopyConstructible and CopyAssignable (until C++11)MoveConstructible and MoveAssignable (since C++11)<br>
 * The object must have a defined != operator for the comparison of cell states.<br>
 *
 * If these requirements are not satisfied then the class will produce undefined behavior.
 *
 * @tparam T : struct/class with a .state property, move operator and assignment operator.<br>
 * @tparam int : when cell states are represented by an integer
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
                           T &new_cell_state, void(custom_rule)(int *, int, T *, int, T &))
    {
        int sum = 0;                         // sum of cells within boundary_radius for Parity rule
        MajorityCounter state_votes_counter; // counter to keep track of votes for Majority rule
        initialize_majority_rule_counter(state_votes_counter, num_states);

        switch (rule_type)
        {
        case CAEnums::Custom:
            if (vector == nullptr && matrix == nullptr && tensor == nullptr)
            {
                new_cell_state.state = -1;
                return CAEnums::CellsAreNull;
            }
            if (custom_rule == nullptr)
            {
                new_cell_state.state = -1;
                return CAEnums::CustomRuleIsNull;
            }
            else
            {
                // custom_rule should set the new_cell_state
                custom_rule(cell_index, index_size, neighborhood_cells, neighborhood_size, new_cell_state);
            }
            break;
        case CAEnums::Parity:
            for (int i = 0; i < neighborhood_size; i++)
            {
                // update sum with current cell value
                sum += neighborhood_cells[i].state;
            }
            new_cell_state.state = sum % num_states; // store the parity state as the new state
            break;
        case CAEnums::Majority:
            for (int i = 0; i < neighborhood_size; i++)
            {
                // increment the cell state's number of votes
                auto it = state_votes_counter.find(neighborhood_cells[i].state);
                if (it != state_votes_counter.end())
                {
                    it->second++;
                }
            }
            // find the max_element in the counter based on the pairs' second variable (number of "votes")
            auto max_elem = max_element(state_votes_counter.begin(), state_votes_counter.end(),
                                        less_than_votes);
            new_cell_state.state = max_elem->first; // set the majority state as the new state
            break;
        }
        return 0;
    }

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
                                       void(custom_rule)(int *, int, T *, int, T &))
    {
        int error_code = 0;         // store error code return by other methods
        int i = cell_index[0];      // get i-th index from array
        int neighborhood_index = 0; // keep track of neighborhood array as we iterate through CA cells

        // allocate memory and check if operation was successful
        T *neighborhood_cells = malloc_neighborhood_array(index_size);
        if (neighborhood_cells == nullptr)
        {
            return CAEnums::NeighborhoodCellsMalloc;
        }

        // VonNeumann and Moore do not differ for 1d (vector) case
        switch (boundary_type)
        {
        case CAEnums::Periodic:
            generate_periodic_neighborhood_1d(cell_index, neighborhood_cells, neighborhood_index);
            error_code = set_new_cell_state(cell_index, index_size, neighborhood_cells, neighborhood_index, new_cell_state, custom_rule);
            break;
        case CAEnums::Walled: // with walled boundaries the edge cells never change
            // check if i is a boundary cell
            if (i == 0 || i == axis1_dim - 1)
            {
                new_cell_state = vector[i]; // keep boundary cell state
                break;
            }
            // fall-through to CutOff case when cell isn't a boundary cell
        case CAEnums::CutOff:
            generate_cutoff_neighborhood_1d(cell_index, neighborhood_cells, neighborhood_index);
            error_code = set_new_cell_state(cell_index, index_size, neighborhood_cells, neighborhood_index, new_cell_state, custom_rule);
            break;
        }
        delete[] neighborhood_cells;
        return error_code;
    }

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
                                       void(custom_rule)(int *, int, T *, int, T &))
    {
        int error_code = 0;         // store error code return by other methods
        int i = cell_index[0];      // get i-th index from array
        int j = cell_index[1];      // get j-th index from array
        int neighborhood_index = 0; // keep track of neighborhood array as we iterate through CA cells

        // allocate memory and check if operation was successful
        T *neighborhood_cells = malloc_neighborhood_array(index_size);
        if (neighborhood_cells == nullptr)
        {
            return CAEnums::NeighborhoodCellsMalloc;
        }

        switch (boundary_type)
        {
        case CAEnums::Periodic:
            generate_periodic_neighborhood_2d(cell_index, neighborhood_cells, neighborhood_index);
            error_code = set_new_cell_state(cell_index, index_size, neighborhood_cells, neighborhood_index, new_cell_state, custom_rule);
            break;
        case CAEnums::Walled: // with walled boundaries the edge cells never change
            // check if i,j is a boundary cell
            if ((i == 0 || i == axis1_dim - 1) || (j == 0 || j == axis2_dim - 1))
            {
                new_cell_state = matrix[i][j]; // keep boundary cell state
                break;
            }
            // fall-through to CutOff case when cell isn't a boundary cell
        case CAEnums::CutOff:
            generate_cutoff_neighborhood_2d(cell_index, neighborhood_cells, neighborhood_index);
            error_code = set_new_cell_state(cell_index, index_size, neighborhood_cells, neighborhood_index, new_cell_state, custom_rule);
            break;
        }
        delete[] neighborhood_cells;
        return error_code;
    }

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
                                       void(custom_rule)(int *, int, T *, int, T &))
    {
        int error_code = 0;         // store error code return by other methods
        int i = cell_index[0];      // get i-th index from array
        int j = cell_index[1];      // get j-th index from array
        int k = cell_index[2];      // get k-th index from array
        int neighborhood_index = 0; // keep track of neighborhood array as we iterate through CA cells

        // allocate memory and check if operation was successful
        T *neighborhood_cells = malloc_neighborhood_array(index_size); // index_size == rank
        if (neighborhood_cells == nullptr)
        {
            return CAEnums::NeighborhoodCellsMalloc;
        }

        switch (boundary_type)
        {
        case CAEnums::Periodic:
            generate_periodic_neighborhood_3d(cell_index, neighborhood_cells, neighborhood_index);
            error_code = set_new_cell_state(cell_index, index_size, neighborhood_cells, neighborhood_index, new_cell_state, custom_rule);
            break;
        case CAEnums::Walled:
            // check if i,j,k is a boundary cell
            if ((i == 0 || i == axis1_dim - 1) || (j == 0 || j == axis2_dim - 1) || (k == 0 || j == axis3_dim - 1))
            {
                // with walled boundaries the edge cells never change
                new_cell_state = tensor[i][j][k]; // keep boundary cell state
                break;
            }
            // fall-through to CutOff case when cell isn't a boundary cell
        case CAEnums::CutOff:
            generate_cutoff_neighborhood_3d(cell_index, neighborhood_cells, neighborhood_index);
            error_code = set_new_cell_state(cell_index, index_size, neighborhood_cells, neighborhood_index, new_cell_state, custom_rule);
            break;
        }
        delete[] neighborhood_cells;
        return error_code;
    }

    /**
     * @brief Adds periodic neighboring cell state to neighborhood_cells array.
     * Handles vector (1d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_periodic_neighborhood_1d(int *cell_index, T *neighborhood_cells, int &neighborhood_index)
    {
        int i = cell_index[0]; // get i-th index from array
        int periodic_i;        // axis1_dim index used by Periodic boundary type

        for (int di = -boundary_radius; di <= boundary_radius; di++)
        {
            periodic_i = get_periodic_index(i, di, axis1_dim);
            // add state to flattened array
            neighborhood_cells[neighborhood_index] = vector[periodic_i];
            neighborhood_index++;
        }
    }

    /**
     * @brief Adds cutoff neighboring cell state to neighborhood_cells array.
     * Handles vector (1d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_cutoff_neighborhood_1d(int *cell_index, T *neighborhood_cells, int &neighborhood_index)
    {
        int i = cell_index[0]; // get i-th index from array
        int neighbor_i;        // neighboring cell's i-th index

        for (int di = -boundary_radius; di <= boundary_radius; di++)
        {
            neighbor_i = di + i;
            // exclude cells that are out of bounds
            if (neighbor_i <= 0 || neighbor_i >= axis1_dim)
            {
                // outside bounds; don't include cell state in the sum/counter
                continue;
            }
            // add state to flattened array
            neighborhood_cells[neighborhood_index] = vector[di];
            neighborhood_index++;
        }
    }

    /**
     * @brief Adds periodic neighboring cell state to neighborhood_cells array.
     * Handles matrix (2d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_periodic_neighborhood_2d(int *cell_index, T *neighborhood_cells, int &neighborhood_index)
    {
        int i = cell_index[0]; // get i-th index from array
        int j = cell_index[1]; // get j-th index from array
        int periodic_i;        // axis1_dim index used by Periodic boundary type
        int periodic_j;        // axis2_dim index used by Periodic boundary type

        if (neighborhood_type == CAEnums::VonNeumann)
        {
            for (int di = -boundary_radius; di <= boundary_radius; di++)
            {
                periodic_i = get_periodic_index(i, di, axis1_dim);
                for (int dj = -boundary_radius; dj <= boundary_radius; dj++)
                {
                    // exclude diagonal cells from neighborhood when VonNeumann is selected
                    if (is_diagonal_neighboring_cell_2d(di, dj))
                    {
                        // current di,dj cell is a diagonal neighbor so exclude it
                        continue;
                    }
                    periodic_j = get_periodic_index(j, dj, axis2_dim);
                    // add state to flattened array
                    neighborhood_cells[neighborhood_index] = matrix[periodic_i][periodic_j];
                    neighborhood_index++;
                }
            }
        }
        else // Moore neighborhood
        {
            for (int di = -boundary_radius; di <= boundary_radius; di++)
            {
                periodic_i = get_periodic_index(i, di, axis1_dim);
                for (int dj = -boundary_radius; dj <= boundary_radius; dj++)
                {
                    periodic_j = get_periodic_index(j, dj, axis2_dim);
                    // add state to flattened array
                    neighborhood_cells[neighborhood_index] = matrix[periodic_i][periodic_j];
                    neighborhood_index++;
                }
            }
        }
    }

    /**
     * @brief Adds cutoff neighboring cell state to neighborhood_cells array.
     * Handles matrix (2d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_cutoff_neighborhood_2d(int *cell_index, T *neighborhood_cells, int &neighborhood_index)
    {
        int i = cell_index[0]; // get i-th index from array
        int j = cell_index[1]; // get j-th index from array
        int neighbor_i;        // neighboring cell's i-th index
        int neighbor_j;        // neighboring cell's j-th index

        if (neighborhood_type == CAEnums::VonNeumann)
        {
            for (int di = -boundary_radius; di <= boundary_radius; di++)
            {
                neighbor_i = di + i;
                for (int dj = -boundary_radius; dj <= boundary_radius; dj++)
                {
                    neighbor_j = dj + j;
                    // exclude diagonal cells from neighborhood when VonNeumann is selected
                    if (is_diagonal_neighboring_cell_2d(di, dj))
                    {
                        // current di,dj cell is a diagonal neighbor so exclude it
                        continue;
                    }
                    // exclude cells that are out of bounds
                    if ((neighbor_i <= 0 || neighbor_i >= axis1_dim) || (neighbor_j <= 0 || neighbor_j >= axis2_dim))
                    {
                        // outside bounds; don't include cell state in the sum/counter
                        continue;
                    }
                    // add state to flattened array
                    neighborhood_cells[neighborhood_index] = matrix[neighbor_i][neighbor_j];
                    neighborhood_index++;
                }
            }
        }
        else // Moore neighborhood
        {
            for (int di = -boundary_radius; di <= boundary_radius; di++)
            {
                neighbor_i = di + i;
                for (int dj = -boundary_radius; dj <= boundary_radius; dj++)
                {
                    neighbor_j = dj + j;
                    // exclude cells that are out of bounds
                    if ((neighbor_i <= 0 || neighbor_i >= axis1_dim) || (neighbor_j <= 0 || neighbor_j >= axis2_dim))
                    {
                        // outside bounds; don't include cell state in the sum/counter
                        continue;
                    }
                    // add state to flattened array
                    neighborhood_cells[neighborhood_index] = matrix[neighbor_i][neighbor_j];
                    neighborhood_index++;
                }
            }
        }
    }

    /**
     * @brief Adds periodic neighboring cell state to neighborhood_cells array.
     * Handles tensor (3d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_periodic_neighborhood_3d(int *cell_index, T *neighborhood_cells, int &neighborhood_index)
    {
        int i = cell_index[0]; // get i-th index from array
        int j = cell_index[1]; // get j-th index from array
        int k = cell_index[2]; // get k-th index from array
        int periodic_i;        // axis1_dim index used by Periodic boundary type
        int periodic_j;        // axis2_dim index used by Periodic boundary type
        int periodic_k;        // axis3_dim index used by Periodic boundary type

        if (neighborhood_type == CAEnums::VonNeumann)
        {
            for (int di = -boundary_radius; di <= boundary_radius; di++)
            {
                periodic_i = get_periodic_index(i, di, axis1_dim);
                for (int dj = -boundary_radius; dj <= boundary_radius; dj++)
                {
                    periodic_j = get_periodic_index(j, dj, axis2_dim);
                    for (int dk = -boundary_radius; dk <= boundary_radius; dk++)
                    {
                        // exclude diagonal cells from neighborhood when VonNeumann is selected
                        if (is_diagonal_neighboring_cell_3d(di, dj, dk))
                        {
                            // current di,dj,dk cell is a diagonal neighbor so exclude it
                            continue;
                        }
                        periodic_k = get_periodic_index(k, dk, axis3_dim);
                        // add state to flattened array
                        neighborhood_cells[neighborhood_index] = tensor[periodic_i][periodic_j][periodic_k];
                        neighborhood_index++;
                    }
                }
            }
        }
        else // Moore neighborhood
        {
            for (int di = -boundary_radius; di <= boundary_radius; di++)
            {
                periodic_i = get_periodic_index(i, di, axis1_dim);
                for (int dj = -boundary_radius; dj <= boundary_radius; dj++)
                {
                    periodic_j = get_periodic_index(j, dj, axis2_dim);
                    for (int dk = -boundary_radius; dk <= boundary_radius; dk++)
                    {
                        periodic_k = get_periodic_index(k, dk, axis3_dim);
                        // add state to flattened array
                        neighborhood_cells[neighborhood_index] = tensor[periodic_i][periodic_j][periodic_k];
                        neighborhood_index++;
                    }
                }
            }
        }
    }

    /**
     * @brief Adds cutoff neighboring cell state to neighborhood_cells array.
     * Handles tensor (3d) case.
     *
     * @param cell_index cell of interest's index
     * @param neighborhood_cells array containing neighboring cell states
     * @param neighborhood_index keep track of number of cell states added
     */
    void generate_cutoff_neighborhood_3d(int *cell_index, T *neighborhood_cells, int &neighborhood_index)
    {
        int i = cell_index[0]; // get i-th index from array
        int j = cell_index[1]; // get j-th index from array
        int k = cell_index[2]; // get k-th index from array
        int neighbor_i;        // neighboring cell's i-th index
        int neighbor_j;        // neighboring cell's j-th index
        int neighbor_k;        // neighboring cell's k-th index

        if (neighborhood_type == CAEnums::VonNeumann)
        {
            for (int di = -boundary_radius; di <= boundary_radius; di++)
            {
                neighbor_i = di + i;
                for (int dj = -boundary_radius; dj <= boundary_radius; dj++)
                {
                    neighbor_j = dj + j;
                    for (int dk = -boundary_radius; dk <= boundary_radius; dk++)
                    {
                        neighbor_k = dk + k;
                        // exclude diagonal cells from neighborhood when VonNeumann is selected
                        if (is_diagonal_neighboring_cell_3d(di, dj, dk))
                        {
                            // current di,dj,dk cell is a diagonal neighbor so exclude it
                            continue;
                        }
                        // exclude cells that are out of bounds
                        if ((neighbor_i <= 0 || neighbor_i >= axis1_dim) ||
                            (neighbor_j <= 0 || neighbor_j >= axis2_dim) ||
                            (neighbor_k <= 0 || neighbor_k >= axis3_dim))
                        {
                            // outside bounds; don't include cell state in the sum/counter
                            continue;
                        }
                        // add state to flattened array
                        neighborhood_cells[neighborhood_index] = tensor[neighbor_i][neighbor_j][neighbor_k];
                        neighborhood_index++;
                    }
                }
            }
        }
        else // Moore neighborhood
        {
            for (int di = -boundary_radius; di <= boundary_radius; di++)
            {
                neighbor_i = di + i;
                for (int dj = -boundary_radius; dj <= boundary_radius; dj++)
                {
                    neighbor_j = dj + j;
                    for (int dk = -boundary_radius; dk <= boundary_radius; dk++)
                    {
                        neighbor_k = dk + k;
                        // exclude cells that are out of bounds
                        if ((neighbor_i <= 0 || neighbor_i >= axis1_dim) ||
                            (neighbor_j <= 0 || neighbor_j >= axis2_dim) ||
                            (neighbor_k <= 0 || neighbor_k >= axis3_dim))
                        {
                            // outside bounds; don't include cell state in the sum/counter
                            continue;
                        }
                        // add state to flattened array
                        neighborhood_cells[neighborhood_index] = tensor[neighbor_i][neighbor_j][neighbor_k];
                        neighborhood_index++;
                    }
                }
            }
        }
    }

    /**
     * @brief Dynamically allocates an array that will contain the neighborhood cell states.
     *
     * @param rank rank of the cells data
     * @return T* pointer to the dynamic array
     */
    T *malloc_neighborhood_array(int rank)
    {
        /*
         * Generate a flatten array of the cell's neighborhood.
         * The neighborhood array can then be utilized for Majority, Parity, or Custom rule
         */
        int max_neighborhood_size = get_neighborhood_size(rank, boundary_radius, neighborhood_type);
        // allocate memory and check if operation was successful
        T *neighborhood_cells = new (std::nothrow) T[max_neighborhood_size];
        return neighborhood_cells;
    }

public:
    /**
     * @brief Get the vector cell grid
     *
     * @return T*
     */
    T *get_vector()
    {
        return vector;
    }

    /**
     * @brief Get the next state vector cell grid
     *
     * @return T*
     */
    T *get_next_vector()
    {
        return next_vector;
    }

    /**
     * @brief Get the matrix cell grid
     *
     * @return T*
     */
    T **get_matrix()
    {
        return matrix;
    }

    /**
     * @brief Get the next state matrix cell grid
     *
     * @return T*
     */
    T **get_next_matrix()
    {
        return next_matrix;
    }

    /**
     * @brief Get the tensor cell grid
     *
     * @return T*
     */
    T ***get_tensor()
    {
        return tensor;
    }

    /**
     * @brief Get the next state tensor cell grid
     *
     * @return T*
     */
    T ***get_next_tensor()
    {
        return next_tensor;
    }

    /**
     * @brief Setup boundary with enum values from boundary and set the boundary radius.
     *
     * @param bound_type enum value for boundary (None, Periodic, Walled, CutOff)
     * @param radius radius for the boundary
     * @return int - error code\n
     * InvalidRadius: radius can't be less than equal to 0\n
     * RadiusLargerThanDimensions: radius must be smaller than half of the dimensions' size.
     * 0: no error
     */
    int setup_boundary(CAEnums::Boundary bound_type, int radius)
    {
        if (radius <= 0)
        {
            return CAEnums::InvalidRadius;
        }
        if (vector != nullptr)
        {
            if (radius > axis1_dim / 2)
            {
                return CAEnums::RadiusLargerThanDimensions;
            }
        }
        else if (matrix != nullptr)
        {
            if (radius > axis1_dim / 2 || radius > axis2_dim / 2)
            {
                return CAEnums::RadiusLargerThanDimensions;
            }
        }
        else if (tensor != nullptr)
        {
            if (radius > axis2_dim / 2 || radius > axis3_dim / 2)
            {
                return CAEnums::RadiusLargerThanDimensions;
            }
        }

        this->boundary_type = bound_type;
        this->boundary_radius = radius;
        return 0;
    }

    /**
     * @brief Construct a new Cellular Automata object.
     * Sets the default value to all class attributes.
     *
     */
    CellularAutomata() : BaseCellularAutomata()
    {
        steps_taken = 0;
        vector = nullptr;
        next_vector = nullptr;
        matrix = nullptr;
        next_matrix = nullptr;
        tensor = nullptr;
        next_tensor = nullptr;
    }

    /**
     * @brief Destroy the Cellular Automata object.
     * Deallocates memory reserved for vector/matrix/tensor.
     *
     */
    ~CellularAutomata()
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
     * @brief Set up 1d array of cell states.
     *
     * @param axis1_dim the size of a one dimensional vector
     * @param fill_value the value to set every cell state to
     * @return int - error code\n
     * CellsAlreadyInitialized: vector was already allocated\n
     * CellsMalloc: couldn't allocate memory for the specified vector size\n
     * 0: no error
     */
    int setup_dimensions_1d(int axis1_dim, int fill_value = 0)
    {
        if (vector != nullptr)
        {
            return CAEnums::CellsAlreadyInitialized;
        }

        this->axis1_dim = axis1_dim;
        vector = new (std::nothrow) T[axis1_dim];
        next_vector = new (std::nothrow) T[axis1_dim];

        if (vector == nullptr || next_vector == nullptr)
        {
            return CAEnums::CellsMalloc;
        }

        // initialize vector filled with zeros
        for (int j = 0; j < axis1_dim; j++)
        {
            vector[j].state = fill_value;
            next_vector[j].state = fill_value;
        }

        return 0;
    }

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
    int setup_dimensions_2d(int axis1_dim, int axis2_dim, int fill_value = 0)
    {
        if (matrix != nullptr)
        {
            return CAEnums::CellsAlreadyInitialized;
        }

        this->axis1_dim = axis1_dim;
        this->axis2_dim = axis2_dim;
        matrix = new (std::nothrow) T *[axis1_dim];
        next_matrix = new (std::nothrow) T *[axis1_dim];

        if (matrix == nullptr || next_matrix == nullptr)
        {
            return CAEnums::CellsMalloc;
        }

        for (int i = 0; i < axis1_dim; i++)
        {
            matrix[i] = new (std::nothrow) T[axis2_dim];
            next_matrix[i] = new (std::nothrow) T[axis2_dim];
        }

        // initialize matrix filled with zeros
        for (int j = 0; j < axis1_dim; j++)
        {
            for (int k = 0; k < axis2_dim; k++)
            {
                matrix[j][k].state = fill_value;
                next_matrix[j][k].state = fill_value;
            }
        }

        return 0;
    }

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
    int setup_dimensions_3d(int axis1_dim, int axis2_dim, int axis3_dim, int fill_value = 0)
    {
        if (tensor != nullptr)
        {
            return CAEnums::CellsAlreadyInitialized;
        }

        this->axis1_dim = axis1_dim;
        this->axis2_dim = axis2_dim;
        this->axis3_dim = axis3_dim;
        tensor = new (std::nothrow) T **[axis1_dim];
        next_tensor = new (std::nothrow) T **[axis1_dim];

        if (tensor == nullptr || next_tensor == nullptr)
        {
            return CAEnums::CellsMalloc;
        }

        for (int i = 0; i < axis1_dim; i++)
        {
            tensor[i] = new (std::nothrow) T *[axis2_dim];
            next_tensor[i] = new (std::nothrow) T *[axis2_dim];

            for (int j = 0; j < axis2_dim; j++)
            {
                tensor[i][j] = new (std::nothrow) T[axis3_dim];
                next_tensor[i][j] = new (std::nothrow) T[axis3_dim];
            }
        }

        // initialize matrix filled with zeros
        for (int i = 0; i < axis1_dim; i++)
        {
            for (int j = 0; j < axis2_dim; j++)
            {
                for (int k = 0; k < axis3_dim; k++)
                {
                    tensor[i][j][k].state = fill_value;
                    next_tensor[i][j][k].state = fill_value;
                }
            }
        }

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
    int init_condition(int x_state, double prob)
    {
        if (!(x_state < num_states))
        {
            return CAEnums::InvalidCellStateCondition;
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
                    vector[i].state = x_state;
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
                        matrix[j][k].state = x_state;
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
                            tensor[i][j][k].state = x_state;
                        }
                    }
                }
            }
        }
        else
        {
            return CAEnums::CellsAreNull;
        }
        return 0;
    }

    /**
     * @brief Simulates a cellular automata step.
     * A new state is generated and stored stored as the new state for subsequent calls to step method.
     *
     * This method supports the use of a custom rule type.
     *
     * If cell states move on the grid, the user is responsible for handling clashes.
     * If two cells move to the same cell position, the grid will retain the most
     * recent cell assignment (new will replace the old).
     *
     * @param custom_rule function that is called when a Custom rule type is specified
     *@return int - error code\n
     * Error codes returned by get_state_from_neighborhood_Xd methods\n
     * 0: no error
     */
    int step(void(custom_rule)(int *, int, T *, int, T &))
    {
        int error_code = 0; // store error code return by other methods
        T new_cell_state;   // stores the cell's new state
        T empty_cell_state; // cell state for zeroing out old states
        int index_size;     // number of indices required to address the cell

        if (vector != nullptr)
        {
            // store the main cell's index in cell_index for custom rule type
            index_size = 1;
#ifdef ENABLE_OMP
#pragma omp parallel for firstprivate(error_code) private(new_cell_state)
#endif
            for (int i = 0; i < axis1_dim; i++)
            {
                int cell_index[index_size] = {i};
                new_cell_state = vector[i];
                error_code = get_state_from_neighborhood_1d(cell_index, index_size, new_cell_state, custom_rule);
                if (error_code < 0)
                {
#ifdef ENABLE_OMP
#pragma omp cancel for
#endif
#ifndef ENABLE_OMP
                    /*
                     * Return error code when omp is not enabled
                     * otherwise the above pragma will exit out of the for loop(s)
                     */
                    return error_code;
#endif
                }
                /*
                 * The update cell if new_cell_state is no empty_state.
                 * Avoids overwriting the motion of cells.
                 */
                if (new_cell_state != empty_cell_state)
                {
                    next_vector[cell_index[0]] = new_cell_state;
                }
            }
            // store next cell state to the current cell state for the next time step
            swap_states<T>(vector, next_vector, axis1_dim);
        }
        else if (matrix != nullptr)
        {
            // store the main cell's index in cell_index for custom rule type
            index_size = 2;
            int cell_index[index_size];
#ifdef ENABLE_OMP
#pragma omp parallel for firstprivate(error_code) private(new_cell_state)
#endif
            for (int i = 0; i < axis1_dim; i++)
            {
                cell_index[0] = i; // store the i-th index
                for (int j = 0; j < axis2_dim; j++)
                {
                    int cell_index[index_size] = {i, j};
                    new_cell_state = matrix[i][j];
                    error_code = get_state_from_neighborhood_2d(cell_index, index_size, new_cell_state, custom_rule);
                    if (error_code < 0)
                    {
#ifdef ENABLE_OMP
#pragma omp cancel for
#endif
#ifndef ENABLE_OMP
                        /*
                         * Return error code when omp is not enabled
                         * otherwise the above pragma will exit out of the for loop(s)
                         */
                        return error_code;
#endif
                    }
                    /*
                     * The update cell if new_cell_state is no empty_state.
                     * Avoids overwriting the motion of cells.
                     */
                    if (new_cell_state != empty_cell_state)
                    {
                        next_matrix[cell_index[0]][cell_index[1]] = new_cell_state;
                    }
                }
            }
            // store next cell state to the current cell state for the next time step
            swap_states<T>(matrix, next_matrix, axis1_dim, axis2_dim);
        }
        else if (tensor != nullptr)
        {
            index_size = 3;
#ifdef ENABLE_OMP
#pragma omp parallel for firstprivate(error_code) private(new_cell_state)
#endif
            for (int i = 0; i < axis1_dim; i++)
            {
                for (int j = 0; j < axis2_dim; j++)
                {
                    for (int k = 0; k < axis3_dim; k++)
                    {
                        int cell_index[index_size] = {i, j, k};
                        new_cell_state = tensor[i][j][k];
                        error_code = get_state_from_neighborhood_3d(cell_index, index_size, new_cell_state, custom_rule);
                        if (error_code < 0)
                        {
#ifdef ENABLE_OMP
#pragma omp cancel for
#endif
#ifndef ENABLE_OMP
                            /*
                             * Return error code when omp is not enabled
                             * otherwise the above pragma will exit out of the for loop(s)
                             */
                            return error_code;
#endif
                        }
                        /*
                         * The update cell if new_cell_state is no empty_state.
                         * Avoids overwriting the motion of cells.
                         */
                        if (new_cell_state != empty_cell_state)
                        {
                            next_tensor[cell_index[0]][cell_index[1]][cell_index[2]] = new_cell_state;
                        }
                    }
                }
            }
            // store next cell state to the current cell state for the next time step
            swap_states<T>(tensor, next_tensor, axis1_dim, axis2_dim, axis3_dim);
        }
        else
        {
            return CAEnums::CellsAreNull;
        }

        steps_taken++;
        return error_code;
    }

    /**
     * @brief Simulates a cellular automata step.
     * A new state is generated and stored as the new state for subsequent calls to step method.
     *
     *@return int - error code\n
     * Error codes returned by get_state_from_neighborhood_Xd methods\n
     * 0: no error
     */
    int step()
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
    int print_grid()
    {
        if (vector != nullptr)
        {
            for (int i = 0; i < axis1_dim; i++)
            {
                std::cout << vector[i].state << " ";
            }
            std::cout << std::endl;
        }
        else if (matrix != nullptr)
        {
            for (int j = 0; j < axis1_dim; j++)
            {
                for (int k = 0; k < axis2_dim; k++)
                {
                    std::cout << matrix[j][k].state << " ";
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
                        std::cout << tensor[i][j][k].state << " ";
                    }
                    std::cout << std::endl;
                }
            }
        }
        else
        {
            return CAEnums::CellsAreNull;
        }
        return 0;
    }

    /**
     * @brief calculates the neighborhood size using the rank, radius and neighborhood type
     *
     * @param rank rank of the cell's tensor
     * @param radius radius of the neighborhood
     * @param neighborhood_type the type of neighborhood
     * @return int
     */
    static int get_neighborhood_size(int rank, int radius, CAEnums::Neighborhood neighborhood_type)
    {
        if (neighborhood_type == CAEnums::VonNeumann)
        {
            return (2 * rank * radius) + 1; // +1 to include cell of interest
        }
        else // CAEnums::Moore neighborhood
        {
            return pow((2 * radius + 1), rank);
        }
    }

    /**
     * @brief initializes a MajorityCounter instance utilized by Majority rule
     *
     * @param counter object to keep number of votes for each particular cell state
     * @param num_states number of different cell states
     */
    static void initialize_majority_rule_counter(MajorityCounter &counter, int num_states)
    {
        // sets the counter for every cell state type to 0
        for (int j = 0; j < num_states; j++)
        {
            counter.insert(std::make_pair(j, 0));
        }
    }
};

// Declare specialized int methods. Methods are defined in cellularautomata.cpp

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
template <>
int CellularAutomata<int>::setup_dimensions_1d(int axis1_dim, int fill_value);

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
template <>
int CellularAutomata<int>::setup_dimensions_2d(int axis1_dim, int axis2_dim, int fill_value);

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
template <>
int CellularAutomata<int>::setup_dimensions_3d(int axis1_dim, int axis2_dim, int axis3_dim, int fill_value);

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
template <>
int CellularAutomata<int>::init_condition(int x_state, double prob);

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
template <>
int CellularAutomata<int>::set_new_cell_state(int *cell_index, int index_size,
                                              int *neighborhood_cells, int neighborhood_size,
                                              int &new_cell_state, void(custom_rule)(int *, int, int *, int, int &));

/**
 * @brief Print the current state of the grid.
 *
 * @return int - error code\n
 * CellsAreNull: neither vector, matrix, nor tensor are initialized\n
 * 0: no error
 */
template <>
int CellularAutomata<int>::print_grid();

/**
 * @brief Simulates a cellular automata step.
 * A new state is generated and stored stored as the new state for subsequent calls to step method.
 *
 * This method supports the use of a custom rule type.
 *
 * If cell states move on the grid, the user is responsible for handling clashes.
 * If two cells move to the same cell position, the grid will retain the most
 * recent cell assignment (new will replace the old).
 *
 * @param custom_rule function that is called when a Custom rule type is specified
 *@return int - error code\n
 * Error codes returned by get_state_from_neighborhood_Xd methods\n
 * 0: no error
 */
template <>
int CellularAutomata<int>::step(void(custom_rule)(int *, int, int *, int, int &));
