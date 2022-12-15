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
#ifdef ENABLE_OMP
#include <omp.h>
#endif

BaseCellularAutomata::BaseCellularAutomata()
{
    axis1_dim = 0;
    axis2_dim = 0;
    axis3_dim = 0;
    num_states = 2;
    boundary_type = CAEnums::Periodic;
    boundary_radius = 1;
    neighborhood_type = CAEnums::Moore;
    rule_type = CAEnums::Majority;
}

int BaseCellularAutomata::setup_neighborhood(CAEnums::Neighborhood neighborhood_type)
{
    this->neighborhood_type = neighborhood_type;
    return 0;
}

int BaseCellularAutomata::setup_cell_states(int num_states)
{
    if (num_states < 2)
    {
        return CAEnums::InvalidNumStates;
    }
    this->num_states = num_states;
    return 0;
}

int BaseCellularAutomata::setup_rule(CAEnums::Rule rule_type)
{
    this->rule_type = rule_type;
    return 0;
}

void BaseCellularAutomata::print_error_status(CAEnums::ErrorCode error)
{
    std::cout << "ERROR [" << error;
    switch (error)
    {
    case CAEnums::CellsAlreadyInitialized:
        std::cout << "]: Can't reinitialize vector, matrix, nor tensor.";
        break;
    case CAEnums::CellsAreNull:
        std::cout << "]: The vector, matrix, and tensor are null.";
        break;
    case CAEnums::CellsMalloc:
        std::cout << "]: Could not allocate memory for either vector, matrix, or tensor.";
        break;
    case CAEnums::InvalidCellState:
        std::cout << "]: Invalid cell state given. Must be greater than equal to 2.";
        break;
    case CAEnums::InvalidCellStateCondition:
        std::cout << "]: Invalid cell state condition given. Must be less than the set number of states.";
        break;
    case CAEnums::InvalidRadius:
        std::cout << "]: Invalid boundary radius given.";
        break;
    case CAEnums::InvalidNumStates:
        std::cout << "]: Invalid number of states given.";
        break;
    case CAEnums::NeighborhoodCellsMalloc:
        std::cout << "]: Could not allocate memory for neighborhood array.";
        break;
    case CAEnums::CustomRuleIsNull:
        std::cout << "]: Custom rule function is null (none given).";
        break;
    case CAEnums::RadiusLargerThanDimensions:
        std::cout << "]: Boundary radius is smaller than one of the following: axis1_dim / 2, axis2_dim / 2, and/or axis3_dim / 2";
    }
    std::cout << "\n";
}

template <>
int CellularAutomata<int>::setup_dimensions_1d(int axis1_dim, int fill_value)
{
    if (vector != nullptr)
    {
        return CAEnums::CellsAlreadyInitialized;
    }

    this->axis1_dim = axis1_dim;
    vector = new (std::nothrow) int[axis1_dim];
    next_vector = new (std::nothrow) int[axis1_dim];

    if (vector == nullptr || next_vector == nullptr)
    {
        return CAEnums::CellsMalloc;
    }

    // initialize vector filled with zeros
    for (int j = 0; j < axis1_dim; j++)
    {
        vector[j] = fill_value;
        next_vector[j] = fill_value;
    }

    return 0;
}

template <>
int CellularAutomata<int>::setup_dimensions_2d(int axis1_dim, int axis2_dim, int fill_value)
{
    if (matrix != nullptr)
    {
        return CAEnums::CellsAlreadyInitialized;
    }

    this->axis1_dim = axis1_dim;
    this->axis2_dim = axis2_dim;
    matrix = new (std::nothrow) int *[axis1_dim];
    next_matrix = new (std::nothrow) int *[axis1_dim];

    if (matrix == nullptr || next_matrix == nullptr)
    {
        return CAEnums::CellsMalloc;
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
            matrix[j][k] = fill_value;
            next_matrix[j][k] = fill_value;
        }
    }

    return 0;
}

template <>
int CellularAutomata<int>::setup_dimensions_3d(int axis1_dim, int axis2_dim, int axis3_dim, int fill_value)
{
    if (tensor != nullptr)
    {
        return CAEnums::CellsAlreadyInitialized;
    }

    this->axis1_dim = axis1_dim;
    this->axis2_dim = axis2_dim;
    this->axis3_dim = axis3_dim;
    tensor = new (std::nothrow) int **[axis1_dim];
    next_tensor = new (std::nothrow) int **[axis1_dim];

    if (tensor == nullptr || next_tensor == nullptr)
    {
        return CAEnums::CellsMalloc;
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
                tensor[i][j][k] = fill_value;
                next_tensor[i][j][k] = fill_value;
            }
        }
    }

    return 0;
}

template <>
int CellularAutomata<int>::init_condition(int x_state, double prob)
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
        return CAEnums::CellsAreNull;
    }
    return 0;
}

template <>
int CellularAutomata<int>::set_new_cell_state(int *cell_index, int index_size,
                                              int *neighborhood_cells, int neighborhood_size,
                                              int &new_cell_state, void(custom_rule)(int *, int, int *, int, int &))
{
    int sum = 0;                         // sum of cells within boundary_radius for Parity rule
    MajorityCounter state_votes_counter; // counter to keep track of votes for Majority rule
    initialize_majority_rule_counter(state_votes_counter, num_states);

    switch (rule_type)
    {
    case CAEnums::Custom:
        if (vector == nullptr && matrix == nullptr && tensor == nullptr)
        {
            new_cell_state = -1;
            return CAEnums::CellsAreNull;
        }
        if (custom_rule == nullptr)
        {
            new_cell_state = -1;
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
            sum += neighborhood_cells[i];
        }
        new_cell_state = sum % num_states; // store the parity state as the new state
        break;
    case CAEnums::Majority:
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
                                    less_than_votes);
        new_cell_state = max_elem->first; // set the majority state as the new state
        break;
    }
    return 0;
}

template <>
int CellularAutomata<int>::print_grid()
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
        return CAEnums::CellsAreNull;
    }
    return 0;
}

template <>
int CellularAutomata<int>::step(void(custom_rule)(int *, int, int *, int, int &))
{
    int error_code = 0;       // store error code return by other methods
    int new_cell_state;       // stores the cell's new state
    int empty_cell_state = 0; // cell state for zeroing out old states
    int index_size;           // number of indices required to address the cell

    if (vector != nullptr)
    {
        index_size = 1;
#ifdef ENABLE_OMP
#pragma omp parallel for firstprivate(error_code) private(new_cell_state)
#endif
        for (int i = 0; i < axis1_dim; i++)
        {
            // store the main cell's index in cell_index for custom rule type
            int cell_index[1] = {i};
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
        swap_states<int>(vector, next_vector, axis1_dim);
    }
    else if (matrix != nullptr)
    {
        index_size = 2;
#ifdef ENABLE_OMP
#pragma omp parallel for firstprivate(error_code) private(new_cell_state)
#endif
        for (int i = 0; i < axis1_dim; i++)
        {
            for (int j = 0; j < axis2_dim; j++)
            {
                // store the main cell's index in cell_index for custom rule type
                int cell_index[2] = {i, j};
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
        swap_states<int>(matrix, next_matrix, axis1_dim, axis2_dim);
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
                    // store the main cell's index in cell_index for custom rule type
                    int cell_index[3] = {i, j, k};
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
        swap_states<int>(tensor, next_tensor, axis1_dim, axis2_dim, axis3_dim);
    }
    else
    {
        return CAEnums::CellsAreNull;
    }

    steps_taken++;
    return error_code;
}
