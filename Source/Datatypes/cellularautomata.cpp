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
#include <random>
#include <ctime>

/**
 * @brief Construct a new Cellular Automata:: Cellular Automata object.
 * Sets the default value to all class attributes
 *
 */
CellularAutomata::CellularAutomata()
{
    this->axis1_dim = 10;
    this->axis2_dim = 10;
    this->axis3_dim = 10;
    num_states = 2;
    boundary boundary_type = None;
    boundary_radius = 1;
    neighborhood neighborhood_type = VonNeumann;
    rule rule_type = Majority;
    this->shortr_weight = 1;
    this->longr_weight = 2;
    steps_taken = 0;
    vector = nullptr;
    next_vector == nullptr;
    matrix = nullptr;
    next_matrix == nullptr;
    tensor = nullptr;
    next_tensor == nullptr;
}

/**
 * @brief Destroy the Cellular Automata:: Cellular Automata object.
 * Deallocates memory reserved for vector/matrix/tensor
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
 * @brief Construct a new CellularAutomata::setup_dimensions object
 *
 * @param axis1_dim the size of a one dimensional vector
 * @returns the error code (-1: vector is already allocated) (0: no errors)
 */
int CellularAutomata::setup_dimensions(int axis1_dim)
{
    if (vector != nullptr)
    {
        return CellsAlreadyInitialized;
    }

    this->axis1_dim = axis1_dim;
    vector = new (nothrow) int[axis1_dim];
    next_vector = new (nothrow) int[axis1_dim];

    // initialize vector filled with zeros
    for (int j = 0; j < axis1_dim; j++)
    {
        vector[j] = 0;
        next_vector[j] = 0;
    }

    return 0;
}

/**
 * @brief Construct a new CellularAutomata::setup_dimensions object
 *
 * @param axis1_dim  The size of the first dimension
 * @param axis2_dim  The size of the second dimension
 * @returns error code (-1: matrix is already allocated) (0: no error)
 */
int CellularAutomata::setup_dimensions(int axis1_dim, int axis2_dim)
{
    if (matrix != nullptr)
    {
        return CellsAlreadyInitialized;
    }

    this->axis1_dim = axis1_dim;
    this->axis2_dim = axis2_dim;
    matrix = new (nothrow) int *[axis1_dim];
    next_matrix = new (nothrow) int *[axis1_dim];
    for (int i = 0; i < axis1_dim; i++)
    {
        matrix[i] = new (nothrow) int[axis2_dim];
        next_matrix[i] = new (nothrow) int[axis2_dim];
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
 * @brief Construct a new CellularAutomata::setup_dimensions object
 *
 * @param axis1_dim The size of the first dimension
 * @param axis2_dim The size of the second dimension
 * @param axis3_dim The size of the third dimension
 * @returns error code (-1: tensor was already allocated) (0: no error)
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
    tensor = new (nothrow) int **[axis1_dim];
    next_tensor = new (nothrow) int **[axis1_dim];
    for (int i = 0; i < axis1_dim; i++)
    {
        tensor[i] = new (nothrow) int *[axis2_dim];
        next_tensor[i] = new (nothrow) int *[axis2_dim];

        for (int j = 0; j < axis2_dim; j++)
        {
            tensor[i][j] = new (nothrow) int[axis3_dim];
            next_tensor[i][j] = new (nothrow) int[axis3_dim];
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
 * @brief setup neighborhood with values from enum neighborhood
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
 * @brief setup boundary with enum values from boundary, include boundary radius
 *
 * @param bound_type enum value for boundary (None, Periodic, Walled, CutOff)
 * @param radius radius for the boundary
 * @return int error code
 */
int CellularAutomata::setup_boundary(boundary bound_type, int radius)
{
    this->boundary_type = bound_type;
    this->boundary_radius = radius;
    return 0;
}

/**
 * @brief describes the range of cell states to be used in the CA object
 *
 * @param n_states describing the range of numbers to use for cell states
 * @return int error code
 */
int CellularAutomata::setup_cell_states(int n_states)
{
    this->num_states = n_states;
}

/**
 * @brief initializes the first state of the grid using random numbers
 *
 * @param x_state choose the cell state to initialize the grid with.
 * @param prob the probability of a cell to turn to state given from x_state
 * @return int error code (-2: tensor not initialized) (0: no error)
 */
int CellularAutomata::init_condition(int x_state, double prob)
{
    // initialize values in matrix
    // initialize matrix filled with random cell values
    srand(time(NULL));
    double random_cell_value;

    if (vector != nullptr)
    {
        for (int i = 0; i < axis1_dim; i++)
        {
            random_cell_value = (double)rand() / RAND_MAX;
            if (random_cell_value < prob)
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
                random_cell_value = (double)rand() / RAND_MAX;
                if (random_cell_value < prob)
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
                    random_cell_value = (double)rand() / RAND_MAX;
                    if (random_cell_value < prob)
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
 * @brief set the short and long weights used in the majority rule
 *
 * @param shortr_weight the short weight
 * @param longr_weight the long weight
 * @return int error code
 */
int CellularAutomata::setup_rule_short_long(double shortr_weight, double longr_weight)
{
    this->shortr_weight = shortr_weight;
    this->longr_weight = longr_weight;
    return 0;
}

/**
 * @brief setup the rule choice
 * used to specify the rule type to be used in CA object
 * @param rule_type enum rule representing the rule type. choose from (Majority, Parity)
 * @return int error code
 */
int CellularAutomata::setup_rule(rule rule_type)
{
    this->rule_type = rule_type;
    return 0;
}

int CellularAutomata::get_state_from_neighborhood(int axis1_dim, int &new_cell_state)
{
    if (neighborhood_type == VonNeumann)
    {
        switch (boundary_type)
        {
        case None:
            break;
        case Periodic:
            break;
        case Walled:
            break;
        case CutOff:
            break;
        }
    }
    else if (neighborhood_type == Moore)
    {
    }
    return 0;
}

int CellularAutomata::get_state_from_neighborhood(int axis1_dim, int axis2_dim, int &new_cell_state)
{
    if (neighborhood_type == VonNeumann)
    {
        switch (boundary_type)
        {
        case None:
            break;
        case Periodic:
            break;
        case Walled:
            break;
        case CutOff:
            break;
        }
    }
    else if (neighborhood_type == Moore)
    {
    }
    return 0;
}

int CellularAutomata::get_state_from_neighborhood(int axis1_dim, int axis2_dim, int axis3_dim, int &new_cell_state)
{
    if (neighborhood_type == VonNeumann)
    {
        switch (boundary_type)
        {
        case None:
            break;
        case Periodic:
            break;
        case Walled:
            break;
        case CutOff:
            break;
        }
    }
    else if (neighborhood_type == Moore)
    {
    }
    return 0;
}

int CellularAutomata::step()
{
    int new_cell_state;
    int error = 0;
    if (vector != nullptr)
    {
        for (int i = 0; i < axis1_dim; i++)
        {
            error = get_state_from_neighborhood(i, new_cell_state);
            if (error < 0)
            {
                return error;
            }
            next_vector[i] = new_cell_state;
        }
        swap_states(vector, next_vector, axis1_dim);
    }
    else if (matrix != nullptr)
    {
        for (int j = 0; j < axis1_dim; j++)
        {
            for (int k = 0; k < axis2_dim; k++)
            {
                error = get_state_from_neighborhood(j, k, new_cell_state);
                if (error < 0)
                {
                    return error;
                }
                next_matrix[j][k] = new_cell_state;
            }
        }
        swap_states(matrix, next_matrix, axis1_dim, axis2_dim);
    }
    else if (tensor != nullptr)
    {
        for (int i = 0; i < axis1_dim; i++)
        {
            for (int j = 0; j < axis2_dim; j++)
            {
                for (int k = 0; k < axis3_dim; k++)
                {
                    error = get_state_from_neighborhood(i, j, k, new_cell_state);
                    if (error < 0)
                    {
                        return error;
                    }
                    next_tensor[i][j][k] = new_cell_state;
                }
            }
        }
        swap_states(tensor, next_tensor, axis1_dim, axis2_dim, axis3_dim);
    }
    else
    {
        return CellsAreNull;
    }

    steps_taken++;
    return 0;
}

/**
 * @brief print the current state of the grid
 *
 * @return int error code (-2: tensor not initialized) (0: no error)
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
