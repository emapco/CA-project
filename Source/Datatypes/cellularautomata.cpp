/**
 * @file cellularautomata.cpp
 * @author Trevor Oldham (trevoldham@berkeley.edu)
 *
 * <b>Contributor(s)</b> <br> &emsp;&emsp; Emmanuel Cortes
 * @brief This file contains the class CellularAutomata functions declared in CAdatatypes.h
 * @date 2022-12-03
 */
#include "CAdatatypes.h"
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
    this->n = 10;
    this->m = 10;
    this->p = 10;
    num_states = 2;
    boundary boundary_type = None;
    boundary_radius = 1;
    neighborhood neighborhood_type = VonNeumann;
    rule rule_type = Majority;
    this->shortr_weight = 1;
    this->longr_weight = 2;
    vector = nullptr;
    matrix = nullptr;
    tensor = nullptr;
}

/**
 * @brief Destroy the Cellular Automata:: Cellular Automata object.
 * Deallocates memory reserved for vector/matrix/tensor
 *
 */
CellularAutomata::~CellularAutomata()
{
    // Free each sub-array
    if (ndims == 1)
    {
        // Free the array of pointers
        delete[] vector;
    }

    if (ndims == 2)
    {
        // Free each sub-array
        for (int i = 0; i < n; ++i)
        {
            delete[] matrix[i];
        }
        // Free the array of pointers
        delete[] matrix;
    }

    if (ndims == 3)
    {
        // Free each sub-array
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; j++)
            {
                delete[] tensor[i][j];
            }

            delete[] tensor[i];
        }

        // Free the array of pointers
        delete[] tensor;
    }
}

/**
 * @brief Construct a new CellularAutomata::setup_dimensions object
 *
 * @param n the size of a one dimensional vector
 * @returns the error code (-1: vector is already allocated) (0: no errors)
 */
int CellularAutomata::setup_dimensions(int n)
{
    if (vector != nullptr)
    {
        return CellsAlreadyInitialized;
    }

    this->ndims = 1;
    this->n = n;
    vector = new (nothrow) int[n];

    // initialize vector filled with zeros
    for (int j = 0; j < n; j++)
    {
        vector[j] = 0;
    }

    return 0;
}

/**
 * @brief Construct a new CellularAutomata::setup_dimensions object
 *
 * @param n  The size of the first dimension
 * @param m  The size of the second dimension
 * @returns error code (-1: matrix is already allocated) (0: no error)
 */
int CellularAutomata::setup_dimensions(int n, int m)
{
    if (matrix != nullptr)
    {
        return CellsAlreadyInitialized;
    }
    this->ndims = 2;
    this->n = n;
    this->m = m;
    matrix = new (nothrow) int *[n];
    for (int i = 0; i < n; i++)
    {
        matrix[i] = new (nothrow) int[m];
    }

    // initialize matrix filled with zeros
    for (int j = 0; j < n; j++)
    {
        for (int k = 0; k < m; k++)
        {
            matrix[j][k] = 0;
        }
    }

    return 0;
}

/**
 * @brief Construct a new CellularAutomata::setup_dimensions object
 *
 * @param n The size of the first dimension
 * @param m The size of the second dimension
 * @param p The size of the third dimension
 * @returns error code (-1: tensor was already allocated) (0: no error)
 */
int CellularAutomata::setup_dimensions(int n, int m, int p)
{
    if (tensor != nullptr)
    {
        return CellsAlreadyInitialized;
    }

    this->ndims = 3;
    this->n = n;
    this->m = m;
    this->p = p;
    tensor = new (nothrow) int **[n];
    for (int i = 0; i < n; i++)
    {
        tensor[i] = new (nothrow) int *[m];

        for (int j = 0; j < m; j++)
        {

            tensor[i][j] = new (nothrow) int[p];
        }
    }

    // initialize matrix filled with zeros
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            for (int k = 0; k < p; k++)
            {
                tensor[i][j][k] = 0;
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
        for (int i = 0; i < n; i++)
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
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < m; k++)
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
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                for (int k = 0; k < p; k++)
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

/**
 * @brief print the current state of the grid
 *
 * @return int error code (-2: tensor not initialized) (0: no error)
 */
int CellularAutomata::print_grid()
{
    if (vector != nullptr)
    {
        for (int i = 0; i < n; i++)
        {
            std::cout << vector[i] << " ";
        }
        std::cout << std::endl;
    }
    else if (matrix != nullptr)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < m; k++)
            {
                std::cout << matrix[j][k] << " ";
            }
            std::cout << std::endl;
        }
    }
    else if (tensor != nullptr)
    {
        for (int i = 0; i < n; i++)
        {
            std::cout << "Printing " << i << "'th slice of Tensor" << std::endl;
            for (int j = 0; j < m; j++)
            {
                for (int k = 0; k < p; k++)
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
