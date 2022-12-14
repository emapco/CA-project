/**
 * @file CAutils.h
 * @author Emmanuel Cortes (ecortes@berkeley.edu)
 *
 * <b>Contributor(s)</b> <br> &emsp;&emsp;
 * @brief This header files contains utility functions used by CellularAutomata class.
 * @date 2022-12-06
 */
#pragma once
#include <utility> // pair
#ifdef ENABLE_OMP
#include <omp.h>
#endif

/**
 * @brief swaps the computed next_vector to the current vector
 *
 * @param vector cellular automata current vector state
 * @param next_vector cellular automata next vector state
 * @param axis1_dim vector dimension
 */
template <typename T>
void swap_states(T *vector, T *next_vector, int axis1_dim)
{
#pragma omp parallel for
    for (int i = 0; i < axis1_dim; i++)
    {
        std::swap(vector[i], next_vector[i]);
    }
}

/**
 * @brief swaps the computed next_matrix to the current matrix
 *
 * @param matrix cellular automata current matrix state
 * @param next_matrix cellular automata next matrix state
 * @param axis1_dim first matrix dimension
 * @param axis2_dim second matrix dimension
 */
template <typename T>
void swap_states(T **matrix, T **next_matrix, int axis1_dim, int axis2_dim)
{
#pragma omp parallel for collapse(2)
    for (int i = 0; i < axis1_dim; i++)
    {
        for (int j = 0; j < axis2_dim; j++)
        {
            std::swap(matrix[i][j], next_matrix[i][j]);
        }
    }
}

/**
 * @brief swaps the computed next_tensor to the current tensor
 *
 * @param tensor cellular automata current tensor state
 * @param next_tensor cellular automata next tensor state
 * @param axis1_dim first tensor dimension
 * @param axis2_dim second tensor dimension
 * @param axis3_dim third tensor dimension
 */
template <typename T>
void swap_states(T ***tensor, T ***next_tensor, int axis1_dim, int axis2_dim, int axis3_dim)
{
#pragma omp parallel for collapse(3)
    for (int i = 0; i < axis1_dim; i++)
    {
        for (int j = 0; j < axis2_dim; j++)
        {
            for (int k = 0; k < axis3_dim; k++)
            {
                std::swap(tensor[i][j][k], next_tensor[i][j][k]);
            }
        }
    }
}

/**
 * @brief determines if a cell is diagonal to the central cell\n
 * diagram of slice: 1 = diagonal cell\n
 * 1 0 1\n
 * 0 0 0\n
 * 1 0 1\n
 *
 * @param i cell's i-th index
 * @param j cell's j-th index
 * @return true: cell is diagonal
 * @return false: cell is not diagonal
 */
bool is_diagonal_neighboring_cell_2d(int i, int j);

/**
 * @brief determines if a cell is diagonal to the central cell\n
 * if k == 0 then both i and j have be non-zero\n
 * diagram of slice: 1 = diagonal cell\n
 * 1 0 1\n
 * 0 0 0\n
 * 1 0 1\n
 *
 * if k != 0 then either i or j can be non-zero\n
 * diagram of slice: 1 = diagonal cell\n
 * 1 1 1\n
 * 1 0 1\n
 * 1 1 1\n
 *
 * @param i cell's i-th index
 * @param j cell's j-th index
 * @param k cell's k-th index
 * @return true: cell is diagonal
 * @return false: cell is not diagonal
 */
bool is_diagonal_neighboring_cell_3d(int i, int j, int k);

/**
 * @brief Determines if a has less votes than b.
 * Used to determine max element in a MajorityCounter instance.
 *
 * @param a pair instance with the votes stored in the second property
 * @param b other pair with the votes stored in the second property
 */
bool less_than_votes(const std::pair<int, int> &a, const std::pair<int, int> &b);

/**
 * @brief get the periodic index used for finding periodic boundary neighbors
 *
 * @param i cell's i-th index
 * @param di neighbor's cell i-th offset
 * @param axis_dim dimension of i-th axis
 * @return int
 */
int get_periodic_index(int i, int di, int axis_dim);

/**
 * @brief Get the periodic Moore neighbor index [x, y, z] from a flattened neighborhood array index.
 *
 * @param rank cell data rank
 * @param radius neighborhood radius
 * @param neighborhood_array_index flat array index to be converted to [x, y, z]
 * @param neighbor_index int array containing x, y and z coordinates
 */
void get_periodic_moore_neighbor_index(int rank, int radius, int neighborhood_array_index, int *neighbor_index);

/**
 * @brief Get the periodic Von  Neumann neighbor index [x, y, z] from a flattened neighborhood array index.
 *
 * @param rank cell data rank
 * @param radius neighborhood radius
 * @param neighborhood_array_index flat array index to be converted to [x, y, z]
 * @param neighbor_index int array containing x, y and z coordinates
 */
void get_periodic_von_neumann_neighbor_index(int rank, int radius, int neighborhood_array_index, int *neighbor_index);
