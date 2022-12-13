/**
 * @file cautils.cpp
 * @author Emmanuel Cortes (ecortes@berkeley.edu)
 *
 * <b>Contributor(s)</b> <br> &emsp;&emsp;
 * @brief Implementation file for the utility functions utilized
 * by CellularAutomata class
 * defined in CAutils.h
 * @date 2022-12-06
 */
#include "CAdatatypes.h"
#include <utility> // swap, pair
#include <cmath>   // pow
#ifdef ENABLE_OMP
#include <omp.h>
#endif

void swap_states(int *vector, int *next_vector, int axis1_dim)
{
#pragma omp parallel for
    for (int i = 0; i < axis1_dim; i++)
    {
        std::swap(vector[i], next_vector[i]);
    }
}

void swap_states(int **matrix, int **next_matrix, int axis1_dim, int axis2_dim)
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

void swap_states(int ***tensor, int ***next_tensor, int axis1_dim, int axis2_dim, int axis3_dim)
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

bool is_diagonal_neighboring_cell_2d(int i, int j)
{
    // diagram of slice: 1 = diagonal cell
    // 1 0 1
    // 0 0 0
    // 1 0 1
    return (i != 0 && j != 0);
}

bool is_diagonal_neighboring_cell_3d(int i, int j, int k)
{
    // assumes that the middle matrix has the index i = 0
    if (i == 0)
    {
        // diagram of slice: 1 = diagonal cell
        // 1 0 1
        // 0 0 0
        // 1 0 1
        return (k != 0 && j != 0);
    }
    else // assumes that non-zero matrices are have i != 0
    {    // (i.e. i = -1; i = 2)
        // diagram of slice: 1 = diagonal cell
        // 1 1 1
        // 1 0 1
        // 1 1 1
        return (k != 0 || j != 0);
    }
}

bool less_than_votes(const std::pair<int, int> &a, const std::pair<int, int> &b)
{
    return a.second < b.second;
}

int get_periodic_index(int i, int di, int axis_dim)
{
    return (i + di + axis_dim) % axis_dim;
}

void get_periodic_moore_neighbor_index(int rank, int radius, int neighborhood_array_index, int *neighbor_index)
{
    int factor = 2 * radius + 1;
    switch (rank)
    {
    case 1: // vector
        neighbor_index[0] = neighborhood_array_index - radius;
        break;
    case 2: // matrix
        /*
         * matrix representation
         * x\y -1 0 1
         *  -1: x x x
         *   0: x x x
         *   1: x x x
         *
         * flatten representation
         * [0, 1, 2, 3, 4, 5, 6, 7, 8]
         * [[-1, -1], [-1, 0], [-1, 1], [0, -1], [0, 0], [0, 1], [1, -1], [1, 0], [1, 1]]
         *
         * factor = 3
         *
         * ex: index = 2 = [-1, 1]
         * index_i = (2 div 3) - 1 = -1
         * index_j = (2 mod 3) - 1 = 1
         * -> [-1, 1]
         *
         * ex: index = 0 = [-1, 1]
         * index_i = (0 div 3) - 1 = -1
         * index_j = (0 mod 3) - 1 = -1
         * -> [-1, -1]
         *
         * ex: index 8 = [1, 1]
         * index_i = (8 div 3) - 1 = 1
         * index_j = (8 mod 3) - 1 = 1
         */
        neighbor_index[0] = (neighborhood_array_index / factor) - radius; // x index
        neighbor_index[1] = (neighborhood_array_index % factor) - radius; // y index
        break;
    case 3: // tensor
        /*
         * Similar as matrix but we have multiple matrices. We first get x with f2.
         * The first matrix (x = -1) has elements [0..8] and second matrix has [9..17] and so forth
         * y and z index as essentially the same as x and y index for the rank 2 (matrix) case.
         * We just have to adjust neighborhood_array_index (adjusted_index) to match the same logic as the rank 2 case
         */
        int f2 = factor * factor;
        int adjusted_index = neighborhood_array_index % f2;
        neighbor_index[0] = (neighborhood_array_index / f2) - radius; // x index
        neighbor_index[1] = (adjusted_index / factor) - radius;       // y index
        neighbor_index[2] = (adjusted_index % factor) - radius;       // z index
        break;
    }
}

void get_periodic_von_neumann_neighbor_index(int rank, int radius, int neighborhood_array_index, int *neighbor_index)
{
    int length = (2 * rank * radius) + 1;
    int middle_cell_index = length / 2;
    switch (rank)
    {
    case 1: // vector
        neighbor_index[0] = neighborhood_array_index - radius;
        break;
    case 2: // matrix
        /*
         * More complex than Moore neighbor implementation
         * due to the exclusion of diagonal neighbors.
         *
         * matrix representation
         * x\y -1 0 1
         *  -1: 0 x 0
         *   0: x x x
         *   1: 0 x 0
         *
         * flatten representation
         * [0, 1, 2, 3, 4, 5, 6, 7, 8]
         * [[-1, 0], [0, -1], [0, 0], [0, 1], [1, 0]]
         *
         * the if block handles north neighbors
         * the else if block handles south neighbors
         * the else handles the middle row
         *
         */
        if (neighborhood_array_index < radius)
        {
            neighbor_index[0] = -radius + neighborhood_array_index;
            neighbor_index[1] = 0;
        }
        else if (length - neighborhood_array_index <= radius)
        {
            neighbor_index[0] = neighborhood_array_index - length + 1 + radius;
            neighbor_index[1] = 0;
        }
        else
        {
            neighbor_index[0] = 0;
            neighbor_index[1] = neighborhood_array_index - middle_cell_index;
        }
        break;
    case 3: // tensor
        /*
         * More complex than Moore neighbor implementation
         * due to the exclusion of diagonal neighbors.
         *
         * the if block handles the case where: x = -radius, -radius+1, ..., -1
         *      these matrices only have one neighbor (-i, 0, 0)
         * the else if block handles the case where: x = 1, 2, ..., radius
         *      these matrices only have one neighbor (i, 0, 0)
         * the else handles the middle matrix
         *
         */
        if (neighborhood_array_index < radius)
        {
            neighbor_index[0] = -radius + neighborhood_array_index;
            neighbor_index[1] = 0;
            neighbor_index[2] = 0;
        }
        else if (length - neighborhood_array_index <= radius)
        {
            neighbor_index[0] = neighborhood_array_index - length + 1 + radius;
            neighbor_index[1] = 0;
            neighbor_index[2] = 0;
        }
        else
        {
            /*
             * Same logic as the matrix case but we need to account
             * for the matrices where x != 0
             *
             * the if block handles north neighbors
             * the else if block handles south neighbors
             * the else handles the middle row
             */
            int adjusted_index = neighborhood_array_index - radius;
            if (adjusted_index < radius)
            {
                neighbor_index[0] = 0;
                neighbor_index[1] = -radius + adjusted_index;
                neighbor_index[2] = 0;
            }
            else if (length - neighborhood_array_index - radius <= radius)
            {
                neighbor_index[0] = 0;
                neighbor_index[1] = neighborhood_array_index - length + 1 + 2 * radius;
                neighbor_index[2] = 0;
            }
            else
            {
                neighbor_index[0] = 0;
                neighbor_index[1] = 0;
                neighbor_index[2] = neighborhood_array_index - middle_cell_index;
            }
        }
        break;
    }
}
