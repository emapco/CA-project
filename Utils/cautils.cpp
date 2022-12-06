/**
 * @file cautils.cpp
 * @author Emmanuel Cortes (ecortes@berkeley.edu)
 *
 * <b>Contributor(s)</b> <br> &emsp;&emsp;
 * @brief Implementation file for the utility functions
 * defined in CAutils.h
 * @date 2022-12-06
 */

#include <utility>

void swap_states(int *vector, int *next_vector, int axis1_dim)
{
    for (int i = 0; i < axis1_dim; i++)
    {
        std::swap(vector[i], next_vector[i]);
    }
}

void swap_states(int **matrix, int **next_matrix, int axis1_dim, int axis2_dim)
{
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
