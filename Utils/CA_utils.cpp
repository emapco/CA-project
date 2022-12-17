/**
 * @file CA_utils.cpp
 * @author Emmanuel Cortes (ecortes@berkeley.edu)
 *
 * <b>Contributor(s)</b> <br> &emsp;&emsp;  Chongye Feng
 * @brief Implementation file for the utility functions utilized
 * by CellularAutomata class
 * defined in CAutils.h
 * @date 2022-12-06
 */
#include "CAdatatypes.h"
#include <utility> // swap, pair
#include <cmath>   // pow
#include <string>
#include <map>
#include <sstream>
#include <fstream>

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

void get_density(std::ifstream &data, std::ofstream &result)
{
    // Starts with dim_3 -> dim_2 -> dim_1
    int n_states, dim[3], n_steps, i, datapoint;
    std::string line, word, rule;
    std::stringstream str;
    std::map<int, int> step;
    std::map<int, int>::iterator itr;

    if (data.is_open())
    {
        // First line: n_states
        std::getline(data, line);
        n_states = std::stoi(line);
        for (i = 0; i < n_states; i++)
        {
            step.insert(std::make_pair(i, 0));
        }
        result << n_states << "\n";

        // Second line: rules and functions
        std::getline(data, line);
        rule = line;

        // Third line: dims
        std::getline(data, line);
        str << line;
        i = 0; // Iterator
        while (std::getline(str, word, ','))
        {
            dim[i] = std::stoi(word);
        }
        for (i = 0; i < 3; i++)
        {
            result << dim[i] << ",";
        }
        result << "\n";

        // Lines after: log
        while (std::getline(data, line))
        {
            // Making the map zeroed
            for (i = 0; i < n_states; i++)
            {
                itr = step.find(i);
                if (itr != step.end())
                {
                    itr->second = 0;
                }
            }

            str << line;
            while (std::getline(str, word, ','))
            {
                datapoint = std::stoi(word);
                itr = step.find(datapoint);
                if (itr != step.end())
                {
                    itr->second += 1;
                }
            }

            for (itr = step.begin(); itr != step.end(); ++itr)
            {
                result << itr->second << ",";
            }
            result << "\n";
            n_steps += 1;
        }
    }
}
