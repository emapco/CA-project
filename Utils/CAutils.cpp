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

int get_neighborhood_size(int rank, int radius, CellularAutomata::neighborhood neighborhood_type)
{
    if (neighborhood_type == CellularAutomata::VonNeumann)
    {
        return (2 * rank * radius) + 1; // +1 to include cell of interest
    }
    else // CellularAutomata::Moore neighborhood
    {
        return pow((2 * radius + 1), rank);
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
    if (k == 0)
    {
        // diagram of slice: 1 = diagonal cell
        // 1 0 1
        // 0 0 0
        // 1 0 1
        return (i != 0 && j != 0);
    }
    else // k != 0
    {
        // diagram of slice: 1 = diagonal cell
        // 1 1 1
        // 1 0 1
        // 1 1 1
        return (i != 0 || j != 0);
    }
}

void initialize_majority_rule_counter(MajorityCounter &counter, int num_states)
{
    // sets the counter for every cell state type to 0
    for (int j = 0; j < num_states; j++)
    {
        counter.insert(std::make_pair(j, 0));
    }
}

bool less_than_votes(const std::pair<int, int> &a, const std::pair<int, int> &b)
{
    return a.second < b.second;
}

int get_periodic_index_axis(int i, int di, int axis_dim)
{
    return (i + di + axis_dim) % axis_dim;
}

/**
 * @brief Get the density object
 * 
 * @param data 
 * @param result Line1 n_states; Line2 dims; Line3 counts of states of each step.
 */
void get_density(std::ifstream& data, std::ofstream& result)
{
    // Starts with dim_3 -> dim_2 -> dim_1
    int n_states, dim[3], n_steps, i, datapoint;
    std::string line, word, rule;
    std::stringstream str;
    std::map<int,int> step;
    std::map<int,int>::iterator itr;

    if (data.isopen())
    {
        // First line: n_states
        getline(data,line);
        n_states = std::stoi(line);
        for (i=0;i<n_states;i++)
        {
            step.insert(i,0);
        }
        result << n_states << "\n";

        // Second line: rules and functions
        getline(data,line);
        rule = line;

        // Third line: dims
        getline(data,line);
        str(line);
        i = 0; // Iterator
        while(getline(str,word,','))
        {
            dim[i] = std::stoi(word);
        }
        for (i=0;i<3;i++)
        {
            result << dim[i] << ",";
        }
        result << "\n";

        // Lines after: log
        while(getline(data,line))
        {
            // Making the map zeroed
            for (i=0;i<n_states;i++)
            {
                itr = step.find(i);
                if (itr != step.end())
                {
                    itr->second = 0;
                }
            }

            str(line);
            while(getline(str,word,','))
            {
                datapoint = std::stoi(word);
                itr = step.find(datapoint);
                if (itr != step.end())
                {
                    itr->second += 1;
                }
            }
            
            for (itr = step.begin();itr!=step.end();++itr)
            {
                result << itr->second << ",";
            }
            result << "\n";
            n_steps += 1;
        }
    }
}