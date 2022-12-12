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

int BaseCellularAutomata::setup_boundary(CAEnums::Boundary bound_type, int radius)
{
    if (radius <= 0)
    {
        return CAEnums::InvalidRadius;
    }
    this->boundary_type = bound_type;
    this->boundary_radius = radius;
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
    std::cout << "ERROR [";
    switch (error)
    {
    case CAEnums::CellsAlreadyInitialized:
        std::cout << error << "]: Can't reinitialize vector, matrix, nor tensor. \n";
        break;
    case CAEnums::CellsAreNull:
        std::cout << error << "]: The vector, matrix, and tensor are null. \n";
        break;
    case CAEnums::CellsMalloc:
        std::cout << error << "]: Could not allocate memory for either vector, matrix, or tensor. \n";
        break;
    case CAEnums::InvalidCellState:
        std::cout << error << "]: Invalid cell state given. Must be greater than equal to 2. \n";
        break;
    case CAEnums::InvalidCellStateCondition:
        std::cout << error << "]: Invalid cell state condition given. Must be less than the set number of states. \n";
        break;
    case CAEnums::InvalidRadius:
        std::cout << error << "]: Invalid boundary radius given. \n";
        break;
    case CAEnums::InvalidNumStates:
        std::cout << error << "]: Invalid number of states given. \n";
        break;
    case CAEnums::NeighborhoodCellsMalloc:
        std::cout << error << "]: Could not allocate memory for neighborhood array. \n";
        break;
    case CAEnums::CustomRuleIsNull:
        std::cout << error << "]: Custom rule function is null (none given). \n";
        break;
    }
}

template <typename T>
CellularAutomata<T>::CellularAutomata() : BaseCellularAutomata()
{
    steps_taken = 0;
    vector = nullptr;
    next_vector = nullptr;
    matrix = nullptr;
    next_matrix = nullptr;
    tensor = nullptr;
    next_tensor = nullptr;
}

template <typename T>
CellularAutomata<T>::~CellularAutomata()
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

template <typename T>
int CellularAutomata<T>::setup_dimensions_1d(int axis1_dim, int fill_value)
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

template <typename T>
int CellularAutomata<T>::setup_dimensions_2d(int axis1_dim, int axis2_dim, int fill_value)
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

template <typename T>
int CellularAutomata<T>::setup_dimensions_3d(int axis1_dim, int axis2_dim, int axis3_dim, int fill_value)
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

template <typename T>
int CellularAutomata<T>::init_condition(int x_state, double prob)
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

template <typename T>
int CellularAutomata<T>::set_new_cell_state(int *cell_index, int index_size,
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
            sum += neighborhood_cells[i];
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

template <typename T>
int CellularAutomata<T>::get_state_from_neighborhood_1d(int *cell_index, int index_size, T &new_cell_state,
                                                        void(custom_rule)(int *, int, T *, int, T &))
{
    int error_code = 0;         // store error code return by other methods
    int i = cell_index[0];      // get i-th index from array
    int neighborhood_index = 0; // keep track of neighborhood array as we iterate through CA cells

    // allocate memory and check if operation was successful
    T *neighborhood_cells = malloc_neighborhood_array(index_size); // index_size == rank
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

template <typename T>
int CellularAutomata<T>::get_state_from_neighborhood_2d(int *cell_index, int index_size, T &new_cell_state,
                                                        void(custom_rule)(int *, int, T *, int, T &))
{
    int error_code = 0;         // store error code return by other methods
    int i = cell_index[0];      // get i-th index from array
    int j = cell_index[1];      // get j-th index from array
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

template <typename T>
int CellularAutomata<T>::get_state_from_neighborhood_3d(int *cell_index, int index_size, T &new_cell_state,
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

template <typename T>
int CellularAutomata<T>::step(void(custom_rule)(int *, int, T *, int, T &))
{
    int error_code = 0;         // store error code return by other methods
    T new_cell_state;           // stores the cell's new state
    T empty_cell_state;         // cell state for zeroing out old states
    empty_cell_state.state = 0; // ensure the state property is set to zero
    int index_size;             // number of indices required to address the cell

    if (vector != nullptr)
    {
        // store the main cell's index in cell_index for custom rule type
        index_size = 1;
        int cell_index[index_size];
#pragma omp parallel for firstprivate(error_code) private(new_cell_state, cell_index)
        for (int i = 0; i < axis1_dim; i++)
        {
            cell_index[0] = i; // store the i-th index
            error_code = get_state_from_neighborhood_1d(cell_index, index_size, new_cell_state, custom_rule);
            if (error_code < 0)
            {
#pragma omp cancel for
#ifndef ENABLE_OMP
                /*
                 * Return error code when omp is not enabled
                 * otherwise the above pragma will exit out of the for loop(s)
                 */
                return error_code;
#endif
            }
            /*
             * The two following statements allow for dynamic systems to be modeled.
             * If the cell does not move then we properly update it in the second statement.
             * If the cell moves then the old cell index is zeroed
             * and the new cell index will contain the newly computed cell state.
             */
            next_vector[i] = empty_cell_state;
            next_vector[cell_index[0]] = new_cell_state;
        }
        // store next cell state to the current cell state for the next time step
        swap_states(vector, next_vector, axis1_dim);
    }
    else if (matrix != nullptr)
    {
        // store the main cell's index in cell_index for custom rule type
        index_size = 2;
        int cell_index[index_size];
#pragma omp parallel for firstprivate(error_code) private(new_cell_state, cell_index)
        for (int i = 0; i < axis1_dim; i++)
        {
            cell_index[0] = i; // store the i-th index
            for (int j = 0; j < axis2_dim; j++)
            {
                cell_index[1] = j; // store the j-th index
                error_code = get_state_from_neighborhood_2d(cell_index, index_size, new_cell_state, custom_rule);
                if (error_code < 0)
                {
#pragma omp cancel for
#ifndef ENABLE_OMP
                    /*
                     * Return error code when omp is not enabled
                     * otherwise the above pragma will exit out of the for loop(s)
                     */
                    return error_code;
#endif
                }
                /*
                 * The two following statements allow for dynamic systems to be modeled.
                 * If the cell does not move then we properly update it in the second statement.
                 * If the cell moves then the old cell index is zeroed
                 * and the new cell index will contain the newly computed cell state.
                 */
                next_matrix[i][j] = empty_cell_state;
                next_matrix[cell_index[0]][cell_index[1]] = new_cell_state;
            }
        }
        // store next cell state to the current cell state for the next time step
        swap_states(matrix, next_matrix, axis1_dim, axis2_dim);
    }
    else if (tensor != nullptr)
    {
        // store the main cell's index in cell_index for custom rule type
        index_size = 3;
        int cell_index[index_size];
#pragma omp parallel for firstprivate(error_code) private(new_cell_state, cell_index)
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
#pragma omp cancel for
#ifndef ENABLE_OMP
                        /*
                         * Return error code when omp is not enabled
                         * otherwise the above pragma will exit out of the for loop(s)
                         */
                        return error_code;
#endif
                    }
                    /*
                     * The two following statements allow for dynamic systems to be modeled.
                     * If the cell does not move then we properly update it in the second statement.
                     * If the cell moves then the old cell index is zeroed
                     * and the new cell index will contain the newly computed cell state.
                     */
                    next_tensor[i][j][k] = empty_cell_state;
                    next_tensor[cell_index[0]][cell_index[1]][cell_index[2]] = new_cell_state;
                }
            }
        }
        // store next cell state to the current cell state for the next time step
        swap_states(tensor, next_tensor, axis1_dim, axis2_dim, axis3_dim);
    }
    else
    {
        return CAEnums::CellsAreNull;
    }

    steps_taken++;
    return error_code;
}

template <typename T>
int CellularAutomata<T>::step()
{
    return step(nullptr); // return step(func) error code
}

template <typename T>
int CellularAutomata<T>::print_grid()
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

template <typename T>
void CellularAutomata<T>::generate_periodic_neighborhood_1d(int *cell_index, T *neighborhood_cells, int &neighborhood_index)
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

template <typename T>
void CellularAutomata<T>::generate_cutoff_neighborhood_1d(int *cell_index, T *neighborhood_cells, int &neighborhood_index)
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

template <typename T>
void CellularAutomata<T>::generate_periodic_neighborhood_2d(int *cell_index, T *neighborhood_cells, int &neighborhood_index)
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

template <typename T>
void CellularAutomata<T>::generate_cutoff_neighborhood_2d(int *cell_index, T *neighborhood_cells, int &neighborhood_index)
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

template <typename T>
void CellularAutomata<T>::generate_periodic_neighborhood_3d(int *cell_index, T *neighborhood_cells, int &neighborhood_index)
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

template <typename T>
void CellularAutomata<T>::generate_cutoff_neighborhood_3d(int *cell_index, T *neighborhood_cells, int &neighborhood_index)
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

template <typename T>
T *CellularAutomata<T>::malloc_neighborhood_array(int rank)
{
    /*
     * Generate a flatten array of the cell's neighborhood.
     * The neighborhood array can then be utilized for Majority, Parity, or Custom rule
     */
    int max_neighborhood_size = get_neighborhood_size(rank, boundary_radius, neighborhood_type);
    // allocate memory and check if operation was successful
    int *neighborhood_cells = new (std::nothrow) T[max_neighborhood_size];
    return neighborhood_cells;
}

CellularAutomata<int>::CellularAutomata() : BaseCellularAutomata()
{
    steps_taken = 0;
    vector = nullptr;
    next_vector = nullptr;
    matrix = nullptr;
    next_matrix = nullptr;
    tensor = nullptr;
    next_tensor = nullptr;
}

CellularAutomata<int>::~CellularAutomata()
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

int CellularAutomata<int>::get_state_from_neighborhood_1d(int *cell_index, int index_size, int &new_cell_state,
                                                          void(custom_rule)(int *, int, int *, int, int &))
{
    int error_code = 0;         // store error code return by other methods
    int i = cell_index[0];      // get i-th index from array
    int neighborhood_index = 0; // keep track of neighborhood array as we iterate through CA cells

    // allocate memory and check if operation was successful
    int *neighborhood_cells = malloc_neighborhood_array(index_size);
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

int CellularAutomata<int>::get_state_from_neighborhood_2d(int *cell_index, int index_size, int &new_cell_state,
                                                          void(custom_rule)(int *, int, int *, int, int &))
{
    int error_code = 0;         // store error code return by other methods
    int i = cell_index[0];      // get i-th index from array
    int j = cell_index[1];      // get j-th index from array
    int neighborhood_index = 0; // keep track of neighborhood array as we iterate through CA cells

    // allocate memory and check if operation was successful
    int *neighborhood_cells = malloc_neighborhood_array(index_size);
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

int CellularAutomata<int>::get_state_from_neighborhood_3d(int *cell_index, int index_size, int &new_cell_state,
                                                          void(custom_rule)(int *, int, int *, int, int &))
{
    int error_code = 0;         // store error code return by other methods
    int i = cell_index[0];      // get i-th index from array
    int j = cell_index[1];      // get j-th index from array
    int k = cell_index[2];      // get k-th index from array
    int neighborhood_index = 0; // keep track of neighborhood array as we iterate through CA cells

    // allocate memory and check if operation was successful
    int *neighborhood_cells = malloc_neighborhood_array(index_size);
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

int CellularAutomata<int>::step(void(custom_rule)(int *, int, int *, int, int &))
{
    int error_code = 0; // store error code return by other methods
    int new_cell_state; // stores the cell's new state
    int index_size;     // number of indices required to address the cell

    if (vector != nullptr)
    {
        // store the main cell's index in cell_index for custom rule type
        index_size = 1;
        int cell_index[index_size];
#pragma omp parallel for firstprivate(error_code) private(new_cell_state, cell_index)
        for (int i = 0; i < axis1_dim; i++)
        {
            cell_index[0] = i; // store the i-th index
            error_code = get_state_from_neighborhood_1d(cell_index, index_size, new_cell_state, custom_rule);
            if (error_code < 0)
            {
#pragma omp cancel for
#ifndef ENABLE_OMP
                /*
                 * Return error code when omp is not enabled
                 * otherwise the above pragma will exit out of the for loop(s)
                 */
                return error_code;
#endif
            }
            /*
             * The two following statements allow for dynamic systems to be modeled.
             * If the cell does not move then we properly update it in the second statement.
             * If the cell moves then the old cell index is zeroed
             * and the new cell index will contain the newly computed cell state.
             */
            next_vector[i] = 0;
            next_vector[cell_index[0]] = new_cell_state;
        }
        // store next cell state to the current cell state for the next time step
        swap_states(vector, next_vector, axis1_dim);
    }
    else if (matrix != nullptr)
    {
        // store the main cell's index in cell_index for custom rule type
        index_size = 2;
        int cell_index[index_size];
#pragma omp parallel for firstprivate(error_code) private(new_cell_state, cell_index)
        for (int i = 0; i < axis1_dim; i++)
        {
            cell_index[0] = i; // store the i-th index
            for (int j = 0; j < axis2_dim; j++)
            {
                cell_index[1] = j; // store the j-th index
                error_code = get_state_from_neighborhood_2d(cell_index, index_size, new_cell_state, custom_rule);
                if (error_code < 0)
                {
#pragma omp cancel for
#ifndef ENABLE_OMP
                    /*
                     * Return error code when omp is not enabled
                     * otherwise the above pragma will exit out of the for loop(s)
                     */
                    return error_code;
#endif
                }
                /*
                 * The two following statements allow for dynamic systems to be modeled.
                 * If the cell does not move then we properly update it in the second statement.
                 * If the cell moves then the old cell index is zeroed
                 * and the new cell index will contain the newly computed cell state.
                 */
                next_matrix[i][j] = 0;
                next_matrix[cell_index[0]][cell_index[1]] = new_cell_state;
            }
        }
        // store next cell state to the current cell state for the next time step
        swap_states(matrix, next_matrix, axis1_dim, axis2_dim);
    }
    else if (tensor != nullptr)
    {
        // store the main cell's index in cell_index for custom rule type
        index_size = 3;
        int cell_index[index_size];
#pragma omp parallel for firstprivate(error_code) private(new_cell_state, cell_index)
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
#pragma omp cancel for
#ifndef ENABLE_OMP
                        /*
                         * Return error code when omp is not enabled
                         * otherwise the above pragma will exit out of the for loop(s)
                         */
                        return error_code;
#endif
                    }
                    /*
                     * The two following statements allow for dynamic systems to be modeled.
                     * If the cell does not move then we properly update it in the second statement.
                     * If the cell moves then the old cell index is zeroed
                     * and the new cell index will contain the newly computed cell state.
                     */
                    next_tensor[i][j][k] = 0;
                    next_tensor[cell_index[0]][cell_index[1]][cell_index[2]] = new_cell_state;
                }
            }
        }
        // store next cell state to the current cell state for the next time step
        swap_states(tensor, next_tensor, axis1_dim, axis2_dim, axis3_dim);
    }
    else
    {
        return CAEnums::CellsAreNull;
    }

    steps_taken++;
    return error_code;
}

int CellularAutomata<int>::step()
{
    return step(nullptr); // return step(func) error code
}

void CellularAutomata<int>::generate_periodic_neighborhood_1d(int *cell_index, int *neighborhood_cells, int &neighborhood_index)
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

void CellularAutomata<int>::generate_cutoff_neighborhood_1d(int *cell_index, int *neighborhood_cells, int &neighborhood_index)
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

void CellularAutomata<int>::generate_periodic_neighborhood_2d(int *cell_index, int *neighborhood_cells, int &neighborhood_index)
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

void CellularAutomata<int>::generate_cutoff_neighborhood_2d(int *cell_index, int *neighborhood_cells, int &neighborhood_index)
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

void CellularAutomata<int>::generate_periodic_neighborhood_3d(int *cell_index, int *neighborhood_cells, int &neighborhood_index)
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

void CellularAutomata<int>::generate_cutoff_neighborhood_3d(int *cell_index, int *neighborhood_cells, int &neighborhood_index)
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

int *CellularAutomata<int>::malloc_neighborhood_array(int rank)
{
    /*
     * Generate a flatten array of the cell's neighborhood.
     * The neighborhood array can then be utilized for Majority, Parity, or Custom rule
     */
    int max_neighborhood_size = get_neighborhood_size(rank, boundary_radius, neighborhood_type);
    // allocate memory and check if operation was successful
    int *neighborhood_cells = new (std::nothrow) int[max_neighborhood_size];
    return neighborhood_cells;
}
