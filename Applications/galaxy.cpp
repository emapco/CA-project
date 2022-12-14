/**
 * @file galaxy.cpp
 * @author Emmanuel Cortes (ecortes@berkeley.edu)
 *
 * <b>Contributor(s)</b> <br> &emsp;&emsp; Trevor Oldham, Chongye Feng
 * @brief The application
 * @date 2022-12-12
 */

#include "CAdatatypes.h"
#include "CAutils.h"
#include <vector>
#include <algorithm> // copy
#include <cmath>

/**
 * @brief Struct that represents a star system.
 * The struct contains a velocity
 *
 */
struct GalaxyCell
{
    int state = 0; // 0: empty space, 1: star
    double velocity[3] = {0, 0, 0};
    double mass = 0;

    bool operator!=(const GalaxyCell &other)
    {
        if (state != other.state || mass != other.mass)
        {
            return true;
        }
        for (int i = 0; i < 3; i++)
        {
            if (velocity[i] != other.velocity[i])
            {
                return true;
            }
        }
        return false;
    }
};

std::vector<double> compute_gravitational_force(const GalaxyCell &cell_of_interest, const GalaxyCell &neighbor_cell,
                                                const std::vector<int> &neighbor_index);
std::vector<double> compute_accel(const std::vector<double> &total_force, double mass);
std::vector<double> compute_velocity(const std::vector<double> &accel, double time_step);
std::vector<double> compute_displacement(const std::vector<double> &velocity, const GalaxyCell &cell_of_interest, double time_step);
double compute_vector_norm(const std::vector<double> &vector);
std::vector<double> compute_vector_difference(const std::vector<double> &vec1, const std::vector<double> &vec2);
void update_velocity_after_collision(GalaxyCell &new_cell_state, GalaxyCell cell);

/**
 * @brief Rounds a double to an int.
 *
 * @param dbl
 * @return int
 */
int round_int(double dbl)
{
    return static_cast<int>(dbl < 0 ? dbl - 0.5 : dbl + 0.5);
}

/**
 * @brief Get the periodic vector.
 *
 * @param cell_index current new_cell_state position
 * @param offset_index used to compute the neighboring cell position
 * @param CA CellularAutomata instance for accesses model parameters
 * @return std::vector<int>
 */
std::vector<int> get_periodic_vector(int *cell_index, const std::vector<int> &offset_index, CellularAutomata<GalaxyCell> &CA)
{
    int x = get_periodic_index(cell_index[0], offset_index[0], CA.axis1_dim);
    int y = get_periodic_index(cell_index[1], offset_index[1], CA.axis2_dim);
    int z = get_periodic_index(cell_index[2], offset_index[2], CA.axis3_dim);
    return {x, y, z};
}

/**
 * @brief Get the cell state object at position cell_index + offset_index
 * from the CellularAutomata instance.
 *
 * @param cell_index current new_cell_state position
 * @param offset_index used to compute the neighboring cell position
 * @param CA CellularAutomata instance for accesses model parameters
 * @return const GalaxyCell&
 */
GalaxyCell &get_cell_state(int *cell_index, const std::vector<int> &offset_index, CellularAutomata<GalaxyCell> &CA)
{
    std::vector<int> vec = get_periodic_vector(cell_index, offset_index, CA);
    GalaxyCell ***next_tensor = CA.get_next_tensor();
    return next_tensor[vec[0]][vec[1]][vec[2]];
}

/**
 * @brief Determines if the new_cell_state collides with a cell at the position at cell_index + offset_index
 * If there is a collision then we set the new_cell_state index to the previous offset (cell_index + prev_offset).
 * Else we copy the current offset to previous offset.
 *
 * @param cell_index new_cell_state position
 * @param offset_index keep track of path the cell takes
 * @param prev_offset used to compute new_cell_state position
 * @param new_cell_state cell of interests
 * @param CA CellularAutomata instance for accesses model parameters
 * @return true : collision occurred
 * @return false : collision didn't occur
 */
bool did_galaxies_collied(int *cell_index, const std::vector<int> &offset_index, std::vector<int> &prev_offset, GalaxyCell &new_cell_state, CellularAutomata<GalaxyCell> &CA)
{
    const GalaxyCell cell = get_cell_state(cell_index, offset_index, CA);
    int cell_velocity_norm; // if cell collied then use this to scale new_cell_state velocity vector

    // check if there was a collision
    if (cell.state) // cell is occupied if state != 0
    {
        // collision occurred with current offset thus update cell_index with prev_offset
        std::vector<int> periodic_vec = get_periodic_vector(cell_index, prev_offset, CA);

        // std::cout << "collision occurred\n";
        // std::cout << "curr offset: " << offset_index[0] << ", " << offset_index[1] << ", " << offset_index[2] << "\n";
        // std::cout << "prev offset: " << prev_offset[0] << ", " << prev_offset[1] << ", " << prev_offset[2] << "\n";
        // std::cout << "new index: " << periodic_vec[0] << ", " << periodic_vec[1] << ", " << periodic_vec[2] << "\n";
        // std::cout << "old index: " << cell_index[0] << ", " << cell_index[1] << ", " << cell_index[2] << "\n";

        std::copy(periodic_vec.begin(), periodic_vec.end(), cell_index);

        // update_velocity_after_collision(new_cell_state, cell);
        return true;
    }
    else
    {
        // copy current offset to prev_offset
        std::copy(offset_index.begin(), offset_index.end(), prev_offset.begin());
        return false;
    }
}

/**
 * @brief Set the new_cell_state's new position.
 * This function uses the Bresenham's line algorithm for computing the path
 * from cell_index to cell_index + displacement_vector.
 * This method accounts for the possibility that cells might collied.
 *
 * Public domain algorithm and open source 3d implementation
 * Reference: https://antofthy.gitlab.io/info/graphics/bresenham.procs
 *
 * @param cell_index new_cell_state position
 * @param new_cell_state cell that contains mass and velocity
 * @param displacement_vector vector to the new desired position
 * @param CA CellularAutomata instance for accesses model parameters
 */
void set_new_position(int *cell_index, GalaxyCell &new_cell_state,
                      const std::vector<double> &displacement_vector,
                      CellularAutomata<GalaxyCell> &CA)
{
    int i, dx, dy, dz, l, m, n, x_inc, y_inc, z_inc,
        err_1, err_2, dx2, dy2, dz2;

    std::vector<int> offset_index(3, 0);
    std::vector<int> previous_offset(3, 0);

    // std::cout << "init old index: " << cell_index[0] << ", " << cell_index[1] << ", " << cell_index[2] << "\n";

    dx = round_int(displacement_vector.at(0));
    dy = round_int(displacement_vector.at(1));
    dz = round_int(displacement_vector.at(2));

    x_inc = (dx < 0) ? -1 : 1;
    l = abs(dx);
    y_inc = (dy < 0) ? -1 : 1;
    m = abs(dy);
    z_inc = (dz < 0) ? -1 : 1;
    n = abs(dz);
    dx2 = l << 1; // left-shift 1 time
    dy2 = m << 1;
    dz2 = n << 1;

    if ((l >= m) && (l >= n))
    {
        err_1 = dy2 - l;
        err_2 = dz2 - l;
        for (i = 0; i < l; i++)
        {
            if (did_galaxies_collied(cell_index, offset_index, previous_offset, new_cell_state, CA))
            {
                // collision occurred thus no need to compute rest of path
                return;
            }

            if (err_1 > 0)
            {
                offset_index[1] += y_inc;
                err_1 -= dx2;
            }
            if (err_2 > 0)
            {
                offset_index[2] += z_inc;
                err_2 -= dx2;
            }
            err_1 += dy2;
            err_2 += dz2;
            offset_index[0] += x_inc;
        }
    }
    else if ((m >= l) && (m >= n))
    {
        err_1 = dx2 - m;
        err_2 = dz2 - m;
        for (i = 0; i < m; i++)
        {
            if (did_galaxies_collied(cell_index, offset_index, previous_offset, new_cell_state, CA))
            {
                // collision occurred thus no need to compute rest of path
                return;
            }

            if (err_1 > 0)
            {
                offset_index[0] += x_inc;
                err_1 -= dy2;
            }
            if (err_2 > 0)
            {
                offset_index[2] += z_inc;
                err_2 -= dy2;
            }
            err_1 += dx2;
            err_2 += dz2;
            offset_index[1] += y_inc;
        }
    }
    else
    {
        err_1 = dy2 - n;
        err_2 = dx2 - n;
        for (i = 0; i < n; i++)
        {
            if (did_galaxies_collied(cell_index, offset_index, previous_offset, new_cell_state, CA))
            {
                // collision occurred thus no need to compute rest of path
                return;
            }

            if (err_1 > 0)
            {
                offset_index[1] += y_inc;
                err_1 -= dz2;
            }
            if (err_2 > 0)
            {
                offset_index[0] += x_inc;
                err_2 -= dz2;
            }
            err_1 += dy2;
            err_2 += dx2;
            offset_index[2] += z_inc;
        }
    }

    if (did_galaxies_collied(cell_index, offset_index, previous_offset, new_cell_state, CA))
    {
        // collision occurred thus no need to compute rest of path
        return;
    }

    // no collision occurred so just update position
    std::vector<int> periodic_vec = get_periodic_vector(cell_index, offset_index, CA);
    std::copy(periodic_vec.begin(), periodic_vec.end(), cell_index);

    // std::cout << "no collision\n";
    // std::cout << "curr offset: " << offset_index[0] << ", " << offset_index[1] << ", " << offset_index[2] << "\n";
    // std::cout << "prev offset: " << previous_offset[0] << ", " << previous_offset[1] << ", " << previous_offset[2] << "\n";
    // std::cout << "new index: " << periodic_vec[0] << ", " << periodic_vec[1] << ", " << periodic_vec[2] << "\n";
    return;
}

/**
 * @brief Custom CA rule for simulating the motion of star systems.
 *
 * @param cell_index new_cell_state position
 * @param index_size size of cell_index
 * @param neighborhood_cells array of neighboring cells for computing gravitational force
 * @param neighborhood_size size of neighborhood_cells array
 * @param new_cell_state reference to cell state that will be inserted into next state
 * @param CA CellularAutomata instance for accesses model parameters
 */
void galaxy_formation_rule(int *cell_index, const int index_size,
                           GalaxyCell *neighborhood_cells, const int neighborhood_size,
                           GalaxyCell &new_cell_state, CellularAutomata<GalaxyCell> &CA)
{
    if (new_cell_state.state == 0)
    {
        return; // empty cell
    }

    int vector_size = 3;    // x, y, z
    double time_step = 1.0; // can allow user III to modify
    std::vector<double> total_force_vector = {3, 0};

    // our neighborhood_cells also include our cell of interest
    int cell_of_interest_index = neighborhood_size / 2;

    // compute total vector
    for (int i = 0; i < neighborhood_size; i++)
    {
        // don't include the cell of interest in our force calculation
        if (cell_of_interest_index == i)
        {
            continue;
        }

        std::vector<int> neighbor_position(vector_size, 0);
        get_periodic_von_neumann_neighbor_index(vector_size, CA.boundary_radius, i, neighbor_position.data());

        std::vector<double> force_vector = compute_gravitational_force(new_cell_state,
                                                                       neighborhood_cells[i],
                                                                       neighbor_position);
        for (int j = 0; j < vector_size; j++)
        {
            total_force_vector.at(j) += force_vector.at(j);
        }
    }

    // std::vector<double> accel_vector = compute_accel(total_force_vector, new_cell_state.mass);
    // std::vector<double> velocity_vector = compute_velocity(accel_vector, time_step);
    // std::vector<double> displacement_vector = compute_displacement(velocity_vector, new_cell_state, time_step);
    std::vector<double> velocity_vector = {0.1, 0.1, 0.1};
    std::vector<double> displacement_vector = {3.1, 1.0, 2.2};

    // update velocity vector; set_new_position will updated new_cell_state.velocity it if a collision occurs
    // std::copy(velocity_vector.begin(), velocity_vector.end(), new_cell_state.velocity);
    new_cell_state.velocity_x = velocity_vector.at(0);
    new_cell_state.velocity_y = velocity_vector.at(1);
    new_cell_state.velocity_z = velocity_vector.at(2);

    set_new_position(cell_index, new_cell_state, displacement_vector, CA);
}

void test_motion_rule(int *cell_index, const int index_size, GalaxyCell *neighborhood_cells, const int neighborhood_size,
                      GalaxyCell &new_cell_state, CellularAutomata<GalaxyCell> &CA)
{
    if (new_cell_state.state != 0)
    {
        new_cell_state.state = 3;
        cell_index[0] = get_periodic_index(cell_index[0], 1, CA.axis1_dim);
        cell_index[1] = get_periodic_index(cell_index[1], 1, CA.axis2_dim);
        cell_index[2] = get_periodic_index(cell_index[2], 1, CA.axis3_dim);
    }
}

int main(int args, char **argv)
{
    CellularAutomata<GalaxyCell> galaxy = CellularAutomata<GalaxyCell>();
    galaxy.setup_boundary(CAEnums::Boundary::Periodic, 6);
    galaxy.setup_rule(CAEnums::Rule::Custom);
    galaxy.setup_dimensions_3d(1, 8, 8);
    galaxy.init_condition(1, 0.10);

    srand(time(NULL));
    int min_mass = 1;   // can allow User III to modify
    int max_mass = 100; // can allow User III to modify

    // set non-empty cell mass to a random value
    GalaxyCell ***all_cells = galaxy.get_tensor();
    for (int i = 0; i < galaxy.axis1_dim; i++)
    {
        for (int j = 0; j < galaxy.axis2_dim; j++)
        {
            for (int k = 0; k < galaxy.axis3_dim; k++)
            {
                if (all_cells[i][j][k].state != 0) // not empty
                {
                    all_cells[i][j][k].mass = (double)(min_mass + rand() % (max_mass - min_mass));
                }
            }
        }
    }

    galaxy.print_grid();
    galaxy.step(galaxy_formation_rule);
    galaxy.step(galaxy_formation_rule);
    galaxy.print_grid();
}
