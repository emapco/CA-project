/**
 * @file galaxy.cpp
 * @author Emmanuel Cortes (ecortes@berkeley.edu)
 *
 * <b>Contributor(s)</b> <br> &emsp;&emsp; Trevor Oldham, Chongye Feng
 * @brief The file contains the implementations for the various
 * classes defined in galaxydatatypes.h
 * @date 2022-12-13
 */
#include "galaxydatatypes.h"
#include "CAdatatypes.h"
#include <algorithm> // copy
#include <iostream>
#include <vector>
#include <cmath> // pow

CellularAutomata<GalaxyCell> Galaxy::CA = CellularAutomata<GalaxyCell>();
double Galaxy::time_step = 0.1;

Galaxy::Galaxy()
{
    time_step = 0.1;
    min_mass = 1;
    max_mass = 100;
    density = 0.3;
    boundary_radius = 3;
    axis1_dim = 1;
    axis2_dim = 6;
    axis3_dim = 6;
}

Galaxy::Galaxy(double time_step, int min_mass, int max_mass, double density, int boundary_radius, int axis1_dim, int axis2_dim, int axis3_dim)
{
    if (time_step <= 0)
    {
        this->time_step = 0.1;
        std::cout << "Invalid time_step. time_step must be > 0. Using default "
                  << this->time_step << "\n";
    }
    else
    {
        this->time_step = time_step;
    }

    if (min_mass < 1)
    {
        this->min_mass = 1;
        std::cout << "Invalid min_mass. min_mass must be >= 1. Using default "
                  << this->min_mass << "\n";
    }
    else
    {
        this->min_mass = min_mass;
    }

    if (max_mass < min_mass)
    {
        this->max_mass = 100;
        std::cout << "Invalid max_mass. max_mass must be >= min_mass. Using default "
                  << this->max_mass << "\n";
    }
    else
    {
        this->max_mass = max_mass;
    }

    if (density <= 0.0 || density > 1.0)
    {
        this->density = 0.3;
        std::cout << "Invalid density. 0 < density must be <= 1. Using default "
                  << this->density << "\n";
    }
    else
    {
        this->density = density;
    }

    if (axis1_dim < 1)
    {
        this->axis1_dim = 1;
        std::cout << "Invalid axis1_dim. axis1_dim must be > 1. Using default "
                  << this->axis1_dim << "\n";
    }
    else
    {
        this->axis1_dim = axis1_dim;
    }

    if (axis2_dim <= 2)
    {
        this->axis2_dim = 6;
        std::cout << "Invalid axis2_dim. axis2_dim must be > 2. Using default "
                  << this->axis2_dim << "\n";
    }
    else
    {
        this->axis2_dim = axis2_dim;
    }

    if (axis3_dim <= 2)
    {
        this->axis3_dim = 6;
        std::cout << "Invalid axis3_dim. axis3_dim must be > 2. Using default "
                  << this->axis3_dim << "\n";
    }
    else
    {
        this->axis3_dim = axis3_dim;
    }

    // find min axis
    int min_axis = this->axis2_dim < this->axis3_dim ? this->axis2_dim : this->axis3_dim;

    if (boundary_radius > min_axis / 2 || boundary_radius <= 0)
    {
        this->boundary_radius = min_axis / 2;
        std::cout << "Invalid boundary_radius. boundary_radius must be <= half the smallest axis dimensions and > 0. "
                  << "Setting to " << this->boundary_radius << "\n";
    }
    else
    {
        this->boundary_radius = boundary_radius;
    }
}

int Galaxy::init_galaxy()
{
    int error = 0; // store error return by CA
    int i, j, k;   // iterators

    error = CA.setup_dimensions_3d(axis1_dim, axis2_dim, axis3_dim);
    if (error == CAEnums::CellsAlreadyInitialized)
    {
        // If already called init_galaxy then remove old instance and deinitialize it.
        // This allows our model to be restarted with in the same application.
        CA = CellularAutomata<GalaxyCell>();
        error = CA.setup_dimensions_3d(axis1_dim, axis2_dim, axis3_dim);
    }

    if (error < 0)
    {
        CA.print_error_status(static_cast<CAEnums::ErrorCode>(error));
        return error;
    }
    error = CA.setup_boundary(CAEnums::Boundary::Periodic, boundary_radius);
    if (error < 0)
    {
        CA.print_error_status(static_cast<CAEnums::ErrorCode>(error));
        return error;
    }
    CA.setup_rule(CAEnums::Rule::Custom);
    CA.init_condition(1, density);

    srand(time(NULL));

    // set non-empty cell mass to a random value
    GalaxyCell ***all_cells = CA.get_tensor();
    for (i = 0; i < CA.axis1_dim; i++)
    {
        for (j = 0; j < CA.axis2_dim; j++)
        {
            for (k = 0; k < CA.axis3_dim; k++)
            {
                if (all_cells[i][j][k].state != 0) // not empty
                {
                    all_cells[i][j][k].mass = (double)(min_mass + rand() % (max_mass - min_mass));
                }
            }
        }
    }

    std::cout << "**** Starting simulation ****\n";
    CA.print_grid();
    return error;
}

int Galaxy::simulation(int steps)
{
    int error, i;
    if (steps < 1)
    {
        std::cout << "Invalid number of steps. steps must be > 0. Setting steps to 1\n";
        steps = 1;
    }

    for (i = 0; i < steps; i++)
    {
        error = CA.step(galaxy_formation_rule);
        if (error < 0)
        {
            CA.print_error_status(static_cast<CAEnums::ErrorCode>(error));
            return error;
        }
    }
    std::cout << "**** Simulation's final state ****\n";
    CA.print_grid();

    return error;
}

void Galaxy::galaxy_formation_rule(int *cell_index, const int index_size,
                                   GalaxyCell *neighborhood_cells, const int neighborhood_size,
                                   GalaxyCell &new_cell_state)
{
    if (new_cell_state.state == 0)
    {
        return; // empty cell
    }

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

        std::vector<int> neighbor_position(index_size, 0);
        get_periodic_von_neumann_neighbor_index(index_size, CA.boundary_radius, i, neighbor_position.data());

        // std::vector<double> force_vector = compute_gravitational_force(new_cell_state,
        //                                                                neighborhood_cells[i],
        //                                                                neighbor_position);
        // for (int j = 0; j < vector_size; j++)
        // {
        //     total_force_vector.at(j) += force_vector.at(j);
        // }
    }

    // std::vector<double> accel_vector = compute_accel(total_force_vector, new_cell_state.mass);
    // std::vector<double> velocity_vector = compute_velocity(accel_vector, time_step);
    // std::vector<double> displacement_vector = compute_displacement(velocity_vector, new_cell_state, time_step);
    std::vector<double> velocity_vector = {0.1, 0.1, 0.1};
    std::vector<double> displacement_vector = {3.1, 1.0, 2.2};

    // update velocity vector
    // set_new_position will updated new_cell_state.velocity if a collision occurs
    std::copy(velocity_vector.begin(), velocity_vector.end(), new_cell_state.velocity);

    set_new_position(cell_index, new_cell_state, displacement_vector);
}

void Galaxy::set_new_position(int *cell_index, GalaxyCell &new_cell_state,
                              const std::vector<double> &displacement_vector)
{
    int i, dx, dy, dz, l, m, n, x_inc, y_inc, z_inc,
        err_1, err_2, dx2, dy2, dz2;

    std::vector<int> offset_index(3, 0);

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
            if (did_galaxies_collide(cell_index, offset_index, new_cell_state))
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
            if (did_galaxies_collide(cell_index, offset_index, new_cell_state))
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
            if (did_galaxies_collide(cell_index, offset_index, new_cell_state))
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

    if (did_galaxies_collide(cell_index, offset_index, new_cell_state))
    {
        // collision occurred thus no need to compute rest of path
        return;
    }

    // no collision occurred so just update position
    std::vector<int> periodic_vec = get_periodic_vector(cell_index, offset_index);
    std::copy(periodic_vec.begin(), periodic_vec.end(), cell_index);

    // std::cout << "no collision\n";
    // std::cout << "curr offset: " << offset_index[0] << ", " << offset_index[1] << ", " << offset_index[2] << "\n";
    // std::cout << "new index: " << periodic_vec[0] << ", " << periodic_vec[1] << ", " << periodic_vec[2] << "\n";
    return;
}

void Galaxy::update_velocity_after_collision(GalaxyCell &new_cell, GalaxyCell collied_cell)
{
    double m1, m2, v1, v2;
    for (int i = 0; i < 3; i++)
    {
        m1 = new_cell.mass;
        m2 = collied_cell.mass;
        v1 = new_cell.velocity[i];
        v2 = collied_cell.velocity[i];
        new_cell.velocity[i] = (m1 * v1 + m2 * v1) / (m1 + m2);
    }
}

int Galaxy::round_int(double dbl)
{
    return static_cast<int>(dbl < 0 ? dbl - 0.5 : dbl + 0.5);
}

std::vector<int> Galaxy::get_periodic_vector(int *cell_index, const std::vector<int> &offset_index)
{
    int x = get_periodic_index(cell_index[0], offset_index[0], CA.axis1_dim);
    int y = get_periodic_index(cell_index[1], offset_index[1], CA.axis2_dim);
    int z = get_periodic_index(cell_index[2], offset_index[2], CA.axis3_dim);
    return {x, y, z};
}

GalaxyCell &Galaxy::get_cell_state(int *cell_index, const std::vector<int> &offset_index)
{
    std::vector<int> vec = get_periodic_vector(cell_index, offset_index);
    GalaxyCell ***next_tensor = CA.get_next_tensor();
    return next_tensor[vec[0]][vec[1]][vec[2]];
}

bool Galaxy::did_galaxies_collide(int *cell_index, const std::vector<int> &offset_index, GalaxyCell &new_cell_state)
{
    const GalaxyCell cell = get_cell_state(cell_index, offset_index);

    // check if there was a collision
    if (cell.state) // cell is occupied if state != 0
    {
        // collision occurred with current offset and collision type is set to inelastic thus merge the cells
        // at the collided index
        std::vector<int> periodic_new_cell_vec = get_periodic_vector(cell_index, offset_index);
        update_velocity_after_collision(new_cell_state, cell);
        std::copy(periodic_new_cell_vec.begin(), periodic_new_cell_vec.end(), cell_index);
        new_cell_state.state += cell.state;
        new_cell_state.mass += cell.mass;
        return true;
    }
    return false;
}

bool GalaxyCell::operator!=(const GalaxyCell &other)
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

GalaxyCell &GalaxyCell::operator=(const GalaxyCell &other)
{
    if (this == &other) // Guard self assignment
        return *this;
    // copy variables
    state = other.state;
    mass = other.mass;
    std::copy(other.velocity, other.velocity + 3, velocity);
    return *this;
}

GalaxyCell::GalaxyCell(const GalaxyCell &other) : state(other.state), mass(other.mass)
{
    std::copy(other.velocity, other.velocity + 3, velocity);
}

GalaxyCell::GalaxyCell() : state(0), mass(0)
{
    // fill array with zeros
    for (int i = 0; i < 3; i++)
    {
        velocity[i] = 0.0;
    }
}
