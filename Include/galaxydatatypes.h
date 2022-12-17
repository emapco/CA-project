/**
 * @file galaxydatatypes.h
 * @author Emmanuel Cortes (ecortes@berkeley.edu)
 *
 * <b>Contributor(s)</b> <br> &emsp;&emsp;
 * @brief Header file with datatypes for creating a galaxy formation
 * cellular automata model.
 * @date 2022-12-13
 */
#pragma once
#include "CAdatatypes.h"
#include <algorithm> // copy
#include <vector>

/**
 * @brief Class that represents a star system.
 * The struct contains a velocity array and mass variable.
 * Satisfies CellularAutomata requirements
 * - state member variable
 * - inequality operator
 * - assignment operator
 * - copy constructor
 * - default constructor (required due to explicit copy constructor)
 */
class GalaxyCell
{
public:
    int state; //!< CellularAutomata requirement; 0: empty space

    // additional application specific member variables and static variables
    double velocity[3]; //!< velocity vector
    double mass;        //!< cell mass

    /**
     * @brief CellularAutomata requirement (inequality operator)
     *
     * @param other instance with which to compare data
     * @return true
     * @return false
     */
    bool operator!=(const GalaxyCell &other);

    /**
     * @brief CellularAutomata requirement (assignment operator)
     *
     * @param other instance from which to copy data
     * @return GalaxyCell&
     */
    GalaxyCell &operator=(const GalaxyCell &other);

    /**
     * @brief CellularAutomata requirement (copy constructor)
     *
     * @param other instance from which to copy data
     */
    GalaxyCell(const GalaxyCell &other);

    /**
     * @brief Construct a new Galaxy Cell object
     *
     */
    GalaxyCell();
};

/**
 * @brief Class for create a galaxy cellular automata model.
 * Only one instance should be created since the CellularAutomata instance is static and will be shared
 * among all instances.
 * CellularAutomata accepts a free function and not a member function thus galaxy_formation_rule needs to be static.
 * In order for galaxy_formation_rule to access the CellularAutomata and other class variable, they need to be static.
 *
 */
class Galaxy
{
public:
    int min_mass;                           //!< minimum cell mass (int required for rand)
    int max_mass;                           //!< maximum cell mass (int required for rand)
    double density;                         //!< density of galaxy
    int boundary_radius;                    //!< cutoff radius above which forces are not considered
    int axis1_dim;                          //!< cellular automata axis1 dimension
    int axis2_dim;                          //!< cellular automata axis2 dimension
    int axis3_dim;                          //!< cellular automata axis3 dimension
    static double time_step;                //!< time_step for computing forces during each simulation step
    static CellularAutomata<GalaxyCell> CA; //!< the CA the simulation utilizes to model the formation of a galaxy

    /**
     * @brief Construct a new Galaxy object using the following default parameters:<br>
     * min_mass = 1;<br>
     * max_mass = 100;<br>
     * density = 0.3;<br>
     * boundary_radius = 3;<br>
     * axis1_dim = 1;<br>
     * axis2_dim = 6;<br>
     * axis3_dim = 6;<br>
     * boundary_radius = 1;
     *
     */
    Galaxy();

    /**
     * @brief Create a Galaxy model instance using the provided values.
     * If given values are invalid then the default values are used.
     *
     * @param time_step simulation time step for computing forces and physical quantities
     * @param min_mass minimum cell mass
     * @param max_mass maximum cell mass
     * @param density density of occupied cell states
     * @param boundary_radius cutoff radius above which forces are not considered
     * @param axis1_dim cellular automata axis1 dimension
     * @param axis2_dim cellular automata axis2 dimension
     * @param axis3_dim cellular automata axis3 dimension
     */
    Galaxy(double time_step, int min_mass, int max_mass, double density, int boundary_radius, int axis1_dim, int axis2_dim, int axis3_dim);

    /**
     * @brief Destroy the Galaxy object
     *
     */
    ~Galaxy() = default;

    /**
     * @brief Sets up the galaxy CA instance.
     * This can be called multiple times if the user wants to restart the simulation
     * using different parameters.
     *
     * @return int - error code\n
     * Error codes returned by CellularAutomata instance\n
     * 0: no error
     */
    int init_galaxy();

    /**
     * @brief Runs the simulation for the specified amount of steps.
     *
     * @param steps
     * @return int - error code\n
     * Error codes returned by CellularAutomata instance\n
     * 0: no error
     */
    int simulation(int steps);

    /**
     * @brief Custom CA rule for simulating the motion of star systems.
     *
     * @param cell_index new_cell_state position
     * @param index_size size of cell_index
     * @param neighborhood_cells array of neighboring cells for computing gravitational force
     * @param neighborhood_size size of neighborhood_cells array
     * @param new_cell_state reference to cell state that will be inserted into next state
     */
    static void galaxy_formation_rule(int *cell_index, const int index_size,
                                      GalaxyCell *neighborhood_cells, const int neighborhood_size,
                                      GalaxyCell &new_cell_state);

    /**
     * @brief Set the new_cell_state's new position.
     * This function uses the Bresenham's line algorithm for computing the path
     * from cell_index to cell_index + displacement_vector.
     * This method accounts for the possibility that cells might collied.
     *
     * Public domain algorithm and open source 3d implementation reference: https://antofthy.gitlab.io/info/graphics/bresenham.procs
     *
     * @param cell_index new_cell_state position
     * @param new_cell_state cell that contains mass and velocity
     * @param displacement_vector vector to the new desired position
     */
    static void set_new_position(int *cell_index, GalaxyCell &new_cell_state,
                                 const std::vector<double> &displacement_vector);

    /**
     * @brief Determines if the new_cell_state collides with a cell at the current position: cell_index + offset_index.
     * If there is a collision then we set the new_cell_state index to the current position. 
     * And also increases the cell state by one to represent the number of collision. The mass and velocity properties are also updated using the laws of physics.
     * Else return false.
     *
     * @param cell_index new_cell_state position
     * @param offset_index keep track of path the cell takes
     * @param new_cell_state cell of interests
     * @return true : collision occurred
     * @return false : collision didn't occur
     */
    static bool did_galaxies_collide(int *cell_index, const std::vector<int> &offset_index, GalaxyCell &new_cell_state);

    /**
     * @brief compute the velocity after an inelastic collision
     * https://www.plasmaphysics.org.uk/collision3d.htm
     *
     * @param new_cell
     * @param collied_cell
     */
    static void update_velocity_after_collision(GalaxyCell &new_cell, GalaxyCell collied_cell);

    /**
     * @brief Get the cell state object at position cell_index + offset_index
     * from the CellularAutomata instance.
     *
     * @param cell_index current new_cell_state position
     * @param offset_index used to compute the neighboring cell position
     * @return const GalaxyCell&
     */
    static GalaxyCell &get_cell_state(int *cell_index, const std::vector<int> &offset_index);

    /**
     * @brief Get the periodic vector.
     *
     * @param cell_index current new_cell_state position
     * @param offset_index used to compute the neighboring cell position
     * @return std::vector<int>
     */
    static std::vector<int> get_periodic_vector(int *cell_index, const std::vector<int> &offset_index);

    /**
     * @brief Rounds a double to an int.
     *
     * @param dbl
     * @return int
     */
    static int round_int(double dbl);

    /**
     * @brief Computes the force between cell of interest and neighbor_index
     * Note: cell of interest is at (0, 0, 0)
     *
     * F = - m_1 * m_2 / | r_12 | * r_hat <br>
     * r_hat = r_12 / | r_12 |
     *
     * @param cell_of_interest
     * @param neighbor_cell
     * @param neighbor_index
     * @return std::vector<double>
     */
    static std::vector<double> compute_gravitational_force(const GalaxyCell &cell_of_interest, const GalaxyCell &neighbor_cell, const std::vector<int> &neighbor_index);

    /**
     * @brief Computes the acceleration vector
     * A = F/M
     *
     * @param total_force Sum of all forces in component form
     * @param mass cell of interet's mass
     * @return std::vector<double>
     */
    static std::vector<double> compute_accel(const std::vector<double> &total_force, double mass);

    /**
     * @brief Computes the velocity vector
     * 
     * V = velocity_i + acceleration_i * time_step
     *
     * @param accel current step acceleration
     * @param cell_of_interest contains initial velocity
     * @param time_step simulation time step
     * @return std::vector<double>
     */
    static std::vector<double> compute_velocity(const std::vector<double> &accel, const GalaxyCell &cell_of_interest, double time_step);

    /**
     * @brief Computes the displacement vector
     * 
     * D = 1/2 * (velocity_i + velocity_f) * time_step
     *
     * @param velocity current step velocity
     * @param cell_of_interest contains initial velocity
     * @param time_step simulation time step
     * @return std::vector<double>
     */
    static std::vector<double> compute_displacement(const std::vector<double> &velocity, const GalaxyCell &cell_of_interest, double time_step);

    /**
     * @brief Computes the norm of a vector. Used by compute_gravitational_force.
     *
     * @param vector
     * @return double
     */
    static double compute_vector_norm(const std::vector<int> &vector);

    /**
     * @brief Computes the difference in vector for determining the direction and magnitude of the force vector.
     *
     * @param vec1 Initial vector
     * @param vec2 Final vector
     * @return std::vector<double>
     */
    static std::vector<double> compute_vector_difference(const std::vector<double> &vec1, const std::vector<double> &vec2);
};
