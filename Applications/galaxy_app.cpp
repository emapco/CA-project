/**
 * @file galaxy_app.cpp
 * @author Emmanuel Cortes (ecortes@berkeley.edu)
 *
 * <b>Contributor(s)</b> <br> &emsp;&emsp;
 * @brief This application runs our galaxy model.
 * The application validates the user input so the application can appropriately.
 * @date 2022-12-14
 */

#include "galaxydatatypes.h"
#include <iostream>
#include <sstream>

/**
 * @brief Clear input buffer if input given is invalid.
 */
void input_failure()
{
    // not a valid number
    std::cout << "Invalid Input! Please input a valid numeric value.\n";
    std::cin.clear();
    while (std::cin.get() != '\n')
        ; // empty loop
}

/**
 * @brief Get a valid int from user in the given range.
 *
 * @param message prompt message
 * @param min minimum value that the user can specify (inclusive)
 * @param max maximum value that the user can specify (exclusive, optional if -1)
 * @return int
 */
int get_numeric_value(std::string message, int min, int max)
{
    int x;
    while (1)
    {
        std::cout << message;
        if (std::cin >> x)
        {
            if (x < min)
            {
                input_failure();
            }
            else if (max != -1 && x >= max)
            {
                input_failure();
            }
            else
            {
                break; // valid number
            }
        }
        else
        {
            input_failure();
        }
    }
    return x;
}

/**
 * @brief Get a valid double from user in the given range.
 *
 * @param message prompt message
 * @param min minimum value that the user can specify (exclusive)
 * @param max maximum value that the user can specify (exclusive, optional if -1)
 * @return double
 */
double get_numeric_value(std::string message, double min, double max)
{
    double x;
    while (1)
    {
        std::cout << message;
        if (std::cin >> x)
        {
            if (x <= min)
            {
                input_failure();
            }
            else if (max != -1.0 && x >= max)
            {
                input_failure();
            }
            else
            {
                break; // valid number
            }
        }
        else
        {
            input_failure();
        }
    }
    return x;
}

/**
 * @brief Gets user input and runs the simulation.
 *
 * @return int
 */
int main()
{
    // variables required for running the model
    int min_mass;
    int max_mass;
    double density;
    int boundary_radius;
    int axis1_dim;
    int axis2_dim;
    int axis3_dim;
    double time_step;
    int steps;
    // stringstreams for generating input prompt messages
    std::stringstream ss_radius;
    std::stringstream ss_min_mass;

    // get grid dimensions
    axis1_dim = get_numeric_value("Input the desired z dimension size (>= 3): ", 3, -1);
    axis2_dim = get_numeric_value("Input the desired x dimension size (>= 3): ", 3, -1);
    axis3_dim = get_numeric_value("Input the desired y dimension size (>= 3): ", 3, -1);

    // get cell mass range
    min_mass = get_numeric_value("Input the minimum mass a may have (>= 1): ", 1, -1);

    // generate message for maximum mass
    ss_min_mass << "Input the maximum mass a cell can have: (> " << min_mass << "): ";
    max_mass = get_numeric_value(ss_min_mass.str(), min_mass + 1, -1);

    // get grid density
    density = get_numeric_value("Input the desired density of the cellular automata grid (0.0 < density <= 1.0): ", 0.0, 1.0);

    // find min axis that user can input
    int min_axis = axis2_dim < axis3_dim ? axis2_dim : axis3_dim;
    int max_radius = min_axis / 2;

    // generate message for maximum boundary_radius
    ss_radius << "Input maximum distance to account for forces (1 <= distance <= " << max_radius << " ): ";

    // get boundary radius
    boundary_radius = get_numeric_value(ss_radius.str(), 1, max_radius + 1);

    // get simulation steps and time_step per step
    time_step = get_numeric_value("Input the desired simulation time_step (>= 0.1): ", 0.1, -1.0);
    steps = get_numeric_value("Input the number of steps the simulation should take (>= 1): ", 1, -1);
    std::cout << "\n";

    // create galaxy using user given parameters
    Galaxy galaxy = Galaxy(time_step, min_mass, max_mass, density,
                           boundary_radius, axis1_dim, axis2_dim, axis3_dim);
    galaxy.init_galaxy();
    // start simulation
    int error = galaxy.simulation(steps);
    return error;
}
