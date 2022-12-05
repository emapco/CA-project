
//  UC Berkeley - MSSE Program
//  Chem 279-B  Introduction to Software Engineering
//  Fall 2022
//
// his file contains the custom datatypes used in the Cellular Automata project
//
#pragma once    // Ensures that this file is only included once
                // during compilation

#include <array>
using namespace std;

enum neighborhood{VonNeumann, Moore};

enum boundary{None, Periodic, Walled, CutOff};

enum rule{Majority, Parity};

class CellularAutomata         // The CA datastructure
{
    public:   

        boundary boundary_type;     //enum code to hold boundary
        int boundary_radius;       //declare a radius for the boundary
        neighborhood neighborhood_type;     //enum code to hold neighboorhood type 
        int num_states;         //integer code for number of states
        
        rule rule_type;          //enum code for rule type 
        double shortr_weight;   //double for short radius
        double longr_weight;        //double for long radius
        int n;              //count of cells in first dimension
        int m;              //count of cellls in second dimension
        int **matrix;        //pointer to 2d array for grid cells holding a state 

        ///////////////////////////////////////////
        CellularAutomata(int n, int m);              // Default constructor 
        ~CellularAutomata();             // Destructor

        //provide a setup function to choose from two neighborhood types
        int setup_neighborhood(neighborhood neighborhood_type);

        //provide a setup function to choose from four boundary types
        int setup_boundary(boundary bound_type, int radius);

        //provide a setup function to choose the number of states
        int setup_cell_states(int n_states);

        //provide a function to initialize the beginning state of the grid
        int init_condition(int x_state, double prob);

        int setup_rule_short_long(double shortr_weight, double longr_weight);

        int setup_rule(rule rule_type);

        //provide a function to print the current grid
        int print_grid();

};
