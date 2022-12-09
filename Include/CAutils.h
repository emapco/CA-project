/**
 * @file CAutils.h
 * @author Emmanuel Cortes (ecortes@berkeley.edu)
 *
 * <b>Contributor(s)</b> <br> &emsp;&emsp;
 * @brief This header files contains utility functions used by CellularAutomata class.
 * @date 2022-12-06
 */

/**
 * @brief swaps the computed next_vector to the current vector
 *
 * @param vector cellular automata current vector state
 * @param next_vector cellular automata next vector state
 * @param axis1_dim vector dimension
 */
void swap_states(int *vector, int *next_vector, int axis1_dim);

/**
 * @brief swaps the computed next_matrix to the current matrix
 *
 * @param matrix cellular automata current matrix state
 * @param next_matrix cellular automata next matrix state
 * @param axis1_dim first matrix dimension
 * @param axis2_dim second matrix dimension
 */
void swap_states(int **matrix, int **next_matrix, int axis1_dim, int axis2_dim);

/**
 * @brief swaps the computed next_tensor to the current tensor
 *
 * @param tensor cellular automata current tensor state
 * @param next_tensor cellular automata next tensor state
 * @param axis1_dim first tensor dimension
 * @param axis2_dim second tensor dimension
 * @param axis3_dim third tensor dimension
 */
void swap_states(int ***tensor, int ***next_tensor, int axis1_dim, int axis2_dim, int axis3_dim);

/**
 * @brief calculates the neighborhood size using the rank, radius and neighborhood type
 *
 * @param rank rank of the cell's tensor
 * @param radius radius of the neighborhood
 * @param neighborhood_type the type of neighborhood
 * @return int
 */
int get_neighborhood_size(int rank, int radius, CellularAutomata::neighborhood neighborhood_type);

/**
 * @brief determines if a cell is diagonal to the central cell\n
 * diagram of slice: 1 = diagonal cell\n
 * 1 0 1\n
 * 0 0 0\n
 * 1 0 1\n
 *
 * @param i cell's i-th index
 * @param j cell's j-th index
 * @return true: cell is diagonal
 * @return false: cell is not diagonal
 */
bool is_diagonal_neighboring_cell_2d(int i, int j);

/**
 * @brief determines if a cell is diagonal to the central cell\n
 * if k == 0 then both i and j have be non-zero\n
 * diagram of slice: 1 = diagonal cell\n
 * 1 0 1\n
 * 0 0 0\n
 * 1 0 1\n
 *
 * if k != 0 then either i or j can be non-zero\n
 * diagram of slice: 1 = diagonal cell\n
 * 1 1 1\n
 * 1 0 1\n
 * 1 1 1\n
 *
 * @param i cell's i-th index
 * @param j cell's j-th index
 * @param k cell's k-th index
 * @return true: cell is diagonal
 * @return false: cell is not diagonal
 */
bool is_diagonal_neighboring_cell_3d(int i, int j, int k);

/**
 * @brief initializes a MajorityCounter instance utilized by Majority rule
 *
 * @param counter object to keep number of votes for each particular cell state
 * @param num_states number of different cell states
 */
void initialize_majority_rule_counter(MajorityCounter &counter, int num_states);

/**
 * @brief get the periodic index used for finding periodic boundary neighbors
 *
 * @param i cell's i-th index
 * @param di neighbor's cell i-th offset
 * @param axis_dim dimension of i-th axis
 * @return int
 */
int get_periodic_index_axis(int i, int di, int axis_dim);
