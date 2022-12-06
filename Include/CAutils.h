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