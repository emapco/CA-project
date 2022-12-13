/**
 * @file unit_test_CA_utils.cpp
 * @author Emmanuel Cortes (ecortes@berkeley.edu)
 *
 * <b>Contributor(s)</b> <br> &emsp;&emsp;
 * @brief The program contains unit tests to test the correctness of
 * the various functions defined in CA_utils.cpp.
 * @date 2022-12-11
 */

#include "CAdatatypes.h"
#include <cassert>
#include <iostream>
#include <vector>
#include <utility> // pair

/**
 * @brief Prints that a specific test passed.
 *
 * @param message the test function name
 */
void print_success(std::string message)
{
    std::cout << "TEST PASSED: " << message << "\n";
}

/**
 * @brief Check returned periodic Moore neighbor index follows the same
 * pattern used to create the flatten neighborhood array.
 *
 * @param radius neighborhood radius
 */
void test_get_periodic_moore_neighbor_index(int radius)
{
    int vector_i = 0;
    int matrix_i = 0;
    int tensor_i = 0;

    int eval[3] = {0, 0, 0};

    // check returned periodic moore neighbor index follows the same
    // pattern used to create the neighborhood
    for (int di = -radius; di <= radius; di++)
    {
        get_periodic_moore_neighbor_index(1, radius, vector_i, eval);

        assert((di == eval[0]));
        vector_i++;
    }

    for (int di = -radius; di <= radius; di++)
    {
        for (int dj = -radius; dj <= radius; dj++)
        {
            get_periodic_moore_neighbor_index(2, radius, matrix_i, eval);
            assert((di == eval[0]));
            assert((dj == eval[1]));
            matrix_i++;
        }
    }

    for (int di = -radius; di <= radius; di++)
    {
        for (int dj = -radius; dj <= radius; dj++)
        {
            for (int dk = -radius; dk <= radius; dk++)
            {
                get_periodic_moore_neighbor_index(3, radius, tensor_i, eval);
                assert((di == eval[0]));
                assert((dj == eval[1]));
                assert((dk == eval[2]));
                tensor_i++;
            }
        }
    }
}

/**
 * @brief Check returned periodic Moore neighbor index follows the same
 * pattern used to create the flatten neighborhood array.
 *
 * @param radius neighborhood radius
 */
void test_get_periodic_von_neumann_neighbor_index(int radius)
{
    int vector_i = 0;
    int matrix_i = 0;
    int tensor_i = 0;
    int eval[3] = {0, 0, 0};

    // vector case
    for (int di = -radius; di <= radius; di++)
    {
        get_periodic_von_neumann_neighbor_index(1, radius, vector_i, eval);
        assert((di == eval[0]));
        vector_i++;
    }
    // matrix case
    for (int di = -radius; di <= radius; di++)
    {
        for (int dj = -radius; dj <= radius; dj++)
        {
            // exclude diagonal cells from neighborhood when VonNeumann is selected
            if (is_diagonal_neighboring_cell_2d(di, dj))
            {
                // current di,dj,dk cell is a diagonal neighbor so exclude it
                continue;
            }
            get_periodic_von_neumann_neighbor_index(2, radius, matrix_i, eval);
            assert((di == eval[0]));
            assert((dj == eval[1]));
            matrix_i++;
        }
    }
    // tensor case
    for (int di = -radius; di <= radius; di++)
    {
        for (int dj = -radius; dj <= radius; dj++)
        {
            for (int dk = -radius; dk <= radius; dk++)
            {
                // exclude diagonal cells from neighborhood when VonNeumann is selected
                if (is_diagonal_neighboring_cell_3d(di, dj, dk))
                {
                    // current di,dj,dk cell is a diagonal neighbor so exclude it
                    continue;
                }
                get_periodic_von_neumann_neighbor_index(3, radius, tensor_i, eval);
                assert((di == eval[0]));
                assert((dj == eval[1]));
                assert((dk == eval[2]));
                tensor_i++;
            }
        }
    }
}

/**
 * @brief Calls get_neighborhood_size and asserts that the returned value
 * equals the expected value.
 *
 * @param tests_and_answer Vector of pairs where first is radius test case and second is the expect value
 * @param rank cell data rank
 * @param nt neighborhood type
 */
void test_get_neighborhood_size_assert(std::vector<std::pair<int, int>> test_arg_and_expected_value, int rank, CAEnums::Neighborhood nt)
{
    for (const auto &pair : test_arg_and_expected_value)
    {
        assert((CellularAutomata<int>::get_neighborhood_size(rank, pair.first, nt) == pair.second));
    }
}

/**
 * @brief  Tests every type of combination possible that get_neighborhood_size is supposed to support.
 */
void test_get_neighborhood_size()
{
    // first is radius; second is expected value for asserting
    std::vector<std::pair<int, int>> vec_moore_tests_expected = {{1, 3}, {2, 5}, {5, 11}, {9, 19}};
    test_get_neighborhood_size_assert(vec_moore_tests_expected, 1, CAEnums::Moore);

    std::vector<std::pair<int, int>> matrix_moore_tests_expected = {{1, 9}, {2, 25}, {5, 121}, {9, 361}};
    test_get_neighborhood_size_assert(matrix_moore_tests_expected, 2, CAEnums::Neighborhood::Moore);

    std::vector<std::pair<int, int>> tensor_moore_tests_expected = {{1, 27}, {2, 125}, {5, 1331}, {9, 6859}};
    test_get_neighborhood_size_assert(tensor_moore_tests_expected, 3, CAEnums::Neighborhood::Moore);

    std::vector<std::pair<int, int>> vec_von_neumann_tests_expected = {{1, 3}, {2, 5}, {5, 11}, {9, 19}};
    test_get_neighborhood_size_assert(vec_von_neumann_tests_expected, 1, CAEnums::Neighborhood::VonNeumann);

    std::vector<std::pair<int, int>> matrix_von_neumann_tests_eval_expected = {{1, 5}, {2, 9}, {5, 21}, {9, 37}};
    test_get_neighborhood_size_assert(matrix_von_neumann_tests_eval_expected, 2, CAEnums::Neighborhood::VonNeumann);

    std::vector<std::pair<int, int>> tensor_von_neumann_tests_eval_expected = {{1, 7}, {2, 13}, {5, 31}, {9, 55}};
    test_get_neighborhood_size_assert(tensor_von_neumann_tests_eval_expected, 3, CAEnums::Neighborhood::VonNeumann);

    print_success("test_get_neighborhood_size");
}

/**
 * @brief Tests the correctness of the is_diagonal_neighboring_cell_2d method.
 */
void test_is_diagonal_neighboring_cell_2d()
{
    std::vector<std::vector<int>> indices = {{0, 0}, {0, 1}, {1, 1}, {2, 1}, {10, 0}, {8, 9}};
    std::vector<bool> expected = {false, false, true, true, false, true};
    for (auto i = 0; i < indices.size(); i++)
    {
        int x = indices.at(i).at(0);
        int y = indices.at(i).at(1);
        assert((is_diagonal_neighboring_cell_2d(x, y) == expected.at(i)));
    }
    print_success("test_is_diagonal_neighboring_cell_2d");
}

/**
 * @brief Tests the correctness of the is_diagonal_neighboring_cell_3d method.
 */
void test_is_diagonal_neighboring_cell_3d()
{
    int x, y, z; // indices
    std::vector<std::vector<int>> indices = {{0, 0, 1}, {0, 1, 1}, {1, 1, -1}, {2, 1, 0}, {-10, 0, -7}, {8, 9, 0}, {-2, 0, 0}};
    std::vector<bool> expected = {false, true, true, true, true, true, false};
    for (auto i = 0; i < indices.size(); i++)
    {
        x = indices.at(i).at(0);
        y = indices.at(i).at(1);
        z = indices.at(i).at(2);
        assert((is_diagonal_neighboring_cell_3d(x, y, z) == expected.at(i)));
    }
    print_success("test_is_diagonal_neighboring_cell_3d");
}

/**
 * @brief Tests the correctness of the get_periodic_index method.
 */
void test_get_periodic_index()
{
    int j, dj, axis_dim; // function args
    std::vector<std::vector<int>> indices = {{0, -1, 2}, {2, 1, 3}, {5, 2, 8}, {0, -2, 6}, {1, -3, 4}};
    std::vector<int> expected = {{1}, {0}, {7}, {4}, {2}};

    for (auto i = 0; i < indices.size(); i++)
    {
        j = indices.at(i).at(0);
        dj = indices.at(i).at(1);
        axis_dim = indices.at(i).at(2);
        assert((get_periodic_index(j, dj, axis_dim) == expected.at(i)));
    }
    print_success("test_get_periodic_index");
}

int main()
{
    // ensure the methods work for various neighborhood radii
    for (int i = 1; i <= 50; i++)
    {
        test_get_periodic_moore_neighbor_index(i);
        test_get_periodic_von_neumann_neighbor_index(i);
    }
    print_success("test_get_periodic_moore_neighbor_index");
    print_success("test_get_periodic_von_neumann_neighbor_index");

    test_get_neighborhood_size();
    test_is_diagonal_neighboring_cell_2d();
    test_is_diagonal_neighboring_cell_3d();
    test_get_periodic_index();

    return 0;
}
