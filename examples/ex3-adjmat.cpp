/**
 * TODO: Display creation of a custom weight function
 */
#include <array>
#include <iostream>

#include <matrixgen/core>

  /**
   * Define a lambda to compute the matrix's values i.e. the weights of the
   * adjacency matrix. The function shall return a scalar convertible to the
   * matrix's scalar type and shall take be invokable with the following
   * argument types:
   *
   * 1. Take no argument.
   *  foo ()
   *
   * 2. Take the matrix entry's coordinates as a pair of (row, column).
   *  foo (std::array<int, 2> coords)
   *
   * 3. Take the geometric positions of the grid node in question and its
   *    neighbor as tuples of (x, y, z). The node's coordinates are passed
   *    first followed by its neighbor's coordinates.
   *  foo (std::array<int, 3> me, std::array<int, 3> neighbor)
   *
   * 4. Combination of (2) and (3). First the matrix entry's coordinates
   *    are passed, then the nodes' geometric positions.
   *  foo (std::array<int, 2> coords, std::array<int, 3> me, std::array<int, 3> neighbor)
   *
   * For this example we choose the most simple weight-function returning a
   * constant value. See the other examples for usages of different weight
   * functions.
   */
int main()
{
  /**
   * Define a weightfunction which shall return `-1` for any diagonal element.
   * The elements in the upper trianglar matrix shall be the sum of the row
   * and column indices (starting at `1` to avoid nulls) while the the lower
   * triangular matrix's values shall be the negatives thereof.
   */
  auto weightfn = [](const std::array<int, 2>& coords) {

    const auto [row, col] = coords;
    if (row == col) { return -1; }            // diagonal
    if (row < col) { return row + col + 2; }  // upper triangular matrix
    return -(row + col + 2);                  // lower triangular matrix
  };

  /**
   * Generate and use. We convert the output to a dense matrix for pretty
   * output.
   */
  const auto mat = matrixgen::adjmat({3, 3, 1}, matrixgen::STENCIL<7>, weightfn);
  std::cout << std::endl << Eigen::MatrixXd(mat) << std::endl;
}
