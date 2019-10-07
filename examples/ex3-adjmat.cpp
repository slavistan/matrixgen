/**
 * Example detailing the definition of custom weight functions.
 */
#include <array>
#include <iostream>

#include <matrixgen/core>

int main()
{
  /**
   * Custom weightfunctions are defined as lambdas with correct signatures.
   * The function shall return a scalar convertible to the matrix's scalar
   * type and shall be invokable with the following argument types:
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
   * During construction of the matrix the grid is traversed in x-then-y-then-z
   * direction and the numerical value for each connection between node and
   * its neighbors is determined according to the weightfunction and the
   * context provided by it.
   *
   * For this example we choose the second signature in order to define a
   * weightfunction which shall return `-1` for any diagonal element.
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

  const auto mat = matrixgen::adjmat({3, 3, 1}, matrixgen::stencil7p(), weightfn);
  std::cout << std::endl << Eigen::MatrixXd(mat) << std::endl;
}
