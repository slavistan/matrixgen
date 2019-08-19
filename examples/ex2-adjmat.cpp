// CONTINUEHERE
// TODO: Show off 
//        different boundary conditions --> Stencil also encodes bcs
//        randomweight
//        convert to MatrixXd for pretty print
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
  auto weightfn = []() { return 1; };
int main()
{
  /**
   * Create the very same matrix from example 1 but use available presets for
   * the stencil and the weightfunction. See 'matrixgen/presets.hpp' for other
   * predefined stencils and weighfunctions.
   */
  const auto mat = matrixgen::adjmat({3, 3, 1}, matrixgen::STENCIL<7>, matrixgen::constweight(1));
  std::cout << std::endl << mat << std::endl;
}
