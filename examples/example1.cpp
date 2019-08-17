/**
 * Detailed and verbosely commented example detailing the generation of basic
 * adjacency matrices.
 */
#include <array>
#include <iostream>

#include <matrixgen/core>

int main()
{

  /**
   * Choose the grid's extent as a `std::array<int, 3>`. The coordinates denote
   * the number of nodes in the x, y, z direction, with the grid's origin being
   * located at (0, 0, 0).
   *
   * For the examples we choose a small grid to keep the output matrix small
   * enough for your poor terminal.
   */
  const auto grid = std::array {3, 3, 1};

  /**
   * Choose a stencil. The stencil encodes which nodes are considered adjacent
   * as a set of 3d offsets. A stencil shall be a range of arbitrary length
   * over elements of type 'std::array<int, 3>' which is traversed once for
   * each node in the grid. We use a symmetric 7p-stencil for this example.
   */
  const auto stencil = std::array<std::array<int, 3>, 7>{{
      {1, 0, 0}, {-1, 0,0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}, {0, 0, 0}
  }};

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

  /**
   * Generate the matrix and use it.
   */
  const auto mat = matrixgen::adjmat(grid, stencil, weightfn);
  std::cout << std::endl << mat << std::endl;
}
