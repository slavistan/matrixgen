/**
 * Example detailing the basic usage of `matrixgen::adjmat` to generate
 * adjacency matrices using predefined weight functions and stencils.
 */
#include <array>
#include <iostream>

#include <matrixgen/core>

int main() {

  /**
   * 1. The Grid
   *
   * Choose the grid's extent as a `std::array<int, 3>`. The coordinates denote
   * the number of nodes in the x, y, z direction, with the grid's origin being
   * located at (0, 0, 0).
   *
   * For the examples we use tiny grids enough to keep the output matrix small
   * enough for your poor terminal.
   */
  const auto grid = std::array {3, 3, 1};

  /**
   * 2. The Stencil.
   *
   * Choose a stencil. The stencil encodes the information about which nodes
   * are to be considered adjacent and is hence responsible for the adjacency
   * matrix's non-zero elements' positions. We'll worry about the technical
   * details in a later example and use a preset symmetric 7-point stencil for
   * now.
   */
  const auto stencilfn = matrixgen::stencil7p();

  /**
   * 3. The Weight Function.
   *
   * The weight function determines the non-zeros' numerical values. We use
   * constant weights for this example implemented in a preset weight
   * function. Implementing a custom weight function will the subject of later
   * examples.
   */
  const auto weightfn = matrixgen::constweight(7.0);

  /**
   * Using the above building blocks we generate the matrix and use it.
   * By default `adjmat` returns an uncompressed row-major
   * `Eigen::SparseMatrix` using a double as its scalar type.
   */
  const auto mat = matrixgen::adjmat(grid, stencilfn, weightfn);
  std::cout << std::endl << mat << std::endl;
}
