/**
 * Example detailing the basic usage of `matrixgen::adjmat` to generate
 * adjacency matrices given a set of boundary conditions using preset weight
 * and stencil functions.
 */
#include <array>
#include <iostream>

#include <matrixgen/core>

int main()
{
  /**
   * Adjacency matrices for structure grids are generated adhering to a simple
   * schema: The stencil function is responsible for the location of all
   * non-zeros while the weight function determines their values.
   *
   * Thus boundary conditions are implemented within the scope of the stencil
   * function. For this example we use the preset stencil function to declare
   * periodic boundary conditions in the X-dimension and Dirichlet bcs for the
   * Y and Z directions. This is a common configuration for the simulation
   * of 2D flows through a pipe, where conservation of flow matter implies
   * periodic boundary conditions along the flow path.
   */
  using matrixgen::BC;
  const auto stencilfn = matrixgen::stencil7p<BC::PERIODIC, BC::DIRICHLET, BC::DIRICHLET>();
  /**                                         ^~~~~ X ~~~^  ^~~~~ Y ~~~~^  ^~~~~ Z ~~~~^
   */

  /**
   * Convert to a `Eigen::Matrix` for pretty-print. Note the inline usage of a
   * different weight function picking random weights from the unit interval.
   *
   * We convert the matrix to a dense `Eigen::Matrix` for its pretty-print
   * capabilities.
   */
  const auto mat = matrixgen::adjmat({4, 4, 1}, stencilfn, matrixgen::randweight(1));
  std::cout << std::endl << Eigen::MatrixXd(mat) << std::endl;
}
