/**
 * Example detailing the basic usage of `matrixgen::adjmat` to generate
 * adjacency matrices given a set of boundary conditions using preset weight
 * and stencil functions.
 */
#include <array>
#include <iostream>

#include <matrixgen/core>

#include <Eigen/Core>

int main()
{
  /**
   * Choose an adjacency function and boundary conditions.
   *
   * Adjacency matrices for grids are generated adhering to a simple
   * schema: The adjacency function is responsible for the location of all
   * non-zeros while the weight function determines their values.
   *
   * Thus boundary conditions are implemented within the scope of the stencil
   * function. For this example we use the preset stencil function to declare
   * periodic boundary conditions in the X-dimension and Dirichlet BCs for the
   * Y and Z directions. This is a common configuration for the simulation
   * of 2D flows through a pipe, where conservation of flow matter implies
   * periodic boundary conditions along the flow path.
   */
  using matrixgen::BC;
  const auto adjfn = matrixgen::stencil7p<BC::PERIODIC, BC::DIRICHLET, BC::DIRICHLET>();
  //                                      ^~~~~ X ~~~^  ^~~~~ Y ~~~~^  ^~~~~ Z ~~~~^

  /**
   * Generate the matrix and, optionally, pick a different output matrix type.
   * Any valid specialization of `Eigen::SparseMatrix` will do.
   *
   * Note the inline usage of a different weight function picking random
   * weights from the unit interval. Also, we pass the grid dimensions as a
   * braced-init-list.
   *
   * Finally, we convert the matrix to a dense `Eigen::Matrix` for its pretty-print
   * capabilities.
   */
  using Scalar_t = float; // Try double here
  using Matrix_t = Eigen::SparseMatrix<Scalar_t, Eigen::ColMajor>; // Try RowMajor here
  const auto mat = matrixgen::adjmat<Matrix_t>({4, 4, 1}, adjfn, matrixgen::randweight());

  using DenseMatrix_t = Eigen::Matrix<Scalar_t, Eigen::Dynamic, Eigen::Dynamic>;
  std::cout << std::endl << DenseMatrix_t(mat) << std::endl;
}
