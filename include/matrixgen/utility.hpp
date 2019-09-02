#pragma once

#include <Eigen/Sparse>

#include <gsl/gsl_assert>

namespace matrixgen
{

/**
 * Return the number of nonzeros in the `ii`-th row (column) for row-major
 * (col-major) sparse matrix `mat`.
 */
template <
  typename Scalar_t,
  int ALIGNMENT,
  typename Index_t
    >
int32_t
num_of_nnz_in_outer(
    const Eigen::SparseMatrix<Scalar_t, ALIGNMENT, Index_t>& mat,
    int32_t ii) {

  Expects( mat.outerSize() > ii );
  Expects( ii >= 0 ) ;

  /**
   * Return the difference of the row-pointers (col-pointers) for row (col)
   * `ii+1` and `ii`, except for the final row (col), which does not have a
   * successor. Use the total count of nonzeros instead.
   */
  if(ii < mat.outerSize() - 1) {
    const auto b = *std::next(mat.outerIndexPtr(), ii + 1);
    const auto a = *std::next(mat.outerIndexPtr(), ii);
    return b - a;
  }
  return mat.nonZeros() - *std::next(mat.outerIndexPtr(), ii);
}


} // namespace matrixgen
