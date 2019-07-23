/**
 * \file
 * \author Stanislaw Hüll
 */
#pragma once

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <gsl/gsl_assert>

#include <initializer_list>

namespace matrixgen::implementation {

template <
  typename OutMatrix_t
    >
struct Create {};

/**
 * \brief Create Eigen::Matrix objects
 *
 * Usage:
 */
template <
  typename EigenScalar_t,
  int32_t ROWS,
  int32_t COLS,
  int32_t OPTIONS,
  int32_t MAXROWS,
  int32_t MAXCOLS
    >
struct Create<Eigen::Matrix<EigenScalar_t, ROWS, COLS, OPTIONS, MAXROWS, MAXCOLS>> {

  using Matrix_t = Eigen::Matrix<EigenScalar_t, ROWS, COLS, OPTIONS, MAXROWS, MAXCOLS>;

  static
  Matrix_t
  create(uint32_t numRows, uint32_t numCols, const std::vector<EigenScalar_t>& list) {

    auto denseMat = Matrix_t(numRows, numCols);
    for(uint32_t row = 0; row < numRows; ++row) {
      for(uint32_t col = 0; col < numCols; ++col) {
         denseMat(row, col) = list[row * numCols + col];
      }
    }
    return denseMat;
  }
};
//  
//template <
//  typename EigenScalar_t,
//  int EigenAlignment,
//  typename EigenIndex_t
//    >
//struct Create<Eigen::SparseMatrix<EigenScalar_t, EigenAlignment, EigenIndex_t>> {
//
//  using Matrix_t = Eigen::SparseMatrix<EigenScalar_t, Eigen::RowMajor, EigenIndex_t>;
//
//  Matrix_t
//  create(uint32_t numRows, uint32_t numCols, std::initializer_list<EigenScalar_t> list) {
//
//    auto denseMat = Eigen::Matrix<EigenScalar_t, Eigen::Dynamic, Eigen::Dynamic>(numRows, numCols);
//    for(uint32_t row = 0; row < numRows; ++row) {
//      for(uint32_t col = 0; col < numCols; ++col) {
//         
//      }
//    }
//  }

} // namespace matrixgen::implementation

namespace matrixgen {

/**
 * Create matrices using a simple syntax:
 *
 * auto myMatrix = matrixgen::create<MyMatrixType>(
 *     numRows, NumCols,
 *     { 3.14, 2.71, .. dense list of all elements .., 0.71 });
 *
 * This function is a simple dispatcher for the chosen matrix type. See the
 * specialization for details.
 */
template <
  typename OutMatrix_t,
  typename Scalar_t = double
    >
OutMatrix_t create(
    uint32_t numRows,
    uint32_t numCols,
    const std::vector<Scalar_t>& list) {

  Expects( numRows * numCols == list.size() );

  return implementation::Create<OutMatrix_t>::create(numRows, numCols, list);
}

}
