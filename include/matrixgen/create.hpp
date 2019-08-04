/**
 * \file
 * \author Stanislaw Hüll
 */
#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <gsl/gsl_assert>

#include <initializer_list>

namespace matrixgen::implementation {

template <
  typename OutMatrix_t,
  typename ListElem_t
    >
struct Create {};

/**
 * \brief Create Eigen::Matrix objects
 *
 * Usage:
 *
 * using DenseMatrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
 * auto myMatrix = matrixgen::create<DenseMatrix_t>(3, 2,
 *     { 3.14,   0,
 *          0, 1.1,
 *        9.2,   0 });
 *
 * std::cout << myMatrix << std::endl;
 */
template <
  typename EigenScalar_t,
  int32_t ROWS,
  int32_t COLS,
  int32_t OPTIONS,
  int32_t MAXROWS,
  int32_t MAXCOLS,
  typename ListElem_t
    >
struct Create<Eigen::Matrix<EigenScalar_t, ROWS, COLS, OPTIONS, MAXROWS, MAXCOLS>, ListElem_t> {

  using Matrix_t = Eigen::Matrix<EigenScalar_t, ROWS, COLS, OPTIONS, MAXROWS, MAXCOLS>;

  static
  Matrix_t
  create(uint32_t numRows, uint32_t numCols, std::initializer_list<ListElem_t> list) {

    auto denseMat = Matrix_t(numRows, numCols);
    for(uint32_t row = 0; row < numRows; ++row) {
      for(uint32_t col = 0; col < numCols; ++col) {
         denseMat(row, col) = static_cast<ListElem_t>(*std::next(list.begin(), row * numCols + col));
      }
    }
    return denseMat;
  }
};

/**
 * \brief Create SparseEigen::Matrix objects
 *
 * Usage:
 *
 * using SparseMatrix_t = Eigen::SparseMatrix<double, Eigen::RowMajor>;
 * auto myMatrix = matrixgen::create<SparseMatrix_t>(3, 2,
 *     { 3.14,   0,
 *          0, 1.1,
 *        9.2,   0 });
 *
 * std::cout << myMatrix << std::endl;
 *
 * Returned sparse matrices are not compressed.
 */
template <
  typename EigenScalar_t,
  int ALIGNMENT,
  typename EigenIndex_t,
  typename ListElem_t
    >
struct Create<Eigen::SparseMatrix<EigenScalar_t, ALIGNMENT, EigenIndex_t>, ListElem_t> {

  using Matrix_t = Eigen::SparseMatrix<EigenScalar_t, ALIGNMENT, EigenIndex_t>;

  static
  Matrix_t
  create(uint32_t numRows, uint32_t numCols, std::initializer_list<ListElem_t> list) {

    // Create a dense matrix, feed it the values and from it create a sparse view.
    auto denseMat = Eigen::Matrix<
        EigenScalar_t,
        Eigen::Dynamic,
        Eigen::Dynamic,
        ALIGNMENT>(numRows, numCols);

    for(uint32_t row = 0; row < numRows; ++row) {
      for(uint32_t col = 0; col < numCols; ++col) {
         denseMat(row, col) = static_cast<EigenScalar_t>(*(std::next(list.begin(), row * numCols + col)));
      }
    }

    Matrix_t outmatrix = denseMat.sparseView();
    return outmatrix;
  }
};

} // namespace matrixgen::implementation

namespace matrixgen {

/**
 * Create matrices using a simple syntax:
 *
 * auto myMatrix = matrixgen::create<MyMatrixType>(
 *     numRows, NumCols,
 *     { 3.14, 2.71, .. dense list of all elements .., 0.71 });
 *
 * The list type must match the matrix's scalar type. Elements are parsed in
 * row-major order (row-by-row starting from the top-left element) irrespective
 * of the output matrix's data layout.
 *
 * This function is a simple dispatcher for the chosen matrix type. See the
 * specialization for details.
 *
 */
template <
  typename OutMatrix_t,
  typename ListElem_t
    >
OutMatrix_t create(
    uint32_t numRows,
    uint32_t numCols,
    std::initializer_list<ListElem_t> list) {

  Expects( numRows * numCols == list.size() );

  return implementation::Create<OutMatrix_t, ListElem_t>::create(numRows, numCols, list);
}

}
