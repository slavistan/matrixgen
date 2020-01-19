#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <gsl/gsl-lite.hpp>

#include <initializer_list>
#include <utility>

namespace matrixgen::implementation {

template <
  typename OutMatrix_t,
  typename InputIter_t
    >
struct Create {};

/**
 * Implementation of 'matrixgen::create' for 'Eigen::Matrix' objects.
 */
template <
  typename EigenScalar_t,
  int32_t ROWS,
  int32_t COLS,
  int32_t OPTIONS,
  int32_t MAXROWS,
  int32_t MAXCOLS,
  typename InputIter_t
    >
struct Create<Eigen::Matrix<EigenScalar_t, ROWS, COLS, OPTIONS, MAXROWS, MAXCOLS>, InputIter_t> {

  using Matrix_t = Eigen::Matrix<EigenScalar_t, ROWS, COLS, OPTIONS, MAXROWS, MAXCOLS>;

  static
  Matrix_t
  create(uint32_t numRows, uint32_t numCols, InputIter_t first, InputIter_t last) {

    Expects( numRows * numCols == std::distance(first, last) );

    auto denseMat = Matrix_t(numRows, numCols);
    for(uint32_t row = 0; row < numRows; ++row) {
      for(uint32_t col = 0; col < numCols; ++col) {
         denseMat(row, col) = static_cast<EigenScalar_t>(*std::next(first, row * numCols + col));
      }
    }
    return denseMat;
  }
};

/**
 * Implementation of 'matrixgen::create' for 'Eigen::SparseMatrix' objects.
 *
 * NOTE: Returned sparse matrices are not compressed. Call method
 *       'makeCompressed' on the object to compress it.
 */
template <
  typename EigenScalar_t,
  int ALIGNMENT,
  typename EigenIndex_t,
  typename InputIter_t
    >
struct Create<Eigen::SparseMatrix<EigenScalar_t, ALIGNMENT, EigenIndex_t>, InputIter_t> {

  using Matrix_t = Eigen::SparseMatrix<EigenScalar_t, ALIGNMENT, EigenIndex_t>;

  static
  Matrix_t
  create(uint32_t numRows, uint32_t numCols, InputIter_t first, InputIter_t last) {

    Expects( numRows * numCols == std::distance(first, last) );

    // Create a dense matrix, feed it the values and from it create a sparse
    // view. This is a bogus implementation suitable for small matrices only.
    // If performance ever becomes an issue here's why.
    auto denseMat = Eigen::Matrix<
        EigenScalar_t,
        Eigen::Dynamic,
        Eigen::Dynamic,
        ALIGNMENT>(numRows, numCols);

    for(uint32_t row = 0; row < numRows; ++row) {
      for(uint32_t col = 0; col < numCols; ++col) {
         denseMat(row, col) = static_cast<EigenScalar_t>(*(std::next(first, row * numCols + col)));
      }
    }

    Matrix_t outmatrix = denseMat.sparseView();
    return outmatrix;
  }
};

} // namespace matrixgen::implementation

namespace matrixgen {

/**
 * Create matrices of various types using a simple syntax:
 *
 * ****************************************************************************
 * #include <matrixgen/core>
 *
 * int main() {
 *   using MyMatrixType_t = Eigen::SparseMatrix<double>;
 *   auto myMatrix = matrixgen::create<MyMatrixType_t>(3, 2,
 *       { 3.14, 2.71 ,
 *         4.01, 0.00 ,
 *         0.00, 1.01 });
 *   // Use 'myMatrix'.
 * }
 * ****************************************************************************
 *
 * The list element's type shall be convertible to the scalar type. Elements are
 * fed into the matrix in row-major order (row-by-row starting from the top-left
 * element) irrespective of the output matrix's data layout.
 *
 * This function is a simple dispatcher for the chosen matrix type. See the
 * specializations above for available matrix types.
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

  using Iter_t = typename std::initializer_list<ListElem_t>::const_iterator;
  return implementation::Create<OutMatrix_t, Iter_t>::create(numRows, numCols, list.begin(), list.end());
}

/**
 * As above for a range of elements.
 */
template <
  typename OutMatrix_t,
  typename InputIter_t
    >
OutMatrix_t create(
    uint32_t numRows,
    uint32_t numCols,
    InputIter_t first,
    InputIter_t last) {

  return implementation::Create<OutMatrix_t, InputIter_t>::create(numRows, numCols, first, last);
}

} // namespace matrixgen::implementation
