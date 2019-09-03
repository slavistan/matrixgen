#pragma once

#include <matrixgen/utility.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>

#include <Eigen/Sparse>

#include <gsl/gsl_assert>

namespace matrixgen::implementation
{

template <
  typename OutMatrix_t,
  typename InMatrixIter_t,
  typename IndexIter_t
    >
struct Assemble;

/**
 * Specialization of `assemble` for `Eigen::SparseMatrix<>` return types.
 */
template <
  typename Scalar_t,
  int ALIGNMENT,
  typename Index_t,
  typename InMatrixIter_t,
  typename IndexIter_t
    >
struct Assemble<Eigen::SparseMatrix<Scalar_t, ALIGNMENT, Index_t>, InMatrixIter_t, IndexIter_t>
{

  using OutMatrix_t = Eigen::SparseMatrix<Scalar_t, ALIGNMENT, Index_t>;

  static
  OutMatrix_t
  invoke(
      InMatrixIter_t matrixFirst,
      InMatrixIter_t matrixLast,
      IndexIter_t indexFirst,
      IndexIter_t indexLast) {

    const auto numOfMatrices = std::distance(matrixFirst, matrixLast);
    const auto numOfIndices = std::distance(indexFirst, indexLast);

    Expects( std::all_of(indexFirst, indexLast,
                [numOfMatrices](auto idx) { return (0 <= idx && idx < numOfMatrices);}) );

    if(numOfMatrices == 0) {
      return OutMatrix_t {};
    }

    // Output matrix's outer size is equal to the number of indices.
    const auto targetMatrixOuterSize = std::distance(indexFirst, indexLast);

    // Infer target matrix's inner size from the source matrices' maximum inner
    // size
    const auto whichMatrix = *std::max_element(matrixFirst, matrixLast,
        [](const auto& a, const auto& b) {
          return a.innerSize() < b.innerSize();});
    const auto targetMatrixInnerSize = whichMatrix.innerSize();

    // Construct a vector of the number of nonzeros for each outer in the
    // output matrix. Used to allocate memory later.
    auto nonzerosInOuter = std::vector<uint32_t>(targetMatrixOuterSize);
    for(auto ii = 0; ii < numOfIndices; ++ii) {
      const auto pMatrix = std::next(matrixFirst, *std::next(indexFirst, ii));
      nonzerosInOuter[ii] = num_of_nnz_in_outer(*pMatrix, ii);
    }

    // Allocate storage space for target matrix.
    OutMatrix_t result;
    if constexpr(ALIGNMENT == Eigen::RowMajor) {
      result.resize(targetMatrixOuterSize, targetMatrixInnerSize);
    }
    else {
      result.resize(targetMatrixInnerSize, targetMatrixOuterSize);
    }
    result.reserve(nonzerosInOuter);

    // Copy outers from source matrices into target matrix.
    // TODO: Improve the implementation. We're just copying segments of memory
    //       from the source matrices into the output matrix. No need for so
    //       many loops and calls to insert.
    for(auto ii = 0; ii < targetMatrixOuterSize; ++ii) {
      const auto pMatrix = std::next(matrixFirst, *std::next(indexFirst, ii));
      const auto innerIndexStart = *std::next(pMatrix->outerIndexPtr(), ii);
      const auto innerIndexSize = nonzerosInOuter[ii];
      for(auto jj = 0; jj < innerIndexSize; ++jj) {
        const auto val = *std::next(pMatrix->valuePtr(), innerIndexStart + jj);
        Index_t innerIndex = *std::next(pMatrix->innerIndexPtr(), innerIndexStart + jj);
        if constexpr(ALIGNMENT == Eigen::RowMajor) {
          result.insert(ii, innerIndex) = val;
        }
        else {
          result.insert(innerIndex, ii) = val;
        }
      }
    };
    return result;
  }
};

} /* asc::matrixgen::implementation */

namespace matrixgen
{
/*
 *  assemble
 *
 *  Returns a new row-major (col-major) matrix built from the row-major
 *  (col-major) matrices in the input range [matrix_first, matrix_last) by
 *  choosing individual rows (columns) according to the indices in the index
 *  range [index_fist, indexLast). The returned matrix's k-th row (column)
 *  is equal to the k-th row (column) of the input matrix pointed by the k-th
 *  index.
 *
 *  The returned matrix is uncompressed. Its height (width) is the total
 *  number of indices, whereas its width (height) is equal to the input
 *  matrices' maximum width (height). Rows (columns) which were drawn from
 *  a matrix whose width (height) is less than the resulting matrix's width
 *  (height) are tail-padded with 0.
 *
 *  Example:
 *  ~~~~~~~~
 *  Consider two matrices A, B (whose elements may or may not be zero)
 *
 *      (a11, a12, a13)
 *      (a21, a22, a23)      (b11, b12)
 *  A = (a31, a32, a23)  B = (b21, b22)
 *      (a41, a42, a33)      (b31, b32)
 *      (a51, a52, a53)
 *
 *  Assuming row-major layout the indices [0,1,1,0,0] would yield the following result C
 *
 *      (a11, a12, a13)
 *      (b21, b22,   0)
 *      (b31, b32,   0)
 *  C = (a41, a42, a43)
 *      (a51, a52, a53)
 */
template <
  typename InMatrixIter_t,
  typename OutMatrix_t = typename std::iterator_traits<InMatrixIter_t>::value_type,
  typename IndexIter_t = void
    >
OutMatrix_t
assemble(
    InMatrixIter_t matFirst, // range over matrices
    InMatrixIter_t matLast,
    IndexIter_t indexFirst,   // range over indices
    IndexIter_t indexLast) {

  using InMatrix_t = typename std::iterator_traits<InMatrixIter_t>::value_type;
  static_assert(std::is_same<OutMatrix_t, InMatrix_t>(),
      "`assemble` does not yet support converting between matrix types. Input "
      "type must match output type.");

  return implementation::Assemble<OutMatrix_t, InMatrixIter_t, IndexIter_t>::
          invoke(matFirst, matLast, indexFirst, indexLast);
}

} // namespace matrixgen
