#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <gsl/gsl-lite.hpp>

#include <matrixgen/utility.hpp>

#include <algorithm>
#include <initializer_list>
#include <random>
#include <utility>

namespace matrixgen::implementation {

template <
  typename Matrix_t,
  typename InputIter_t
    >
struct Perturb {};

/**
 * Implementation of 'matrixgen::perturb' for 'Eigen::Matrix' objects. ??
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
struct Perturb<Eigen::Matrix<EigenScalar_t, ROWS, COLS, OPTIONS, MAXROWS, MAXCOLS>, InputIter_t>;

/**
 * Implementation of 'matrixgen::perturb' for 'Eigen::SparseMatrix' objects.
 *
 * NOTE: Returned sparse matrices are not compressed. Call method
 *       'makeCompressed' on the object to compress it.
 */
template <
  typename Scalar_t,
  int ALIGNMENT,
  typename Index_t,
  typename InputIter_t
    >
struct Perturb<Eigen::SparseMatrix<Scalar_t, ALIGNMENT, Index_t>, InputIter_t> {

  using Matrix_t = Eigen::SparseMatrix<Scalar_t, ALIGNMENT, Index_t>;

  static
  Matrix_t
  perturb(
      const Matrix_t& matrix,
      InputIter_t outerIndicesFirst,
      InputIter_t outerIndicesLast,
      uint64_t seed) {

    Expects(std::distance(outerIndicesFirst, outerIndicesLast) >= 0);
    Expects(std::all_of(outerIndicesFirst, outerIndicesLast, [matrix] (auto index) {
          return 0 <= index && index < matrix.outerSize();}));

    if constexpr (ALIGNMENT == Eigen::RowMajor) {

      Matrix_t result = matrix;
      if (!result.isCompressed()) {
        result.makeCompressed(); // Might be redundant; Copy-ctor seems to create a compressed matrix
      }

      auto engine = std::default_random_engine(seed);

      // Generate all possible inner indices. We shuffle this container and
      // draw the first `outerSize' to randomize the inner indices as we cannot
      // draw random numbers due to possible duplicates.
      auto innerIndices = std::vector<Index_t> {};
      innerIndices.resize(matrix.innerSize());
      std::generate(innerIndices.begin(), innerIndices.end(), [ii = (Index_t)0] () mutable { return ii++; });
      auto innerIndicesIter = innerIndices.begin();

      // Distribution from which values will be drawn
      auto dist = std::uniform_real_distribution<Scalar_t>((Scalar_t)1.0, (Scalar_t)2.0);
      auto randval = [&]() { return dist(engine); };

      for (auto it = outerIndicesFirst; it != outerIndicesLast; ++it) {
        const auto outerIndex = *it;
        const auto outerOffset = *std::next(result.outerIndexPtr(), outerIndex); // Offset into V and CI for this row
        const auto nnzInOuter = *std::next(result.outerIndexPtr(), outerIndex + 1)
                                - *std::next(result.outerIndexPtr(), outerIndex);

        // Randomize inner indices. We write out the inner indices in ascending
        // order which makes the output more predictable.
        if (std::distance(innerIndicesIter, innerIndices.end()) <= nnzInOuter) {
          std::shuffle(innerIndices.begin(), innerIndices.end(), engine);
          innerIndicesIter = innerIndices.begin();
        }
        std::sort(innerIndicesIter, std::next(innerIndicesIter, nnzInOuter));
        std::copy_n(innerIndicesIter, nnzInOuter, result.innerIndexPtr() + outerOffset);
        std::advance(innerIndicesIter, nnzInOuter);

        // Randomize values.
        std::generate_n(result.valuePtr() + outerOffset, nnzInOuter, randval);
      }
      return result;
    } else {
      static_assert(!std::is_same<Scalar_t, Scalar_t>(), "ColMajor not yet implemented.");
    }
  }

};

} // namespace matrixgen::implementation

namespace matrixgen {

template <
  typename Matrix_t,
  typename ListElem_t
    >
Matrix_t perturb(
    const Matrix_t& matrix,
    std::initializer_list<ListElem_t> list,
    uint64_t seed) {

  using Iter_t = typename std::initializer_list<ListElem_t>::const_iterator;
  return implementation::Perturb<Matrix_t, decltype(list.begin())>::perturb(matrix, list.begin(), list.end(), seed);
}

/**
 * Perturb selected rows of a matrix.
 *
 * Same as above with a range of 0-indexed row numbers.
 */
template <
  typename Matrix_t,
  typename InputIter_t
    >
Matrix_t perturb(
    const Matrix_t& matrix,
    InputIter_t outerIndicesFirst,
    InputIter_t outerIndicesLast,
    uint64_t seed) {

  return implementation::Perturb<Matrix_t, InputIter_t>::perturb(matrix, outerIndicesFirst, outerIndicesLast, seed);
}

} // namespace matrixgen::implementation
