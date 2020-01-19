#pragma once

#include <matrixgen/assemble.hpp>

#include <iterator>
#include <random>

#include <Eigen/Sparse>

#include <gsl/gsl-lite.hpp>

namespace matrixgen::implementation {

template <
  typename OutMatrix_t,
  typename InMatrixIter_t,
  typename PropIter_t
    >
struct Interleave;

/**
 * Specialization of `interleave` for `Eigen::SparseMatrix<>` return types.
 */
template <
  typename Scalar_t,
  int ALIGNMENT,
  typename Index_t,
  typename InMatrixIter_t,
  typename PropIter_t
    >
struct Interleave<Eigen::SparseMatrix<Scalar_t, ALIGNMENT, Index_t>, InMatrixIter_t, PropIter_t>
{

  using OutMatrix_t = Eigen::SparseMatrix<Scalar_t, ALIGNMENT, Index_t>;

  static
  OutMatrix_t
  invoke(
      InMatrixIter_t matFirst,
      InMatrixIter_t matLast,
      PropIter_t propFirst,
      PropIter_t propLast,
      int32_t coupling,
      int64_t seed) {

  const auto numOfMatrices = std::distance(matFirst, matLast);
  const auto numOfProportions = std::distance(propFirst, propLast);

  // Require one proportion per matrix and that all matrices have the same outer size.
  Expects( numOfMatrices > 0);
  Expects( numOfMatrices == numOfProportions );
  Expects( std::all_of(matFirst, matLast, [&matFirst](auto mat) {return matFirst->outerSize() == mat.outerSize();}) );

  std::size_t outerSize = matFirst->outerSize();

  //
  // (1) Generate indices
  //
  // Generate random uniformly distributed numbers in [0, 1]
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> distribution(0, 1);
  std::vector<double> runif(outerSize);
  std::transform(runif.begin(), runif.end(), runif.begin(), // TODO: Replace by `std::generate`
    [&runif, &distribution, &generator](double)
    {
      return distribution(generator);
    }
  );

  // Apply closed-loop moving mean to runifs
  closed_loop_moving_mean(runif.begin(), runif.end(), runif.begin(), 0, 1, coupling);

  // Generate indices from runifs and proportions
  std::vector<std::size_t> indices(outerSize);
  darts_sampling(propFirst, propLast, runif.begin(), runif.end(), indices.begin());

  // (2) Contruct matrix from indexed rows
  return assemble(matFirst, matLast, indices.cbegin(), indices.cend());
}
};
}

namespace matrixgen {

template <
  typename InMatrixIter_t,
  typename OutMatrix_t = typename std::iterator_traits<InMatrixIter_t>::value_type,
  typename PropIter_t = void
    >
typename std::iterator_traits<InMatrixIter_t>::value_type
interleave(
    InMatrixIter_t matrixFirst,
    InMatrixIter_t matrixLast,
    PropIter_t propFirst,
    PropIter_t propLast,
    int32_t coupling = 0,
    int64_t seed = 42) {

  using InMatrix_t = typename std::iterator_traits<InMatrixIter_t>::value_type;
  static_assert(std::is_same<OutMatrix_t, InMatrix_t>(),
      "`interleave` does not yet support converting between matrix types. Input "
      "type must match output type.");

  return implementation::Interleave<OutMatrix_t, InMatrixIter_t, PropIter_t>::
          invoke(matrixFirst, matrixLast, propFirst, propLast, coupling, seed);
}
} // namespace matrixgen
