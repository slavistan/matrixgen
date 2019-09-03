#pragma once

#include <Eigen/Sparse>

#include <gsl/gsl_assert>

#include <execution>
#include <numeric>

namespace matrixgen
{

/**
 * pi
 *
 * Compile-time constant for Pi.
 */
template <
  typename Scalar_t = double
    >
constexpr
Scalar_t
pi() {

  return std::atan(1) * 4;
}

/**
 * num_of_nnz_in_outer
 *
 * Returns the number of nonzeros in the `ii`-th row (column) for row-major
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


/**
 * central_moving_sum
 *
 * Computes the central moving sum of radius 'radius' over a range of
 * elements [first, last) and writes the outputs to the range beginning at
 * 'outFirst'. Border elements are calculated using the available neighbors
 * only (e.g. using radius 2 the value for element at index 1 is the sum of
 * elements {0, 1, 2, 3}, whereas the value for the element at index 2 is the
 * accumulate of the elements {0, 1, 2, 3, 4}.
 *
 * Can be used in-place.
 */
template <
  typename RandAccIter_t,
  typename OutputIter_t
    >
void
central_moving_sum(
    RandAccIter_t first,
    RandAccIter_t last,
    OutputIter_t outFirst,
    int32_t radius) {

  const auto numOfElements = std::distance(first, last);

  Expects( numOfElements > 0);
  Expects( radius >= 0 );

  /**
   * As opposed to a naive implementation requiring the summation of
   * `2 * radius + 1` elements per single input value a much more
   * performant evaluation can be achieved using the exclusive scan `X` of the
   * inputs. The CMS for a radius 'r' can be computed using the relation
   *
   *  CMS_i = a_{i-r} + ... + a_{i + r} = X_{i+r+1} - X_{i-r}
   *
   * Thus the evaluation is performed by
   * (1) computing the scan of the input values
   * (2) taking the difference of the scan's elements which lie one radius
   *     to the right and left of the element in question, requiring a single
   *     addition per input value.
   */
  using Input_t = typename std::iterator_traits<RandAccIter_t>::value_type;

  // Generate the exclusive scan
  // TODO: Can be parallelized
  std::vector<Input_t> xscan(numOfElements + 1);
  std::exclusive_scan(std::execution::seq, first, last, std::begin(xscan), static_cast<Input_t>(0));
  xscan.back() = xscan[xscan.size() - 2] + *std::prev(last);

  // Compute the moving sum via the elements to the left and right. Stop at
  // the array's boundaries.
  // TODO: Can be parallelized
  for(auto ii = 0; ii < numOfElements; ++ii) {
    const int64_t left = std::max<int64_t>(0, ii - radius);
    const int64_t right = std::min<int64_t>(ii + radius + 1, xscan.size() - 1);
    *std::next(outFirst, ii) = *std::next(std::cbegin(xscan), right) - *std::next(std::cbegin(xscan), left);
  };
}

/**
 * closed_loop_moving_mean
 *
 * TODO: Clarify
 *
 * Computes the central moving mean at radius 'radius' over a range of
 * elements [first, last) and writes the outputs to the range pointed by
 * 'outFirst'. Border elements of the input range are treated according to
 * 'central_moving_sum'.
 *
 * The closed-loop moving mean requires the specification of the closed loop
 * by means of 'loopMin' and 'loopMax'. The numeric value will roll over to
 * 'loopMin' if it is greater than or equal to 'loopMax' while values less
 * than 'loopMin' will roll over to 'loopMax'.
 *
 * Closed-loop moving mean is computed by vectorial addition of the values'
 * unit vectors of a window which correspond to the unit vectors of the sphere
 * whose angle spans [loopMin, loop_mmax) and then normalizing the input.
 *
 * Can be used in-place.
 */
template <
  typename InputIter_t,
  typename OutputIter_t
    >
void
closed_loop_moving_mean(
    InputIter_t first,
    InputIter_t last,
    OutputIter_t outFirst,
    typename std::iterator_traits<InputIter_t>::value_type loopMin,
    typename std::iterator_traits<InputIter_t>::value_type loopMax,
    std::size_t radius) {

  const auto numOfElements = std::distance(first, last);

  Expects( numOfElements > 0 );
  Expects( std::all_of(first, last, [loopMin, loopMax](auto x) {return (loopMin <= x) && (x <= loopMax);}) );

  //
  // Working mechanism
  //
  // (1) Create angles from inputs
  // (2) Perform "moving-sum" vector addition on unit vectors corresponding to angles
  // (3) Convert results of vector addition back into angles
  // (4) Map angles back into user-specified input interval
  //

  // (1.1) Create phase angles in [0; 2*pi) from inputs
  auto angles = std::vector<double>(numOfElements);
  const auto rangeWidth = loopMax - loopMin;
  std::transform(
      first,
      last,
      std::begin(angles), [rangeWidth, loopMin](auto x) {
        return 2 * pi() * (x - loopMin) / rangeWidth;
      });

  // (1.2) Generate complex doubles from phase angles to perform vector arithmetic on
  auto cplx = std::vector<std::complex<double>>(angles.size());
  std::transform(
      angles.cbegin(),
      angles.cend(),
      cplx.begin(),
      [](double d) {
        return std::polar(1.0, d);
      });

  // (2) Compute central moving mean
  //
  // In order to efficiently compute the central moving sum for large radii the evaluation scheme is optimized
  // (see central_movin_sum) and requires the summation of potentially very many elements. In order to avoid
  // the accumulation of rounding errors the summation is performed on integers. It is assumed that the rounding errors
  // of the accumulation of double values outweigh the errors obtained from converting the double values to scaled
  // intergers.
  //
  // TODO: Investigate mechanism in case of overflow
  //
  // (2.1) Set up integer complex
  const auto scale = static_cast<double>(std::pow(2, 50));
  auto cplxInt = std::vector<std::complex<int64_t>>(cplx.size());
  std::transform(
      std::cbegin(cplx),
      std::cend(cplx),
      std::begin(cplxInt),
      [scale](const std::complex<double>& dcplx) {
        return static_cast<std::complex<int64_t>>(scale * dcplx);
      });

  // (2.2) Compute central moving sum (vector addition)
  auto smoothedCplxInt = std::vector<std::complex<int64_t>>(cplxInt.size());
  central_moving_sum(std::cbegin(cplxInt), std::cend(cplxInt), std::begin(smoothedCplxInt), radius);

  // (3) Revert to phase angles
  std::transform(
      std::cbegin(smoothedCplxInt),
      std::cend(smoothedCplxInt),
      std::begin(angles),
      [](std::complex<int64_t> icplx) {
        if(icplx.real() == 0 && icplx.imag() == 0) {
          return static_cast<double>(0); // If, by chance, 0 + 0i is generated return angle '0'.
        }
        const std::complex<double> dcplx(icplx.real(), icplx.imag());
        double arg = std::arg(dcplx);
        if(arg <= 0) {
          arg += 2 * pi();
        }
        return arg;
      });

  // (4) Map smoothed angles into input range
  std::transform(
      std::cbegin(angles),
      std::cend(angles),
      outFirst,
      [rangeWidth, loopMin](double ang) {
        return loopMin + rangeWidth * (ang)/(2 * pi());
      });
}

/**
 * darts_sampling
 *
 * TODO: Documentation
 */
template <
  typename InputIter1_t,
  typename InputIter2_t,
  typename OutputIter_t
    >
void
darts_sampling(
    InputIter1_t quotaFirst,
    InputIter1_t quotaLast,
    InputIter2_t bulletsFirst,
    InputIter2_t bulletsLast,
    OutputIter_t outFirst) {

  // (1) Create target bins (set up 'dartboard')
  std::vector<double> ratios(std::distance(quotaFirst, quotaLast));
  const double sum = std::accumulate(quotaFirst, quotaLast, static_cast<double>(0));
  std::transform(quotaFirst, quotaLast, std::begin(ratios), [sum](auto val) -> double {return val/sum;});

  // (2) Generate indices ('Evaluate hits')
  std::vector<double> ratiosIncScan(ratios.size());
  std::inclusive_scan(std::execution::seq, std::cbegin(ratios), std::cend(ratios), std::begin(ratiosIncScan));
  const auto nbullet = std::distance(bulletsFirst, bulletsLast);
  for(auto loopIndex = 0; loopIndex < nbullet; ++loopIndex) {
    auto it = std::next(bulletsFirst, loopIndex);
    auto resultIter = std::lower_bound(std::cbegin(ratiosIncScan), std::cend(ratiosIncScan), *it);
    auto index = std::distance(std::cbegin(ratiosIncScan), resultIter);

    auto offset = std::distance(bulletsFirst, it);
    *std::next(outFirst, offset) = index;
  };
}

} // namespace matrixgen
