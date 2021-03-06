#pragma once

#include <Eigen/Sparse>

#include <gsl/gsl-lite.hpp>

#include <chrono>
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
  // TODO :Comment format
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

template <
  typename InMatrix_t,
  typename CoordIndex_t,
  typename Filler_t
    >
struct Insert;


/**
 * Specialization of `insert` for scalar fillers.
 *
 * Calls `Eigen::SparseMatrix::insert(row, col) = value`.
 */
template <
  typename Scalar_t,
  int ALIGNMENT,
  typename Index_t,
  typename CoordIndex_t,
  typename Filler_t
    >
struct Insert<Eigen::SparseMatrix<Scalar_t, ALIGNMENT, Index_t>, CoordIndex_t, Filler_t>
{
  using InMatrix_t = Eigen::SparseMatrix<Scalar_t, ALIGNMENT, Index_t>;

  static void invoke(
      InMatrix_t& smat,
      CoordIndex_t row,
      CoordIndex_t col,
      Filler_t object) {

    static_assert(std::is_convertible<Filler_t, Scalar_t>(), "Incompatible filler type.");
    smat.insert(row, col) = object;
  }
};

/**
 * Specialization of `insert` for fillers of type `Eigen::Matrix`
*/
template <
  typename Scalar_t,
  int ALIGNMENT,
  typename Index_t,
  typename CoordIndex_t,
  int NUMROWS,
  int NUMCOLS,
  int ALIGNMENT2,
  int MAXNUMROWS,
  int MAXNUMCOLS
    >
struct Insert<Eigen::SparseMatrix<Scalar_t, ALIGNMENT, Index_t>,
              CoordIndex_t,
              Eigen::Matrix<Scalar_t, NUMROWS, NUMCOLS, ALIGNMENT2, MAXNUMROWS, MAXNUMCOLS>>
{

  using InMatrix_t = Eigen::SparseMatrix<Scalar_t, ALIGNMENT, Index_t>;
  using DenseMatrix_t = Eigen::Matrix<Scalar_t, NUMROWS, NUMCOLS, ALIGNMENT2, MAXNUMROWS, MAXNUMCOLS>;

  static void invoke(
      InMatrix_t& smat,
      CoordIndex_t row,
      CoordIndex_t col,
      const DenseMatrix_t& dmat) {

    for(auto ii = 0; ii < dmat.rows(); ++ii) {
      for(auto jj = 0; jj < dmat.cols(); ++jj) {
        smat.insert(ii + row, jj + col) = dmat(ii, jj);
      }
    }
  }
};


/**
 * insert
 *
 * Insert dense Eigen::matrix objects into sparse Eigen::SparseMatrix objects.
 *
 * This function generalizes the functionality provided by the 
 * `Eigen::SparseMatrix::insert()` member function in that
 * it allows for dense matrices to be inserted using the same syntax.
 *
 * Assuming some 'Eigen::Sparsematrix smat' object calling
 * 'smat.insert(3, 4) = 1' will insert a single value '1' at coordinates
 * (3, 4). In contrast, calling 'insert(smat, 3, 4, dmat)' for some 2x3
 * Eigen::Matrix dense matrix object 'dmat' will yield the calls
 *
 * smat.insert(3, 4) = dmat(0, 0);
 * smat.insert(3, 5) = dmat(0, 1);
 * smat.insert(3, 6) = dmat(0, 2);
 * smat.insert(4, 4) = dmat(1, 0);
 * smat.insert(4, 5) = dmat(1, 1);
 * smat.insert(4, 6) = dmat(1, 2);
 */
template <
  typename InMatrix_t,   // matrix type
  typename CoordIndex_t, // type of the row and column indices
  typename Filler_t      // type of object that's inserted
    >
void insert(
    InMatrix_t& smat,
    CoordIndex_t row,
    CoordIndex_t col,
    const Filler_t& object) {

  Insert<InMatrix_t, CoordIndex_t, Filler_t>::invoke(smat, row, col, object);
}

// TODO: Unify into Coords3d_t<T>. No need to have two types here.
template <typename Scalar_t>
using Coords3d_t = std::array<Scalar_t, 3>;

/**
 * Element-wise addition of arrays.
 */
template <typename Scalar_t>
Coords3d_t<Scalar_t> 
operator+(
    const Coords3d_t<Scalar_t>& a,
    const Coords3d_t<Scalar_t>& b) {

  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

/**
 * Element-wise subtraction of arrays.
 */
template <typename Scalar_t>
Coords3d_t<Scalar_t> 
operator-(
    const Coords3d_t<Scalar_t>& a,
    const Coords3d_t<Scalar_t>& b) {

  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

/**
 * Geometric midpoint
 */
template <typename OutScalar_t, typename InScalar_t = void>
Coords3d_t<OutScalar_t>
midpoint(const Coords3d_t<InScalar_t>& a, const Coords3d_t<InScalar_t>& b) {
  static_assert(std::is_floating_point<OutScalar_t>(), "Output type must be real.");

  return Coords3d_t<OutScalar_t> {
    a[0] + (static_cast<OutScalar_t>((b[0] - a[0])) / 2),
    a[1] + (static_cast<OutScalar_t>((b[1] - a[1])) / 2),
    a[2] + (static_cast<OutScalar_t>((b[2] - a[2])) / 2)};
}

/**
 * Modulus which wraps around at 0. Scalar `b` is mapped into the interval
 * [0; `mod`) if `mod > 0` or (`mod`; 0] if `mod < 0`.
 *
 * Used to implement periodic boundary conditions.
 */
template <typename Scalar_t>
Scalar_t
mod(Scalar_t n, Scalar_t mod) {

  Expects( mod != 0 );
  static_assert(std::numeric_limits<Scalar_t>::is_signed);

  const auto m = n % mod;
  if (m >= 0) {
    return m;
  }

  return mod + m;
}

template <typename Index_t>
Coords3d_t<Index_t>
modplus(
    const Coords3d_t<Index_t>& a,
    const Coords3d_t<Index_t>& b,
    const Coords3d_t<Index_t>& modulus) {

  return {mod((a[0] + b[0]), modulus[0]),
          mod((a[1] + b[1]), modulus[1]),
          mod((a[2] + b[2]), modulus[2])};
}

/**
 * Check whether node at `coords` is an inner node with respect to an extent
 * `EXTENT`, i.e. whether it does not reside on the outer `EXTENT` layers.
 */
template <
  typename Index_t = int,
  uint32_t EXTENT = 1
    >
bool is_inner_node(
    const Coords3d_t<Index_t>& coords,
    const Coords3d_t<Index_t>& gridDimensions) {

  Expects( coords[0] >= 0                );
  Expects( coords[1] >= 0                );
  Expects( coords[2] >= 0                );
  Expects( coords[0] < gridDimensions[0] );
  Expects( coords[1] < gridDimensions[1] );
  Expects( coords[2] < gridDimensions[2] );

  return (coords[0] >= EXTENT && coords[0] < gridDimensions[0] - EXTENT &&
          coords[1] >= EXTENT && coords[1] < gridDimensions[1] - EXTENT &&
          coords[2] >= EXTENT && coords[2] < gridDimensions[2] - EXTENT);
}

/**
 * Return true if node at 'myCoords' is contained by the grid.
 */
template <typename Index_t = int>
bool
is_inside_grid(
  const Coords3d_t<Index_t>& coords,
  const Coords3d_t<Index_t>& gridDimensions) {

  Expects( 0 < gridDimensions[0] );
  Expects( 0 < gridDimensions[1] );
  Expects( 0 < gridDimensions[2] );

  return (0 <= coords[0] && coords[0] < gridDimensions[0] &&
          0 <= coords[1] && coords[1] < gridDimensions[1] &&
          0 <= coords[2] && coords[2] < gridDimensions[2]);
}

/**
 * Compile-time boundary conditions
 *
 * Used in the implementation of `matrixgen::stencil7p`.
 */
enum class BC {
  DIRICHLET,
  NEUMANN,
  PERIODIC
};

/**
 * Static symmetric 7p-stencil
 *
 * Used in the implementation of `matrixgen::stencil7p`.
 */
template <std::size_t SIZE, typename Index_t = int>
const std::array<Coords3d_t<Index_t>, SIZE> STENCIL;

template <typename Index_t>
const std::array<Coords3d_t<Index_t>, 7> STENCIL<7, Index_t> =
  {{{ 0,  0,  0},
    {-1,  0,  0}, {1, 0, 0},   // X .. don't change the order.
    { 0, -1,  0}, {0, 1, 0},   // Y
    { 0,  0, -1}, {0, 0, 1}}}; // Z

/**
 * Generates a uin64_t seed from the system time.
 */
template <
  typename T = void
    >
uint64_t seed_from_time() {
  std::chrono::system_clock::time_point tp = std::chrono::system_clock::now();
  std::chrono::system_clock::duration dtn = tp.time_since_epoch();
  return static_cast<uint64_t>(dtn.count());
}

} // namespace matrixgen
