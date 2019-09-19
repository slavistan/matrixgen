/**
 * Collection of predefined adjacency functions and weight functions.
 * Dependencies reside in `matrixgen/utility.hpp` for the sake of legibility.
 */
#pragma once

#include <matrixgen/adjmat.hpp>
#include <matrixgen/utility.hpp>

#include <gsl/gsl_assert>

#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <utility>
#include <vector>

namespace matrixgen {

/**
 * matrixgen::stencil7p()
 *
 * Returns a lambda which implements the symmetric 7p stencil. Boundary
 * conditions may be chosen independently for x, y and z dimensions
 * independent of each other.
 */
template <
  auto XBC = BC::DIRICHLET, // Boundary conditions for the X, ..
  auto YBC = BC::DIRICHLET, // .. the Y,
  auto ZBC = BC::DIRICHLET, // .. and the Z-dimension, respectively.
  typename Index_t = int
    >
auto stencil7p() {

  return [offsets = std::array<Coords3d_t<Index_t>, 7> {}] (
    const Coords3d_t<Index_t>& coords,
    const std::array<int, 3>& gridDimensions) mutable {

    Expects( coords[0] >= 0                );
    Expects( coords[1] >= 0                );
    Expects( coords[2] >= 0                );
    Expects( coords[0] < gridDimensions[0] );
    Expects( coords[1] < gridDimensions[1] );
    Expects( coords[2] < gridDimensions[2] );
    static_assert( std::is_same<decltype(XBC), BC>() );

    /**
     * Check if we're being called on an inner node in which case the node
     * exhibits its full adjacency pattern. This will apply to the overwhelming
     * majority of nodes for most matrices thus we assume this scenario right
     * away.
     */
    if (is_inner_node(coords, gridDimensions)) {
      return std::pair {STENCIL<7>.cbegin(), STENCIL<7>.cend()};
    }

    /**
     * For all the outer nodes we start out with the null-offset and add other
     * offsets according to our position in the grid and the chosen boundary
     * conditions.
     */
    offsets.front() = STENCIL<7>.front();     // add the null offset and ...
    auto end = std::next(offsets.begin(), 1); // .. keep track of our end ptr.
    if constexpr (XBC == BC::DIRICHLET) {
      /**
       * Add any offset in {(-1, 0, 0), (1, 0, 0)} which points inside the
       * grid and discard the rest. This implementation is repeated for the
       * x and y dimensions, respecitvely.
       */
      end = std::copy_if(
          std::next(STENCIL<7>.cbegin(), 1), // Range over the x-offsets in
          std::next(STENCIL<7>.cbegin(), 3), // `STENICIL<7>`. See declaration.
          end,
          [&coords, &gridDimensions](const auto& offset) {
            return is_inside_grid(coords + offset, gridDimensions);
          });
    }
    if constexpr (XBC == BC::PERIODIC) {
      /**
       * Periodic boundary conditions are implemented by performing a simple
       * (multi-dimensional) modulus addition of the offset to the node's
       * coordinates which maps to the other side of the grid, in case the
       * offset points outside the grid.
       */
      end = std::transform(
          std::next(STENCIL<7>.cbegin(), 1),
          std::next(STENCIL<7>.cbegin(), 3),
          end,
          [&coords, &gridDimensions](const auto& offset) {
            // The stencil is supposed to return offsets, thus we subtract
            // the coordinates after the modplus to obtain what we need.
            return modplus(coords, offset, gridDimensions) - coords;
          });
    }
    if constexpr (XBC == BC::NEUMANN) {
      static_assert(!std::is_same<Index_t, Index_t>(), "TODO: Not yet implemented.");
    }
    if constexpr (YBC == BC::DIRICHLET) {
      end = std::copy_if(
          std::next(STENCIL<7>.cbegin(), 3),
          std::next(STENCIL<7>.cbegin(), 5),
          end,
          [&coords, &gridDimensions](const auto& offset) {
            return is_inside_grid(coords + offset, gridDimensions);
          });
    }
    if constexpr (YBC == BC::PERIODIC) {
      end = std::transform(
          std::next(STENCIL<7>.cbegin(), 3),
          std::next(STENCIL<7>.cbegin(), 5),
          end,
          [&coords, &gridDimensions](const auto& offset) {
            return modplus(coords, offset, gridDimensions) - coords;
          });
    }
    if constexpr (YBC == BC::NEUMANN) {
      static_assert(!std::is_same<Index_t, Index_t>(), "TODO: Not yet implemented.");
    }
    if constexpr (ZBC == BC::DIRICHLET) {
      end = std::copy_if(
          std::next(STENCIL<7>.cbegin(), 5),
          std::next(STENCIL<7>.cbegin(), 7),
          end,
          [&coords, &gridDimensions](const auto& offset) {
            return is_inside_grid(coords + offset, gridDimensions);
          });
    }
    if constexpr (ZBC == BC::PERIODIC) {
      end = std::transform(
          std::next(STENCIL<7>.cbegin(), 5),
          std::next(STENCIL<7>.cbegin(), 7),
          end,
          [&coords, &gridDimensions](const auto& offset) {
            return modplus(coords, offset, gridDimensions) - coords;
          });
    }
    if constexpr (ZBC == BC::NEUMANN) {
      static_assert(!std::is_same<Index_t, Index_t>(), "TODO: Not yet implemented.");
    }

    using Iter_t = typename decltype(offsets)::const_iterator;
    return std::pair {offsets.cbegin(), static_cast<Iter_t>(end)};
  };
}

/*************************************
 ********* Weight functions **********
 *************************************/

 //TODO: Is this the correct way to implement preset lambdas? Pass them as return
 //      "value"? Do multiple calls return different lambdas or is the data 
 //      shared between call sites in the case of mutable lambdas (as in the
 //      randweight below)?

/**
 * matrixgen::constweight()
 *
 * Weighfunction returning a constant.
 */
template <
  typename Scalar_t = double
    >
auto constweight(Scalar_t val = 1) {

  return [val]() { return static_cast<Scalar_t>(val); };
}

/**
 * matrixgen::randweight()
 *
 * Weightfunction drawing values from a uniform real distribution over [0; 1]
 */
template <
  typename Scalar_t = double,
  typename Seed_t = uint64_t
    >
auto randweight(Seed_t seed = 1) {

  using Engine_t = std::default_random_engine;
  using Dist_t = std::uniform_real_distribution<Scalar_t>;

  return
    [engine = Engine_t(seed), dist = Dist_t(static_cast<Scalar_t>(0.0), static_cast<Scalar_t>(1.0))]
    () mutable {
      return static_cast<Scalar_t>(dist(engine));
    };
}

/**
 * matrixgen::sinusoid_add
 *
 * Additive sinusoidal.
 * TODO: Add a little text.
 */
template <
  typename Scalar_t = double,
  typename Index_t = int32_t
    >
auto sinusoid_add(Scalar_t nx, Scalar_t ny, Scalar_t nz) {
  return
    [nx=static_cast<Scalar_t>(nx), ny=static_cast<Scalar_t>(ny), nz=static_cast<Scalar_t>(nz)]
    (Coords3d_t<Index_t> coords, Coords3d_t<Index_t> neighborCoords, Coords3d_t<Index_t> gridDimensions) {
      const auto midpt = midpoint<Scalar_t>(coords, neighborCoords);
      const Scalar_t xrel = (midpt[0]) / gridDimensions[0];
      const Scalar_t yrel = (midpt[1]) / gridDimensions[1];
      const Scalar_t zrel = (midpt[2]) / gridDimensions[2];
      return ((std::sin(matrixgen::pi<Scalar_t>() * nx * xrel)) +
              (std::sin(matrixgen::pi<Scalar_t>() * ny * yrel)) +
              (std::sin(matrixgen::pi<Scalar_t>() * nz * zrel)) / 3);
    };
}

/**
 * matrixgen::sinusoid_add_bias
 *
 * Biased additive sinusoidal.
 */
template <
  typename Scalar_t = double,
  typename Index_t = int32_t
    >
auto sinusoid_add_bias(Scalar_t nx, Scalar_t ny, Scalar_t nz) {
  return
    [nx=static_cast<Scalar_t>(nx), ny=static_cast<Scalar_t>(ny), nz=static_cast<Scalar_t>(nz)]
    (Coords3d_t<Index_t> coords, Coords3d_t<Index_t> neighborCoords, Coords3d_t<Index_t> gridDimensions) {
      const auto midpt = midpoint<Scalar_t>(coords, neighborCoords);
      const Scalar_t xrel = (midpt[0]) / gridDimensions[0];
      const Scalar_t yrel = (midpt[1]) / gridDimensions[1];
      const Scalar_t zrel = (midpt[2]) / gridDimensions[2];
      return ((std::sin(matrixgen::pi<Scalar_t>() * nx * xrel)) +
              (std::sin(matrixgen::pi<Scalar_t>() * ny * yrel)) +
              (std::sin(matrixgen::pi<Scalar_t>() * nz * zrel)) / 3) + 1;
    };
}

/**
 * matrixgen::sinusoid_mul
 *
 * Multiplicative sinusoidal.
 */
template <
  typename Scalar_t = double,
  typename Index_t = int32_t
    >
auto sinusoid_mul(Scalar_t nx, Scalar_t ny, Scalar_t nz) {
  return
    [nx=static_cast<Scalar_t>(nx), ny=static_cast<Scalar_t>(ny), nz=static_cast<Scalar_t>(nz)]
    (Coords3d_t<Index_t> coords, Coords3d_t<Index_t> neighborCoords, Coords3d_t<Index_t> gridDimensions) {
      const auto midpt = midpoint<Scalar_t>(coords, neighborCoords);
      const Scalar_t xrel = (midpt[0]) / gridDimensions[0];
      const Scalar_t yrel = (midpt[1]) / gridDimensions[1];
      const Scalar_t zrel = (midpt[2]) / gridDimensions[2];
      return ((std::sin(matrixgen::pi<Scalar_t>() * nx * xrel)) *
              (std::sin(matrixgen::pi<Scalar_t>() * ny * yrel)) *
              (std::sin(matrixgen::pi<Scalar_t>() * nz * zrel)));
    };
}

/**
 * matrixgen::sinusoid_mul_bias
 *
 * Biased multiplicative sinusoidal.
 */
template <
  typename Scalar_t = double,
  typename Index_t = int32_t
    >
auto sinusoid_mul_bias(Scalar_t nx, Scalar_t ny, Scalar_t nz) {
  return
    [nx=static_cast<Scalar_t>(nx), ny=static_cast<Scalar_t>(ny), nz=static_cast<Scalar_t>(nz)]
    (Coords3d_t<Index_t> coords, Coords3d_t<Index_t> neighborCoords, Coords3d_t<Index_t> gridDimensions) {
      const auto midpt = midpoint<Scalar_t>(coords, neighborCoords);
      const Scalar_t xrel = (midpt[0]) / gridDimensions[0];
      const Scalar_t yrel = (midpt[1]) / gridDimensions[1];
      const Scalar_t zrel = (midpt[2]) / gridDimensions[2];
      return ((std::sin(matrixgen::pi<Scalar_t>() * nx * xrel)) *
              (std::sin(matrixgen::pi<Scalar_t>() * ny * yrel)) *
              (std::sin(matrixgen::pi<Scalar_t>() * nz * zrel))) + 1;
    };
}

/*************************************
 *********** Full Wrappers ***********
 *************************************/

/**
 * structured_grid_sinusoidal
 *
 * Returns a diagonally dominant matrix whose non-diagonal values are
 * determinated according to the biased multiplicative sinusoid.
 */
template <
  typename AdjFn_t,
  int ALIGNMENT = Eigen::RowMajor,
  typename Scalar_t = double,
  typename Index_t = int
    >
Eigen::SparseMatrix<Scalar_t, ALIGNMENT, Index_t>
structured_grid_sinusoidal(
    const Coords3d_t<Index_t>& gridDimensions,
    AdjFn_t adjfn,
    Scalar_t nx,
    Scalar_t ny,
    Scalar_t nz){

  // Generate the baseline matrix.
  using OutMatrix_t = Eigen::SparseMatrix<Scalar_t, ALIGNMENT, Index_t>;
  auto matrix = adjmat<OutMatrix_t>(gridDimensions, adjfn, matrixgen::sinusoid_add_bias(nx, ny, nz));

  // Make matrix diagonally dominant by setting the diagonal elements to their
  // row's sum over the other elements' absolute values.
  for(auto ii = 0; ii < matrix.rows(); ++ii) {
    auto begin = matrix.valuePtr() + (matrix.outerIndexPtr())[ii];
    auto end = std::next(begin, matrix.row(ii).nonZeros());
    Scalar_t rowSum = std::accumulate(
        begin,
        end,
        static_cast<Scalar_t>(1), [](auto sum, auto elem){ // Start acc. at 1 to be strictly ddom
          return sum + std::abs(elem);
        });
    matrix.coeffRef(ii, ii) = -(std::abs(rowSum) + 1);
  }

  return matrix;
}

} // namespace matrixgen

// TODO: Implement 19p & 27p stencils
