#pragma once

#include <gsl/gsl_assert>

#include <array>
#include <limits>
#include <utility>
#include <vector>

namespace matrixgen {

/**
 * Collection of predefined stencil function.
 */

enum class BC {
  DIRICHLET,
  NEUMANN,
  PERIODIC
};

template <typename Index_t = int>
using DiscreteCoords3d_t = std::array<Index_t, 3>;

/**
 * Element-wise addition of arrays.
 */
template <typename Index_t>
DiscreteCoords3d_t<Index_t> 
operator+(
    const DiscreteCoords3d_t<Index_t>& a,
    const DiscreteCoords3d_t<Index_t>& b) {

  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

/**
 * Element-wise addition of arrays.
 */
template <typename Index_t>
DiscreteCoords3d_t<Index_t> 
operator-(
    const DiscreteCoords3d_t<Index_t>& a,
    const DiscreteCoords3d_t<Index_t>& b) {

  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
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
DiscreteCoords3d_t<Index_t>
modplus(
    const DiscreteCoords3d_t<Index_t>& a,
    const DiscreteCoords3d_t<Index_t>& b,
    const DiscreteCoords3d_t<Index_t>& modulus) {

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
    const DiscreteCoords3d_t<Index_t>& coords,
    const DiscreteCoords3d_t<Index_t>& gridDimensions) {

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
  const DiscreteCoords3d_t<Index_t>& coords,
  const DiscreteCoords3d_t<Index_t>& gridDimensions) {

  Expects( 0 < gridDimensions[0] );
  Expects( 0 < gridDimensions[1] );
  Expects( 0 < gridDimensions[2] );

  return (0 <= coords[0] && coords[0] < gridDimensions[0] &&
          0 <= coords[1] && coords[1] < gridDimensions[1] &&
          0 <= coords[2] && coords[2] < gridDimensions[2]);
}

/**
 * Symmetric 7p-stencil
 */
template <std::size_t SIZE, typename Index_t = int>
const std::array<DiscreteCoords3d_t<Index_t>, SIZE> STENCIL;

template <typename Index_t>
const std::array<DiscreteCoords3d_t<Index_t>, 7> STENCIL<7, Index_t> =
  {{{ 0,  0,  0},
    {-1,  0,  0}, {1, 0, 0},   // X .. don't change the order.
    { 0, -1,  0}, {0, 1, 0},   // Y
    { 0,  0, -1}, {0, 0, 1}}}; // Z

/**
 * Returns a lambda which implements the symmetric 7p stencil given a set of
 * boundary conditions.
 *
 * TODO: Implement Neumann boundary conditions.
 */
template <
  auto XBC = BC::DIRICHLET,
  auto YBC = BC::DIRICHLET,
  auto ZBC = BC::DIRICHLET,
  typename Index_t = int
    >
auto stencil7p() {

  return [offsets = std::array<DiscreteCoords3d_t<Index_t>, 7>{}] (
    const DiscreteCoords3d_t<Index_t>& coords,
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

/**
 * Collection of predefined weight functions.
 *
 * TODO: Is this the correct way to implement preset lambdas? Pass them as return
 *       "value"?
 */
/**
 * Weighfunction returning a constant.
 */
template <
  typename Scalar_t = double
    >
auto constweight(Scalar_t val = 1) {

  return ([val]() { return static_cast<Scalar_t>(val); });
}

} // namespace matrixgen
