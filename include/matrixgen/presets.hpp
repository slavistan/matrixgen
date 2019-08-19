#pragma once

#include <gsl/gsl_assert>

#include <array>
#include <utility>
#include <vector>

namespace matrixgen {

/**
 * Collection of predefined stencil function.
 */

enum class BOUNDCOND {
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
  {{{  0,  0,  0},
    { -1,  0,  0}, {  0,  -1,  0}, {  0,  0, -1},
    {  1,  0,  0}, {  0,   1,  0}, {  0,  0,  1}}};

/**
 * Returns a lambda which implements the symmetric 7p stencil.
 */
template <
  auto BC = BOUNDCOND::DIRICHLET,
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
    static_assert( std::is_same<decltype(BC), BOUNDCOND>() );

    /**
     * Check if we're being called on an inner node in which case the node
     * a predefined static stencil. This will apply to the overwhelming
     * majority of nodes for most matrices irrespective of boundary
     * conditions.
     */
    if (is_inner_node(coords, gridDimensions)) {
      return std::pair {STENCIL<7>.cbegin(), STENCIL<7>.cend()};
    }

    /**
     * For the outer nodes modify the set of offsets according to our boundary
     * conditions.
     */
    std::copy(STENCIL<7>.cbegin(), STENCIL<7>.cend(), offsets.begin());
    if constexpr (BC == BOUNDCOND::DIRICHLET) {
      /**
       * For any outer node start out at the full 7p stencil and remove any
       * offsets which are not compatible with our node's position.
       * This is a generic implementation suitable for all kinds of grids.
       */
      const auto end = std::remove_if(offsets.begin(), offsets.end(),
                  // "Remove" any offsets which point outside the grid.
                  // `remove_if` swaps the bad elements to the end of the array.
                  [&coords, &gridDimensions](const auto& offset) {
                    const auto neighborCoords = coords + offset;
                    return !is_inside_grid(neighborCoords, gridDimensions);
                  });

      using Iter_t = typename decltype(offsets)::const_iterator;
      return std::pair {offsets.cbegin(), static_cast<Iter_t>(end)};
    }
    if constexpr (BC == BOUNDCOND::PERIODIC) {
      // TODO
    }
    if constexpr (BC == BOUNDCOND::NEUMANN) {
      // TODO
    }
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
