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

template <typename Index_t>
DiscreteCoords3d_t<Index_t> 
operator+(
    const DiscreteCoords3d_t<Index_t>& a,
    const DiscreteCoords3d_t<Index_t>& b) {

  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

/**
 * Check whether node at `coords` is an inner node with respect to a symmetric
 * 7-point stencil, i.e. whether it does not reside on any of the outermost
 * sides of the grid.
 */
template <
  typename Index_t = int
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

  if (coords[0] >= 1 && coords[0] < gridDimensions[0] - 1 &&
      coords[1] >= 1 && coords[1] < gridDimensions[1] - 1 &&
      coords[2] >= 1 && coords[2] < gridDimensions[2] - 1) {
    return true;

  }

  return false;
}

/**
 * Return true if node at 'myCoords' is contained inside the grid.
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

template <
  auto BoundCond = BOUNDCOND::DIRICHLET, // TODO: Implement other boundary conds
  typename Index_t = int
    >
auto stencil7p = [offsets = std::array<DiscreteCoords3d_t<Index_t>, 7>{}] (
    const DiscreteCoords3d_t<Index_t>& coords,
    const std::array<int, 3>& gridDimensions) mutable {

    Expects( coords[0] >= 0 && coords[0] < gridDimensions[0] );
    Expects( coords[1] >= 0 && coords[1] < gridDimensions[1] );
    Expects( coords[2] >= 0 && coords[2] < gridDimensions[2] );

    /**
     * Check if we're being called on an inner node. This will be the case
     * most of the times for any non-trivial matrix. Return iterators to
     * a predefined static stencil.
     */
    if (is_inner_node(coords, gridDimensions)) {
      return std::pair {STENCIL<7>.cbegin(), STENCIL<7>.cend()};
    }
    /**
     * For any outer node start out at the full 7p stencil and remove any
     * offsets which are not compatible with our node's position.
     *
     * TODO: std::array does not have a erase method. Thus we stuff a dummy
     *       vector with the full stencil and erase individual offsets from
     *       that vector. Finally we copy the content of the vector into
     *       out static `offsets` array.
     *       This is a bogus implementation. Work on the static array directly.
     *
     */
    auto offsets_dummy = std::vector (STENCIL<7>.cbegin(), STENCIL<7>.cend());
    for( const auto& offset: STENCIL<7> ){
      const auto neighborCoords = coords + offset;
       if(!is_inside_grid(neighborCoords, gridDimensions)) {
         const auto where = std::find(offsets_dummy.begin(), offsets_dummy.end(), offset);
         offsets_dummy.erase(where);
       }
    }
    std::copy(offsets_dummy.begin(), offsets_dummy.end(), offsets.begin());
    return std::pair {offsets.cbegin(), offsets.cbegin() + offsets_dummy.size()};
};

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
