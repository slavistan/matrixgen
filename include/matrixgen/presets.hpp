#pragma once

namespace matrixgen {

/**
 * Collection of predefined stencils.
 */


/**
 * Collection of predefined weight functions.
 *
 * TODO: Is this the correct way to implement preset lambdas? Pass them as return
 *       "value"?
 */
template <
  typename Scalar_t = double
    >
auto constweight(Scalar_t val = 1) {

  return ([val]() { return static_cast<Scalar_t>(val); });
}

} // namespace matrixgen
