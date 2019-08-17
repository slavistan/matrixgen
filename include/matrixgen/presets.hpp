#pragma once

#include <array>

namespace matrixgen {

/**
 * Collection of predefined stencils.
 */
template <std::size_t SIZE, typename Scalar_t = int>
const std::array<std::array<Scalar_t, 3>, SIZE> STENCIL;

/**
 * Symmetric 7p-stencil
 */
template <typename Scalar_t>
const std::array<std::array<Scalar_t, 3>, 7> STENCIL<7, Scalar_t> =
  {{{  0,  0,  0},
    { -1,  0,  0}, {  0,  -1,  0}, {  0,  0, -1},
    {  1,  0,  0}, {  0,   1,  0}, {  0,  0,  1}}};

/**
 * Symmetric 19p-stencil
 */
template <typename Scalar_t>
const std::array<std::array<Scalar_t, 3>, 19> STENCIL<19, Scalar_t> =
  {{{  0,  0,  0},
    { -1,  0,  0}, {  0, -1,  0}, {  0,  0, -1},
    {  1,  0,  0}, {  0,  1,  0}, {  0,  0,  1},
    { -1, -1,  0}, { -1,  1,  0}, {  1, -1,  0}, {  1,  1,  0},
    { -1,  0, -1}, { -1,  0,  1}, {  1,  0, -1}, {  1,  0,  1},
    {  0, -1, -1}, {  0, -1,  1}, {  0,  1, -1}, {  0,  1,  1}}};

/**
 * Symmetric 27p-stencil
 */
template <typename Scalar_t>
const std::array<std::array<Scalar_t, 3>, 27> STENCIL<27, Scalar_t> =
  {{{  0,  0,  0},
    { -1,  0,  0}, {  0, -1,  0}, {  0,  0, -1},
    {  1,  0,  0}, {  0,  1,  0}, {  0,  0,  1},
    { -1, -1,  0}, { -1,  1,  0}, {  1, -1,  0}, {  1,  1,  0},
    { -1,  0, -1}, { -1,  0,  1}, {  1,  0, -1}, {  1,  0,  1},
    {  0, -1, -1}, {  0, -1,  1}, {  0,  1, -1}, {  0,  1,  1},
    { -1, -1, -1}, { -1, -1,  1}, { -1,  1, -1}, { -1,  1,  1},
    {  1, -1, -1}, {  1, -1,  1}, {  1,  1, -1}, {  1,  1,  1}}};

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
