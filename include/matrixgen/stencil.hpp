#pragma once

#include <array>

namespace matrixgen {

template <std::size_t SIZE, typename Scalar_t = int>
const std::array<std::array<Scalar_t, 3>, SIZE> STENCIL;

/*
 *  7p - STENCIL
 */
template <typename Scalar_t>
const std::array<std::array<Scalar_t, 3>, 7> STENCIL<7, Scalar_t> =
  {{{  0,  0,  0},
    { -1,  0,  0}, {  0,  -1,  0}, {  0,  0, -1},
    {  1,  0,  0}, {  0,   1,  0}, {  0,  0,  1}}};

/*
 *  19p - STENCIL
 */
template <typename Scalar_t>
const std::array<std::array<Scalar_t, 3>, 19> STENCIL<19, Scalar_t> =
  {{{  0,  0,  0},
    { -1,  0,  0}, {  0, -1,  0}, {  0,  0, -1},
    {  1,  0,  0}, {  0,  1,  0}, {  0,  0,  1},
    { -1, -1,  0}, { -1,  1,  0}, {  1, -1,  0}, {  1,  1,  0},
    { -1,  0, -1}, { -1,  0,  1}, {  1,  0, -1}, {  1,  0,  1},
    {  0, -1, -1}, {  0, -1,  1}, {  0,  1, -1}, {  0,  1,  1}}};
/*
 *  27p - STENCIL
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

} // matrixgen
