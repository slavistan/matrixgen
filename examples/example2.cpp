/**
 * Example displaying the generation of adjacency matrices using a compact
 * syntax. Equivalent to example 1.
 */
#include <array>
#include <iostream>

#include <matrixgen/core>

int main()
{
  /**
   * Create the very same matrix from example 1 but use available presets for
   * the stencil and the weightfunction. See 'matrixgen/presets.hpp' for other
   * predefined stencils and weighfunctions.
   */
  const auto mat = matrixgen::adjmat({3, 3, 1}, matrixgen::STENCIL<7>, matrixgen::constweight(1));
  std::cout << std::endl << mat << std::endl;
}
