/**
 * Example displaying the inner workings of the stencil abstraction by
 * implementing a 1-p stencil.
 */
#include <array>
#include <iostream>
#include <utility>

#include <matrixgen/core>

int main() {

  using DiscreteCoords3d_t = std::array<int32_t, 3>;
  const auto grid = std::array {4, 4, 1};

  /**
   * TODO: Explain
   */
  const auto offsets = std::array<DiscreteCoords3d_t, 1> {{0, 0, 0}};
  auto adjfn = [&offsets](
      const DiscreteCoords3d_t& /* node */,
      const DiscreteCoords3d_t& /* gridDimensions */) {
    return std::pair {offsets.cbegin(), offsets.cend()};
  };

  const auto mat = matrixgen::adjmat(grid, adjfn, matrixgen::constweight(5));
  std::cout << std::endl << mat << std::endl;
}
