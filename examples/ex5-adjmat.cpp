/**
 * Example displaying the inner workings of the adjacency function abstraction
 * by implementing a more complex example.
 */
#include <array>
#include <iostream>
#include <utility>

#include <matrixgen/core>

#include <Eigen/Sparse>

int main() {

  using DiscreteCoords3d_t = std::array<int32_t, 3>;
  const auto grid = std::array {5, 5, 1};

  /**
   * TODO: Explain
   * TODO: Note about second version of adjfn + hardcoding information about the grid
   * TODO: "Check out the preset stencil7p"
   */

  auto adjfn = [offsets = std::array<DiscreteCoords3d_t, 1> {{0, 0, 0}},
                gridCenter = DiscreteCoords3d_t {2, 2, 0}]
    (const DiscreteCoords3d_t& nodeCoords) mutable {

    offsets.front()[0] = gridCenter[0] - nodeCoords[0];
    offsets.front()[1] = gridCenter[1] - nodeCoords[1];
    offsets.front()[2] = gridCenter[2] - nodeCoords[2];
    return std::pair {offsets.begin(), offsets.end()};
  };

  const auto mat = matrixgen::adjmat(grid, adjfn, matrixgen::constweight(7));
  std::cout << std::endl << Eigen::MatrixXd(mat) << std::endl;
}
