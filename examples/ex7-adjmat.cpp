/**
 * Example demonstrating the usage of full wrapper presets (wrappers around adjmat).
 */
#include <matrixgen/core>

#include <Eigen/Dense>

int main() {

  /**
   * Some matrices cannot be generated efficiently using the sematics of
   * `matrixgen::adjmat` and thus we provide convenience wrappers which
   * hide steps of pre- and post-processing.
   *
   * This example utilizes `structured_grid_sinusoidal` which evaluates
   * a sinusoidal at the midpoint between the node and its neighbor.
   * In order to ensure diagonal dominance the matrix's diagonal values are
   * set greater than the  accumulate of the row's sum of absolute values.
   *
   * You're encouraged to have a look at the implementation of the wrapper
   * in 'presets.hpp'.
   */
  const auto mat = matrixgen::structured_grid_sinusoidal({4, 4, 4}, matrixgen::stencil7p(), 1.1, 1.2, 1.3);
  std::cout << std::endl << Eigen::MatrixXd(mat) << std::endl;
}
