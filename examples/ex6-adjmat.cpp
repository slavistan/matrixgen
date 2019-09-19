/**
 * Display of useful preset weight functions.
 * TODO: Add a little text
 */

#include <matrixgen/core>

int main() {
  {
  auto weightfn = matrixgen::sinusoid_add(1.1, 0.3, 0.2);
  const auto mat = matrixgen::adjmat({3, 3, 1}, matrixgen::stencil7p(), weightfn);
  std::cout << std::endl << Eigen::MatrixXd(mat) << std::endl;
  }
  {
  auto weightfn = matrixgen::sinusoid_add_bias(1.1, 0.3, 0.2);
  const auto mat = matrixgen::adjmat({3, 3, 1}, matrixgen::stencil7p(), weightfn);
  std::cout << std::endl << Eigen::MatrixXd(mat) << std::endl;
  }
  {
  auto weightfn = matrixgen::sinusoid_mul(1.1, 0.3, 0.2);
  const auto mat = matrixgen::adjmat({3, 3, 2}, matrixgen::stencil7p(), weightfn);
  std::cout << std::endl << Eigen::MatrixXd(mat) << std::endl;
  }
  {
  auto weightfn = matrixgen::sinusoid_mul_bias(1.1, 0.3, 0.2);
  const auto mat = matrixgen::adjmat({3, 3, 2}, matrixgen::stencil7p(), weightfn);
  std::cout << std::endl << Eigen::MatrixXd(mat) << std::endl;
  }
}
