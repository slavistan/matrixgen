#include <array>
#include <iostream>

#include <matrixgen/core>

int main()
{
  /**
   * Generate the matrix and use it.
   */
  using matrixgen::BOUNDCOND;
  const auto mat = matrixgen::adjmat(
      {2, 2, 2},
      matrixgen::stencil7p<BOUNDCOND::PERIODIC>(),
      matrixgen::constweight(1));

  std::cout << std::endl << Eigen::MatrixXd {mat} << std::endl;
}
