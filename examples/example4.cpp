#include <array>
#include <iostream>

#include <matrixgen/core>

int main()
{
  /**
   * Choose the grid's extent as a `std::array<int, 3>`. The coordinates denote
   * the number of nodes in the x, y, z direction, with the grid's origin being
   * located at (0, 0, 0).
   *
   * For the examples we choose a small grid to keep the output matrix small
   * enough for your poor terminal.
   */
  auto weightfn = []() { return 1; };

  /**
   * Generate the matrix and use it.
   */
  using matrixgen::BOUNDCOND;
  const auto mat = matrixgen::adjmat(
      {2, 2, 2},
      matrixgen::stencil7p<BOUNDCOND::PERIODIC>(),
      weightfn);

  std::cout << std::endl << Eigen::MatrixXd {mat} << std::endl;
}
