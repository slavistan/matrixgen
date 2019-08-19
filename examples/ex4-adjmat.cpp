#include <array>
#include <iostream>

#include <matrixgen/core>

int main()
{
  /**
   * Generate the matrix and use it.
   */
  using matrixgen::BC;
  auto stencil = matrixgen::stencil7p<
                  BC::DIRICHLET,
                  BC::PERIODIC,
                  BC::DIRICHLET>();
   const auto mat = matrixgen::adjmat(
       {1, 1, 1},
       stencil,
       matrixgen::constweight(1));

   std::cout << std::endl << Eigen::MatrixXd {mat} << std::endl;
}
