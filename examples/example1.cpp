#include <iostream>
#include <array>

#include <matrixgen/adjacency-matrix.hpp>

auto lambda = [](int a) { return a; };

int main()
{
  const auto stencil = std::array<matrixgen::DiscretePoint3d, 7>{{
      {1, 0, 0}, {-1, 0,0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}, {0, 0, 0}
  }};

  const auto mat = matrixgen::adjmat1(3, 3, 3, stencil);
  std::cout << std::endl << Eigen::MatrixXd(mat) << std::endl;
}
