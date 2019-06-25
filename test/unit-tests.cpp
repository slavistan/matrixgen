#include <iostream>
#include <array>

#include <matrixgen/adjacency-matrix.hpp>
#include <Eigen/Dense>

auto lambda = [](int a) { return a; };

int main()
{
  const auto stencil = std::array<std::array<int, 3>, 7>{{
      {1, 0, 0}, {-1, 0,0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}, {0, 0, 0}}};
}
