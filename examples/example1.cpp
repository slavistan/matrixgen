#include <iostream>
#include <array>

#include <matrixgen/adjacency-matrix.hpp>

int main()
{
  // Choose a stencil. We use a symmetric 7p-stencil for this example. 
  const auto stencil = std::array<matrixgen::DiscretePoint3d, 7>{{
      {1, 0, 0}, {-1, 0,0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}, {0, 0, 0}
  }};

  // Generate the matrix from the grid's dimensions and the stencil. We use a 3x3x1 grid to keep the output small
  // enough for the poor terminal.
  const auto mat = matrixgen::adjmat1(3, 3, 1, stencil);

  // Print the matrix. Converting to a dense matrix beforehand make the output look cleaner.
  std::cout << std::endl << Eigen::MatrixXd(mat) << std::endl;
}
