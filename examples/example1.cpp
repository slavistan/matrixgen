#include <iostream>
#include <array>

#include <matrixgen/adjacency-matrix.hpp>

int main()
{
  // Choose a stencil. The stencil encodes which nodes are considered adjacent as a list of 3d offsets.
  // A stencil shall be a 'std::array<T, n>' of arbitrary length 'n' where T is 'std::array<int, 3>'.
  // We use a symmetric 7p-stencil for this example.
  const auto stencil = std::array<std::array<int, 3>, 7>{{
      {1, 0, 0}, {-1, 0,0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}, {0, 0, 0}
  }};

  // Choose a function to compute the matrix's values aka. the weights of the adjacency matrix. The function shall
  // accept three arguments:
  //
  // 1. a 'std::array<int, 2>' corresp. to the matrix entry's coordinates (i, j)
  // 2. a 'std::array<int, 3>' corresp. to the node's coordinates for which the adjacency pattern is generated
  // 3. a 'std::array<int, 3>' corresp. to the neighboring node's coordinates
  //
  // We choose a function which sets all diagonal elements stemming from the {{0, 0, 0}} offset in the stencil
  // to -1 while all non-diagonals are set to their row's index (starting at 1).
  auto computeValue = [](
      const std::array<int, 2>& matrixEntryCoords,
      const std::array<int, 3>& nodeCoords,
      const std::array<int, 3>& neighborCoords){
    if(matrixEntryCoords[0] == matrixEntryCoords[1]) return -1;
    return matrixEntryCoords[0] + 1;
  };

  // Now generate the matrix from the grid's dimensions, the stencil and the weight function. We use a 3x3x1 grid to
  // keep the output small enough for your poor terminal. The grid shall be specified as a 'std::array<int, 3>'.
  const auto mat = matrixgen::adjmat({{3, 3, 1}}, stencil, computeValue);

  // Print the matrix. Convert to a dense matrix to make the output look cleaner.
  std::cout << std::endl << Eigen::MatrixXd(mat) << std::endl;
}
