#include <array>
#include <iostream>

#include <matrixgen/core>

int main()
{
  // Choose a stencil. The stencil encodes which nodes are considered adjacent as a list of 3d offsets.
  // A stencil shall be a range of arbitrary length over elements of type 'std::array<int, 3>'.
  // Note that the stencil is traversed once for each node of the grid. Choose it's type carefully.
  // We use a symmetric 7p-stencil for this example.
  const auto stencil = std::array<std::array<int, 3>, 7>{{
      {1, 0, 0}, {-1, 0,0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}, {0, 0, 0}
  }};

  // Choose a function to compute the matrix's values aka. the weights of the adjacency matrix. The function shall
  // be compatible with either of the following calls:
  //
  // 1. Take no argument.
  //  foo ()
  //
  // 2. Take the matrix entry's coordinates as a tuple of (row, column).
  //  foo (std::array<int, 2>)
  //
  // 3. Take the geometric positions of the node in question and its neighbor as tuples of (x, y, z). The node's
  //    coordinates are passed first.
  //  foo (std::array<int, 3>, std::array<int, 3>)
  //
  // 4. Combination of (2) and (3). First the matrix entry's coordinates are passed, then the nodes' geometric
  //    positions.
  //  foo (std::array<int, 2>, std::array<int, 3>, std::array<int, 3>)
  //
  // For this example we choose the most simple weight-function returning a constant value.
  auto weightfn = []() { return 1; };

  // Now generate the matrix from the grid's dimensions, the stencil and the weight function. We use a 3x3x1 grid to
  // keep the output small enough for your poor terminal. The grid shall be specified as a 'std::array<int, 3>'.
  const auto mat = matrixgen::adjmat({{3, 3, 1}}, stencil, weightfn);

  // Print the matrix. Convert to a dense matrix to make the output look cleaner.
  std::cout << std::endl << Eigen::MatrixXd(mat) << std::endl;
}
