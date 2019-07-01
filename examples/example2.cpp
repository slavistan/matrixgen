#include <array>
#include <iostream>

#include <matrixgen/adjacency-matrix.hpp>
#include <matrixgen/stencil.hpp>

int main()
{
  // In this example we choose a weight-function which returns -1 for the diagonal elements of the matrix. The elements
  // in the upper triangular matrix are set to the sum of the row and column indices. The values in the lower
  // triangular matrix are the negatives thereof.
  auto weightfn = [](const std::array<int, 2>& coords) {

    const auto [row, col] = coords;

    // diagonal
    if (row == col) {
      return -1;
    }
    // upper triangular matrix
    if (row < col) {
      return row + col + 2;
    }
    // lower triangular matrix
    return -(row + col + 2);
  };

  // As in the previous example we use a symmetric 7p stencil, but instead of defining it ourselves, we use a
  // predefined variable template. See the header file for more pre-defined stencils.
  const auto mat = matrixgen::adjmat({{3, 3, 1}}, matrixgen::STENCIL<7>, weightfn);

  // Print the matrix. Convert to a dense matrix to make the output look cleaner.
  std::cout << std::endl << Eigen::MatrixXd(mat) << std::endl;
}
