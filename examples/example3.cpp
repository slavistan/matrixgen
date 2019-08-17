#include <array>
#include <iostream>

#include <matrixgen/core>

int main()
{
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
  const auto mat = matrixgen::adjmat({{3, 3, 1}}, matrixgen::STENCIL<19>, weightfn);

  // Print the matrix. Convert to a dense matrix to make the output look cleaner.
  std::cout << std::endl << Eigen::MatrixXd(mat) << std::endl;
}
