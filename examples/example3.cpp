/**
 * Example generating an adjacency matrix using a more complex weight function.
 */
#include <array>
#include <iostream>

#include <matrixgen/core>

int main()
{
  /**
   * Define a weightfunction which shall return `-1` for any diagonal element.
   * The elements in the upper trianglar matrix shall be the sum of the row
   * and column indices (starting at `1` to avoid nulls) while the the lower
   * triangular matrix's values shall be the negatives thereof.
   */
  auto weightfn = [](const std::array<int, 2>& coords) {

    const auto [row, col] = coords;
    if (row == col) { return -1; }            // diagonal
    if (row < col) { return row + col + 2; }  // upper triangular matrix
    return -(row + col + 2);                  // lower triangular matrix
  };

  /**
   * Generate and use. We convert the output to a dense matrix for pretty
   * output.
   */
  const auto mat = matrixgen::adjmat({3, 3, 1}, matrixgen::STENCIL<7>, weightfn);
  std::cout << std::endl << Eigen::MatrixXd(mat) << std::endl;
}
