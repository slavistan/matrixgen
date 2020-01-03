/**
 * Example illustrating the basic usage of `matrixgen::create'.
 */
#include <iostream>

#include <matrixgen/core>

int main() {

  using SparseMat_t = Eigen::SparseMatrix<double, Eigen::RowMajor>;

  /**
   * `create' is used to construct matrix objects from a matrix's dimensions
   * and a dense list of its values. The output format is specified as the
   * single function template parameter.
   *
   * Supported output formats are all specializations of `Eigen::Matrix' and
   * `Eigen::SparseMatrix'.`
   */
  auto matrix =
    matrixgen::create<SparseMat_t>(3, 4,
      {3.14, 0, 0, 5.3,
       1, 0, 1, 0,
       1, 1, 0, 0});

  /**
   * Convert to an `Eigen::Matrix' for pretty-printing.
   */
  std::cout << std::endl << Eigen::MatrixXd(matrix) << std::endl;
}
