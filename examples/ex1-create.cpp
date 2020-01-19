/**
 * Example illustrating the basic usage of `matrixgen::create'.
 */
#include <iostream>

#include <matrixgen/core>

int main() {

  /**
   * `create' is used to construct matrix objects from a matrix's dimensions
   * and a dense list of its values. The output format is specified as the
   * single function template parameter.
   *
   * Supported output formats are all specializations of `Eigen::Matrix' and
   * `Eigen::SparseMatrix'. For this example we choose a dense matrix type.
   */
  using Matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

  auto matrix =
    matrixgen::create<Matrix_t>(3 /* rows */, 4 /* cols */,
      {3.1, 0.0, 0.0, 5.2,
       1.0, 0.0, 1.0, 0.0,
       1.0, 7.7, 0.0, 0.0});

  /**
   * Print
   */
  std::cout << Eigen::MatrixXd(matrix) << std::endl;
}
