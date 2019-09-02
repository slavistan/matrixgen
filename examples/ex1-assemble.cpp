/**
 * Example illustrating the basic usage of `assemble`.
 */
#include <iostream>

#include <matrixgen/core>

int main() {

  using SparseMat_t = Eigen::SparseMatrix<double, Eigen::RowMajor>;
  using DenseMat_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

  /**
   * `assemble` traverses a range of index-pointers pointing matrices from
   * which rows (columns) are drawn row by row (column by column).
   *
   * Depending on the source matrices' layout the output matrix is assembled
   * from the source matrices' rows or columns. Col-major source matrices
   * will produce a col-major output matrix from their columns and v.v for
   * row-major source matrices.
   *
   * For this example we pick three rows from three row-major matrices.
   */
  const auto indices = std::vector {
    2 /* 1st row is drawn from matrix 2 */,
    1 /* 2nd row is drawn from matrix 1 */,
    0 /* 3rd row is drawn from matrix 0 */};
  const auto matrices = std::vector {
    matrixgen::create<SparseMat_t>(3, 4,
      {1, 0, 0, 1,
       1, 0, 1, 0,
       1, 1, 0, 0}), // <- 3rd row
    matrixgen::create<SparseMat_t>(3, 4,
      {0, 0, 2, 2,
       2, 0, 0, 2, // <- 2nd row
       0, 2, 0, 0}),
    matrixgen::create<SparseMat_t>(3, 4,
      {0, 0, 3, 3, // <- 1st row
       3, 0, 0, 3,
       0, 3, 0, 0})
  };

  const auto result = matrixgen::assemble(
      matrices.begin(),
      matrices.end(),
      indices.begin(),
      indices.end());

  /**
   * Convert to a dense matrix for pretty-printing.
   */
  std::cout << std::endl << Eigen::MatrixXd(result) << std::endl;
}
