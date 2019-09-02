/**
 * Example illustrating displaying 0-padding for `assemble`.
 */
#include <iostream>

#include <matrixgen/core>

int main() {

  using SparseMat_t = Eigen::SparseMatrix<double, Eigen::ColMajor>;
  using DenseMat_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

  /**
   * If the source matrices' inner dimesions [*] do not match, the output
   * matrix's inner dimension is determined by the maximum of the source
   * matrices maximum. Missing values in the smaller matrices are 0-padded.
   *
   * For this example we choose a col-major layout and create a matrix from
   * the source matrices' columns. The first source matrix's height is less
   * than that of the other matrices which leads us to pad the missing values
   * with 0.
   */
  const auto indices = std::vector {
    0 /* 1st column is drawn from matrix 0 */,
    1 /* 2nd column is drawn from matrix 1 */,
    2 /* 3rd column is drawn from matrix 2 */,
    1 /* 4th column is drawn from matrix 0 */};

  const auto matrices = std::vector {
    matrixgen::create<SparseMat_t>(2, 4,
      {1, 0, 0, 1,
       1, 0, 1, 0}),
    // ^~~ 1st column
    matrixgen::create<SparseMat_t>(3, 4,
      {0, 0, 2, 2,
       2, 0, 0, 2,
       0, 2, 0, 0}),
    //    ^~~~~~^~~ 2nd and 4th columns
    matrixgen::create<SparseMat_t>(3, 4,
      {0, 0, 3, 3,
       3, 0, 0, 3,
       0, 3, 0, 0})
    //       ^~~ 3rd column
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

/**
 * [*] `inner` and `outer` dimension are Eigen-lingo referring to the matrix
 *     dimension corresponding to the matrices data layout. Hence for a row-
 *     major matrix `inner` refers to the columns and `outer` to the rows and
 *     v.v.
 */
