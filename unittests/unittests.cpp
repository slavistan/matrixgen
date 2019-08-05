/**
 * \file
 * \author Stanislaw Hüll
 * \brief Unit tests
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <matrixgen/core>

#include <Eigen/Dense>
#include <Eigen/Sparse>

TEMPLATE_TEST_CASE("Matrices are created", "[create]", double, float, int32_t, int64_t) {

  using DenseColMajMat_t = Eigen::Matrix<TestType, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
  using DenseRowMajMat_t = Eigen::Matrix<TestType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using SparseRowMajMat_t = Eigen::SparseMatrix<TestType, Eigen::RowMajor>;
  using SparseColMajMat_t = Eigen::SparseMatrix<TestType, Eigen::ColMajor>;

  /**
   * Create a reference matrix by hand against we'll compare the test objects.
   * Use the row-major layout used by 'matrixgen::create' and extract the values
   * into a vector for later usage.
   */
  const auto rows = 3;
  const auto cols = 2;
  auto referenceMat = DenseRowMajMat_t(rows, cols);
  referenceMat <<
    3, 0,
    0, 1,
    9, 4;

  const auto elems = std::vector<TestType>(
      referenceMat.data(),
      std::next(referenceMat.data(), rows * cols));

  SECTION("Dense matrices work for row-major and col-major layouts.")
  {
  const auto myMatrix = matrixgen::create<DenseColMajMat_t>(rows, cols, elems.begin(), elems.end());
  const auto myMatrix2 = matrixgen::create<DenseRowMajMat_t>(rows, cols, elems.begin(), elems.end());

  REQUIRE(referenceMat == myMatrix);
  REQUIRE(referenceMat == myMatrix2);
  }

  SECTION("Sparse matrices work for row-major and col-major layouts.")
  {
  const auto myMatrix = matrixgen::create<SparseRowMajMat_t>(rows, cols, elems.begin(), elems.end());
  const auto myMatrix2 = matrixgen::create<SparseColMajMat_t>(rows, cols, elems.begin(), elems.end());

  // Convert to dense matrix so that we may use operator==
  REQUIRE(referenceMat == DenseRowMajMat_t(myMatrix));
  REQUIRE(referenceMat == DenseRowMajMat_t(myMatrix2));
  }
}
