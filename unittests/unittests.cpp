/**
 * \file
 * \author Stanislaw Hüll
 * \brief Unit tests
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <matrixgen/core>

#include <Eigen/Core>

TEMPLATE_TEST_CASE("Matrices are created", "[create]", double, float, int32_t, int64_t) {

  using DenseColMajMat_t = Eigen::Matrix<TestType, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
  using DenseRowMajMat_t = Eigen::Matrix<TestType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

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

  const auto myMatrix = matrixgen::create<DenseColMajMat_t>(rows, cols, elems);
  const auto myMatrix2 = matrixgen::create<DenseRowMajMat_t>(rows, cols, elems);

  REQUIRE(referenceMat == myMatrix);
  REQUIRE(referenceMat == myMatrix2);
  }
//  std::cout << myMatrix << std::endl;
//  using DenseMatrix2_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
//  auto myMatrix2 = matrixgen::create<DenseMatrix2_t>(3, 2,
//      { 3.14,   0,
//           0, 1.1,
//         9.2,   0 });
//  std::cout << myMatrix2 << std::endl;
//  (void)myMatrix;
}
