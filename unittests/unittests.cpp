/**
 * \file
 * \author Stanislaw Hüll
 * \brief Unit tests
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <matrixgen/core>

#include <Eigen/Core>

TEST_CASE("Creation of matrices works") {

  using DenseMatrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  auto myMatrix = matrixgen::create<DenseMatrix_t>(3, 2,
      { 3.14,   0,
           0, 1.1,
         9.2,   0 });

  std::cout << myMatrix << std::endl;
  using DenseMatrix2_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  auto myMatrix2 = matrixgen::create<DenseMatrix2_t>(3, 2,
      { 3.14,   0,
           0, 1.1,
         9.2,   0 });
  std::cout << myMatrix2 << std::endl;
  (void)myMatrix;
}
