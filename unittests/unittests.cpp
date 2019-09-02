// TODO: Unittests adjmat for a few small matrices
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <matrixgen/core>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using Scalar_t = double;
using DenseRowMajMat_t  = Eigen::Matrix<Scalar_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

using matrixgen::create;

// Check a subset of all possible output matrix types we want to support.
#define Types                                                               \
  Eigen::Matrix<Scalar_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>, \
  Eigen::Matrix<Scalar_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, \
  Eigen::SparseMatrix<Scalar_t, Eigen::RowMajor>,                           \
  Eigen::SparseMatrix<Scalar_t, Eigen::ColMajor>,                           \
  Eigen::SparseMatrix<Scalar_t, Eigen::ColMajor, int64_t>
TEST_CASE_TEMPLATE("create", MatType_t, Types) {

  /**
   * Create a reference matrix by hand with which to compare the test objects.
   */
  const auto rows = 3;
  const auto cols = 2;
  auto m = DenseRowMajMat_t(rows, cols);
  m << 3, 0,
       0, 1,
       9, 4;

  MatType_t subject;

  SUBCASE("Elements by range") {
    subject = create<MatType_t>(
        rows, cols,
        m.data(),
        std::next(m.data(), rows * cols));
  }

  SUBCASE("Elements by list") {
    subject = create<MatType_t>(rows, cols,
        {3, 0,
         0, 1,
         9, 4});
  }

  CHECK(m == DenseRowMajMat_t(subject));
}

TEST_CASE("assemble") {

  using SparseMatRowMaj_t = Eigen::SparseMatrix<double, Eigen::RowMajor>;
  using SparseMatColMaj_t = Eigen::SparseMatrix<double, Eigen::ColMajor>;
  using DenseMat_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

  SUBCASE("trivial-rowmaj") {
    const auto matrices = std::vector {
      matrixgen::create<SparseMatRowMaj_t>(3, 4,
        {1, 0, 0, 1,
         1, 0, 1, 0,
         1, 1, 0, 0}),
      matrixgen::create<SparseMatRowMaj_t>(3, 4,
        {0, 0, 2, 2,
         2, 0, 0, 2,
         0, 2, 0, 0}),
      matrixgen::create<SparseMatRowMaj_t>(3, 4,
        {0, 0, 3, 3,
         3, 0, 0, 3,
         0, 3, 0, 0})
    };
    const auto indices = std::vector {0, 1, 2};

    const auto result = matrixgen::assemble(
        matrices.begin(),
        matrices.end(),
        indices.begin(),
        indices.end());

    const auto target = matrixgen::create<DenseMat_t>(3, 4,
        {1, 0, 0, 1,
         2, 0, 0, 2,
         0, 3, 0, 0});

    REQUIRE(DenseMat_t(result) == target);
  }

  SUBCASE("different-dimensions-rowmaj") {
    const auto matrices = std::vector {
      matrixgen::create<SparseMatRowMaj_t>(3, 3,
        {1, 0, 1,
         1, 1, 0,
         0, 0, 0}),
      matrixgen::create<SparseMatRowMaj_t>(3, 4,
        {0, 0, 2, 2,
         2, 0, 0, 2,
         0, 2, 0, 0}),
      matrixgen::create<SparseMatRowMaj_t>(3, 4,
        {0, 0, 3, 3,
         3, 0, 0, 3,
         0, 3, 0, 0})
    };
    const auto indices = std::vector {0, 1, 2};

    const auto result = matrixgen::assemble(
        matrices.begin(),
        matrices.end(),
        indices.begin(),
        indices.end());

    const auto target = matrixgen::create<DenseMat_t>(3, 4,
        {1, 0, 1, 0,
         2, 0, 0, 2,
         0, 3, 0, 0});

    REQUIRE(DenseMat_t(result) == target);
  }

  SUBCASE("trivial-colmaj") {
    const auto matrices = std::vector {
      matrixgen::create<SparseMatColMaj_t>(3, 4,
        {1, 0, 0, 1,
         1, 0, 1, 0,
         1, 1, 0, 0}),
      matrixgen::create<SparseMatColMaj_t>(3, 4,
        {0, 0, 2, 2,
         2, 0, 0, 2,
         0, 2, 0, 0}),
      matrixgen::create<SparseMatColMaj_t>(3, 4,
        {0, 0, 3, 3,
         3, 0, 0, 3,
         0, 3, 0, 0})
    };
    const auto indices = std::vector {0, 1, 2, 0};

    const auto result = matrixgen::assemble(
        matrices.begin(),
        matrices.end(),
        indices.begin(),
        indices.end());

    const auto target = matrixgen::create<DenseMat_t>(3, 4,
        {1, 0, 3, 1,
         1, 0, 0, 0,
         1, 2, 0, 0});

    REQUIRE(DenseMat_t(result) == target);
  }

  SUBCASE("different-dimensions-colmaj") {
    const auto matrices = std::vector {
      matrixgen::create<SparseMatColMaj_t>(2, 4,
        {1, 0, 0, 1,
         1, 0, 1, 0}),
      matrixgen::create<SparseMatColMaj_t>(3, 4,
        {0, 0, 2, 2,
         2, 0, 0, 2,
         0, 2, 0, 0}),
      matrixgen::create<SparseMatColMaj_t>(3, 4,
        {0, 0, 3, 3,
         3, 0, 0, 3,
         0, 3, 0, 0})
    };
    const auto indices = std::vector {0, 1, 2, 2};

    const auto result = matrixgen::assemble(
        matrices.begin(),
        matrices.end(),
        indices.begin(),
        indices.end());

    const auto target = matrixgen::create<DenseMat_t>(3, 4,
        {1, 0, 3, 3,
         1, 0, 0, 3,
         0, 2, 0, 0});

    REQUIRE(DenseMat_t(result) == target);
  }
}

