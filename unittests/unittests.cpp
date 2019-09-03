// TODO: Unittests for `adjmat` for a few small matrices
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <matrixgen/core>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using Scalar_t = double;
using DenseRowMajMat_t  = Eigen::Matrix<Scalar_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

using matrixgen::create;

TEST_CASE_TEMPLATE("create", OutMatrix_t,
  Eigen::Matrix<Scalar_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>,
  Eigen::Matrix<Scalar_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>,
  Eigen::SparseMatrix<Scalar_t, Eigen::RowMajor>,
  Eigen::SparseMatrix<Scalar_t, Eigen::ColMajor>,
  Eigen::SparseMatrix<Scalar_t, Eigen::ColMajor, int64_t>
    ) {

  // Create a reference matrix by hand with which to compare the test objects.
  const auto rows = 3;
  const auto cols = 2;
  auto m = DenseRowMajMat_t(rows, cols);
  m << 3, 0,
       0, 1,
       9, 4;

  OutMatrix_t subject;

  SUBCASE("Elements by range") {
    subject = create<OutMatrix_t>(
        rows, cols,
        m.data(),
        std::next(m.data(), rows * cols));
  }

  SUBCASE("Elements by list") {
    subject = create<OutMatrix_t>(rows, cols,
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

TEST_CASE("utility") {

  SUBCASE("Central Moving Sum") {
    SUBCASE("Null-radius does not change elements") {
      const auto input = std::vector<double> {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

      auto result = std::vector<double>(input.size());
      matrixgen::central_moving_sum(input.begin(), input.end(), result.begin(), 0);

      const auto target = input;

      REQUIRE(result == target);
    }

    SUBCASE("Radius == 1") {
      const auto input = std::vector<double> {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

      auto result = std::vector<double>(input.size());
      matrixgen::central_moving_sum(input.begin(), input.end(), result.begin(), 1);

      const auto target = std::vector {3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 19.0};

      REQUIRE(result == target);
    }

    SUBCASE("Radius capturing every element leads to every output being the same") {
      const auto input = std::vector {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

      auto result = std::vector<double>(input.size());
      matrixgen::central_moving_sum(input.begin(), input.end(), result.begin(), 10);

      const auto target = std::vector<double>(10, 55.0);

      REQUIRE(result == target);
    }

    SUBCASE("Radius capturing every element leads to every output being the same. "
            "Radii larger than actual data array do not cause trouble.") {
      auto input = std::vector {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

      auto result = std::vector<double>(input.size());
      matrixgen::central_moving_sum(input.begin(), input.end(), result.begin(), 100);

      const auto target = std::vector<double>(10.0, 55.0);

      REQUIRE(result == target);
    }

    SUBCASE("Test a large input array") {
      const auto input = std::vector<int32_t>(100000, 1);

      auto result = std::vector<int32_t>(input.size());
      matrixgen::central_moving_sum(input.begin(), input.end(), result.begin(), 1);

      auto target = std::vector<int32_t>(input.size(), 3);
      target.front() = 2;
      target.back() = 2;

      REQUIRE(result == target);
    }
  }

  // CONTINUEHERE: Migrate unit tests
  // SUBCASE("Closed-Loop Moving Mean")
  // {
  //   {
  //   const std::vector<double> numbers = {0, 0.25, 0.5, 0.75};
  //   {
  //   //
  //   // Null-radius leads to untouched values unless they're at the boundary.
  //   //
  //   std::vector<double> result(numbers.size());
  //   matrixgen::closed_loop_moving_mean(numbers.begin(), numbers.end(), result.begin(), 0, 1, 0);
  //   REQUIRE(result == std::vector<double>({1.0, 0.25, 0.5, 0.75}));
  //   }
  //   {
  //   //
  //   // Symmetric inner values are untouched + boundary values use all available neighbors only
  //   //
  //   std::vector<double> result(numbers.size());
  //   matrixgen::closed_loop_moving_mean(numbers.begin(), numbers.end(), result.begin(), 0, 1, 1);
  //   REQUIRE(result == std::vector<double>({0.125, 0.25, 0.5, 0.625}));
  //   }
  //   }
  //   {
  //   //
  //   // Non-symmetric inner values get correct treatment
  //   //
  //   const std::vector<double> numbers = {0, 1, 6, 8};
  //   std::vector<double> result(numbers.size());
  //   matrixgen::closed_loop_moving_mean(numbers.begin(), numbers.end(), result.begin(), 0, 10, 1);

  //   std::vector<double> expected = {0.5, 10.0, 8.0, 7.0};
  //   for(auto ii = 0u; ii < result.size(); ++ii)
  //     REQUIRE(result[ii] == expected[ii]);
  //   }
  //   {
  //   //
  //   // Non-symmetric inner values get correct treatment
  //   //
  //   const std::vector<double> numbers = {0, 6, 1, 9};
  //   std::vector<double> result(numbers.size());
  //   matrixgen::closed_loop_moving_mean(numbers.begin(), numbers.end(), result.begin(), 0, 10, 1);

  //   std::vector<double> expected = {8.0, 10.0, 9.0, 10.0};
  //   for(auto ii = 0u; ii < result.size(); ++ii)
  //     REQUIRE(result[ii] ==expected[ii]);
  //   }
  // }
  // SUBCASE("Darts Sampling")
  // {
  //   std::vector<double> quota = {1, 1, 2};
  //   std::vector<double> bullets(10);
  //   std::generate(bullets.begin(), bullets.end(), [index = 0]() mutable -> double { return index++/10.0; });

  //   std::vector<int> indices(bullets.size());
  //   matrixgen::darts_sampling(quota.begin(), quota.end(), bullets.begin(), bullets.end(), indices.begin());
  //   REQUIRE(indices == std::vector<int>({0, 0, 0, 1, 1, 1, 2, 2, 2, 2}));
  // }
}
