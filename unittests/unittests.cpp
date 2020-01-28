// TODO: Unittests for `adjmat` for a few small matrices
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <matrixgen/core>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <iostream>

using Scalar_t = double;
using DenseRowMajMat_t  = Eigen::Matrix<Scalar_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

using matrixgen::create;

/* Provide ostream insertion for `std::vector` for stringification by doctest. */
namespace std {
template <typename Elem_t>
std::ostream& operator<<(std::ostream& os, const std::vector<Elem_t>& vec) {
  os << "{";
  for(auto ii = 0u; ii < vec.size() - 1; ++ii) {
    os << vec[ii] << ",";
  }
  if (!vec.empty()) {
    os << vec.back();
  }
  os << "}";
  return os;
}
}

/* Floating-point comparison on `std::vector` to be used by doctest. */
template <typename Elem_t>
bool approx_eq(
    const std::vector<Elem_t>& a,
    const std::vector<Elem_t>& b,
    double eps = 0.02) { // maximum allowed relative deviation

  if (a.size() != b.size()) {
    return false;
  }

  for(auto ii = 0u; ii < a.size(); ++ii) {
    if(a[ii] != doctest::Approx(b[ii]).epsilon(eps)) {
      return false;
    }
  }

  return true;
}


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

TEST_CASE("perturb") {

  using Matrix_t = Eigen::SparseMatrix<double, Eigen::RowMajor>;
  auto matrix = matrixgen::create<Matrix_t>(3, 4,
                  {0, 0, 0, 3,
                   0, 0, 2, 9,
                   0, 3, 4, 8});

  auto perturbed_matrix = matrixgen::perturb(matrix, {0, 1, 2}, matrixgen::seed_from_time());

  std::cout << Eigen::MatrixXd(perturbed_matrix) << std::endl;
  // ??
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

  SUBCASE("Closed-Loop Moving Mean") {

    SUBCASE("Null-radius does nothing but roll over loop min to loop max") {
      const auto input = std::vector {0.0, 0.25, 0.5, 0.75};
      const auto loop_min = 0.0;
      const auto loop_max = 1.0;
      const auto radius = 0;

      auto result = std::vector<double>(input.size());
      matrixgen::closed_loop_moving_mean(
          input.begin(), input.end(), result.begin(),
          loop_min, loop_max, radius);

      const auto targeta = std::vector {1.0, 0.25, 0.5, 0.75};
      const auto targetb = std::vector {0.0, 0.25, 0.5, 0.75};
      REQUIRE((result == targeta || result == targetb));
    }

    SUBCASE("Boundary values use available neighbors only. Symmetric inner values are untouched") {
      const auto input = std::vector {0.0, 0.25, 0.5, 0.75};
      const auto loop_min = 0.0;
      const auto loop_max = 1.0;
      const auto radius = 1;

      auto result = std::vector<double>(input.size());
      matrixgen::closed_loop_moving_mean(
          input.begin(), input.end(), result.begin(),
          loop_min, loop_max, radius);

      const auto target = std::vector {0.125, 0.25, 0.5, 0.625};
      REQUIRE(approx_eq(result, target));
    }

    SUBCASE("Boundary values use available neighbors only. "
        "Asymmetric inner values are treated correctly") {
      const auto input = std::vector {0, 1, 6, 8};
      const auto loop_min = 0.0;
      const auto loop_max = 10.0;
      const auto radius = 1;

      auto result = std::vector<double>(input.size());
      matrixgen::closed_loop_moving_mean(
          input.begin(), input.end(), result.begin(),
          loop_min, loop_max, radius);

      const auto target = std::vector {0.5, 10.0, 8.0, 7.0};
      REQUIRE(approx_eq(result, target));
    }

    SUBCASE("Asymmetric inner values are treated correctly") {
      const auto input = std::vector {0, 6, 1, 9};
      const auto loop_min = 0.0;
      const auto loop_max = 10.0;
      const auto radius = 1;

      auto result = std::vector<double>(input.size());
      matrixgen::closed_loop_moving_mean(
          input.begin(), input.end(), result.begin(),
          loop_min, loop_max, radius);

      const auto target = std::vector {8.0, 10.0, 9.0, 10.0};
      REQUIRE(approx_eq(result, target));
    }
  }

  SUBCASE("Darts Sampling") {
    std::vector<double> quota = {1, 1, 2};
    std::vector<double> bullets(10);
    std::generate(bullets.begin(), bullets.end(), [index = 0]() mutable -> double { return index++/10.0; });

    std::vector<int> indices(bullets.size());
    matrixgen::darts_sampling(quota.begin(), quota.end(), bullets.begin(), bullets.end(), indices.begin());
    REQUIRE(indices == std::vector<int>({0, 0, 0, 1, 1, 1, 2, 2, 2, 2}));
  }

  SUBCASE("Insert") {

    using DenseMatrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

    SUBCASE("Can insert scalars into colmajor matrices") {
      using SparseMatrix_t = Eigen::SparseMatrix<double, Eigen::ColMajor>;

      auto result = SparseMatrix_t(3, 2);
      matrixgen::insert(result, 0, 1, 1.0);
      matrixgen::insert(result, 2, 1, 3.0);

      const auto target = matrixgen::create<SparseMatrix_t>(3, 2,
          {0.0, 1.0,
           0.0, 0.0,
           0.0, 3.0});

      REQUIRE(Eigen::MatrixXd(result) == Eigen::MatrixXd(target));
    }

    SUBCASE("Can insert scalars into rowmajor matrices") {
      using SparseMatrix_t = Eigen::SparseMatrix<double, Eigen::RowMajor>;

      auto result = SparseMatrix_t(3, 2);
      matrixgen::insert(result, 0, 1, 1.0);
      matrixgen::insert(result, 2, 1, 3.0);

      const auto target = matrixgen::create<SparseMatrix_t>(3, 2,
          {0.0, 1.0,
           0.0, 0.0,
           0.0, 3.0});

      REQUIRE(Eigen::MatrixXd(result) == Eigen::MatrixXd(target));
    }

    SUBCASE("Can insert dense matrices into colmajor matrices") {
      using SparseMatrix_t = Eigen::SparseMatrix<double, Eigen::ColMajor>;

      auto filler = DenseMatrix_t(2, 2);
      filler << 1.0, 1.0, 1.0, 1.0;

      auto result = SparseMatrix_t(4, 4);
      matrixgen::insert(result, 2, 2, filler);

      const auto target = matrixgen::create<SparseMatrix_t>(4, 4,
          {0.0, 0.0, 0.0, 0.0,
           0.0, 0.0, 0.0, 0.0,
           0.0, 0.0, 1.0, 1.0,
           0.0, 0.0, 1.0, 1.0});

      REQUIRE(Eigen::MatrixXd(result) == Eigen::MatrixXd(target));
    }

    SUBCASE("Can insert dense matrices into rowmajor matrices") {
      using SparseMatrix_t = Eigen::SparseMatrix<double, Eigen::RowMajor>;

      auto filler = DenseMatrix_t(2, 2);
      filler << 1.0, 1.0, 1.0, 1.0;

      auto result = SparseMatrix_t(4, 4);
      matrixgen::insert(result, 2, 2, filler);

      const auto target = matrixgen::create<SparseMatrix_t>(4, 4,
          {0.0, 0.0, 0.0, 0.0,
           0.0, 0.0, 0.0, 0.0,
           0.0, 0.0, 1.0, 1.0,
           0.0, 0.0, 1.0, 1.0});

      REQUIRE(Eigen::MatrixXd(result) == Eigen::MatrixXd(target));
    }
  }
}
