/**
 * \file
 * \author Stanislaw Hüll
 */
#pragma once

#include <Eigen/Sparse>
#include <gsl/gsl_assert>

#include <initializer_list>
namespace matrixgen {

template <
  typename OutMatrix_t,
  typename ListType_t = double
    >
OutMatrix_t create(
    uint32_t numRows,
    uint32_t numCols,
    std::initializer_list<ListType_t> list) {


}

}
