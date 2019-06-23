#pragma once

#include <iostream>

#include <Eigen/Sparse>

namespace matrixgen
{

using DiscretePoint3d = std::array<int, 3>;

// Return a node's index within its grid. Counting is always performed in x-then-y-then-z direction.
int get_node_index(const DiscretePoint3d& gridDimensions, const DiscretePoint3d& myCoords){
  return myCoords[0] +
         myCoords[1] * gridDimensions[0] +
         myCoords[2] * gridDimensions[0] * gridDimensions[1];
}

// Return coordinates of an adjacency matrix's entry corresponding to a node at 'myCoords' and its neighbor at
// 'neighborCoords'.
std::array<int, 2> get_matrix_entry_coordinates(
  const DiscretePoint3d& gridDimensions,
  const DiscretePoint3d& myCoords,
  const DiscretePoint3d& neighborCoords)
{
  const int ii = get_node_index(gridDimensions, myCoords);
  const int jj = get_node_index(gridDimensions, neighborCoords);
  return {{ii, jj}};
}

template <
  typename Stencil_t, // std::vector/array over 3d integral points
  typename Scalar_t = double, 
  int StorageOrder = Eigen::RowMajor,
  typename Index_t = int
    >
Eigen::SparseMatrix<Scalar_t, StorageOrder, Index_t>
adjmat1(const int dimx, // number of grid nodes in the x-dimension
        const int dimy,
        const int dimz,
        const Stencil_t& stencil,
        const Scalar_t value = 1 ) // the nonzeros' value 
{
  const auto matrixDim = dimx * dimy * dimz; // height/width of adj. matrix

  // Now insert elements according to the chosen stencil. Traverse the grid in x-then-y-then-z direction.
  // We use initialization via Triplets as described here:
  // http://eigen.tuxfamily.org/dox/group__TutorialSparse.html#TutorialSparseFilling
  const auto upperLimToCountOfNnz = matrixDim * std::size(stencil);
  std::vector<Eigen::Triplet<Scalar_t>> triplets;
  triplets.reserve(upperLimToCountOfNnz);

  for(auto zz = 0; zz < dimz; ++zz){
    for(auto yy = 0; yy < dimy; ++yy){
      for(auto xx = 0; xx < dimx; ++xx){
        // for each neighboring node, if its inside the grid compute the entry's coordinates (i,j) and
        // store it away as triplet.
        for(const auto& point: stencil){
          const auto neighborX = xx + point[0];
          const auto neighborY = yy + point[1];
          const auto neighborZ = zz + point[2];
          if(0 <= neighborX && neighborX < dimx &&
             0 <= neighborY && neighborY < dimy &&
             0 <= neighborZ && neighborZ < dimz)
          {
            const auto [ii, jj] = get_matrix_entry_coordinates(
                {{dimx, dimy, dimz}},
                {{xx, yy, zz}},
                {{xx + point[0], yy + point[1], zz + point[2]}});
            triplets.push_back(Eigen::Triplet<Scalar_t>{ii, jj, value});
          }
        }
      }
    }
  }

  // create the sparse matrix from the triplet list
  auto result = Eigen::SparseMatrix<Scalar_t, StorageOrder, Index_t>(matrixDim, matrixDim);
  result.setFromTriplets(triplets.begin(), triplets.end());
  return result;
}

} // namespace matrixgen
