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

// Return true if node at 'myCoords' is contained inside the grid.
bool is_inside_grid(const DiscretePoint3d& gridDimensions,
                    const DiscretePoint3d& myCoords){
  return (0 <= myCoords[0] && myCoords[0] < gridDimensions[0] &&
          0 <= myCoords[1] && myCoords[1] < gridDimensions[1] &&
          0 <= myCoords[2] && myCoords[2] < gridDimensions[2]);
}

template <
  typename Stencil_t, // std::vector/array over 3d integral points
  typename ValueFunc_t, // lambda computing the values
  typename Scalar_t = double, 
  int StorageOrder = Eigen::RowMajor, // TODO: implement ColMajor
  typename Index_t = int
    >
Eigen::SparseMatrix<Scalar_t, StorageOrder, Index_t>
adjmat(const DiscretePoint3d gridDimensions,
       const Stencil_t& stencil,
       ValueFunc_t computeValue)
{
  // Note that the adjacency matrix is a square matrix. Wherever its width is required we use its height.
  const auto matrixHeight = gridDimensions[0] *
                            gridDimensions[1] *
                            gridDimensions[2];

  // Now insert elements according to the chosen stencil. Traverse the grid in x-then-y-then-z direction.
  // We use initialization via Triplets as described here:
  // http://eigen.tuxfamily.org/dox/group__TutorialSparse.html#TutorialSparseFilling
  const auto upperLimToCountOfNnz = matrixHeight * std::size(stencil);
  std::vector<Eigen::Triplet<Scalar_t>> triplets;
  triplets.reserve(upperLimToCountOfNnz);

  for(auto zz = 0; zz < gridDimensions[2]; ++zz){
    for(auto yy = 0; yy < gridDimensions[1]; ++yy){
      for(auto xx = 0; xx < gridDimensions[0]; ++xx){
        // for each neighboring node, if it's inside the grid compute the entry's coordinates (i,j) and
        // store it away as triplet.
        for(const auto& offset: stencil){
          const auto myCoords = DiscretePoint3d{{xx, yy, zz}};
          const auto neighborCoords = DiscretePoint3d{{
              myCoords[0] + offset[0],
              myCoords[1] + offset[1],
              myCoords[2] + offset[2]
          }};
          if(is_inside_grid(gridDimensions, neighborCoords)){
            const auto [ii, jj] = get_matrix_entry_coordinates(
                gridDimensions,
                myCoords,
                neighborCoords
            );

            const auto value = static_cast<Scalar_t>(computeValue({{ii, jj}}, myCoords, neighborCoords));
            triplets.push_back(Eigen::Triplet<Scalar_t>{ii, jj, value});
          }
        }
      }
    }
  }

  // create the sparse matrix from the triplets
  auto result = Eigen::SparseMatrix<Scalar_t, StorageOrder, Index_t>(matrixHeight, matrixHeight);
  result.setFromTriplets(triplets.begin(), triplets.end());
  return result;
}

} // namespace matrixgen
