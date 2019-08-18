#pragma once

#include <array>
#include <iostream>
#include <type_traits>

#include <Eigen/Sparse>

#include <gsl/gsl_assert>

namespace matrixgen
{

template <typename Index_t = int>
using DiscreteCoords2d_t = std::array<Index_t, 2>;

template <typename Index_t = int>
using DiscreteCoords3d_t = std::array<Index_t, 3>;

/**
 * Return a node's index within its grid. Indexing is always performed in
 * x-then-y-then-z direction. Indices start at 0.
 */
template <typename Index_t = int>
Index_t
get_node_index(
  const DiscreteCoords3d_t<Index_t>& gridDimensions,
  const DiscreteCoords3d_t<Index_t>& nodeCoords) {

  Expects ( nodeCoords[0] >= 0                );
  Expects ( nodeCoords[1] >= 0                );
  Expects ( nodeCoords[2] >= 0                );
  Expects ( nodeCoords[0] < gridDimensions[0] );
  Expects ( nodeCoords[1] < gridDimensions[1] );
  Expects ( nodeCoords[2] < gridDimensions[2] );

  // TODO: Reduce multiplications by factorization
  const auto index = nodeCoords[0] +
                     nodeCoords[1] * gridDimensions[0] +
                     nodeCoords[2] * gridDimensions[0] * gridDimensions[1];

  Ensures ( index >= 0 );
  Ensures ( index < gridDimensions[0] * gridDimensions[1] * gridDimensions[2] );

  return index;
}

/**
 * Return coordinates (i, j) of an adjacency matrix's entry corresponding to a
 * node at 'myCoords' and its neighbor at 'neighborCoords'.
 */
template <typename Index_t = int>
std::array<Index_t, 2>
get_matrix_entry_coordinates(
  const DiscreteCoords3d_t<Index_t>& gridDimensions,
  const DiscreteCoords3d_t<Index_t>& myCoords,
  const DiscreteCoords3d_t<Index_t>& neighborCoords) {

  Expects ( myCoords[0] >= 0                      );
  Expects ( myCoords[1] >= 0                      );
  Expects ( myCoords[2] >= 0                      );
  Expects ( myCoords[0] < gridDimensions[0]       );
  Expects ( myCoords[1] < gridDimensions[1]       );
  Expects ( myCoords[2] < gridDimensions[2]       );
  Expects ( neighborCoords[0] >= 0                );
  Expects ( neighborCoords[1] >= 0                );
  Expects ( neighborCoords[2] >= 0                );
  Expects ( neighborCoords[0] < gridDimensions[0] );
  Expects ( neighborCoords[1] < gridDimensions[1] );
  Expects ( neighborCoords[2] < gridDimensions[2] );

  const Index_t ii = get_node_index(gridDimensions, myCoords);
  const Index_t jj = get_node_index(gridDimensions, neighborCoords);

  return {{ii, jj}};
}

/**
 * Given the dimensions of a grid and a stencil generate an adjacency matrix.
 * The numeric values of the matrix entries are determined by the weight 
 * function lambda. See examples for different use cases.
 *
 * Returns an Eigen::SparseMatrix whose template parameters may be freely
 * chosed according to the signature of this function template.
 */
template <
  typename StencilFn_t,
  typename WeightFn_t, // TODO: Provide a default weightfn
  typename Scalar_t = double,
  int STORAGE_ORDER = Eigen::RowMajor,
  typename Index_t  = int
    >
Eigen::SparseMatrix<Scalar_t, STORAGE_ORDER, Index_t>
adjmat(
    // TODO: Assert stencil's signature and return type
  const DiscreteCoords3d_t<Index_t> gridDimensions,
  StencilFn_t stencil,
  WeightFn_t computeValue) {

  // The adjacency matrix is a square matrix. Store its height.
  const auto matrixHeight = gridDimensions[0] *
                            gridDimensions[1] *
                            gridDimensions[2];

  /*
   * Insert elements according to the chosen stencil. Traverse the grid in
   * x-then-y-then-z direction and compute the element's numeric value acc. to
   * the weight function.
   * We use initialization via Triplets as described here:
   * http://eigen.tuxfamily.org/dox/group__TutorialSparse.html#TutorialSparseFilling
   *
   * TODO: Parallelize via tbb. Might have to use different initialization
   *       scheme or let each thread create its vector of triplets and then
   *       create a (custom?) range by connecting all the vectors.
   *       Note that parallelization requires that the lambdas don't carry any
   *       state.
   */
  // TODO: I removed the preallocation of memory to avoid passing another
  //       parameter (the number of nonzeros is unknown, as the stencil is now
  //       dynamic).
  // const auto upperLimToCountOfNnz = matrixHeight * std::size(stencil);
  auto triplets = std::vector<Eigen::Triplet<Scalar_t>> {};
  //triplets.reserve(upperLimToCountOfNnz);

  for(auto zz = 0; zz < gridDimensions[2]; ++zz){
    for(auto yy = 0; yy < gridDimensions[1]; ++yy){
      for(auto xx = 0; xx < gridDimensions[0]; ++xx){

        /**
         * For each neighboring node, if it's inside the grid compute the entry's
         * coordinates (i,j) and store it away as triplet.
         */
        const auto myCoords = DiscreteCoords3d_t<Index_t> {xx, yy, zz};
        // TODO: Comment new implementation
        const auto offset_range = stencil(myCoords, gridDimensions);
        for(auto offset_it = offset_range.first; offset_it != offset_range.second; ++offset_it){

          const auto neighborCoords = DiscreteCoords3d_t<Index_t> {
            myCoords[0] + (*offset_it)[0],
            myCoords[1] + (*offset_it)[1],
            myCoords[2] + (*offset_it)[2]
          };

          const auto [ii, jj] = get_matrix_entry_coordinates(
            gridDimensions,
            myCoords,
            neighborCoords
          );

          /**
           * Select the correct implementation depending on the
           * weight-function's signature.
           */
          Scalar_t value;

          // A. WeightFn takes no arguments (e.g. constant weights)
          if constexpr (std::is_invocable_r<Scalar_t, WeightFn_t>()) {
            value = static_cast<Scalar_t>(computeValue());
          }
          // B. WeighFn computes values from the matrix element's positions
          //    (row, column).
          else if constexpr (std::is_invocable_r<Scalar_t, WeightFn_t, std::array<int, 2>>()) {
            value = static_cast<Scalar_t>(computeValue({{ii, jj}}));
          }
          // C. WeightFn computes values from the geometric position of the
          //    and its neighbor node.
          else if constexpr (std::is_invocable_r<
                              Scalar_t,
                              WeightFn_t,
                              DiscreteCoords3d_t<int>,
                              DiscreteCoords3d_t<int>
                             >()) {
            value = static_cast<Scalar_t>(computeValue(myCoords, neighborCoords));
          }
          // D. WeightFn computes values from the matrix element's position
          //    and geometric positions of the node and its neighbor.
          else if constexpr (std::is_invocable_r<
                              Scalar_t,
                              WeightFn_t,
                              DiscreteCoords2d_t<int>,
                              DiscreteCoords3d_t<int>,
                              DiscreteCoords3d_t<int>
                               >()) {
            value = static_cast<Scalar_t>(computeValue({{ii, jj}}, myCoords, neighborCoords));
          }
          // E. WeightFn had too much to drink again.
          else {
            static_assert(
              // TODO: Is there a non-hacky solution to this? Something like
              //       std::abort_compilation("Error: ... ");
              !std::is_same<Scalar_t, Scalar_t>(),
              "Function computing the weights has incompatible signature.");
          }
          triplets.push_back(Eigen::Triplet<Scalar_t>{ii, jj, value});
        }
      }
    }
  }

  // create the sparse matrix from the triplets
  auto result = Eigen::SparseMatrix<Scalar_t, STORAGE_ORDER, Index_t>(matrixHeight, matrixHeight);
  result.setFromTriplets(triplets.begin(), triplets.end());
  return result;
}

} // namespace matrixgen
