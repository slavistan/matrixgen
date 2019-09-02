#pragma once

#include <array>
#include <iostream>
#include <type_traits>

#include <Eigen/Sparse>

#include <gsl/gsl_assert>

namespace matrixgen::implementation
{

template <typename Index_t = int>
using DiscreteCoords2d_t = std::array<Index_t, 2>;

template <typename Index_t = int>
using DiscreteCoords3d_t = std::array<Index_t, 3>;

/**
 * Return a node's index within its grid. Indexing is always performed in
 * x-then-y-then-z direction. Indices start at 0.
 *
 * Used as a helper within `get_matrix_entry_coordinates`.
 */
template <typename Index_t = int>
Index_t
get_node_index(
  const DiscreteCoords3d_t<Index_t>& nodeCoords,
  const DiscreteCoords3d_t<Index_t>& gridDimensions) {

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
 * Given a node at `coords` and its neighbor at `neighborCoords` return the
 * corresponding adjacency matrix's non-zero's position as (row, column).
 *
 * Used to determine non-zero coordinates' within the adjacency matrices. A
 * connection between a node with index `n` and its neighbor with index `m`
 * will produce a non-zero matrix entry at (n, m).
 */
template <typename Index_t = int>
std::array<Index_t, 2>
get_matrix_entry_coordinates(
  const DiscreteCoords3d_t<Index_t>& coords,
  const DiscreteCoords3d_t<Index_t>& neighborCoords,
  const DiscreteCoords3d_t<Index_t>& gridDimensions) {

  Expects ( coords[0] >= 0                        );
  Expects ( coords[1] >= 0                        );
  Expects ( coords[2] >= 0                        );
  Expects ( coords[0] < gridDimensions[0]         );
  Expects ( coords[1] < gridDimensions[1]         );
  Expects ( coords[2] < gridDimensions[2]         );
  Expects ( neighborCoords[0] >= 0                );
  Expects ( neighborCoords[1] >= 0                );
  Expects ( neighborCoords[2] >= 0                );
  Expects ( neighborCoords[0] < gridDimensions[0] );
  Expects ( neighborCoords[1] < gridDimensions[1] );
  Expects ( neighborCoords[2] < gridDimensions[2] );

  const Index_t ii = get_node_index(coords, gridDimensions);
  const Index_t jj = get_node_index(neighborCoords, gridDimensions);

  return {{ii, jj}};
}

template <
  typename OutMatrix_t,
  typename StencilFn_t,
  typename WeightFn_t,
  typename Index_t
    >
struct Adjmat {};

/**
 * Specialization of `adjmat` for `Eigen::SparseMatrix<>` return types.
 * This is the only output matrix type we'll support for now due to
 * the fact that adjacency matrices are inherently sparse.
 * TODO: static_assert the above in the dispatcher.
 */
template <
  typename Scalar_t,
  int ALIGNMENT,
  typename EigenIndex_t,
  typename StencilFn_t,
  typename WeightFn_t,
  typename Index_t
    >
struct Adjmat<
  Eigen::SparseMatrix<Scalar_t, ALIGNMENT, EigenIndex_t>,
  StencilFn_t,
  WeightFn_t,
  Index_t
    >
{

  using Matrix_t = Eigen::SparseMatrix<Scalar_t, ALIGNMENT, Index_t>;

  /**
   * Given the dimensions of a grid and a stencil generate an adjacency matrix.
   * The numeric values of the matrix entries are determined by the weight 
   * function lambda. See examples for different use cases.
   *
   * Returns an Eigen::SparseMatrix whose template parameters may be freely
   * chosed according to the signature of this function template.
   */
  static
  Matrix_t
  invoke(
    const DiscreteCoords3d_t<Index_t> gridDimensions,
    StencilFn_t stencilfn,
    WeightFn_t weightfn) {

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
     */
    auto triplets = std::vector<Eigen::Triplet<Scalar_t>> {};
    // TODO: As the stencil is now dynamically computed for each node there's no
    //       such thing as a static size of the stencil for me to assume here.
    //       Hence I removed the call to `triplets.reserve()` and use dynamic
    //       reallocation via `push_back()`. Implement some option for the user
    //       to specify an upper limit to the the triplets vector's size.
    //       Below is the previously used allocation code.
    // const auto upperLimToCountOfNnz = matrixHeight * std::size(stencilfn);
    // triplets.reserve(upperLimToCountOfNnz);

    for(auto zz = 0; zz < gridDimensions[2]; ++zz){
      for(auto yy = 0; yy < gridDimensions[1]; ++yy){
        for(auto xx = 0; xx < gridDimensions[0]; ++xx){

          /**
           * For each neighboring node, if it's inside the grid compute the entry's
           * coordinates (i,j) and store it away as triplet.
           */
          const auto myCoords = DiscreteCoords3d_t<Index_t> {xx, yy, zz};
          // TODO: Comment new implementation
          const auto offsetRange = stencilfn(myCoords, gridDimensions);
          for(auto offsetIt = offsetRange.first; offsetIt != offsetRange.second; ++offsetIt){

            const auto neighborCoords = DiscreteCoords3d_t<Index_t> {
              myCoords[0] + (*offsetIt)[0],
              myCoords[1] + (*offsetIt)[1],
              myCoords[2] + (*offsetIt)[2]
            };

            const auto [ii, jj] = implementation::get_matrix_entry_coordinates(
              myCoords,
              neighborCoords,
              gridDimensions
            );

            /**
             * Select the correct implementation depending on the
             * weight-function's signature.
             */
            Scalar_t value;

            // A. WeightFn takes no arguments (e.g. constant weights)
            if constexpr (std::is_invocable_r<Scalar_t, WeightFn_t>()) {
              value = static_cast<Scalar_t>(weightfn());
            }
            // B. WeighFn computes values from the matrix element's positions
            //    (row, column).
            else if constexpr (std::is_invocable_r<Scalar_t, WeightFn_t, std::array<int, 2>>()) {
              value = static_cast<Scalar_t>(weightfn({{ii, jj}}));
            }
            // C. WeightFn computes values from the geometric position of the
            //    and its neighbor node.
            else if constexpr (std::is_invocable_r<
                                Scalar_t,
                                WeightFn_t,
                                DiscreteCoords3d_t<int>,
                                DiscreteCoords3d_t<int>
                               >()) {
              value = static_cast<Scalar_t>(weightfn(myCoords, neighborCoords));
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
              value = static_cast<Scalar_t>(weightfn({{ii, jj}}, myCoords, neighborCoords));
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

    /**
     * Generate the matrix from the set of triplets. Note that `setFromTriplets`
     * does respect multiple triplets for the same matrix entry in which case
     * the values are simply added. This circumstance is relevant for boundary
     * conditions which map onto the same neighboring node from disparate
     * offsets (e.g. Neumann which mirrors any offset outside the grid ad the 
     * grid boundary). All is well.
     */
    auto result = Eigen::SparseMatrix<Scalar_t, ALIGNMENT, Index_t>(matrixHeight, matrixHeight);
    // TODO: What about makeCompressed()? Is it required when we construct via triplets?
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
  }
};

} // namespace matrixgen::implementation

namespace matrixgen
{

/**
 * Dispatcher for `adjmat` (workaround for a function template partial
 * specialization). See above implementation for details.
 */
template <
  typename OutMatrix_t = Eigen::SparseMatrix<double, Eigen::RowMajor>,
  typename StencilFn_t = void,
  typename WeightFn_t = void,
  typename Index_t = int
    >
OutMatrix_t
adjmat(
    const implementation::DiscreteCoords3d_t<Index_t>& gridDimensions,
    StencilFn_t stencilfn,
    WeightFn_t weightfn) {

  // TODO: Assert stencilfn's signature to avoid funny compilation errors.
  // TODO: Assert weightfn's signature to avoid funny compilation errors.

  return implementation::Adjmat<OutMatrix_t, StencilFn_t, WeightFn_t, Index_t>::
          invoke(gridDimensions, stencilfn, weightfn);
}

} // namespace matrixgen
