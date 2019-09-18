/**
 * Example displaying the inner workings of the adjacency function abstraction
 * by implementing a more complex example.
 */
#include <array>
#include <iostream>
#include <utility>

#include <matrixgen/core>

#include <Eigen/Sparse>

int main() {

  using DiscreteCoords3d_t = std::array<int32_t, 3>;
  const auto grid = std::array {5, 5, 1};

  /**
   * This example's adjacency function encodes a connectivity pattern in which
   * every node has a connection to the grid's center node. Thus
   * each node has a single "neighbor" whose offset is different for each
   * node, as opposed to the previous example in which the only offset was
   * (0, 0, 0).
   *
   * We are hence required to implement logic determining the node we're
   * processing in order to generate the correct offset. This is exemplary to
   * the general case in which there are several classes of nodes whose
   * adjacency patterns differ (e.g. the symmetric 7p stencil, which exhibits
   * its full adjacency pattern only for inner grid nodes, while any nodes
   * at the surfaces, edges and corners must be treated separately).
   *
   * By choice of the adjacency pattern for this example the logic is kept
   * very simple. The single offset we'll return is determined by determining
   * the distance between the node's central node and our node in question. No
   * further dispatching logic is required to distinguish between different
   * classes of nodes.
   *
   * We implement the adjacency function as a mutable lambda storing the range
   * of offset as a data member, aswell as the grid's central node's coords.
   */
  auto adjfn = [offsets = std::array<DiscreteCoords3d_t, 1> {{0, 0, 0}},
                gridCenter = DiscreteCoords3d_t {2, 2, 0}]
    (const DiscreteCoords3d_t& nodeCoords) mutable {

    offsets.front()[0] = gridCenter[0] - nodeCoords[0];
    offsets.front()[1] = gridCenter[1] - nodeCoords[1];
    offsets.front()[2] = gridCenter[2] - nodeCoords[2];
    return std::pair {offsets.begin(), offsets.end()};
  };

  const auto mat = matrixgen::adjmat(grid, adjfn, matrixgen::constweight(7));
  std::cout << std::endl << Eigen::MatrixXd(mat) << std::endl;

  /**
   * You are encouraged to have a look at the implementation for the symmetric
   * 7p stencil preset `matrixgen::stencil7p`, which was used in a previous
   * example. It details the usage of a second possible signature of the
   * adjacency function which passes the grid's size a second argument in
   * order to create generic adjacency functions for which it's not possible
   * to hard-code information about the grid into the implementation.
   */
}
