/**
 * Example displaying the inner workings of the stencil abstraction by
 * implementing a simple 1-p stencil.
 */
#include <array>
#include <iostream>
#include <utility>

#include <matrixgen/core>

int main() {

  using DiscreteCoords3d_t = std::array<int32_t, 3>;
  const auto grid = std::array {4, 4, 1};

  /**
   * The grid is traversed one node at a time in x-then-y-then-z direction,
   * starting at the grid's origin (0, 0, 0). For every node in the grid the
   * adjacency function is called on the node's  coordinate triplet (x, y, z).
   * The adjancency function shall return a range of discrete offsets, which
   * encode the neighbors for the node in question. The range shall be
   * represented by a std::pair of iterators.
   *
   * For every offset in the range the neighbor's node index j is computed,
   * which, in addition to the index i of our node being processed, comprises
   * the coordinates (i, j), where the adjacency matrix will contain a nonzero
   * entry. The numerical value of that entry is determined by the weight
   * function, which is completely independent of the adjacency function (see
   * the previous example).
   *
   * In this trivial example we simply return the begin and end iterators of
   * a single-element array which contains the null-offset, hence the matrix
   * will be diagonal. There's no computation or logic required inside the
   * adjacency function since the range of offsets is identical for each and
   * every node in the grid. This is quite different to adjacency functions
   * where the adjacency pattern changes at different regions of the grid, as
   * is the case for a stencil, which has to obey boundary conditions. A more
   * complex example will be the topic of the next example.
   */
  const auto offsets = std::array<DiscreteCoords3d_t, 1> {{0, 0, 0}};
  auto adjfn = [&offsets](
      const DiscreteCoords3d_t& /* node */,
      const DiscreteCoords3d_t& /* gridDimensions */) { // TODO: add a second signature without the grid
    return std::pair {offsets.cbegin(), offsets.cend()};
  };

  const auto mat = matrixgen::adjmat(grid, adjfn, matrixgen::constweight(5));
  std::cout << std::endl << mat << std::endl;
}
