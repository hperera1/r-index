#ifndef R_INDEX_R_A_TREE_
#define R_INDEX_R_A_TREE_

#include <definitions.hpp>
#include <r_index.hpp>

using namespace sdsl;

namespace ri {
class rads_tree {
public:
  std::vector<std::pair<long long int, long long int>> tree; // nodes are pairs that represent: (edge cost, edge threshold).

  rads_tree(){};
  rads_tree(std::vector<uint> &cycle, std::vector<std::pair<ulint, ulint>> &bounds){
    size_t n = cycle.size(); // path size
    tree = std::vector<std::pair<long long int, long long int>>((size_t)1<<(size_t)(ceil(log2(n)) + 1), std::make_pair(0,0)); // tree size initialization.
    constructor_helper(cycle, bounds, 0, 0, n-1);
  }

  // recursively constructs a balanced binary tree from the leaves up to the root.
  void constructor_helper(std::vector<uint> &cycle, std::vector<std::pair<ulint, ulint>> &bounds, size_t node, size_t begin, size_t end) {
    size_t mid = (begin + end)/2;
    if(begin >= end) { // condition which creates the leaf nodes.
      tree[node].first = (long long int)bounds[begin].first; // edge cost
      tree[node].second = (long long int)bounds[begin].second; // edge threshold
      return;
    }

    constructor_helper(cycle, bounds, 2*node+1, begin, mid); // left child
    constructor_helper(cycle, bounds, 2*node+2, mid+1, end); // right child
    tree[node].first = tree[2*node+1].first + tree[2*node+2].first; // edge sum of left and right child
    tree[node].second = std::min(tree[2*node+1].second, tree[2*node+2].second) - tree[2*node+1].first; // edge threshold 
  }
};
}

#endif
