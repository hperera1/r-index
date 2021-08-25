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
  rads_tree(std::vector<uint> &cycle, std::vector<std::pair<ulint, ulint>> &bounds, uint tree_num, std::vector<std::tuple<ulint, ulint, uint>> &tree_pointers){
    size_t n = cycle.size(); // path size
    tree = std::vector<std::pair<long long int, long long int>>((size_t)1<<(size_t)(ceil(log2(n)) + 1), std::make_pair(0,0)); // tree size initialization.
    constructor_helper(cycle, bounds, 0, 0, n-1, tree_num, tree_pointers);
  }

  // so many arguments...
  // recursively constructs a balanced binary tree from the leaves up to the root.
  void constructor_helper(std::vector<uint> &cycle, std::vector<std::pair<ulint, ulint>> &bounds, size_t node, size_t begin, size_t end, uint tree_num, std::vector<std::tuple<ulint, ulint, uint>> &tree_pointers) {
    size_t mid = (begin + end)/2;
    if(begin >= end) { // condition which creates the leaf nodes.
      tree[node].first = (long long int)bounds[begin].first; // edge cost
      tree[node].second = (long long int)bounds[begin].second; // edge threshold
      tree_pointers.push_back(std::make_tuple(cycle[begin], tree_num, node));
      return;
    }

    constructor_helper(cycle, bounds, 2*node+1, begin, mid, tree_num, tree_pointers); // left child
    constructor_helper(cycle, bounds, 2*node+2, mid+1, end, tree_num, tree_pointers); // right child
    tree[node].first = tree[2*node+1].first + tree[2*node+2].first; // edge sum of left and right child
    tree[node].second = std::min(tree[2*node+1].second, tree[2*node+2].second) - tree[2*node+1].first; // edge threshold
  }

  // args: leaf_pos, cost, d | returns: should return a pair methinks
  //                           (sa', d)
  std::pair<ulint, ulint> query(ulint leaf_pos, uint cumm_cost, uint d) {
    // in short we go up ^
    // check if we can keep going ->
    // if not, we go down v
    // first to the left and then to the right
    // if were coming from the left go right
    // if were at the end of the left tree from the root, we can skip
    // to the next node with some clever bit tricks :)

    // this is to go up until we can't no more
    // cond 1: cumm_cost must be within current threshold
    // cond 2: parent threshold must be positive
    // cond 3: cumm_cost + cost must be within parent threshold
    // only then can we go up
    while((cumm_cost <= tree[leaf_pos].second) && (tree[leaf_pos>>1].second >= 0) && (cumm_cost + tree[leaf_pos].first <= tree[leaf_pos>>1].second)) { // do we make this a function? probably
      cumm_cost += tree[leaf_pos].first;
      leaf_pos = leaf_pos>>1;
    }



    // at this point we've gone up, now we need to check if we can go down
    // going down means what? we start going to the right and asking different questions
    // while we go down

    return (std::make_pair(leaf_pos, d));
  }
};
}

#endif
