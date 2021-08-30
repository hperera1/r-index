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
    tree = std::vector<std::pair<long long int, long long int>>((size_t)1<<(size_t)(ceil(log2(n)) + 1), std::make_pair(-1,0)); // tree size initialization.
    constructor_helper(cycle, bounds, 1, 0, n-1, tree_num, tree_pointers);
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

    constructor_helper(cycle, bounds, 2*node, begin, mid, tree_num, tree_pointers); // left child
    constructor_helper(cycle, bounds, 2*node+1, mid+1, end, tree_num, tree_pointers); // right child
    tree[node].first = tree[2*node].first + tree[2*node+1].first; // edge sum of left and right child
    tree[node].second = std::min(tree[2*node].second, tree[2*node+1].second) - tree[2*node].first; // edge threshold
  }

  // args: leaf_pos, cost, d | returns: should return a pair methinks
  //                           (sa', d)
  std::pair<ulint, ulint> query(ulint leaf_pos, uint cost, uint d) {
    // begin climbing the tree and only collect when we go to a sibling.
    // when you go to a parent do not collect the money, that is if the parent has the same left most node.
    // if youre moving to a subtree where the leftmost node is not the same then we have to collect the money.
    climb(leaf_pos, cost, d);

    return (std::make_pair(leaf_pos, d));
  }

  void climb(ulint node_pos, uint cost, uint d) {
    // climbing conditions
    // 1. node is non-negative
    // 2. cost does not exceed nodes threshold
    if(node_pos>>(__builtin_ctzll(~node_pos)) == 2) {
      node_pos = 3;
      descend(node_pos, cost, d);
      return;
    }

    while((tree[node_pos].second > 0) && (cost <= tree[node_pos].second)) {
      node_pos = node_pos >> 1;
    }

    cost += tree[node_pos].first;
    node_pos = (node_pos << 1) + 1;
    descend(node_pos, cost, d);
  }

  void descend(ulint node_pos, uint cost, uint d) {
    // 1. check node_pos that is passed in
    // 2. if ✓ then check your immediate sibling
    // 3. if ✕ then check left child until you get a green light
    // while we are not at a leaf node
    while((node_pos < tree.size()) && (tree[node_pos].first >= 0)) {
      if(node_pos & 1 == 0) { // left node on go -> go to right sibling
        if((tree[node_pos].second > 0) && (cost <= tree[node_pos].second)) {
          node_pos += 1;
        }
        else {
          node_pos = node_pos << 1;
        }
      }
      else { // right node on go -> go to right child
        if((tree[node_pos].second > 0) && (cost <= tree[node_pos].second)) {
          node_pos = (node_pos << 1) + 1;
        }
        else {
          node_pos = node_pos << 1;
        }
      }
    }

    return;
  }
};
}

#endif
