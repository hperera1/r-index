#ifndef R_INDEX_R_A_TREE_
#define R_INDEX_R_A_TREE_

#include <definitions.hpp>
#include <r_index.hpp>
#include "sparse_sd_vector.hpp"
#include "sparse_hyb_vector.hpp"

using namespace sdsl;
namespace ri {
  template<class sparse_bv_type = sparse_sd_vector,
           class rle_string_t = rle_string_sd>

class rads_tree {
public:
  std::vector<std::pair<long long int, long long int>> tree; // nodes are pairs that represent: (edge cost, edge threshold).
  std::vector<ulint> leaf_samples; // actual sa sample values represented by run index
  sparse_bv_type leaf_node_bv; // bv telling us which node is a leaf node
  uint left_most_i;
  ulint height;

  rads_tree(){};
  rads_tree(std::vector<ulint> &cycle, std::vector<std::pair<ulint, ulint>> &bounds, uint tree_num, std::vector<std::tuple<ulint, ulint, uint>> &tree_pointers) {
    size_t n = cycle.size(); // path size
    height = ceil(log2(n)); // height of the tree
    leaf_samples = cycle; // store the samples in the cycle
    tree = std::vector<std::pair<long long int, long long int>>((size_t)1<<(size_t)(height + 1), std::make_pair(-1, 0)); // tree size initialization, default the costs to -1
    auto temp_leaf_bv = vector<bool>(tree.size(), false);

    constructor_helper(bounds, 1, 0, n - 1, tree_num, tree_pointers, temp_leaf_bv);
    leaf_node_bv = sparse_bv_type(temp_leaf_bv);
    left_most_i = (tree.size() >> 1);
  }

  rads_tree(const rads_tree &rads_tree_) {
    this->tree = rads_tree_.tree;
    this->leaf_samples = rads_tree_.leaf_samples;
    this->leaf_node_bv = rads_tree_.leaf_node_bv;
    this->left_most_i = rads_tree_.left_most_i;
    this->height = rads_tree_.height;
  }

  rads_tree(rads_tree &&rads_tree_) {
    this.tree = rads_tree_.tree;
    this.leaf_samples = rads_tree_.leaf_samples;
    this.leaf_node_bv = rads_tree_.leaf_node_bv;
    this.left_most_i = rads_tree_.left_most_i;
    this.height = rads_tree_.height;
    rads_tree_.tree = NULL;
    rads_tree_.leaf_samples = NULL;
    rads_tree_.leaf_node_bv = NULL;
    rads_tree_.left_most_i = NULL;
    rads_tree_.height = NULL;
  }

  // so many arguments...
  // recursively constructs a balanced binary tree from the leaves up to the root.
  void constructor_helper(std::vector<std::pair<ulint, ulint>> &bounds, size_t node, size_t begin, size_t end, uint tree_num, std::vector<std::tuple<ulint, ulint, uint>> &tree_pointers, std::vector<bool> &leaf_bv) {
    size_t mid = (begin + end)/2;
    if(begin >= end) { // condition which creates the leaf nodes.
      tree[node].first = (long long int)bounds[begin].first; // edge cost
      tree[node].second = (long long int)bounds[begin].second; // edge threshold
      tree_pointers.push_back(std::make_tuple(leaf_samples[begin], tree_num, node));
      leaf_bv[node] = true;
      return;
    }

    constructor_helper(bounds, 2*node, begin, mid, tree_num, tree_pointers, leaf_bv); // construct left child
    constructor_helper(bounds, 2*node+1, mid+1, end, tree_num, tree_pointers, leaf_bv); // construct right child
    tree[node].first = tree[2*node].first + tree[2*node+1].first; // edge sum of left and right child
    tree[node].second = std::min(tree[2*node].second, tree[2*node+1].second - tree[2*node].first); // edge threshold
  }

  // args: leaf_pos, cost, d | returns: returns new sample and distance left
  //                           (sa', d)
  std::pair<ulint, ulint> query(ulint leaf_pos, uint cost, uint d) {
    // begin climbing the tree and only collect when we go to a sibling.
    // when you go to a parent do not collect the money, that is if the parent has the same left most node.
    // if youre moving to a subtree where the leftmost node is not the same then we have to collect the money.
    cout << "Querying tree ..." << endl;
    cout << "leaf_pos: " << leaf_pos << ", cost: " << cost << ", d: " << d << endl;
    return (climb(leaf_pos, cost, d));
  }

  // climbing conditions
  // 1. node is non-negative
  // 2. cost does not exceed nodes threshold
  std::pair<ulint, ulint> climb(ulint &node_pos, uint &cost, uint d) {
    // this provides us with the index of the sample that is being used to query
    // left_most_i is 0 and we get the distance from our start_pos
    // at the end we do leaf_index + the distance that was travelled to get the new sample
    // eg. start_pos = 21 -> leaf_index = 0; distance travelled = 2; leaf_index + distance travelled = 2; new_sample = leaf_samples[2] (11)
    cout << "Climbing ..." << endl;
    ulint current_height = height; // we always enter at a leaf node so we are at level with the height of the tree
    ulint start_pos = node_pos;

    // check if when we shift up we get to the root of the left subtree
    // meaning we can skip forward to the beginning of the right subtree and descend
    if(node_pos>>(__builtin_ctzl(~node_pos)) == 2) {
      cout << "\nSkipping to the right subtree and desceding ..." << endl;
      node_pos = 3;
      current_height = 1;
      descend(start_pos, node_pos, cost, d, current_height);
      ulint distance = calculate_d(start_pos, node_pos);

      cout << "leaf sample: " << leaf_samples[start_pos + distance] + cost << endl;
      cout << "distance left: " << d - distance << endl;

      return (std::make_pair(leaf_samples[start_pos + distance] + cost, d - distance)); // return new sample and new distance
    }

    ulint shifts = 0;
    int max_d_travelled = 0;
    cout << "node pos: " << node_pos << endl;
    while((node_pos != 1) && (tree[node_pos].second > 0) && (cost <= tree[node_pos].second)) { // climb up while the upper_bounds let us
      int distance_diff = (int) d - max_d_travelled;
      cout << "distance_diff: " << distance_diff << endl;
      if(distance_diff <= 0) {
        // this means that the sample we are looking for is in the current interval that we are covering
        // if this is true then we should start descending from this node
        cout << "max_d check true" << endl;
        cost += tree[node_pos << 1].first; // add cost of node that denies us
        node_pos = (node_pos << 1) + 1; // move node_pos to the right child of the current node
        current_height += 1;
        descend(start_pos, node_pos, cost, d, current_height);
        ulint distance = calculate_d(start_pos, node_pos); // calculate the distance from start to end
        return (std::make_pair(leaf_samples[start_pos + distance] + cost, d - distance)); // return new sample and new distance
      }

      node_pos = node_pos >> 1;
      cout << "\nnew node pos: " << node_pos << endl;
      current_height -= 1;
      max_d_travelled = (int) calculate_d(start_pos, ((node_pos << (height - current_height)) + ((1 << (height - current_height)) - 1))); // (node_pos gets shifted by the distance to the leaves) + (1 shifted that many times left - 1)
      cout << "max_d_travelled: " << max_d_travelled << endl;
    }

    cout << "\nwe're at the root or we've been denied by a bound." << endl;
    current_height += 1;
    cost += tree[node_pos << 1].first; // add cost of node that denies us
    node_pos = (node_pos << 1) + 1; // move node_pos to the right child of the current node
    descend(start_pos, node_pos, cost, d, current_height); // descend using new node_pos, cost, and d
    ulint distance = calculate_d(start_pos, node_pos); // calculate the distance from start to end

    // this doesn't work. probably simple fix cause i'm tired and can't tell
    cout << "leaf sample: " << leaf_samples[start_pos + distance] + cost << endl;
    cout << "distance left: " << d - distance << endl;

    return (std::make_pair(leaf_samples[start_pos + distance] + cost, d - distance)); // return new sample and new distance
  }

  // during the descent the node_pos gets shifted around so that calculate_d can do its job
  // args: node_pos: current node_pos, cost: cumulated cost, d: distance of where we came in from to the sample we're looking for
  void descend(ulint start_pos, ulint &node_pos, uint &cost, uint d, ulint current_height) {
    // 1. check node_pos that is passed in
    // 2. if YAY then check your immediate sibling
    // 3. if NAY then check left child until you get a green light
    // while((node_pos < tree.size()) && (tree[node_pos].first >= 0)) { // 1. check if we are out of bounds 2. check if there is even a node here
    cout << "\nDescending ..." << endl;
    cout << "start pos: " << start_pos << endl;
    cout << "node pos: " << node_pos << endl;
    cout << "d: " << d << endl;
    cout << "current height: " << current_height << endl;

    while((node_pos < leaf_node_bv.size()) && !leaf_node_bv[node_pos]) { // while we are not at a leaf node
      // check min_d_travelled before descending to the next node
      // if its < 0 then that means we have moved past the interval our sample is contained in. if this is the case go to our left sibling
      cout << "testing end_pos: " << (node_pos << (height - current_height)) << endl;
      cout << "new node pos: " << node_pos << endl;
      int min_d_travelled;
      int last_bit = (node_pos & 1);
      cout << "last_bit: " << last_bit << endl;

      if(last_bit == 0) { // this means we are at a left child node
        cout << "left node." << endl;
        if((tree[node_pos].second > 0) && (cost <= tree[node_pos].second)) { // if good -> go to right sibling
          cout << "go to right sibling." << endl;
          node_pos += 1;
          min_d_travelled = (int) calculate_d(start_pos, (node_pos << (height - current_height)));
          if(((int) d - min_d_travelled) < 0) {
            cout << "we've gone too far. go to left sib right child." << endl;
            node_pos -= 1;
            node_pos = node_pos << 1;
            node_pos += 1;
            current_height += 1;
          }
          else {
            cost += tree[node_pos].first;
            node_pos += 1;
          }
        }
        else {
          cout << "go to left child." << endl;
          current_height += 1;
          node_pos = node_pos << 1; // if not good then go to left child and check on next iteration
        }
      }
      else { // right node
        cout << "right node." << endl;
        if((tree[node_pos].second > 0) && (cost <= tree[node_pos].second)) { // if good -> go to right child
          cout << "go to right child." << endl;
          current_height += 1;
          node_pos = (node_pos << 1) + 1;
          min_d_travelled = (int) calculate_d(start_pos, (node_pos << (height - current_height)));
          if(((int) d - min_d_travelled) < 0) {
            cout << "right child doesn't have sample. go to left child." << endl;
            node_pos -= 1;
          }
        }
        else {
          cout << "go to left child." << endl;
          current_height += 1;
          node_pos = node_pos << 1; // if not good then go to left child and check that guy on the next iterations
        }
      }
    }

    int min_d_travelled = (int) calculate_d(start_pos, (node_pos << (height - current_height)));
    if(((int) d - min_d_travelled) < 0) {
      node_pos -= 1;
    }

    cout << "node pos after shifts: " << node_pos << endl;
  }

  // calculates distance between leaf nodes start_pos and end_pos
  ulint calculate_d(ulint start_pos, ulint end_pos) {
    ulint d = 0;

    if(start_pos >= left_most_i) { // starting from bottom layer
      if(end_pos >= left_most_i) { // going to bottom layer bot->bot
        d += leaf_node_bv.rank(end_pos);
        d -= leaf_node_bv.rank(start_pos);
        d -= leaf_node_bv.rank((start_pos >> 1));
        d += leaf_node_bv.rank((end_pos >> 1));
      }
      else { // going to top layer bot->top
        d += leaf_node_bv.rank(end_pos);
        d -= leaf_node_bv.rank(start_pos);
        d -= leaf_node_bv.rank((start_pos >> 1));
        d += leaf_node_bv.rank((end_pos << 1));
      }
    }
    else { // starting from top layer
      if(end_pos >= left_most_i) {// going to bottom layer top->bot
        d += leaf_node_bv.rank(end_pos);
        d -= leaf_node_bv.rank(start_pos);
        d -= leaf_node_bv.rank((start_pos << 1));
        d += leaf_node_bv.rank((end_pos >> 1));
      }
      else { // going to top layer top->top
        d += leaf_node_bv.rank(end_pos);
        d -= leaf_node_bv.rank(start_pos);
        d -= leaf_node_bv.rank((start_pos << 1));
        d += leaf_node_bv.rank((end_pos << 1));
      }
    }

    cout << "calc_d d: " << d << endl;
    return d;
  }

  void print_array() {
    cout << "Printing array ..." << endl;
    for (size_t i = 0; i < tree.size(); i++) {
      cout << "(" << tree[i].first << ", " << tree[i].second << ") ";
    }

    cout << endl;
    int leaf_counter = 0;
    for (size_t i = 0; i < leaf_node_bv.size(); i++) {
      if(leaf_node_bv[i]) {
        cout << " 1,";
        leaf_counter++;
      }
      else {
        cout << " -,";
      }
    }
    cout << endl;
  }
};
}

#endif
