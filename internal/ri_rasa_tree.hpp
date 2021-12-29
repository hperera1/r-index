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
//! Tree data structure that stores the samples contained within a cycle. Performs queries.
class rads_tree {
public:
  std::vector<std::pair<long long int, long long int>> tree; // nodes are pairs that represent: (edge cost, edge threshold).
  std::vector<ulint> leaf_samples; // actual sa sample values represented by run index
  sparse_bv_type leaf_node_bv; // bv telling us which node is a leaf node
  uint left_most_i; // index at which the left most leaf is stored || sdsl::serialize
  ulint height; // height of the tree starting at 0 || sdsl::serialize

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

  rads_tree(const rads_tree &other_tree) {
    this->tree = other_tree.tree;
    this->leaf_samples = other_tree.leaf_samples;
    this->leaf_node_bv = other_tree.leaf_node_bv;
    this->left_most_i = other_tree.left_most_i;
    this->height = other_tree.height;
  }

  rads_tree(rads_tree &&other_tree) noexcept
  : tree(move(other_tree.tree))
  , leaf_samples(move(other_tree.leaf_samples))
  , leaf_node_bv(move(other_tree.leaf_node_bv))
  , left_most_i(move(other_tree.left_most_i))
  , height(move(other_tree.height))
  {}

  // we need a copy and swap operator?
  rads_tree& operator=(const rads_tree &other_tree) {
    return *this = rads_tree(other_tree);
  }

  rads_tree& operator=(rads_tree &&other_tree) {
    swap(tree, other_tree.tree);
    swap(leaf_samples, other_tree.leaf_samples);
    swap(leaf_node_bv, other_tree.leaf_node_bv);
    swap(left_most_i, other_tree.left_most_i);
    swap(height, other_tree.height);
    return *this;
  }

  // recursively constructs a balanced binary tree from the leaves up to the root.
  void constructor_helper(std::vector<std::pair<ulint, ulint>> &bounds, size_t node, size_t begin, size_t end, uint tree_num, std::vector<std::tuple<ulint, ulint, uint>> &tree_pointers, std::vector<bool> &leaf_bv) {
    size_t mid = (begin + end)/2;
    if(begin >= end) { // condition which creates the leaf nodes.
      tree[node].first = (long long int)bounds[leaf_samples[begin]].first; // edge cost
      tree[node].second = (long long int)bounds[leaf_samples[begin]].second; // edge threshold
      if(begin != leaf_samples.size()-1) { // if it is not the last sample in the tree, add a tree pointer to be used when querying
        tree_pointers.push_back(std::make_tuple(leaf_samples[begin], tree_num, node));
      }
      leaf_bv[node] = true;
      return;
    }

    constructor_helper(bounds, 2*node, begin, mid, tree_num, tree_pointers, leaf_bv); // construct left child
    constructor_helper(bounds, 2*node+1, mid+1, end, tree_num, tree_pointers, leaf_bv); // construct right child
    tree[node].first = tree[2*node].first + tree[2*node+1].first; // edge sum of left and right child
    tree[node].second = std::min(tree[2*node].second, tree[2*node+1].second - tree[2*node].first); // edge threshold
  }

  // args: leaf_pos, cost, d | returns: returns new sample index and distance travelled (sa', d)
  std::tuple<ulint, ulint, ulint> query(ulint leaf_pos, uint cost, uint d) {
    // begin climbing the tree and only collect when we go to a sibling.
    // when you go to a parent do not collect the cost, that is if the parent has the same left most node.
    // if youre moving to a subtree where the leftmost node is not the same then we have to collect the money.
    // cout << "\nleaf_pos: " << leaf_pos << ", cost: " << cost << ", d: " << d << endl;
    return (climb(leaf_pos, cost, d));
  }

  // climbing conditions
  // 1. node is non-negative
  // 2. cost does not exceed nodes threshold
  std::tuple<ulint, ulint, ulint> climb(ulint node_pos, uint &cost, uint d) {
    // this provides us with the index of the sample that is being used to query
    // left_most_i is 0 and we get the distance from our start_pos
    // at the end we do leaf_index + the distance that was travelled to get the new sample
    // eg. start_pos = 21 -> leaf_index = 0; distance travelled = 2; leaf_index + distance travelled = 2; new_sample = leaf_samples[2] (11)
    // cout << "Climbing ..." << endl;
    ulint current_height = height; // we always enter at a leaf node so we are at level with the height of the tree
    ulint start_pos = node_pos;
    ulint prev_pos = node_pos;

    if((node_pos << 1) < tree.size()) { // this means that we are actually at height - 1 because we entered from a leaf node one height up.
      // cout << "height adjusted." << endl;
      current_height -= 1;
    }

    ulint current_leftmost_node = (node_pos << (height - current_height));
    int max_d_travelled = 0;
    int last_bit = 0;
    while(((node_pos != 1) && (tree[node_pos].second > 0) && (cost < tree[node_pos].second)) || (start_pos == node_pos)) { // climb up while the upper_bounds let us
      int distance_diff = (int) d - max_d_travelled;
      last_bit = (node_pos & 1);
      if(distance_diff <= 0) {
        // this means that the sample we are looking for is in the current interval that we are covering
        // if this is true then we should start descending from this node
        descend(start_pos, node_pos, cost, d, current_height);
        ulint distance = calculate_d(start_pos, node_pos); // calculate the distance from start to end
        ulint leaf_sample_index = calculate_d(left_most_i, node_pos);
        return (std::make_tuple(leaf_samples[leaf_sample_index], distance, cost)); // return new sample and new distance
      }

      // check if when we shift up we get to the root of the left subtree
      // meaning we can skip forward to the beginning of the right subtree and descend
      prev_pos = node_pos;
      if(node_pos>>(__builtin_ctzl(~node_pos)) == 2) {
        // cout << "\nskipping to right subtree." << endl;
        node_pos = 3;
        current_height = 1;
        cost += tree[prev_pos].first;
        descend(start_pos, node_pos, cost, d, current_height);
        ulint distance = calculate_d(start_pos, node_pos);
        ulint leaf_sample_index = calculate_d(left_most_i, node_pos);
        return (std::make_tuple(leaf_samples[leaf_sample_index], distance, cost)); // return new sample and new distance
      }

      if(last_bit == 0) {
        // cout << "left node up." << endl;
        node_pos = node_pos >> 1;
      }
      else {
        // cout << "right node up and right." << endl;
        if(node_pos != start_pos) { // check parent node if we are no longer at our start position
          if((tree[node_pos >> 1].second > 0) && (cost < tree[node_pos >> 1].second)) {
            node_pos = (node_pos >> 1) + 1;
          }
          else {
            break; // get out cause youre gonna make an illegal jump!
          }
        }
        else {
          node_pos = (node_pos >> 1) + 1;
        }
      }

      current_height -= 1;
      if((node_pos << (height - current_height)) != current_leftmost_node) {
        cost += tree[prev_pos].first;
        current_leftmost_node = (node_pos << (height - current_height));
      }

      ulint max_node_pos = ((node_pos << (height - current_height)) + ((1 << (height - current_height)) - 1));
      ulint max_d = calculate_d(start_pos, max_node_pos); // (node_pos gets shifted by the distance to the leaves) + (1 shifted that many times left - 1)
      max_d_travelled = (int) max_d;
    }

    // cout << "done climbing." << endl;
    // if we climbed up from a right node, we might encounter a node that denies us,
    // and the nodes below that might also deny us, therefore we shouldn't collect.
    // this descends until we have found a node that lets us pass and we collect that.
    last_bit = (prev_pos & 1);
    if(last_bit == 1) {
      while(!leaf_node_bv[node_pos]) {
        if(!(tree[node_pos].second > 0) || !(cost < tree[node_pos].second)) { // this means we keep getting blocked
          node_pos = (node_pos << 1);
          current_height += 1;
        }
        else {
          descend(start_pos, node_pos, cost, d, current_height); // descend using new node_pos, cost, and d
          break;
        }
      }
    }
    else {
      cost += tree[node_pos << 1].first; // add cost of node that denies us
      node_pos = (node_pos << 1) + 1; // move node_pos to the right child of the current node // or should it be the left child of the node
      current_height += 1;
      descend(start_pos, node_pos, cost, d, current_height); // descend using new node_pos, cost, and d
    }

    ulint distance = calculate_d(start_pos, node_pos); // calculate the distance from start to end
    ulint leaf_sample_index = calculate_d(left_most_i, node_pos);
    return (std::make_tuple(leaf_samples[leaf_sample_index], distance, cost)); // return new sample and new distance
  }

  // during the descent the node_pos gets shifted around so that calculate_d can do its job
  // args: node_pos: current node_pos, cost: cumulated cost, d: distance of where we came in from to the sample we're looking for.
  void descend(ulint start_pos, ulint &node_pos, uint &cost, uint d, ulint current_height) {
    ulint prev_pos = node_pos;
    while((node_pos < leaf_node_bv.size()) && !leaf_node_bv[node_pos]) { // while we are in the bounds of the tree and not at a leaf node
      // check min_d_travelled before descending to the next node
      // if its < 0 then that means we have moved past the interval our sample is contained in. if this is the case go to our left sibling
      // cout << "testing end_pos: " << (node_pos << (height - current_height)) << endl;
      // cout << "new node pos: " << node_pos << endl;
      int min_d_travelled;
      int last_bit = (node_pos & 1);
      prev_pos = node_pos;
      if(last_bit == 0) { // this means we are at a left child node
        // cout << "left node." << endl;
        if((tree[node_pos].second > 0) && (cost < tree[node_pos].second)) { // if good, go to right sibling
          // cout << "go to right sibling." << endl;
          node_pos += 1;
          ulint min_node = (node_pos << (height - current_height));
          min_d_travelled = (int) calculate_d(start_pos, min_node);
          // this means we went too far when going to right sibling and we go to the right child instead
          if(((int) d - min_d_travelled) < 0) {
            node_pos -= 1;
            node_pos = node_pos << 1;
            node_pos += 1;
            current_height += 1;
            cost += tree[node_pos - 1].first;
          }
          else {
            cost += tree[prev_pos].first;
          }
        }
        else {
          // cout << "go to left child." << endl;
          current_height += 1;
          node_pos = node_pos << 1; // if not good then go to left child and check on next iteration
        }
      }
      else { // right node
        // cout << "right node." << endl;
        if((tree[node_pos].second > 0) && (cost < tree[node_pos].second)) { // if good -> go to right child
          // cout << "go to right child." << endl;
          ulint min_node = (node_pos << (height - current_height)); // this is how you can get the left mode node of the current subtree
          min_d_travelled = (int) calculate_d(start_pos, min_node);
          if(((int) d - min_d_travelled) < 0) {
            // right node doesn't have sample. go to left child
            node_pos -= 1;
            continue;
          }

          current_height += 1;
          cost += tree[node_pos << 1].first;
          node_pos = (node_pos << 1) + 1;
        }
        else {
          // cout << "go to left child." << endl;
          current_height += 1;
          node_pos = node_pos << 1; // if not good then go to left child and check that guy on the next iterations
        }
      }
    }

    ulint min_node = (node_pos << (height - current_height));
    int min_d_travelled = (int) calculate_d(start_pos, min_node);
    if(((int) d - min_d_travelled) < 0) {
      // cout << "we go back one again." << endl;
      node_pos -= 1;
      cost -= tree[node_pos].first;
      descend(start_pos, node_pos, cost, d, current_height);
    }
  }

  void left_climb() {}

  void left_descend() {}

  void right_climb() {}

  void right_descend() {}

  // calculates the distance between leaf nodes: start_pos and end_pos. distance meaning the number of leaf nodes between these two indices.
  ulint calculate_d(ulint start_pos, ulint end_pos) {
    ulint d = 0;
    // cout << "start: " << start_pos << endl;
    // cout << "end: " << end_pos << endl;

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

    // cout << "calc_d d: " << d << endl;
    return d;
  }

  void print_tree_info() {
    cout << "tree size: " << tree.size() << endl;
    cout << "# of samples: " << leaf_samples.size() << endl;
    print_bounds();
    print_samples();
    print_nodebv();
    cout << "left_most_i: " << left_most_i << endl;
    cout << "height: " << height << endl;
  }

  void print_bounds() {
    cout << "tree bounds: ";
    for (size_t i = 0; i < tree.size(); i++) {
      cout << "[" << tree[i].first << "," << tree[i].second << "] ";
    }
    cout << endl;
  }

  void print_samples() {
    cout << "tree samples: ";
    for(size_t i = 0; i < leaf_samples.size(); i++) {
      cout << "[" << leaf_samples[i] << "] ";
    }
    cout << endl;
  }

  void print_nodebv() {
    cout << "node bv: ";
    for(size_t i = 0; i < leaf_node_bv.size(); i++) {
      cout << "[" << leaf_node_bv[i] << "] ";
    }
    cout << endl;
  }

  // serialize tree, leaf_samples, leaf_node_bv, left_most_i, height
  ulint serialize(std::ostream& out) {
    ulint w_bytes = 0;

    w_bytes += sdsl::serialize(tree.size(), out);
    out.write((char*)tree.data(), tree.size()*sizeof(tree[0]));
    w_bytes += sizeof(tree[0])*tree.size();

    w_bytes += sdsl::serialize(leaf_samples.size(), out);
    out.write((char*)leaf_samples.data(), leaf_samples.size()*sizeof(leaf_samples[0]));
    w_bytes += sizeof(leaf_samples[0])*leaf_samples.size();

    w_bytes += leaf_node_bv.serialize(out);
    w_bytes += sdsl::serialize(left_most_i, out);
    w_bytes += sdsl::serialize(height, out);

    return w_bytes;
  }

  // load tree, leaf_samples, leaf_node_bv, left_most_i, height
  void load(std::istream& in) {
    size_t temp_size;

    in.read((char*)&temp_size, sizeof(temp_size));
    tree.resize(temp_size);
    in.read((char*)tree.data(), temp_size*sizeof(tree[0]));

    in.read((char*)&temp_size, sizeof(temp_size));
    leaf_samples.resize(temp_size);
    in.read((char*)leaf_samples.data(), temp_size*sizeof(leaf_samples[0]));

    leaf_node_bv.load(in);

    in.read((char*)&left_most_i, sizeof(left_most_i));
    in.read((char*)&height, sizeof(height));
  }
};
}

#endif
