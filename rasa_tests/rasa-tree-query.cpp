#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <cstring>
#include <string>

#include "../internal/rle_string.hpp"
#include "../internal/r_index.hpp"
#include "../internal/ri_rasa_tree.hpp"
#include "../internal/utils.hpp"
#include "../internal/sparse_sd_vector.hpp"
#include "../internal/sparse_hyb_vector.hpp"

using namespace sdsl;
using namespace ri;
using namespace std;

// testing the query function for the rads_tree class.
// create several trees.
int main() {
  std::vector<rads_tree<>> trees;
  std::vector<std::tuple<ulint, ulint, uint>> tree_pointers;

  std::vector<ulint> tree_0 = {0, 1, 2, 3, 4, 5, 6, 7}; // basic tree with 8 elements to test traversal unhindered by bounds.
  std::vector<std::pair<ulint, ulint>> tree_0_bounds = {{0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}};

  std::vector<ulint> tree_1 = {26, 3, 11, 20, 0, 17}; // example with actual bounds.
  std::vector<std::pair<ulint, ulint>> tree_1_bounds = {{5,1}, {0,6}, {0,3}, {2,1}, {0,3}, {4,1}};

  std::vector<ulint> tree_2 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}; // larger example with 14 samples to test traversal again.
  std::vector<std::pair<ulint, ulint>> tree_2_bounds = {{0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}};

  std::vector<ulint> tree_3 = {0, 1, 2, 3, 4, 5, 6, 7, 8}; // 9 samples to test the new rung of leaves that get made.
  std::vector<std::pair<ulint, ulint>> tree_3_bounds = {{0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}};

  std::vector<ulint> tree_4; // a lot of samples to push it.
  std::vector<std::pair<ulint, ulint>> tree_4_bounds;
  for(size_t i = 0; i <= 123; i++) {
    tree_4.push_back(i);
    tree_4_bounds.push_back(std::make_pair(0,1));
  }

  std::vector<ulint> tree_5 = {26, 13, 3, 11, 16, 20, 0, 17, 1, 2, 19, 24};
  std::vector<std::pair<ulint, ulint>> tree_5_bounds = {{5,1}, {0,3}, {0,6}, {0,3}, {3,0}, {2,1}, {0,3}, {4,1}, {0,3}, {0,3}, {0,3}, {0,3}};

  trees.push_back(rads_tree(tree_0, tree_0_bounds, 0, tree_pointers));
  trees.push_back(rads_tree(tree_1, tree_1_bounds, 1, tree_pointers));
  trees.push_back(rads_tree(tree_2, tree_2_bounds, 2, tree_pointers));
  trees.push_back(rads_tree(tree_3, tree_3_bounds, 3, tree_pointers));
  trees.push_back(rads_tree(tree_4, tree_4_bounds, 4, tree_pointers));
  trees.push_back(rads_tree(tree_5, tree_5_bounds, 5, tree_pointers));

  std::tuple<ulint, ulint, ulint> sample_and_delta;
  // retrieve every sample using the left most sample.
  for(size_t i = 1; i <= 6; i++) {
    sample_and_delta = trees[0].query(trees[0].left_most_i, 0, i);
    assert(std::get<0>(sample_and_delta) == i);
  }

  // couple of tests to see if we get the samples that we expect from the tree with bounds.
  sample_and_delta = trees[1].query(8, 0, 1);
  assert((std::get<0>(sample_and_delta) + std::get<2>(sample_and_delta)) == 3);

  sample_and_delta = trees[1].query(8, 0, 2);
  assert((std::get<0>(sample_and_delta) + std::get<2>(sample_and_delta)) == 13);

  sample_and_delta = trees[1].query(9, 0, 1);
  assert((std::get<0>(sample_and_delta) + std::get<2>(sample_and_delta)) == 13);

  sample_and_delta = trees[1].query(5, 0, 1);
  assert(std::get<0>(sample_and_delta) == 20);

  sample_and_delta = trees[1].query(5, 0, 2);
  assert(std::get<0>(sample_and_delta) == 0);

  sample_and_delta = trees[1].query(5, 0, 3);
  assert(std::get<0>(sample_and_delta) == 17);

  for(size_t i = 1; i <= 13; i++) {
    sample_and_delta = trees[2].query(trees[2].left_most_i, 0, i);
    assert(std::get<0>(sample_and_delta) == i);
  }

  for(size_t i = 1; i <= 8; i++) {
    sample_and_delta = trees[3].query(trees[3].left_most_i, 0, i);
    assert(std::get<0>(sample_and_delta) == i);
  }

  for(size_t i = 1; i <= 123; i++) {
    sample_and_delta = trees[4].query(trees[4].left_most_i, 0, i);
    assert(std::get<0>(sample_and_delta) == i);
  }

  // retrieve every sample from every possible leaf.
  for(size_t i = 0; i < trees[4].leaf_node_bv.size(); i++) {
    if(trees[4].leaf_node_bv[i] == 1) {
      ulint node_distance = trees[4].calculate_d(trees[4].left_most_i, i);
      ulint samples_left = 123 - node_distance;
      for(size_t j = 1; j <= samples_left; j++) {
        sample_and_delta = trees[4].query(i, 0, j);
        assert(std::get<0>(sample_and_delta) == (node_distance + j));
      }
    }
  }

  sample_and_delta = trees[5].query(9, 0, 3);
  assert(std::get<0>(sample_and_delta) == 11);

  cout << "all good!" << endl;

  return 0;
}
