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

// Testing rads_tree query function.
int main() {
  std::vector<rads_tree<>> trees;
  std::vector<std::tuple<ulint, ulint, uint>> tree_pointers;

  std::vector<ulint> tree_0 = {0, 1, 2, 3, 4, 5, 6, 7};
  std::vector<std::pair<ulint, ulint>> tree_0_bounds = {{0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}};

  std::vector<ulint> tree_1 = {26, 3, 11, 20, 0, 17};
  std::vector<std::pair<ulint, ulint>> tree_1_bounds = {{5,1}, {0,6}, {0,3}, {2,1}, {0,3}, {4,1}};

  trees.push_back(rads_tree(tree_0, tree_0_bounds, 0, tree_pointers));
  trees.push_back(rads_tree(tree_1, tree_1_bounds, 1, tree_pointers));

  std::tuple<ulint, ulint, ulint> sample_and_delta;
  // traversal test without bounds
  for(size_t i = 1; i <= 7; i++) {
    sample_and_delta = trees[0].query(8, 0, i);
    assert(std::get<0>(sample_and_delta) == i);
  }

  // query tests with bounds
  sample_and_delta = trees[1].query(8, 0, 1);
  assert(std::get<0>(sample_and_delta) == 3);

  sample_and_delta = trees[1].query(8, 0, 2);
  assert(std::get<0>(sample_and_delta) == 11);

  sample_and_delta = trees[1].query(5, 0, 1);
  assert(std::get<0>(sample_and_delta) == 20);

  sample_and_delta = trees[1].query(5, 0, 2);
  assert(std::get<0>(sample_and_delta) == 0);

  sample_and_delta = trees[1].query(5, 0, 3);
  assert(std::get<0>(sample_and_delta) == 17);

  cout << "all good!" << endl;

  return 0;
}
