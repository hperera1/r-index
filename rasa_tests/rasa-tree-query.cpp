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
  std::vector<ulint> tree_1 = {26, 3, 11, 20, 0, 17};
  std::vector<std::pair<ulint, ulint>> tree_1_bounds = {{5,1}, {0,6}, {0,3}, {2,1}, {0,3}, {4,1}};

  trees.push_back(rads_tree(tree_1, tree_1_bounds, 0, tree_pointers));
  trees[0].print_tree_info();
  cout << endl;

  // query tests
  std::tuple<ulint, ulint, ulint> sample_and_delta = trees[0].query(8, 0, 1);
  cout << "sample: " << std::get<0>(sample_and_delta) << endl;
  cout << "delta: " << std::get<1>(sample_and_delta) << endl;

  sample_and_delta = trees[0].query(8, 0, 2);
  cout << "sample: " << std::get<0>(sample_and_delta) << endl;
  cout << "delta: " << std::get<1>(sample_and_delta) << endl;

  ofstream out("../rasa_tests/test_rads", ios::binary);
  trees[0].serialize(out);

  return 0;
}
