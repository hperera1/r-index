#include <iostream>
#include <random>
#include <set>
#include <cstring>
#include <string>

#include "internal/rle_string.hpp"
#include "internal/r_index.hpp"
#include "internal/ri_rasa_tree.hpp"
#include "internal/utils.hpp"

using namespace ri;
using namespace std;

void build_tree(std::vector<rads_tree<>> &trees, std::vector<ulint> samples_cycle, std::vector<std::pair<ulint, ulint>> sample_bounds, ulint tree_num, std::vector<std::tuple<ulint, ulint, uint>> &tree_pointers) {
  rads_tree test_tree = rads_tree(samples_cycle, sample_bounds, tree_num, tree_pointers);
  trees.push_back(test_tree);
}

int main(int argc, char** argv) {
  std::vector<rads_tree<>> trees;
  std::vector<ulint> samples_cycle = {1,2,3,4,5,6,7,8};
  std::vector<std::pair<ulint, ulint>> sample_bounds = {{0,1}, {0,2}, {0,3}, {0,4}, {0,5}, {0,6}, {0,7}, {0,8}};
  std::vector<std::tuple<ulint, ulint, uint>> tree_pointers;

  cout << "Building trees ..." << endl;
  ulint tree_num = 0;
  for(size_t i = 0; i < 10; i++) {
    //trees.push_back(example_tree);
    build_tree(trees, samples_cycle, sample_bounds, tree_num, tree_pointers);
    tree_num++;
  }

  cout << "\nTrying a regular query ..." << endl;
  trees[1].query(8, 0, 7);
}
