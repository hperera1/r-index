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

// we need to take in "samples" which can just be a text file.
// as well as bounds for each sample in the cycle.

int main(int argc, char** argv) {
  std::vector<rads_tree<>> trees;
  std::vector<ulint> samples_cycle = {1,2,3,4,5,6,7,8};
  // ,9,10,11,12,13,14 || , {0,9}, {0,10}, {0,11}, {0,12}, {0,13}, {0,14}
  std::vector<std::pair<ulint, ulint>> sample_bounds = {{0,1}, {0,2}, {0,3}, {0,4}, {0,5}, {0,6}, {0,7}, {0,8}};

  //std::vector<ulint> samples_cycle = {14, 22, 13, 4, 1, 6, 10, 12, 20, 18, 30};
  //std::vector<std::pair<ulint, ulint>> sample_bounds = {{0,3}, {4,1}, {3,2}, {0,1}, {2,1}, {0,1}, {1,1}, {0,2}, {3,3}, {2,2}, {5,5}};
  std::vector<std::tuple<ulint, ulint, uint>> tree_pointers;

  rads_tree example_tree = rads_tree(samples_cycle, sample_bounds, 1, tree_pointers);
  example_tree.print_array();

  for(size_t i = 0; i < 10; i++) {
    trees.push_back(example_tree);
  }

  cout << "\nPrinting tree pointers ..." << endl;
  for(int i = 0; i < tree_pointers.size(); i++) {
    cout << std::get<0>(tree_pointers[i]) << " " << std::get<1>(tree_pointers[i]) << " " << std::get<2>(tree_pointers[i]) << endl;
  }

  cout << "\nTrying rank queries ..." << endl;
  for(ulint i = 1; i < example_tree.leaf_node_bv.size(); i++) {
    cout << trees[0].leaf_node_bv.rank(i) << endl;
  }

  cout << "\nTrying a regular query ..." << endl;
  trees[1].query(11, 0, 4);
}
