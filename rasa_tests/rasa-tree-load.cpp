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

int main() {
  ifstream in("../rasa_tests/test_rads", ios::binary);
  if(!in) {
    cout << "invalid instream." << endl;
  }
  rads_tree<> test_rads;
  test_rads.load(in);

  test_rads.print_tree_info();
}
