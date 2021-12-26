#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <cstring>
#include <string>

#include "../internal/rle_string.hpp"
#include "../internal/r_index.hpp"
#include "../internal/ri_rasa.hpp"
#include "../internal/ri_rasa_tree.hpp"
#include "../internal/utils.hpp"
#include "../internal/sparse_sd_vector.hpp"
#include "../internal/sparse_hyb_vector.hpp"

using namespace sdsl;
using namespace ri;
using std::cout, std::endl;

namespace ri {
  template    <class sparse_bv_type = sparse_sd_vector,
               class rle_string_t = rle_string_sd>
    void start_test() {
      std::vector<std::pair<ulint, ulint>> unsorted_ssa = {{26,0}, {8,0}, {6,0}, {23,0}, {5,0}, {9,0}, {0,0}, {17,0}, {7,0}, {3,0}, {11,0}, {20,0}, {2,0}};
      std::vector<std::pair<ulint, ulint>> ssa = {{26,0}, {8,0}, {6,0}, {23,0}, {5,0}, {9,0}, {0,0}, {17,0}, {7,0}, {3,0}, {11,0}, {20,0}, {2,0}};
      std::vector<ulint> esa = {26, 21, 14, 18, 22, 9, 0, 17, 24, 3, 11, 20, 19};
      sparse_bv_type pred;
      auto pred_bv = std::vector<bool>(27, false);

      std::sort(ssa.begin(), ssa.end());
      for(auto ele : ssa) {
        assert(ele.first < pred_bv.size());
        pred_bv[ele.first] = true;
      }

      pred = pred_bv;
      rads test_rads = rads(unsorted_ssa, ssa, esa, pred);
    }
}

int main() {
  start_test();
  return 0;
}
