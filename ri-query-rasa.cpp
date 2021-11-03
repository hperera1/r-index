#include <iostream>
#include <random>
#include <set>
#include <cstring>
#include <string>

#include "internal/rle_string.hpp"
#include "internal/r_index.hpp"
#include "internal/ri_rasa_tree.hpp"
#include "internal/utils.hpp"
#include "sparse_sd_vector.hpp"
#include "sparse_hyb_vector.hpp"

using namespace std;
using namespace ri;
namespace ri {
  template<class sparse_bv_type = sparse_sd_vector,
           class rle_string_t = rle_string_sd>

class bv_test {
public:
  bv_test() {
    auto test_bv = vector<bool>(10, false);
    sparse_bv_type test_sparse_bv;
    for(int i = 1; i < test_bv.size(); i+=2) {
      test_bv[i] = true;
    }

    test_sparse_bv = sparse_bv_type(test_bv);
    cout << "rank query." << endl;
    cout << test_sparse_bv.rank(4) << endl;
  }
};
}

int main(int argc, char** argv) {
  cout << "testing bv." << endl;
  bv_test ds = bv_test();

  return 0;
}
