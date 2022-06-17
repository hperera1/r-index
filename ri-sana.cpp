#include <iostream>
#include <random>
#include <set>
#include <cstring>
#include <string>
#include "internal/rle_string.hpp"
#include "internal/r_index.hpp"
#include "internal/utils.hpp"

using namespace ri;
using namespace std;

string bwt;
rle_string rle_bwt;

void help() {
  exit(0);
}

void parse_args(char** argv, int argc, int &ptr) {
  assert(ptr < argc);
  string s(argv[ptr]);
  ptr++;
  if(argc < 2) {
    help();
  }
}

vector<ulint>& get_sa(std::string fname, ulint n, vector<ulint> &sa) {
    sa.clear();
    std::ifstream ifs(fname);
    uint64_t x = 0;
    uint64_t y = 0;
    while (ifs.read((char*) &x, 5) && ifs.read((char*) &y, 5)) {
        sa.push_back(y ? y-1 : n-1);
    }

    return sa;
}

void run(r_index<> idx, std::vector<ulint> &ssa, std::vector<ulint> &bwt_pos, std::vector<ulint> &heads, std::vector<ulint> &tails) {
    idx.get_table_vals(ssa, bwt_pos, tails);
    cout << "done" << endl;
}

int main(int argc, char** argv) {
  int ptr = 1;
  if(argc < 2) {
    help();
  }

  while(ptr < (argc - 1)) {
    parse_args(argv, argc, ptr);
  }

  string idx_name = argv[1];
  string ssa_filename = argv[2];
  std::vector<ulint> bwt_pos;
  std::vector<ulint> heads;
  std::vector<ulint> tails;
  r_index<> idx;
  std::vector<ulint> ssa;
  idx.load_from_file(idx_name.c_str());
  bwt = idx.get_bwt();
  rle_bwt = rle_string(bwt);
  get_sa(ssa_filename, rle_bwt.size(), heads);
  bwt_pos.reserve(heads.size());
  tails.reserve(heads.size());
  
  cout << "index, bwt, and rle loaded." << endl;
  run(idx, ssa, bwt_pos, heads, tails);
}