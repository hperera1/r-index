// driver file to retrieve samples from an r-index using a combination of the tree
// data structure and phi.

#include <iostream>
#include <random>
#include <set>
#include <cstring>
#include <string>
#include "../internal/rle_string.hpp"
#include "../internal/r_index.hpp"
#include "../internal/utils.hpp"

using namespace ri;
using namespace std;

bool hyb = false;
string bwt;
rle_string rle_bwt;
vector<ulint> esa;

void help() {
  cout << "rasa-query: random access to the csa data structure." << endl;
  cout << "Usage:       ri-query <index> <end samples> <# of queries>" << endl;
  cout << "   <index>   index file (with extension .ri)";
  cout << "\n   <# of queries>   get all samples up to 'n'" << endl;
  exit(0);
}

void parse_args(char** argv, int argc, int &ptr) {
  assert(ptr < argc);
  string s(argv[ptr]);
  ptr++;
  if(argc < 3) {
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

ulint query_n(r_index<> &idx, ulint sa_n, std::vector<ulint> &ssa) {
  ulint query_val = 0;
  ulint n = rle_bwt.size();
  ulint r = idx.number_of_runs();
  ulint j;
  ulint run;
  r = esa.size();

  query_val = idx.query_csa(sa_n, ssa);

  return query_val;
}

template<typename ri_t>
void run(string filename, string ssa_filename, ulint num_samples) {
  ri_t idx;
  std::vector<ulint> ssa;
  ifstream in_samples("../rasa_tests/samples.txt");
  idx.load_from_file(filename.c_str());
  bwt = idx.get_bwt();
  rle_bwt = rle_string(bwt);
  get_sa(ssa_filename, rle_bwt.size(), ssa);
  cout << "index, bwt, and rle loaded." << endl;

  string string_sample;
  ulint ulint_sample;
  ulint query_result;
  for(size_t i = 0; i <= num_samples; i++) {
    if(i == 542) {
      cout << "break" << endl;
    }

    getline(in_samples, string_sample);
    ulint_sample = strtoul(string_sample.c_str(), NULL, 0);
    query_result = query_n(idx, i, ssa);

    // assert(query_result == ulint_sample);

    if(query_result != ulint_sample) {
      cout << "index: " << i << endl;
      cout << "correct sample: " << ulint_sample << endl;
      cout << "retrieved sample: " << query_result << endl;
    }
  }

  cout << "all good!" << endl;
  // query_result = query_n(idx, num_samples);
  // cout << "result: " << query_result << endl;
}

int main(int argc, char** argv) {
  int ptr = 1;
  if(argc < 3) {
    help();
  }

  while(ptr < (argc - 1)) {
    parse_args(argv, argc, ptr);
  }

  run<r_index<>>(argv[1], argv[2], atoi(argv[3]));
}
