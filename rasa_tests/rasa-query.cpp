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
  cout << "\n   <start samples>   .ssa file";
  cout << "\n   (optional) <-u | -n>   query up to the 'n'-th sample starting at zero, or query the 'n'-th sample";
  cout << "\n   (optional) <n>   n for the previous argument" << endl;
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

ulint query_n(r_index<> &idx, std::vector<ulint> &ssa, ulint sa_n) {
  ulint query_val = 0;
  ulint n = rle_bwt.size();
  ulint r = idx.number_of_runs();
  ulint j;
  ulint run;
  r = esa.size();

  query_val = idx.query_csa(sa_n, ssa);

  return query_val;
}

void run(r_index<> idx, std::vector<ulint> &ssa, ulint num_samples) {\
  string string_sample;
  ulint ulint_sample;
  ulint query_result;
  ifstream in_samples("../rasa_tests/chr19_1_samples.txt");
  for(size_t i = 0; i <= num_samples; i++) {
    getline(in_samples, string_sample);
    ulint_sample = strtoul(string_sample.c_str(), NULL, 0);
    query_result = query_n(idx, ssa, i);

    // assert(query_result == ulint_sample);

    if(query_result != ulint_sample) {
      cout << "index: " << i << endl;
      cout << "correct sample: " << ulint_sample << endl;
      cout << "retrieved sample: " << query_result << endl;
    }
  }
}

int main(int argc, char** argv) {
  int ptr = 1;
  if(argc < 3) {
    help();
  }

  while(ptr < (argc - 1)) {
    parse_args(argv, argc, ptr);
  }

  string idx_name = argv[1];
  string ssa_filename = argv[2];

  r_index<> idx;
  std::vector<ulint> ssa;
  idx.load_from_file(idx_name.c_str());
  bwt = idx.get_bwt();
  rle_bwt = rle_string(bwt);
  get_sa(ssa_filename, rle_bwt.size(), ssa);
  cout << "index, bwt, and rle loaded." << endl;

  if(strcmp(argv[3], "-u") == 0) {
    run(idx, ssa, atoi(argv[4]));
  }
  else if(strcmp(argv[3], "-n") == 0) {
    query_n(idx, ssa, stoul(argv[4]));
  }
}
