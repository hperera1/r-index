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

ulint query_n(r_index<> &idx, ulint sa_n) {
  ulint query_val = 0;
  ulint n = rle_bwt.size();
  ulint r = idx.number_of_runs();
  ulint j;
  ulint run;
  r = esa.size();

  query_val = idx.query_csa(sa_n);

  return query_val;
}

template<typename ri_t>
void run(string filename,ulint num_samples) {
  ri_t idx;
  ifstream in_samples("../rasa_tests/samples.txt");
  idx.load_from_file(filename.c_str());
  bwt = idx.get_bwt();
  rle_bwt = rle_string(bwt);
  cout << "index, bwt, and rle loaded." << endl;

  string string_sample;
  ulint ulint_sample;
  ulint query_result;
  for(size_t i = 0; i <= num_samples; i++) {
    getline(in_samples, string_sample);
    ulint_sample = strtoul(string_sample.c_str(), NULL, 0);
    query_result = query_n(idx, i);

    // if(query_result != ulint_sample) {
    //   cout << "correct sample: " << ulint_sample << endl;
    //   cout << "retrieved sample: " << query_result << endl;
    // }
    // else {
      cout << i << ": " << query_result << endl;
    // }
  }

  // query_result = query_n(idx, num_samples);
  // cout << "result: " << query_result << endl;
}

int main(int argc, char** argv) {
  int ptr = 1;
  if(argc < 2) {
    help();
  }

  while(ptr < (argc - 1)) {
    parse_args(argv, argc, ptr);
  }

  run<r_index<>>(argv[1], atoi(argv[2]));
}
