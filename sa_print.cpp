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

void help() {
  cout << "sa-print: Prints the range of values of the SA samples of the file provided." << endl;
  cout << "Usage:       sa_print <filename> <lower bound> <upper bound>" << endl;
  cout << "\n   <filename>   .xsa file";
  cout << "\n   <lower bound>   lower bound of values";
  cout << "\n   <upper bound>   upper bound of values" << endl;
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

void run(string filename, ulint lower, ulint upper) {
  std::vector<ulint> sa;
  get_sa(filename, 59128984, sa);

  for(size_t i = lower; i <= upper; i++) {
    cout << i << ": " << sa[i] << endl;
    // if(sa[i] == 17671980) {
    //   cout << "found it, " << i << endl;
    // }
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

  run(argv[1], atoi(argv[2]), atoi(argv[3]));
}
