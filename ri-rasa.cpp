#include <cstring>
#include <iostream>
#include <string>

#include "internal/r_index.hpp"
#include "internal/ri_rasa.hpp"

using namespace ri;
using namespace std;

void help() {
  cout << "ri-rasa: main file for query data structure for now." << endl;
  exit(0);
}

void parse_args(char** argv, int argc, int &ptr) {
  assert(ptr<argc);
  string s(argv[ptr]);
  ptr++;
  if(argc < 2)
    help();
}

int main(int argc, char** argv) {
  int ptr = 1;
  if(argc<2) help();

  while(ptr<argc-1)
    parse_args(argv, argc, ptr);
}
