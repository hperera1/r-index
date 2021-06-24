#include <iostream>
#include <random>
#include <set>
#include <cstring>
#include <string>

#include "internal/r_index.hpp"
#include "internal/utils.hpp"

using namespace ri;
using namespace std;

bool hyb=false;

void help() {
  cout << "ri-query: RA to SA using phi" << endl;
  cout << "Usage:       ri-query <index> <# of samples>" << endl;
  //cout << "   -h        use hybrid bitvectors instead of elias-fano in both RLBWT and predecessor structures. -h is required "<<endl;
  //cout << "             if the index was built with -h options enabled."<<endl;
  cout << "   <index>   index file (with extension .ri)" << endl;
  exit(0);
}


void parse_args(char** argv, int argc, int &ptr) {
  assert(ptr<argc);
  string s(argv[ptr]);
  ptr++;

  /*if(s.compare("-h")==0){
  hyb=true;
}else*/
  {
    cout << "Error: unknown option " << s << endl;
    help();
  }
}

// void generateSamples(std::set<ulint> &samples, ulint amount) {
//   ulint min = 1;
//   ulint max = samples_last.size() - 1;
//   ulint seed = 2706265417;
//   std::mt19937_64 generator(seed);
//   std::uniform_int_distribution<unsigned long long> dist{min, max};
//
//   while(samples.size() != amount)
//   samples.insert(dist(generator));
// }
//
// ulint SA(ulint i, rle_string_t &bwt) {
//   ulint run = bwt.run_of_position(i);
//   ulint j = bwt.run_range(run).second;
//   ulint phiVal = samples_last[run];
//
//   for (size_t iter = 0; iter < (j-i); iter++)
//     phiVal = Phi(phiVal);
//
//   return phiVal;
// }

// , std::set<ulint> &querySamples, ulint numSamples
template<typename ri_t>
void run(string filename) {
  cout << "do we get here" << endl;
  ri_t idx;
  idx.load_from_file(filename.c_str());
  auto space = idx.print_space();
  cout << "\nTOT space: " << space << " Bytes" << endl;

  // auto t1 = std::chrono::high_resolution_clock::now();
  // cout << "\nPerforming queries." << endl;
  // for (std::set<ulint>::iterator iter = querySamples.begin(); iter != querySamples.end(); iter++) {
  //   ulint sample = *iter;
  //   SA(sample);
  // }

  // this is to query the whole thing.
  // for (size_t z = 1; z < samples_last.size(); z++) {
  //   queryPhi(z);
  // }

  // auto t2 = std::chrono::high_resolution_clock::now();
  // ulint total1 = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1>>>(t2 - t1).count();
  // ulint total2 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
  // cout << "time taken | " << total1 << " seconds, " << total2 << " ms.\naverage query time | " << (static_cast<float>(total2) / static_cast<float>(numSamples)) << "ms.\n" <<  endl;
}

int main(int argc, char** argv) {
  int ptr = 1;
  if(argc<2) help();

  while(ptr<argc-1)
    parse_args(argv, argc, ptr);

  // generate random samples to be queryed
  std::set<ulint> querySamples;
  //generateSamples(querySamples, argv[ptr+1]);

  if(hyb)
    run<r_index<sparse_hyb_vector,rle_string_hyb>>(argv[ptr]);
  else
    run<r_index<>>(argv[ptr]);

  // if(hyb)
  // {
  // 	r_index<sparse_hyb_vector,rle_string_hyb> idx;
  // 	idx.load_from_file(argv[ptr]);
  // 	auto space = idx.print_space();
  // 	cout << "\nTOT space: " << space << " Bytes" <<endl;
  // }
  // else
  // {
  // 	r_index<> idx;
  // 	idx.load_from_file(argv[ptr]);
  // 	auto space = idx.print_space();
  // 	cout << "\nTOT space: " << space << " Bytes" <<endl;
  // }
}
