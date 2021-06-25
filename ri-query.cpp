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

bool hyb=false;

// do we need the esa file? 
void help() {
  cout << "ri-query: RA to SA using phi" << endl;
  cout << "Usage:       ri-query <index> <end samples> <# of queries>" << endl;
  //cout << "   -h        use hybrid bitvectors instead of elias-fano in both RLBWT and predecessor structures. -h is required "<<endl;
  //cout << "             if the index was built with -h options enabled."<<endl;
  cout << "   <index>   index file (with extension .ri)\n   <end samples>   .esa file\n   <# of queries>   inputting '0', queries the entire index" << endl;
  exit(0);
}

void parse_args(char** argv, int argc, int &ptr) {
  assert(ptr<argc);
  string s(argv[ptr]);
  ptr++;
  if(argc < 4)
    help();
}

vector<ulint>& get_esa(std::string fname, ulint n, vector<ulint> &esa) {
    esa.clear();
    std::ifstream ifs(fname);
    uint64_t x = 0;
    uint64_t y = 0;
    while (ifs.read((char*) &x, 5) && ifs.read((char*) &y, 5)) {
        esa.push_back(y ? y-1 : n-1);
    }
    return esa;
}

void SA(std::set<ulint> &samples, string esaFilename, r_index<> idx) {
  string bwt = idx.get_bwt(); // see if theres a way to avoid this. because if there is you didnt see it.
  rle_string rleBwt = rle_string(bwt); // same with this.
  ulint n = rleBwt.size();
  ulint r = idx.number_of_runs();
  ulint j;
  ulint run;

  std::vector<ulint> esa;
  get_esa(esaFilename, n, esa);
  r = esa.size();

  // here we find out how long it takes to query all the samples being asked for.
  cout << "\nperforming " << samples.size() << " queries:" << endl;
  auto t1 = std::chrono::high_resolution_clock::now();
  for (std::set<ulint>::iterator iter = samples.begin(); iter != samples.end(); iter++) {
    ulint i = *iter;
    run = rleBwt.run_of_position(i);
    j = rleBwt.run_range(run).second;
    ulint phiVal = esa[run];

    for (size_t iter = 0; iter < (j-i); iter++)
       phiVal = idx.Phi(phiVal);

    //cout << "i: " << i << ", run: " << run << ", j: " << j << ", phiVal: " << phiVal << endl;
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  ulint total1 = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1>>>(t2 - t1).count();
  ulint total2 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
  cout << "time taken | " << total1 << " seconds, " << total2 << " ms.\naverage query time | " << (static_cast<float>(total2) / static_cast<float>(samples.size())) << "ms.\n" <<  endl;
  cout << "# of queries: " << samples.size() << endl;
  cout << "# of runs (r): " << r << endl;
  cout << "n: " << n << endl;
  cout << "n/r: " << (static_cast<float>(n)/static_cast<float>(r)) << endl;
}

void SA_all(string esaFilename, r_index<> idx) {
  string bwt = idx.get_bwt();
  rle_string rleBwt = rle_string(bwt);
  ulint n = rleBwt.size();
  ulint r = idx.number_of_runs();
  ulint j;
  ulint run;

  std::vector<ulint> esa;
  get_esa(esaFilename, n, esa);
  r = esa.size();

  // here we find out how long it takes to query all the samples being asked for.
  cout << "\nperforming all queries:" << endl;
  auto t1 = std::chrono::high_resolution_clock::now();

  for (size_t z = 1; z <= n; z++)
    idx.Phi(z);

  //cout << "i: " << i << ", run: " << run << ", j: " << j << ", phiVal: " << phiVal << endl;
  auto t2 = std::chrono::high_resolution_clock::now();
  ulint total1 = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1>>>(t2 - t1).count();
  ulint total2 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
  cout << "time taken | " << total1 << " seconds, " << total2 << " ms.\naverage query time | " << (static_cast<float>(total2) / static_cast<float>(r)) << "ms.\n" <<  endl;
  cout << "# of queries: " << n << endl;
  cout << "# of runs (r): " << r << endl;
  cout << "n: " << n << endl;
  cout << "n/r: " << (static_cast<float>(n)/static_cast<float>(r)) << endl;
}

template<typename ri_t>
void generateSamples(std::set<ulint> &samples, ulint amount, ri_t idx) {
  cout << "generating samples:" << endl;
  ulint min = 1;
  ulint max = idx.number_of_runs() - 1;
  cout << "min: " << min << ", max: " << max << ", amount: " << amount << endl;

  ulint seed = 2706265417;
  std::mt19937_64 generator(seed);
  std::uniform_int_distribution<unsigned long long> dist{min, max};

  while(samples.size() != amount)
    samples.insert(dist(generator));
  cout << "done." << endl;
}

template<typename ri_t>
void run(string filename, string esaFilename, ulint numSamples) {
  ri_t idx;
  idx.load_from_file(filename.c_str());
  std::set<ulint> querySamples;
  generateSamples(querySamples, numSamples, idx);

  if(numSamples == 0)
    SA_all(esaFilename, idx);
  else {
    SA(querySamples, esaFilename, idx);
  }
}

int main(int argc, char** argv) {
  int ptr = 1;
  if(argc<2) help();

  while(ptr<argc-1)
    parse_args(argv, argc, ptr);

  // if(hyb)
  //   run<r_index<sparse_hyb_vector,rle_string_hyb>>(argv[1], argv[2],atoi(argv[3]));
  // else

  run<r_index<>>(argv[1], argv[2],atoi(argv[3]));
}
