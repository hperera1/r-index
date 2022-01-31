// driver file to query an r-index using only phi.

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

bool hyb = false;
string bwt;
rle_string rle_bwt;
std::vector<ulint> esa;

void help() {
  cout << "ri-query: RA to SA using phi." << endl;
  cout << "Usage:       ri-query <index> <end samples> <# of queries>" << endl;
  cout << "   <index>   index file (with extension .ri)";
  cout << "\n   <end samples>   .esa file";
  cout << "\n   <# of queries>   inputting '0', queries the entire index";
  cout << "\n   (optional) <index of query>   provide specific index to query";
  cout << "\n   eg. ri-query abc.fasta.ri abc.fasta.esa 0 100"
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

// queries a particular SA[n].
ulint query_n(string esa_filename, r_index<> &idx, uint sa_n) {
  ulint r = idx.number_of_runs();
  ulint j;
  ulint run;
  r = esa.size();

  run = rle_bwt.run_of_position(sa_n); // run number
  j = rle_bwt.run_range(run).second; // run ends at this value
  ulint phi_val = esa[run]; // end sample at specified run

  for (size_t iter = 0; iter < (j - sa_n); iter++) {
    // cout << phi_val << endl;
    phi_val = idx.Phi(phi_val); // sa-1 until we find n (j-i times)
  }

  return phi_val;
}

// queries every generated sample.
void query_random(std::set<ulint> &samples, string esa_filename, r_index<> idx) {
  string bwt = idx.get_bwt();
  rle_string rle_bwt = rle_string(bwt);
  ulint n = rle_bwt.size();
  ulint r = idx.number_of_runs();
  ulint j;
  ulint run;
  std::vector<ulint> esa;
  get_esa(esa_filename, n, esa);
  r = esa.size();

  // beginning of querying all the samples.
  cout << "\nperforming " << samples.size() << " queries:" << endl;
  auto t1 = std::chrono::high_resolution_clock::now();
  for (std::set<ulint>::iterator sIter = samples.begin(); sIter != samples.end(); sIter++) {
    ulint i = *sIter; // sample.
    run = rle_bwt.run_of_position(i); // run sample belongs to.
    j = rle_bwt.run_range(run).second; // size of run.
    ulint phi_val = esa[run]; // sa value.

    for (size_t iter = 0; iter < (j-i); iter++)
       phi_val = idx.Phi(phi_val); // iterate backwards j-i times to find value of i.
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

// queries the entire r-index.
void query_all(string esa_filename, r_index<> idx) {
  string bwt = idx.get_bwt();
  rle_string rle_bwt = rle_string(bwt);
  ulint n = rle_bwt.size();
  ulint r = idx.number_of_runs();
  ulint j;
  ulint run;

  std::vector<ulint> esa;
  get_esa(esa_filename, n, esa);
  r = esa.size();

  // here we find out how long it takes to query all the samples being asked for.
  cout << "\nperforming all queries:" << endl;
  auto t1 = std::chrono::high_resolution_clock::now();
  for (size_t i = 1; i <= r; i++) {
    run = rle_bwt.run_of_position(i);
    j = rle_bwt.run_range(run).second;
    ulint phi_val = esa[run];

    for (size_t iter = 0; iter < (j-i); iter++)
      phi_val = idx.Phi(phi_val);
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  ulint total1 = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1>>>(t2 - t1).count();
  ulint total2 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
  cout << "time taken | " << total1 << " seconds, " << total2 << " ms.\naverage query time | " << (static_cast<float>(total2) / static_cast<float>(r)) << "ms.\n" <<  endl;
  cout << "# of queries: " << n << endl;
  cout << "# of runs (r): " << r << endl;
  cout << "n: " << n << endl;
  cout << "n/r: " << (static_cast<float>(n)/static_cast<float>(r)) << endl;
}

// generates however many samples we are looking to query.
template<typename ri_t>
void generate_samples(std::set<ulint> &samples, ulint amount, ri_t idx) {
  // cout << "generating samples:" << endl;
  ulint min = 1;
  ulint max = idx.number_of_runs() - 1;
  // cout << "min: " << min << ", max: " << max << ", amount: " << amount << endl;

  ulint seed = 2706265417;
  std::mt19937_64 generator(seed);
  std::uniform_int_distribution<unsigned long long> dist{min, max};

  while(samples.size() != amount)
    samples.insert(dist(generator));

  // cout << "done." << endl;
}

// filename: r-index | esa_filename: end samples | num_samples: how many we want to query
template<typename ri_t>
void run(string filename, string esa_filename, ulint num_samples, uint sa_n) {
  ri_t idx;
  idx.load_from_file(filename.c_str());
  ulint n = rle_bwt.size();
  std::set<ulint> query_samples;
  generate_samples(query_samples, num_samples, idx);

  // 0 means we query the whole index but if followed by a nonzero integer we query that value.
  if(num_samples == 0) {
    if(sa_n == 0) {
      query_all(esa_filename, idx);
    }
    else {
      bwt = idx.get_bwt();
      rle_bwt = rle_string(bwt);
      ulint n = rle_bwt.size();
      get_esa(esa_filename, n, esa);

      ofstream samples;
      samples.open("../rasa_tests/chr19_2_samples.txt");
      for(size_t i = 0; i < sa_n; i++) {
        ulint query_result = query_n(esa_filename, idx, i);
        samples << query_result;
        samples << "\n";
      }
      samples.close();

      // ulint result = query_n(esa_filename, idx, sa_n);
      // cout << "sa[i]: " << result << endl;
    }
  }
  else {
    query_random(query_samples, esa_filename, idx);
  }
}

// the issue with using a templated run in order to work with a hybrid r-index
// is that they didn't have the same functions.
// if(hyb)
//   run<r_index<sparse_hyb_vector,rle_string_hyb>>(argv[1], argv[2],atoi(argv[3]));
int main(int argc, char** argv) {
  int ptr = 1;
  if(argc < 3) {
    help();
  }

  while(ptr < (argc - 1)) {
    parse_args(argv, argc, ptr);
  }

  run<r_index<>>(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
}
