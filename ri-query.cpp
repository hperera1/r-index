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
  cout << "\n   (optional) <-r | -u | -n>   randomly query 'n' samples, query up to the 'n'-th sample starting at zero, or query the 'n'-th sample";
  cout << "\n   (optional) <n>   n for the previous argument";
  cout << "\n   eg. ri-query abc.fasta.ri abc.fasta.esa -r 100 (queries 100 random samples)" << endl;
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

vector<ulint>& get_sa(string fname, ulint n, vector<ulint> &sa) {
    sa.clear();
    std::ifstream ifs(fname);
    uint64_t x = 0;
    uint64_t y = 0;
    while (ifs.read((char*) &x, 5) && ifs.read((char*) &y, 5)) {
        sa.push_back(y ? y-1 : n-1);
    }

    return sa;
}

// generates however many samples we are looking to query.
void generate_samples(std::set<ulint> &samples, ulint amount, r_index<> idx) {
  ulint min = 1;
  ulint max = idx.number_of_runs() - 1;

  ulint seed = 2706265417;
  std::mt19937_64 generator(seed);
  std::uniform_int_distribution<unsigned long long> dist{min, max};

  while(samples.size() != amount) {
    samples.insert(dist(generator));
  }

  // cout << "done." << endl;
}

// queries a particular SA[n].
ulint query_n(r_index<> idx, string sa_filename, uint sa_n) {
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
void query_random(r_index<> idx, string sa_filename, int n_samples) {
  std::set<ulint> query_samples;
  generate_samples(query_samples, n_samples, idx);

  string bwt = idx.get_bwt();
  rle_string rle_bwt = rle_string(bwt);
  ulint n = rle_bwt.size();
  ulint r = idx.number_of_runs();
  ulint j;
  ulint run;
  std::vector<ulint> sa;
  get_sa(sa_filename, n, sa);
  r = sa.size();

  // beginning of querying all the samples.
  cout << "\nperforming " << query_samples.size() << " queries:" << endl;
  auto t1 = std::chrono::high_resolution_clock::now();
  for (std::set<ulint>::iterator sIter = query_samples.begin(); sIter != query_samples.end(); sIter++) {
    ulint i = *sIter; // sample.
    run = rle_bwt.run_of_position(i); // run sample belongs to.
    j = rle_bwt.run_range(run).second; // size of run.
    ulint phi_val = sa[run]; // sa value.

    for (size_t iter = 0; iter < (j-i); iter++) {
      phi_val = idx.Phi(phi_val); // iterate backwards j-i times to find value of i.
    }
  }

  auto t2 = std::chrono::high_resolution_clock::now();
  ulint total1 = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1>>>(t2 - t1).count();
  ulint total2 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
  cout << "time taken | " << total1 << " seconds, " << total2 << " ms.\naverage query time | " << (static_cast<float>(total2) / static_cast<float>(query_samples.size())) << "ms.\n" <<  endl;
  cout << "# of queries: " << query_samples.size() << endl;
  cout << "# of runs (r): " << r << endl;
  cout << "n: " << n << endl;
  cout << "n/r: " << (static_cast<float>(n)/static_cast<float>(r)) << endl;
}

// queries the entire r-index.
void query_all(r_index<> idx, string sa_name) {
  string bwt = idx.get_bwt();
  rle_string rle_bwt = rle_string(bwt);
  ulint n = rle_bwt.size();
  ulint r = idx.number_of_runs();
  ulint j;
  ulint run;

  std::vector<ulint> sa;
  get_sa(sa_name, n, sa);
  r = sa.size();

  // here we find out how long it takes to query all the samples being asked for.
  cout << "\nperforming all queries:" << endl;
  auto t1 = std::chrono::high_resolution_clock::now();
  for (size_t i = 1; i <= r; i++) {
    run = rle_bwt.run_of_position(i);
    j = rle_bwt.run_range(run).second;
    ulint phi_val = sa[run];

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

int main(int argc, char** argv) {
  int ptr = 1;
  if(argc < 3) {
    help();
  }

  while(ptr < (argc - 1)) {
    parse_args(argv, argc, ptr);
  }

  string idx_name = argv[1];
  r_index<> idx;
  idx.load_from_file(idx_name.c_str());

  if(argc == 3) { // query whole index (idx, sa)
    query_all(idx, argv[2]);
  }
  else if(argc == 5) {
    cout << argv[3] << endl;
    if(strcmp(argv[3], "-r") == 0) {
      query_random(idx, argv[2], atoi(argv[4]));
    }
    else if(strcmp(argv[3], "-u") == 0) {
      int sample_limit = atoi(argv[4]);
      // ofstream samples;
      // samples.open("../rasa_tests/chr19_1_samples.txt");
      for(size_t i = 0; i < sample_limit; i++) {
        ulint query_result = query_n(idx, argv[2], i);
        // samples << query_result;
        // samples << "\n";
      }
      // samples.close();
    }
    else if(strcmp(argv[3], "-n") == 0) {
      query_n(idx, argv[2], stoul(argv[4]));
    }
  }

  // run<r_index<>>(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
}
