#ifndef R_INDEX_R_A_
#define R_INDEX_R_A_

#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <definitions.hpp>
#include <r_index.hpp>
#include <ri_rasa_tree.hpp>
#include "sparse_sd_vector.hpp"
#include "sparse_hyb_vector.hpp"

using namespace sdsl;

// temp notes for me:
// - pred_to_run is basically sa_map
// - fix the way the samples get built

namespace ri {
  template<class sparse_bv_type = sparse_sd_vector,
           class rle_string_t = rle_string_sd>
class rads {
public:
  rads(){};
  rads(std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa) {
    init_by_value_phi(ssa, esa);
    cout << "Done. Now listing paths." << endl;
    list_paths(ssa);
    cout << "Done. Optional debug info." << endl<<endl;
    // for (size_t i = 0; i < 30; i++) {
    //   cout << i << ": " << std::get<0>(tree_pointers[i]) << " " << std::get<1>(tree_pointers[i]) << " " << std::get<2>(tree_pointers[i]) << endl;
    // }
  }

  // TO-DO: change name of both of these function, think about the sorting that you do here
  // construction for phi_inverse
  void init_by_value(std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa) {
    assert(ssa.size() == esa.size());
    std::vector<ulint> esa_sorted = esa;
    sa_map.reserve(esa.size());
    bounds.resize(esa.size());
    phi_inv_sa.resize(esa.size());

    // initialization of map variables
    for(ulint i = 0; i < esa.size(); i++) {
      sa_map[esa[i]] = i; // ! this will be replaced !
      if(i < ssa.size() - 1) {
        pis_inv[ssa[i+1].first] = i;
      }
    }

    // esa needs to be sorted for everything past this.
    std::sort(ssa.begin(), ssa.end()); // this wont be needed once phi_inv_sa is sorted.
    std::sort(esa_sorted.begin(), esa_sorted.end());

    ulint i = 0; // ssa iterator
    ulint j = 0; // esa iterator
    ulint node = sa_map[esa_sorted.back()]; // init node as biggest value in pred ds because its circular.

    // these bounds are not getting calculated properly...
    // computing the predecessor values using a merge-sort like combine phase
    while((i < ssa.size()) && (j < esa.size())) {
      if(ssa[i].first < esa_sorted[j]) {
        assert(node <= ssa[i].first);
        phi_inv_sa[pis_inv[ssa[i].first]] = node; // this has to be improved, the hash is unnecessary
        bounds[pis_inv[ssa[i].first]].first = ssa[i].first - esa_sorted[j-1];
        bounds[pis_inv[ssa[i].first]].second = esa_sorted[j] - ssa[i].first;
        i += 1;
      }
      else {
        node = sa_map[esa_sorted[j]];
        j += 1;
      }
    }

    while(i < ssa.size()) {
      phi_inv_sa[pis_inv[ssa[i].first]] = node;
      bounds[pis_inv[ssa[i].first]].first = ssa[i].first - esa_sorted[j-1];
      bounds[pis_inv[ssa[i].first]].second = esa_sorted[j] - ssa[i].first;
      i += 1;
    }
  }

  // TO-DO: change name of both of these function, think about the sorting that you do here
  // construction for phi
  void init_by_value_phi(std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa) {
    assert(ssa.size() == esa.size());
    std::vector<std::pair<ulint, ulint>> ssa_sorted = ssa;
    std::vector<ulint> esa_sorted = esa;
    sa_map.reserve(esa.size());
    phi_inv_sa.resize(esa.size());
    bounds.resize(esa.size());

    cout << "building SA map." << endl;
    for(ulint i = 0; i < esa.size(); i++) {
      sa_map[ssa[i].first] = i;
      if(i > 0) {
        pis_inv[esa[i-1]] = i;
      }
    }

    std::sort(ssa_sorted.begin(), ssa_sorted.end());
    std::sort(esa_sorted.begin(), esa_sorted.end());

    ulint i = 0;
    ulint j = 0;
    ulint node = sa_map[ssa_sorted.back().first];

    cout << "building graph." << endl;
    while((i < ssa.size()) && (j < esa.size())) {
      if(esa_sorted[j] < ssa_sorted[i].first) {
        assert(node <= esa_sorted[j]);
        phi_inv_sa[pis_inv[esa_sorted[j]]] = node;
        bounds[pis_inv[esa_sorted[j]]].first = esa_sorted[j] - ssa_sorted[i-1].first;
        bounds[pis_inv[esa_sorted[j]]].second = ssa_sorted[i].first - esa_sorted[j];
        j += 1;
      }
      else {
        node = sa_map[ssa_sorted[i].first];
        i += 1;
      }
    }

    cout << "final graph loop." << endl;
    while(j < esa.size()) {
      phi_inv_sa[pis_inv[esa_sorted[j]]] = node;
      bounds[pis_inv[esa_sorted[j]]].first = esa_sorted[j] - ssa_sorted[i-1].first;
      bounds[pis_inv[esa_sorted[j]]].second = ssa_sorted[i].first - esa_sorted[j];
      j += 1;
    }
  }

  // list_paths() finds all cycles in our graph.
  void list_paths(std::vector<std::pair<ulint, ulint>> &ssa) {
    std::vector<uint> indegrees(phi_inv_sa.size(), 0); // indegrees of the nodes
    std::vector<bool> visited(phi_inv_sa.size(), false); // visited nodes so far
    auto temp_trees_bv = vector<bool>(phi_inv_sa.size(), false);

    // counting indegrees of the nodes
    for(size_t i = 0; i < phi_inv_sa.size(); i++) {
      if(phi_inv_sa[i] >= 0)
        indegrees[phi_inv_sa[i]] += 1;
    }

    // for all nodes with indegree 0, we check if they are a cycle.
    for(size_t i = 0; i < phi_inv_sa.size(); i++) {
      if(indegrees[i] == 0) {
        std::vector<ulint> current_path;
        int u = i;
        int v = phi_inv_sa[u];
        visited[u] = true;
        current_path.push_back(u);

        while(visited[v] == false) {
          current_path.push_back(v);
          visited[v] = true;
          v = phi_inv_sa[v];
        }

        // scan the current path and see if v is in it
        bool is_cycle = false;
        for (size_t i = 0; !is_cycle && i < current_path.size(); i++) {
          if(v == current_path[i])
            is_cycle = true;
        }

        // we can implement a min. path length threshold to include paths, that
        // may not be cycles but can still be used to traverse samples.
        if(is_cycle) { // if the path is a cycle we construct a tree
          rads_tree branch = rads_tree(current_path, bounds, trees.size()+1, tree_pointers);
          trees.push_back(branch);
          for(size_t i = 0; i < current_path.size()-1; i++) // -1 because the the last leaf node in our cycle is useless for queries
            temp_trees_bv[current_path[i]] = true; // set the nodes that are in the cycle to true in our bitvector
        }
      }
    }

    std::sort(tree_pointers.begin(), tree_pointers.end());
    trees_bv = sparse_bv_type(temp_trees_bv);
    bounds.clear();
  }

  // sa_i is what we are querying
  // we can change std::vector<int> esa to int_vector<> esa (r_index copies samples_last_vec to samples_last)
  void query(ulint sa_i, rle_string_t &bwt, int_vector<> &pred_to_run, sparse_bv_type &pred, std::vector<ulint> &esa) {
    // should we do some assertions like phi does?
    // pass pred into the query helper
    // is pred in a cycle -> bv check
    // if so use tree, else use phi

    // return SA value at position i (sa_i)
    ulint run = bwt.run_of_position(sa_i);
    ulint run_l = bwt.run_range(run).second;
    ulint sa_j = esa[run]; // sa value at position i

    cout << "run: " << run << endl;
    cout << "run l: " << run_l << endl;
    cout << "sample: " << sa_j << endl;

    helper_query(sa_j, run_l-sa_i, run, bwt, pred_to_run, pred, esa);
  }

  // args: sa_j & d (j-i) | returns: what do we need to return? just d?
  void helper_query(ulint &sa_j, ulint d, ulint run, rle_string_t &bwt, int_vector<> &pred_to_run, sparse_bv_type &pred, std::vector<ulint> &esa) {
    cout << "distance: " << d << endl;

    while(d > 0) {
      ulint sa_jr = pred.predecessor_rank_circular(sa_j); // start sample predecessor
      ulint sa_prime = pred.select(sa_jr); // actual predecessor value
      ulint cost = sa_j - sa_prime; // difference between end sample and predecessor

      cout << "sa_jr: " << sa_jr << endl;
      cout << "sa_prime: " << sa_prime << endl;
      cout << "cost: " << cost << endl;

      if(in_cycle(sa_prime)) { // dont need to use sa_prime, you can use sa_jr because the rank tells you what index to query
        cout << "\nwere in a cycle!" << endl;

        std::tuple<ulint, ulint, uint> tree_info = tree_pointers[trees_bv.rank(run)];

        cout << "sample run: " << std::get<0>(tree_info) << ", sample tree: " << std::get<1>(tree_info) << ", sample leaf: " << std::get<2>(tree_info) << endl;

        std::pair<ulint, ulint> sa_prime_and_d = trees[std::get<1>(tree_info)].query(std::get<2>(tree_info), cost, d);

        sa_j = pred_to_run[sa_prime_and_d.first];
        d = d - sa_prime_and_d.second;
      }
      else { // continue the iteration using phi
        ulint delta = sa_prime<sa_j ? sa_j-sa_prime : sa_j+1;
        ulint prev_sample = esa[ pred_to_run[sa_jr]-1 ]; // we dont have samples_last, need to pass it in.
        sa_j = (prev_sample + delta) % bwt.size();
      }
    }
  }

  bool in_cycle(ulint sa) {
    if(trees_bv[sa_map[sa]])
      return true;

    return false;
  }

  inline int get_size() {
    return phi_inv_sa.size();
  }

  inline int get_num_treeptr() {
    return tree_pointers.size();
  }

  inline int get_num_paths() {
    return trees.size();
  }

  void get_largest_tree() {

  }

  void print_tree_runs(std::vector<ulint> &esa) {
    for(size_t i = 0; i < 25; i++) {
      cout << std::get<0>(tree_pointers[i]) << " : " << esa[std::get<0>(tree_pointers[i])] << " : " << std::get<1>(tree_pointers[i]) << endl;
    }

    cout << "\n";
    cout << "trees[1] size: " << trees[1].leaf_samples.size() << endl<<endl;
  }

protected:
  std::unordered_map<ulint, ulint> sa_map; // this is just pred_to_run
  std::unordered_map<ulint, ulint> pis_inv;

  std::vector<std::tuple<ulint, ulint, uint>> tree_pointers; // pointers to the corresponding run & tree & leaf node.
  std::vector<std::pair<ulint,ulint>> bounds; // lower and upper bounds of each node in the sa graph. // can be deleted at some point
  std::vector<ulint> phi_inv_sa; // adj. list representing the sa graph.
  std::vector<rads_tree<>> trees; // list of cycle trees.
  sparse_bv_type trees_bv; // bitvector that tells us which samples are in trees.
};
}

#endif
