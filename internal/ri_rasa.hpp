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
  rads(std::vector<std::pair<ulint, ulint>> &unsorted_ssa, std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa, sparse_bv_type &pred, int_vector<> &pred_to_run) {
    build_rads_phi(ssa, esa, pred, pred_to_run);
    cout << "Done. Now listing paths." << endl;
    list_paths(ssa);
  }

  // TO-DO: change name of both of these function, think about the sorting that you do here
  // construction for phi_inverse
  void build_rads_phi_inv(std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa, sparse_bv_type &pred) {
    assert(ssa.size() == esa.size());
    std::vector<ulint> esa_sorted = esa;
    sa_map.reserve(esa.size());
    bounds.resize(esa.size());
    sa_graph.resize(esa.size());

    // initialization of map variables
    for(ulint i = 0; i < esa.size(); i++) {
      sa_map[esa[i]] = i; // ! this will be replaced !
      if(i < ssa.size() - 1) {
        phi_map[ssa[i+1].first] = i;
      }
    }

    // esa needs to be sorted for everything past this.
    std::sort(ssa.begin(), ssa.end()); // this wont be needed once sa_graph is sorted.
    std::sort(esa_sorted.begin(), esa_sorted.end());

    ulint i = 0; // ssa iterator
    ulint j = 0; // esa iterator
    ulint node = sa_map[esa_sorted.back()]; // init node as biggest value in pred ds because its circular.

    // computing the predecessor values using a merge-sort like combine phase
    while((i < ssa.size()) && (j < esa.size())) {
      if(ssa[i].first < esa_sorted[j]) {
        assert(node <= ssa[i].first);
        ulint successor_rank = pred.predecessor_rank_circular(phi_map[ssa[i].first]) + 1;
        ulint successor = pred.select(successor_rank);
        sa_graph[phi_map[ssa[i].first]] = node; // this has to be improved, the hash is unnecessary
        bounds[phi_map[ssa[i].first]].first = ssa[i].first - esa_sorted[j-1]; // start_sample - pred(start_sample)
        bounds[phi_map[ssa[i].first]].second = successor - esa[phi_map[ssa[i].first]]; // successor of sample - pred of sample
        i += 1;
      }
      else {
        node = sa_map[esa_sorted[j]];
        j += 1;
      }
    }

    while(i < ssa.size()) {
      ulint successor_rank = pred.predecessor_rank_circular(esa[phi_map[ssa[i].first]]) + 1;
      ulint successor = pred.select(successor_rank);
      ulint predecessor = pred.select(pred_rank - 1);
      sa_graph[phi_map[ssa[i].first]] = node;
      bounds[phi_map[ssa[i].first]].first = ssa[i].first - esa_sorted[j-1];
      bounds[phi_map[ssa[i].first]].second = successor - predecessor;
      i += 1;
    }
  }

  // TO-DO: change name of both of these function, think about the sorting that you do here
  // construction for phi
  void build_rads_phi(std::vector<std::pair<ulint, ulint>> &unsorted_ssa, std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa, sparse_bv_type &pred, int_vector<> &pred_to_run) {
    cout << "Building the rads for Phi() ..." << endl;

    assert(ssa.size() == esa.size());
    std::vector<ulint> esa_sorted = esa; // we currently need to sort the end samples for later
    sa_graph.resize(esa.size());
    bounds.resize(esa.size());

    // we need to build the inverse map for pred_to_run
    cout << "Building Inv SA map ..." << endl;
    for(ulint i = 1; i < esa.size(); i++) {
      phi_map[esa[i-1]] = i;
    }

    ulint i = 0;
    ulint j = 0;
    ulint node = pred_to_run[ssa_sorted.back().first];
    std::sort(esa_sorted.begin(), esa_sorted.end());

    // building edges between samples (u,v)
    cout << "Building graph ..." << endl;
    while((i < ssa.size()) && (j < esa.size())) {
      if(esa_sorted[j] < ssa_sorted[i].first) {
        assert(node <= esa_sorted[j]);
        ulint successor_rank = pred.predecessor_rank_circular(unsorted_ssa[phi_map[esa_sorted[j]]].first) + 1;
        ulint successor = pred.select(successor_rank);
        ulint predecessor = pred.select(successor_rank - 1);
        sa_graph[phi_map[esa_sorted[j]]] = node;
        bounds[phi_map[esa_sorted[j]]].first = esa_sorted[j] - ssa_sorted[i - 1].first; // v - pred(v)
        bounds[phi_map[esa_sorted[j]]].second = successor - predecessor; // succ(u) - pred(u)
        j += 1;
      }
      else {
        node = pred_to_run[ssa_sorted[i].first];
        i += 1;
      }
    }

    cout << "Final graph loop ..." << endl;
    while(j < esa_sorted.size()) {
      ulint successor_rank = pred.predecessor_rank_circular(unsorted_ssa[phi_map[esa_sorted[j]]].first) + 1;
      ulint successor = pred.select(successor_rank);
      ulint predecessor = pred.select(successor_rank - 1);
      sa_graph[phi_map[esa_sorted[j]]] = node;
      bounds[phi_map[esa_sorted[j]]].first = esa_sorted[j] - ssa_sorted[i - 1].first; // v - pred(v)
      bounds[phi_map[esa_sorted[j]]].second = successor - predecessor; // succ(u) - pred(u)
      j += 1;
    }
  }

  // list_paths() finds all cycles in our graph.
  void list_paths(std::vector<std::pair<ulint, ulint>> &ssa) {
    std::vector<uint> indegrees(sa_graph.size(), 0); // indegrees of the nodes
    std::vector<bool> visited(sa_graph.size(), false); // visited nodes so far
    auto temp_trees_bv = vector<bool>(sa_graph.size(), false);

    // counting indegrees of the nodes
    for(size_t i = 0; i < sa_graph.size(); i++) {
      if(sa_graph[i] >= 0)
        indegrees[sa_graph[i]] += 1;
    }

    // for all nodes with indegree 0, we check if they are a cycle.
    for(size_t i = 0; i < sa_graph.size(); i++) {
      if(indegrees[i] == 0) {
        std::vector<ulint> current_path;
        int u = i;
        int v = sa_graph[u];
        visited[u] = true;
        current_path.push_back(u);

        while(visited[v] == false) {
          current_path.push_back(v);
          visited[v] = true;
          v = sa_graph[v];
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
    if(trees_bv[pred_to_run[sa]])
      return true;

    return false;
  }

  inline int get_size() {
    return sa_graph.size();
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
  std::unordered_map<ulint, ulint> phi_map;
  // int_vector<> *pred_to_run_ptr; // might want to add pointers to these things so that we dont need to pass them around everywhere
  // sparse_bv_type *predecessors_ptr;

  std::vector<std::tuple<ulint, ulint, uint>> tree_pointers; // pointers to the corresponding run & tree & leaf node.
  std::vector<std::pair<ulint,ulint>> bounds; // lower and upper bounds of each node in the sa graph. // can be deleted at some point
  std::vector<ulint> sa_graph; // adj. list representing the sa graph.
  std::vector<rads_tree<>> trees; // list of cycle trees.
  sparse_bv_type trees_bv; // bitvector that tells us which samples are in trees.
};
}

#endif
