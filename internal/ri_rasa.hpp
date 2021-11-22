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
// - pred_to_run is basically sa_map so that is to be removed
// - do we need to do < or <= when checking the costs

namespace ri {
  template<class sparse_bv_type = sparse_sd_vector,
           class rle_string_t = rle_string_sd>
class rads {
public:
  rads(){};
  rads(std::vector<std::pair<ulint, ulint>> &unsorted_ssa, std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa, sparse_bv_type &pred) {
    build_rads_phi(unsorted_ssa, ssa, esa, pred);
    cout << "Searching paths ..." << endl;
    find_cycles(ssa);
    cout << "--------------------------" << endl;
    cout << "Debug info:" << endl;
    cout << "Largest tree: " << this->get_largest_tree() << endl;
  }

  // TO-DO: how to make efficient?
  // construction for phi_inverse
  // unsorted_ssa: samples_first before sort | ssa: sorted samples_first (post phi) | esa: samples_last | pred: predecessor data structure used to calculate successors and preds
  void build_rads_phi_inv(std::vector<std::pair<ulint, ulint>> &unsorted_ssa, std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa, sparse_bv_type &pred) {
    assert(ssa.size() == esa.size());
    std::vector<ulint> esa_sorted = esa;
    sa_map.reserve(esa.size());
    bounds.resize(esa.size());
    sa_graph.resize(esa.size());

    // initialization of map variables
    for(ulint i = 0; i < esa.size(); i++) {
      sa_map[esa[i]] = i; // sa_map can be removed I think but for now lets leave it
      if(i < ssa.size() - 1) {
        phi_x_inv[ssa[i+1].first] = i;
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
        ulint successor_rank = pred.predecessor_rank_circular(phi_x_inv[ssa[i].first]) + 1;
        ulint successor = pred.select(successor_rank);
        ulint predecessor = pred.select(successor_rank - 1);
        sa_graph[phi_x_inv[ssa[i].first]] = node; // this has to be improved, the hash is unnecessary
        bounds[phi_x_inv[ssa[i].first]].first = ssa[i].first - esa_sorted[j - 1]; // start_sample - pred(start_sample)
        bounds[phi_x_inv[ssa[i].first]].second = successor - esa[phi_x_inv[ssa[i].first]]; // successor of sample - pred of sample
        i += 1;
      }
      else {
        node = sa_map[esa_sorted[j]];
        j += 1;
      }
    }

    while(i < ssa.size()) {
      ulint successor_rank = pred.predecessor_rank_circular(esa[phi_x_inv[ssa[i].first]]) + 1;
      ulint successor = pred.select(successor_rank);
      ulint predecessor = pred.select(successor_rank - 1);
      sa_graph[phi_x_inv[ssa[i].first]] = node;
      bounds[phi_x_inv[ssa[i].first]].first = ssa[i].first - esa_sorted[j - 1];
      bounds[phi_x_inv[ssa[i].first]].second = successor - predecessor;
      i += 1;
    }
  }

  // TO-DO: how to make efficient?
  // construction for phi
  // unsorted_ssa: samples_first before sort | ssa: sorted samples_first (post phi) | esa: samples_last | pred: predecessor data structure used to calculate successors and preds
  void build_rads_phi(std::vector<std::pair<ulint, ulint>> &unsorted_ssa, std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa, sparse_bv_type &pred) {
    cout << "Building rads for phi ..." << endl;
    assert(ssa.size() == esa.size());
    std::vector<ulint> esa_sorted = esa; // we sort a temp array to not mess with esa
    sa_map.reserve(esa.size());
    sa_graph.resize(esa.size());
    bounds.resize(esa.size());

    cout << "Building SA map ..." << endl;
    for(ulint i = 0; i < esa.size(); i++) {
      sa_map[unsorted_ssa[i].first] = i; // sa_map just assignes run values to the ssa samples
      if(i > 0) {
        phi_x_inv[esa[i - 1]] = i; // array to tell you which samples phi takes you to
      }
    }

    ulint i = 0;
    ulint j = 0;
    ulint curr_pred = sa_map[ssa.back().first]; // assign pred to the last start sample
    std::sort(esa_sorted.begin(), esa_sorted.end()); // sort our temporary esa array

    // calculating predecessors and setting them as edges
    // sa_graph[...] = curr_pred; this means that the sample at the i-th run points to curr_pred
    cout << "Building graph ..." << endl;
    while((i < ssa.size()) && (j < esa.size())) {
      if(esa_sorted[j] < ssa[i].first) { // this will always be false on the first iteration. this means so long as the end sample is smaller than the start sample
        // assert(i == 0 || curr_pred <= esa_sorted[j]); | this assertion does not work because we would need to get esa_sorted[j]'s run position, not sample position
        ulint successor_rank = pred.predecessor_rank_circular(ssa[phi_x_inv[esa_sorted[i]]].first) + 1;
        ulint predecessor = pred.select(successor_rank - 1);
        ulint successor = 0;
        if(successor_rank >= pred.number_of_1()) { // this case happens when we are getting the successor of the largest sample
          successor = predecessor + 1; // if true, then successor is just one bigger because the cost will be 1
        }
        else {
          successor = pred.select(successor_rank);
        }

        sa_graph[phi_x_inv[esa_sorted[j]]] = curr_pred; // sample -> curr_pred
        bounds[phi_x_inv[esa_sorted[j]]].first = esa_sorted[j] - ssa[i - 1].first; // lower_bound
        bounds[phi_x_inv[esa_sorted[j]]].second = successor - predecessor; // upper_bound
        j += 1;
      }
      else {
        curr_pred = sa_map[ssa[i].first];
        i += 1;
      }
    }

    // if there are any leftover samples that need to point over to the last curr_pred that gets set, this is necessary
    while(j < esa.size()) {
      ulint successor_rank = pred.predecessor_rank_circular(ssa[phi_x_inv[esa_sorted[i]]].first) + 1;
      ulint predecessor = pred.select(successor_rank - 1);
      ulint successor = 0;
      if(successor_rank >= pred.number_of_1()) {
        successor = predecessor + 1;
      }
      else {
        successor = pred.select(successor_rank);
      }

      sa_graph[phi_x_inv[esa_sorted[j]]] = curr_pred;
      bounds[phi_x_inv[esa_sorted[j]]].first = esa_sorted[j] - ssa[i - 1].first;
      bounds[phi_x_inv[esa_sorted[j]]].second = successor - predecessor;
      j += 1;
    }
  }

  // find_cycles() finds all cycles in our graph.
  void find_cycles(std::vector<std::pair<ulint, ulint>> &ssa) {
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

        // we can implement a min. path length threshold to include paths that
        // may not be cycles but can still be used to traverse samples.
        if(is_cycle) { // if the path is a cycle we construct a tree
          rads_tree branch = rads_tree(current_path, bounds, trees.size(), tree_pointers);
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

  // sa_i is what we want to find
  // we can change std::vector<int> esa to int_vector<> esa (r_index copies samples_last_vec to samples_last)
  void query(ulint sa_i, rle_string_t &bwt, int_vector<> &pred_to_run, sparse_bv_type &pred, std::vector<ulint> &sa) {
    // pass pred into the query helper
    // is pred in a cycle -> bv check
    // if so use tree, else use phi

    // return SA value at position i (sa_i)
    ulint run = bwt.run_of_position(sa_i);
    ulint run_l = bwt.run_range(run).second;
    ulint sa_j = sa[run]; // sa value at position i

    cout << "run: " << run << endl;
    cout << "run l: " << run_l << endl;
    cout << "sample: " << sa_j << endl;

    helper_query(sa_j, run_l-sa_i, run, bwt, pred_to_run, pred, sa); // after helper_query, sa_j will be the
  }

  // args: sa_j & d (j-i) | returns: what do we need to return? just d?
  void helper_query(ulint &sa_j, ulint d, ulint run, rle_string_t &bwt, int_vector<> &pred_to_run, sparse_bv_type &pred, std::vector<ulint> &sa) {
    ulint sa_jr;
    ulint sa_prime;
    ulint cost;

    while(d > 0) {
      sa_jr = pred.predecessor_rank_circular(sa_j); // sample predecessor run
      sa_prime = pred.select(sa_jr); // actual predecessor of sa_j
      cost = sa_j - sa_prime; // distance between sample and predecessor

      cout << "sa_jr: " << sa_jr << endl;
      cout << "sa_prime: " << sa_prime << endl;
      cout << "cost: " << cost << endl;

      // check if pred is in a cycle
      if(in_cycle(sa_prime)) { // dont need to use sa_prime, you can use sa_jr because the rank tells you what index to query
        cout << "\nwere in a cycle ..." << endl;
        std::tuple<ulint, ulint, uint> tree_info = tree_pointers[trees_bv.rank(run)];
        cout << "sample run: " << std::get<0>(tree_info) << ", sample tree: " << std::get<1>(tree_info) << ", sample leaf: " << std::get<2>(tree_info) << endl;
        std::pair<ulint, ulint> sa_prime_and_d = trees[std::get<1>(tree_info)].query(std::get<2>(tree_info), cost, d);
        // sa_prime is the new sample (leaf) that we got from the tree.

        sa_j = sa_prime_and_d.first; // review these two // sa_j is sa_primes run.
        d = d - sa_prime_and_d.second; // this is the distance left over.
      }
      else { // continue the iteration using phi
        ulint delta = sa_prime < sa_j ? sa_j - sa_prime : sa_j + 1;
        ulint prev_sample = sa[pred_to_run[sa_jr] - 1]; // we dont have samples_last, need to pass it in.
        sa_j = (prev_sample + delta) % bwt.size();
        d -= 1;
      }
    }
  }

  bool in_cycle(ulint sa) {
    if(trees_bv[sa_map[sa]])
      return true;

    return false;
  }

  inline int get_num_trees() {
    return tree_pointers.size();
  }

  int get_largest_tree() {
    int max_tree_size = 0;
    for(int i = 0; i < trees.size(); i++) {
      if(trees[i].leaf_samples.size() > max_tree_size)
        max_tree_size = trees[i].leaf_samples.size();
    }

    return max_tree_size;
  }

  void print_tree_runs(std::vector<std::pair<ulint, ulint>> &ssa) {
    for(size_t i = 0; i < 10; i++) {
      cout << "i: " << i << ", s_sample: " << ssa[std::get<0>(tree_pointers[i])].first << endl;
      cout << "run #: " << std::get<0>(tree_pointers[i]) << ", tree #: " << std::get<1>(tree_pointers[i]) << ", leaf node: " << std::get<2>(tree_pointers[i]) << endl;
    }

    // for(size_t i = 0; i < tree_pointers.size(); i++) {
    //   if(std::get<1>(tree_pointers[i]) == 0) {
    //     cout << "i: " << i << ", s_sample: " << ssa[std::get<0>(tree_pointers[i])].first << endl;
    //     cout << "run #: " << std::get<0>(tree_pointers[i]) << ", tree #: " << std::get<1>(tree_pointers[i]) << ", leaf node: " << std::get<2>(tree_pointers[i]) << endl;
    //   }
    // }

    // cout << endl;
    // for(size_t i = 0; i < trees[0].leaf_samples.size(); i++) {
    //   cout << trees[0].leaf_samples[i] << endl;
    // }

    cout << endl;
    for(size_t i = 0; i < 8; i++) {
      cout << trees[0].leaf_samples[i] << endl;
    }
    cout << endl;

    // cout << "\ntree[0] size: " << trees[0].leaf_samples.size() << endl;
    // cout << "tree[1] size: " << trees[1].leaf_samples.size() << endl;
    // cout << "tree[2] size: " << trees[2].leaf_samples.size() << endl;
  }

protected:
  std::unordered_map<ulint, ulint> sa_map; // this is just pred_to_run but pred_to_run wasn't working?
  std::unordered_map<ulint, ulint> phi_x_inv; // map of phi or phi_inverse values
  std::vector<rads_tree<>> trees; // list of trees.
  std::vector<std::tuple<ulint, ulint, uint>> tree_pointers; // pointers to the corresponding run & tree & leaf node.
  std::vector<std::pair<ulint,ulint>> bounds; // lower and upper bounds of each node in the sa graph. // can be deleted at some point
  std::vector<ulint> sa_graph; // adj. list representing the sa graph.
  sparse_bv_type trees_bv; // bitvector that tells us which samples are in trees.
};
}

#endif
