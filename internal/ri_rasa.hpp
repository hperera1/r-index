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

// temp notes for me: pred_to_run is basically esa_map

namespace ri {
  template<class sparse_bv_type = sparse_sd_vector>
class rads {
public:
  rads(){};
  rads(std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa) {
    init_by_value(ssa, esa);
    cout << "Done. Now listing paths." << endl;
    list_paths(ssa);
  }

  void init_by_value(std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa) {
    assert(ssa.size() == esa.size());
    std::vector<ulint> esa_sorted = esa;
    esa_map.reserve(esa.size());
    phi_inv_sa.resize(esa.size());
    bounds.resize(esa.size());

    // initialization of map variables
    for(ulint i = 0; i < esa.size(); i++) {
      esa_map[esa[i]] = i;
      if(i < ssa.size() - 1) {
        pis_inv[ssa[i+1].first] = i;
      }
    }

    // esa needs to be sorted for everything past this.
    std::sort(ssa.begin(), ssa.end()); // this wont be needed once phi_inv_sa is sorted.
    std::sort(esa_sorted.begin(), esa_sorted.end());

    ulint i = 0; // ssa iterator
    ulint j = 0; // esa iterator
    ulint node = esa_map[esa_sorted.back()]; // init node as biggest value in pred ds because its circular.

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
        node = esa_map[esa_sorted[j]];
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
        std::vector<uint> current_path;
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
          rads_tree branch = rads_tree(esa_map, current_path, bounds, trees.size()+1, tree_pointers);
          trees.push_back(branch);
          for(size_t i = 0; i < current_path.size(); i++)
            temp_trees_bv[current_path[i]] = true; // set the nodes that are in the cycle to true in our bitvector
        }
      }
    }
    std::sort(tree_pointers.begin(), tree_pointers.end());
    trees_bv = sparse_bv_type(temp_trees_bv);
  }

  void query(ulint sa_i, sparse_bv_type &pred, int_vector<> &pred_to_run) {
    // do we need some assertions like phi does?
    // find run & pred
    // is pred in a cycle -> bv check
    // if so use tree, else use phi
    ulint sa_i_r = pred.predecessor_rank_circular(sa_i); // rank of the pred of sa_i
    ulint sa_i_j = pred.select(sa_i_r); // select the actual pred
    if(trees_bv[pred_to_run[sa_i_j]]) { // if true, this sample is in a tree
      // rank query the index of the run of the sample to get which tree pointer to use
      // this gives us tree # and leaf node in the tree array
      trees[std::get<1>(tree_pointers[trees_bv.rank(pred_to_run[sa_i_j])])].query_helper();
    }
    else {

    }
  }

  inline int get_size() {
    return phi_inv_sa.size();
  }

  inline int get_num_paths() {
    return trees.size();
  }

protected:
  std::unordered_map<ulint, ulint> esa_map;
  std::unordered_map<ulint, ulint> pis_inv;

  std::vector<std::tuple<ulint, uint, uint>> tree_pointers; // pointers to the corresponding run & tree & leaf node.
  std::vector<std::pair<ulint,ulint>> bounds; // lower and upper bounds of each node in the sa graph.
  std::vector<ulint> phi_inv_sa; // adj. list representing the sa graph.
  std::vector<rads_tree> trees; // list of cycle trees.
  sparse_bv_type trees_bv; // bitvector that tells us which samples are in trees.
};
}

#endif
