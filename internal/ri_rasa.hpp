#ifndef R_INDEX_R_A_
#define R_INDEX_R_A_

#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <definitions.hpp>
#include <r_index.hpp>
#include <ri_rasa_tree.hpp>

using namespace sdsl;

namespace ri {
class rads {
public:
  rads(){};
  rads(std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa) {
    init_by_value(ssa, esa);
  }

  /* rads(r_index &ridx) {
   *   init_by_reference(ridx);
   * }
   *
   * void init_by_reference(r_index &ridx) {
   *   end_map.reserve(ridx.samples_last.size());
   *   phi_inv_sa.resize(ridx.samples_last.size());
   *   for(ulint i = 0; i < ridx.samples_last.size(); i++) {
   *     esa_map[ridx.samples_last[i]] = i;
   *   }
   * }
   */

  void init_by_value(std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa) {
    assert(ssa.size() == esa.size());
    esa_map.reserve(esa.size());
    phi_inv_sa.resize(esa.size()); // becomes adj. list.
    bounds.resize(esa.size());

    // initialization of map variables
    for(ulint i = 0; i < esa.size(); i++) {
      esa_map[esa[i]] = i;
      if(i < ssa.size() - 1) {
        pis_inv[ssa[i+1].first] = i;
      }
    }

    std::sort(ssa.begin(), ssa.end()); // this wont be needed once phi_inv_sa is sorted.
    std::sort(esa.begin(), esa.end());

    ulint i = 0; // ssa iterator
    ulint j = 0; // esa iterator
    ulint node = esa_map[esa.back()]; // init node as biggest value in pred ds because its circular.

    // computing the predecessor values using a merge-sort like combine phase
    while((i < ssa.size()) && (j < esa.size())) {
      if(ssa[i].first < esa[j]) {
        assert(node <= ssa[i].first);
        phi_inv_sa[pis_inv[ssa[i].first]] = node; // this has to be improved, the hash is unnecessary
        bounds[pis_inv[ssa[i].first]].first = ssa[i].first - esa[j-1];
        bounds[pis_inv[ssa[i].first]].second = esa[j] - ssa[i].first;
        i += 1;
      }
      else {
        node = esa_map[esa[j]];
        j += 1;
      }
    }

    while(i < ssa.size()) {
      phi_inv_sa[pis_inv[ssa[i].first]] = node;
      bounds[pis_inv[ssa[i].first]].first = ssa[i].first - esa[j-1];
      bounds[pis_inv[ssa[i].first]].second = esa[j] - ssa[i].first;
      i += 1;
    }
  }

  /*
   * get comments in on everything
   * with comments you get the pseudocode written, the paper written!
   * easier to revisit
   */

  // list_paths() finds all cycles in our graph.
  void list_paths(std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa) {
    std::vector<uint> indegrees(phi_inv_sa.size(), 0); // indegrees of the nodes
    std::vector<bool> visited(phi_inv_sa.size(), false); // visited nodes so far

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

        // min path length; this can be used to store a path that is long enough
        // if(is_cycle || currentPath.size() >= min_path_length) {
        // ds with bools that mark whether a node is in a path or not
        if(is_cycle) {
          rads_tree branch = rads_tree(current_path);
          paths.push_back(current_path);
        }
      }
    }
  }

  int get_size() {
    return phi_inv_sa.size();
  }

  int get_num_paths() {
    return paths.size();
  }

  int get_avg_path_l() {
    int path_length = 0;
    for (size_t i = 0; i < paths.size(); i++) {
      path_length += paths[i].size();
    }

    return (path_length / paths.size());
  }

protected:
  std::unordered_map<ulint, ulint> esa_map;
  std::unordered_map<ulint, ulint> pis_inv;

  std::vector<ulint> phi_inv_sa;
  std::vector<std::pair<ulint,ulint>> bounds;
  std::vector<vector<uint>> paths;
  std::vector<rads_tree> trees;
};
}

#endif
