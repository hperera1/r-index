#ifndef R_INDEX_R_A_
#define R_INDEX_R_A_

#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <definitions.hpp>
#include <r_index.hpp>

using namespace sdsl;

namespace ri {
class rads {
public:
  rads(){};
  rads(std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa) {
    init_by_value(ssa, esa);
  }
  // rads(r_index &ridx) {
  //   init_by_reference(ridx);
  // }

  // void init_by_reference(r_index &ridx) {
  //   end_map.reserve(ridx.samples_last.size());
  //   phi_inv_sa.resize(ridx.samples_last.size());
  //   for(ulint i = 0; i < ridx.samples_last.size(); i++) {
  //     esa_map[ridx.samples_last[i]] = i;
  //   }
  // }

  void init_by_value(std::vector<std::pair<ulint, ulint>> &ssa, std::vector<ulint> &esa) {
    assert(ssa.size() == esa.size());
    esa_map.reserve(esa.size());
    phi_inv_sa.resize(esa.size());
    for(ulint i = 0; i < esa.size(); i++) {
      esa_map[esa[i]] = i;
      phi_inv_sa[i] = ssa[i].second;
      pis_inv[ssa[i].second] = i;
    }
    std::sort(ssa.begin(), ssa.end());
    std::sort(esa.begin(), esa.end());
    esa.push_back(esa.back() + 1);
    ulint i = 0;
    ulint j = 0;
    ulint node = esa_map[esa.back()];
    while(i < ssa.size()) {
      while((i < ssa.size()) && (j < esa.size()) && (ssa[i].second < esa[j])) {
        phi_inv_sa[pis_inv[ssa[i].second]] = node;
        i += 1;
      }

      if(j < esa.size() - 1) {
        node = esa_map[esa[j]];
        j += 1;
      }
    }
  }

  void list_paths() {
    std::vector<int> indegrees(phi_inv_sa.size(), 0);
    std::vector<bool> visited(phi_inv_sa.size(), false);
    std::stack<int> sources;

    for(int i = 0; i < phi_inv_sa.size(); i++) {
      if(phi_inv_sa[i] >= 0)
        indegrees[phi_inv_sa[i]] += 1;
    }

    for(int i = 0; i < indegrees.size(); i++) {
      if(indegrees[i] == 0)
        sources.push(i);
    }

    while(sources.size() > 0) {
      std::vector<int> currentPath;
      int u = sources.top();
      int v = phi_inv_sa[u];
      sources.pop();
      visited[u] = true;
      currentPath.push_back(u); // pushes 7 onto path

      while(v >= 0 && (visited[v] == false)) {
        currentPath.push_back(v);
        visited[v] = true;
        v = phi_inv_sa[v];
      }

      if(v >= 0 && (visited[v] == false)) {
        sources.push(v);
      }

      paths.push_back(currentPath);
    }
  }

  // find, break cycles, and mark which set used to be part of a cycle
  void find_cycles() {
    for(int i = 0; i < paths.size(); i++) {
      std::unordered_set<int> curr_visited;
      for(int j = 0; j < paths[i].size(); j++) {
        int u = paths[i][j];
        if(curr_visited.find(u) == curr_visited.end()) // if not found
          curr_visited.insert(u);
        else {
          cycles.insert(i);
          paths[i].erase(paths[i].begin() + j);
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

  int get_num_cycles() {
    return cycles.size();
  }
protected:
  std::unordered_set<int> cycles;
  std::unordered_map<ulint, ulint> esa_map;
  std::unordered_map<ulint, ulint> pis_inv;
  std::vector<ulint> phi_inv_sa;
  std::vector<vector<int>> paths;
};
}

#endif
