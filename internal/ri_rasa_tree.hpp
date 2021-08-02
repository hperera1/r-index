#ifndef R_INDEX_R_A_TREE_
#define R_INDEX_R_A_TREE_

#include <definitions.hpp>
#include <r_index.hpp>

using namespace sdsl;

namespace ri {
class rads_tree {
public:
  std::vector<uint> path;
  std::vector<std::pair<uint,uint>> bounds;

  rads_tree(){};
  rads_tree(std::vector<uint> &cycle){
    
  }
};
}

#endif
