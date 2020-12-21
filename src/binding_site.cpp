#include "binding_site.hpp"
#include "protofilament.hpp"

BindingSite *BindingSite::GetNeighbor(int dir) {

  if (dir != 1 and dir != -1) {
    Sys::ErrorExit("BindingSite::GetNeighb()");
  }
  for (auto const &neighb : neighbors_) {
    if (index_ + dir == neighb->index_) {
      return neighb;
    }
  }
  return nullptr;
}