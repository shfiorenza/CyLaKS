#include "cylaks/binding_site.hpp"
#include "cylaks/binding_head.hpp"
#include "cylaks/protofilament.hpp"

int BindingSite::GetNumNeighborsOccupied() {

  if (_n_neighbs_max == 0) {
    return 0;
  }
  int n_neighbs{0};
  for (auto const &site : neighbors_) {
    if (site->occupant_ != nullptr) {
      if (occupant_ != nullptr) {
        if (occupant_->parent_ != site->occupant_->parent_) {
          n_neighbs++;
        }
      } else {
        n_neighbs++;
      }
    }
  }
  return n_neighbs;
}

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

void BindingSite::AddForce(Vec<double> f) { filament_->AddForce(this, f); }
void BindingSite::AddTorque(double tq) { filament_->AddTorque(tq); }