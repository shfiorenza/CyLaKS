#include "cylaks/binding_site.hpp"
#include "cylaks/protofilament.hpp"

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