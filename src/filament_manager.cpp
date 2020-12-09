#include "filament_manager.hpp"
#include "curator.hpp"

void FilamentManager::GenerateFilaments() {

  list_.emplace_back(wally_, params_, 0, 0, 5.0);
}

void FilamentManager::UpdateUnoccupied() {

  if (up_to_date_) {
    return;
  }
  up_to_date_ = true;
  for (auto &&pop : unocc_) {
    pop.second.ZeroOut();
  }
  for (auto &&site : sites_) {
    if (site->occupant_ != nullptr) {
      continue;
    }
    unocc_["motors"].AddEntry(site);
    unocc_["xlinks"].AddEntry(site, site->GetNeighborCount());
  }
}

void FilamentManager::RunBD() {

  if (immobile_) {
    return;
  }
  // FIXME take from add_continous_mt_movement branch
}