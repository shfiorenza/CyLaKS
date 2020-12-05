#include "filament_manager.hpp"

void FilamentManager::GenerateFilaments() {}

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

void FilamentManager::RunBD() {}