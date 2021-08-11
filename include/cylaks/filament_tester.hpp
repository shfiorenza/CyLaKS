#ifndef _CYLAKS_FILAMENT_TESTER_HPP_
#define _CYLAKS_FILAMENT_TESTER_HPP_
#include "filament_manager.hpp"

class ProteinTester;

class FilamentTester : public FilamentManager {
protected:
  ProteinTester *proteins_{nullptr};

public:
protected:
public:
  FilamentTester() {}
  ~FilamentTester() {}
  void Initialize(ProteinTester *proteins);
  void UpdateUnoccupied() {
    if (up_to_date_) {
      return;
    }
    up_to_date_ = true;
    for (auto &&pop : unoccupied_) {
      pop.second.ZeroOut();
    }
    // Add sites to unoccupied_ and update weights
    for (auto &&site : sites_) {
      for (auto &&pop : unoccupied_) {
        pop.second.Sort(site);
      }
      // if (Sys::test_mode_ != "motor_lattice_step") {
      int n_neighbs{site->GetNumNeighborsOccupied()};
      site->SetWeight_Bind(Sys::weight_neighb_bind_[n_neighbs]);
      site->SetWeight_Unbind(Sys::weight_neighb_unbind_[n_neighbs]);
      // }
    }
    // if (Sys::test_mode_ != "motor_lattice_step") {
    UpdateLattice();
    // }
  }
};
#endif