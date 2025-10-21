#ifndef _CYLAKS_FILAMENT_TESTER_HPP_
#define _CYLAKS_FILAMENT_TESTER_HPP_
#include "filament_manager.hpp"

class ProteinTester;

class FilamentTester : public FilamentManager {
protected:
  ProteinTester *proteins_{nullptr};

  Vec<double> recorded_force_;

public:
protected:
  void UpdateForces();

public:
  FilamentTester() {}
  ~FilamentTester() {
    if (Sys::test_mode_.empty()) {
      return;
    }
    double avg_force{0.0};
    for (int i_step{0}; i_step < recorded_force_.size(); i_step++) {
      avg_force += recorded_force_[i_step];
    }
    avg_force = avg_force / recorded_force_.size();
    double var_force{0.0};
    for (int i_step{0}; i_step < recorded_force_.size(); i_step++) {
      double var{recorded_force_[i_step] - avg_force};
      var_force += var * var;
    }
    var_force = std::sqrt(var_force) / recorded_force_.size();
    Sys::Log("Avg force on mobile MT: %g +/- %g pN\n", avg_force, var_force);
  }
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
      int n_neighbs{site->GetNumNeighborsOccupied_Motor()};
      site->SetWeight_Bind(Sys::weight_neighb_bind_[n_neighbs]);
      site->SetWeight_Unbind(Sys::weight_neighb_unbind_[n_neighbs]);
    }
    UpdateLattice();
  }
};
#endif