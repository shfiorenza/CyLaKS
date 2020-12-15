#ifndef _CYLAKS_FILAMENT_MANAGER_HPP_
#define _CYLAKS_FILAMENT_MANAGER_HPP_
#include "population.hpp"
#include "protofilament.hpp"
#include "system_namespace.hpp"
#include "system_parameters.hpp"
#include "system_rng.hpp"

class ProteinManager;

class FilamentManager {
private:
  bool up_to_date_{false};

  size_t n_bd_iterations_{0};
  double dt_eff_{0.0};

  Vec<BindingSite *> sites_;

  ProteinManager *proteins_{nullptr};

public:
  bool mobile_{false};
  Vec<Protofilament> proto_;

  UMap<Str, Population<BindingSite>> unoccupied_;

private:
  void SetParameters();
  void GenerateFilaments();
  void InitializeTestEnvironment();

  bool NoMobileFilamentsYet();

  void UpdateProteins();
  void UpdateLattice();

public:
  FilamentManager() {}
  void Initialize(ProteinManager *proteins) {
    proteins_ = proteins;
    if (!Sys::test_mode_.empty()) {
      return;
    }
    SetParameters();
    GenerateFilaments();
  }
  void FlagForUpdate() { up_to_date_ = false; }
  void UpdateUnoccupied() {
    if (up_to_date_) {
      return;
    }
    up_to_date_ = true;
    for (auto &&pop : unoccupied_) {
      pop.second.ZeroOut();
    }
    for (auto &&site : sites_) {
      if (site->occupant_ != nullptr) {
        continue;
      }
      unoccupied_["motors"].AddEntry(site);
      unoccupied_["xlinks"].AddEntry(site, site->GetNeighborCount());
    }
    UpdateLattice();
  }
  void RunBD() {
    if (NoMobileFilamentsYet()) {
      return;
    }
    for (int i_itr{0}; i_itr < n_bd_iterations_; i_itr++) {
      UpdateProteins();
      for (auto &&filament : proto_) {
        filament.UpdatePosition();
      }
    }
  }
};
#endif