#ifndef _CYLAKS_FILAMENT_MANAGER_HPP_
#define _CYLAKS_FILAMENT_MANAGER_HPP_
#include "population.hpp"
#include "protofilament.hpp"
#include "system_namespace.hpp"
#include "system_parameters.hpp"
#include "system_rng.hpp"

class ProteinManager;

// FilamentManager: Manages the brownian dynamics and binding sites of filaments
class FilamentManager {
protected:
  bool up_to_date_{false};

  // Some temporary hacky stuff for WCA potential
  double sigma_{25.0};    // nm
  double epsilon_{1.0};   // kbT
  double threshold_{0.0}; // nm

  size_t n_bd_iterations_{0};
  double dt_eff_{0.0};

  ProteinManager *proteins_{nullptr};

public:
  bool mobile_{false};
  Vec<Protofilament> protofilaments_;
  Vec<BindingSite *> sites_;

  Map<Str, Population<Object>> unoccupied_;

protected:
  void SetParameters();
  void GenerateFilaments();

  bool AllFilamentsImmobile();

  virtual void UpdateForces();
  virtual void UpdateLattice();

public:
  FilamentManager() {}
  virtual ~FilamentManager() {}
  virtual void Initialize(ProteinManager *proteins) {
    proteins_ = proteins;
    SetParameters();
    GenerateFilaments();
  }
  void AddPop(Str name, Fn<Vec<Object *>(Object *)> sort) {
    unoccupied_.emplace(name, Population<Object>(name, sort, sites_.size()));
  }
  void AddPop(Str name, Fn<Vec<Object *>(Object *)> sort, Vec<size_t> i_size,
              Vec<int> i_min, Fn<Vec<int>(Object *)> get_i) {
    Vec<size_t> sz{i_size[0], i_size[1], i_size[2], sites_.size()};
    unoccupied_.emplace(name, Population<Object>(name, sort, sz, i_min, get_i));
  }
  void FlagForUpdate() { up_to_date_ = false; }
  virtual void UpdateUnoccupied() {
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
      int n_neighbs{site->GetNumNeighborsOccupied()};
      site->SetWeight_Bind(Sys::weight_neighb_bind_[n_neighbs]);
      site->SetWeight_Unbind(Sys::weight_neighb_unbind_[n_neighbs]);
    }
    UpdateLattice();
  }
  void RunBD() {
    if (AllFilamentsImmobile()) {
      return;
    }
    for (int i_itr{0}; i_itr < n_bd_iterations_; i_itr++) {
      UpdateForces();
      for (auto &&entry : protofilaments_) {
        entry.UpdatePosition();
      }
    }
  }
};
#endif