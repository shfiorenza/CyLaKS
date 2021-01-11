#ifndef _CYLAKS_PROTEIN_MANAGER_HPP_
#define _CYLAKS_PROTEIN_MANAGER_HPP_
#include "event_manager.hpp"
#include "motor.hpp"
#include "protein.hpp"
#include "reservoir.hpp"
#include "system_namespace.hpp"
#include "system_parameters.hpp"
#include "system_rng.hpp"

class FilamentManager;

class ProteinManager {
private:
  Map<Str, Vec<double>> test_ref_;
  Map<Str, Vec<Pair<size_t, size_t>>> test_stats_;

  FilamentManager *filaments_{nullptr};

public:
  Reservoir<Motor> motors_;
  Reservoir<Protein> xlinks_;
  EventManager kmc_;

private:
  void GenerateReservoirs();
  void InitializeWeights();
  void SetParameters();
  void InitializeTestEnvironment();
  void InitializeTestEvents();
  void InitializeEvents();

  void FlagFilamentsForUpdate();
  void UpdateFilaments();

public:
  ProteinManager() {}
  ~ProteinManager() {
    for (auto const &entry : test_stats_) {
      printf("For event %s:\n", entry.first.c_str());
      for (int index{0}; index < entry.second.size(); index++) {
        auto stats = entry.second[index];
        double p{double(stats.first) / stats.second};
        double ref{test_ref_.at(entry.first)[index]};
        printf("  p[%i] = %.3g (%.3g expected) [%zu events]\n", index, p, ref,
               stats.first);
      }
    }
  }
  void Initialize(FilamentManager *filaments) {
    filaments_ = filaments;
    GenerateReservoirs();
    InitializeWeights();
    SetParameters();
    if (!Sys::test_mode_.empty()) {
      InitializeTestEnvironment();
      InitializeTestEvents();
      return;
    }
    InitializeEvents();
  }
  void UpdateLatticeDeformation() {}
  void UpdateExtensions() {
    bool forced_unbind{xlinks_.UpdateExtensions()};
    if (forced_unbind) {
      xlinks_.FlagForUpdate();
      FlagFilamentsForUpdate();
    }
  }
  void RunKMC() {
    UpdateFilaments();
    motors_.PrepForKMC();
    xlinks_.PrepForKMC();
    kmc_.ExecuteEvents();
  }
};

#endif