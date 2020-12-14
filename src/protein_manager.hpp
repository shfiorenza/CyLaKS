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
  UMap<Str, Vec<double>> test_ref_;
  UMap<Str, Vec<Pair<size_t, size_t>>> test_stats_;

  SysRNG *gsl_{nullptr};
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

public:
  ProteinManager() {}
  void Initialize(SysRNG *gsl, FilamentManager *filaments) {
    gsl_ = gsl;
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
  void UpdateExtensions() {}
  void RunKMC() {
    motors_.PrepForKMC();
    xlinks_.PrepForKMC();
    kmc_.ExecuteEvents();
  }
};

#endif