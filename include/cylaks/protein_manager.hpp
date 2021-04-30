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

// ProteinManager: Initializes proteins and the kMC events that target them
//                 Keeps track of different protein populations as sim evolves
class ProteinManager {
private:
  Map<Str, Vec<Pair<size_t, size_t>>> test_stats_;
  Map<Str, Vec<double>> test_ref_;

  FilamentManager *filaments_{nullptr};

public:
  Reservoir<Motor> motors_;
  Reservoir<Protein> xlinks_;
  EventManager kmc_;

private:
  void GenerateReservoirs();
  void InitializeWeights();
  void SetParameters();
  void InitializeEvents();

  void FlagFilamentsForUpdate();
  void UpdateFilaments();

  void ReportTestStatistics();

  void SetTestMode();
  void InitializeTest_Filament_Ablation();
  void InitializeTest_Filament_Separation();
  void InitializeTest_Filament_HeteroTubulin();
  void InitializeTest_Motor_Heterodimer();
  void InitializeTest_Motor_LatticeStep();
  void InitializeTest_Motor_LatticeBind();
  void InitializeTest_Xlink_Diffusion();
  void InitializeTest_Xlink_Bind_II();

public:
  ProteinManager() {}
  ~ProteinManager() { ReportTestStatistics(); }
  void Initialize(FilamentManager *filaments) {
    filaments_ = filaments;
    GenerateReservoirs();
    InitializeWeights();
    SetParameters();
    InitializeEvents();
  }
  void InitializeTest(FilamentManager *filaments) {
    filaments_ = filaments;
    SetTestMode();
  }
  void UpdateLatticeDeformation() { motors_.UpdateLatticeDeformation(); }
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