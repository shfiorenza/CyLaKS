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
protected:
  FilamentManager *filaments_{nullptr};

public:
  Reservoir<Motor> motors_;
  Reservoir<Protein> xlinks_;
  EventManager kmc_;

protected:
  void FlagFilamentsForUpdate();
  virtual void UpdateFilaments();

  void GenerateReservoirs();
  void InitializeWeights();
  void SetParameters();
  void InitializeEvents();

public:
  ProteinManager() {}
  virtual ~ProteinManager() {}
  virtual void Initialize(FilamentManager *filaments) {
    filaments_ = filaments;
    GenerateReservoirs();
    InitializeWeights();
    SetParameters();
    InitializeEvents();
  }
  void UpdateLatticeDeformation() { motors_.UpdateLatticeDeformation(); }
  void UpdateExtensions() {
    bool motors_forced_unbind{motors_.UpdateExtensions()};
    bool xlinks_forced_unbind{xlinks_.UpdateExtensions()};
    if (motors_forced_unbind or xlinks_forced_unbind) {
      motors_.FlagForUpdate();
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
  void RunBD() {}
};
#endif