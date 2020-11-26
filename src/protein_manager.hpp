#ifndef _CYLAKS_PROTEIN_MANAGER_HPP_
#define _CYLAKS_PROTEIN_MANAGER_HPP_
#include "event_manager.hpp"
#include "motor.hpp"
#include "protein.hpp"
#include "reservoir.hpp"

class Curator;
class FilamentManager;
class RandomNumberManager;
struct SysParameters;

class ProteinManager {
private:
  UMap<Str, Vec<double>> test_ref_;
  UMap<Str, Vec<Pair<size_t, size_t>>> test_stats_;
  // WALLACE, MISTA
  Curator *wally_{nullptr};
  // Pointer to class that manages GSL functions (RNG, sampling, etc.)
  RandomNumberManager *gsl_{nullptr};
  // Pointers to global system params & props; same for all classes
  SysParameters *params_{nullptr};

public:
  EventManager kmc_;
  Reservoir<Motor> motors_;
  Reservoir<Protein> xlinks_;
  FilamentManager *filaments_;

private:
  void GenerateReservoirs();
  void InitializeWeights();
  void SetParameters();
  void InitializeTestEnvironment();
  void InitializeTestEvents();
  void InitializeEvents();

public:
  ProteinManager();
  void Initialize(Curator *wallace);
  void UpdateLatticeDeformation();
  void UpdateExtensions();
  void RunKMC();
};

#endif