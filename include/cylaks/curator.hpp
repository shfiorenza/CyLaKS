#ifndef _CYLAKS_CURATOR_HPP_
#define _CYLAKS_CURATOR_HPP_
#include "definitions.hpp"
#include "filament_manager.hpp"
#include "filament_tester.hpp"
#include "protein_manager.hpp"
#include "protein_tester.hpp"
#include "system_namespace.hpp"
#include "system_parameters.hpp"
#include "system_rng.hpp"

// Curator: Initializes and runs the simulation; controls data output
class Curator {
private:
  // Available test modes; initialize custom scenarios and kMC events
  Vec<Str> test_modes_{"xlink_bind_ii",       "xlink_diffusion",
                       "motor_lattice_bind",  "motor_lattice_step",
                       "filament_separation", "filament_ablation",
                       "hetero_tubulin",      "kinesin_mutant"};

  UMap<Str, Sys::DataFile> data_files_;

  size_t n_steps_per_snapshot_{0};
  size_t n_sites_max_{0};

  SysTimepoint start_time_;

public:
  ProteinManager proteins_;
  FilamentManager filaments_;

  ProteinTester test_proteins_;
  FilamentTester test_filaments_;

private:
  void CheckArgs(int argc, char *agrv[]);
  void GenerateLog();
  void ParseParameters();
  void InitializeSimulation();
  void GenerateDataFiles();

  void UpdateObjects();
  void CheckPrintProgress();
  void OutputData();

public:
  Curator(int argc, char *argv[]) {
    CheckArgs(argc, argv);
    GenerateLog();
    ParseParameters();
    InitializeSimulation();
    GenerateDataFiles();
  }
  void EvolveSimulation() {
    UpdateObjects();
    CheckPrintProgress();
    OutputData();
  }
};
#endif