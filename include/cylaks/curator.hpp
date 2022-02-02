#ifndef _CYLAKS_CURATOR_HPP_
#define _CYLAKS_CURATOR_HPP_
#include "filament_tester.hpp"
#include "protein_tester.hpp"
#include "system_definitions.hpp"
#include "system_namespace.hpp"
#include "system_parameters.hpp"
#include "system_rng.hpp"

// Curator: Initializes and runs the simulation; controls data output
class Curator {
private:
  int n_tests_{4};
  // Available test modes; initialize custom scenarios and kMC events
  Vec<Str> test_modes_{"xlink_bind_ii",      "xlink_diffusion",
                       "motor_lattice_bind", "motor_lattice_step",
                       "hetero_tubulin",     "kinesin_mutant",
                       "filament_ablation",  "filament_separation"};
  Vec<Str> test_param_files_{"overlap", "prc1",  "kif4a",  "kif4a",
                             "k401",    "kif4a", "endtag", "overlap"};
  size_t n_sites_max_{0}; // Largest microtubule length (for padding data)
  size_t n_steps_per_snapshot_{0}; // kMC steps per data output
  SysTimepoint start_time_;
  UMap<Str, Sys::DataFile> data_files_;

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