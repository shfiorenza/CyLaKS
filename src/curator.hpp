#ifndef _CYLAKS_CURATOR_HPP_
#define _CYLAKS_CURATOR_HPP_
#include "filament_manager.hpp"
#include "protein_manager.hpp"
#include "system_files.hpp"
#include "system_parameters.hpp"

class Curator {
private:
  char *param_file_{nullptr};
  char *sim_name_{nullptr};

  size_t n_steps_pre_equil_{0};
  size_t n_steps_equil_{0};
  size_t n_steps_run_{0};
  size_t n_steps_snapshot_{0};
  size_t verbosity_{0};

  ProteinManager proteins_;
  FilamentManager filaments_;
  FileManager files_;
  SysTimepoint start_time_;

public:
  char *test_mode_{nullptr};

  bool sim_running_{true};
  bool sim_equilibrating_{true};

  size_t n_sim_objs_{0};
  size_t n_sim_species_{0};

  size_t i_step_{0};
  size_t i_datapoint_{0};

  SysParameters params_;
  RandomNumberGenerator gsl_;

private:
  void CheckArgs(char *agrv[]);
  void ParseParameters();
  void InitializeSimulation();
  void GenerateDataFiles();
  void CheckPrintProgress();
  void OutputData();

public:
  Curator(char *argv[]) {
    CheckArgs(argv);
    files_.GenerateLog(sim_name_);
    ParseParameters();
    InitializeSimulation();
    GenerateDataFiles();
  }
  void ErrorExit(const char *function_name) {
    Log("\nFatal error in %s\n", sim_name_);
    Log("   Function name: %s\n", function_name);
    Log("   N_STEPS = %zu\n", i_step_ - n_steps_equil_);
    Log("   N_DATAPOINTS = %zu\n", i_datapoint_);
    exit(1);
  }
  // Overloaded log function for verbosity-based alerts
  template <typename... Args>
  void Log(size_t alert_lvl, const char *msg, const Args... args) {
    if (verbosity_ < alert_lvl) {
      return;
    }
    Log(msg, args...);
  }
  // Log function; writes to both terminal and .log file
  template <typename... Args> void Log(const char *msg, const Args... args) {
    // This is technically a horrendous vulnerability, but we don't care about
    // 'hackers' in our sim; also should never be linked to input
    int chars_printed{printf(msg, args..., "MISSING STRING")};
    int chars_written{fprintf(files_.log_, msg, args..., "MISSING STRING")};
    if (chars_printed < 0 or chars_written < 0) {
      printf("Fatal error in Curator::Log()\n *** EXITING ***\n");
      fprintf(files_.log_, "Fatal error in Curator::Log()\n *** EXITING ***\n");
      exit(1);
    }
  }
  void EvolveSimulation();
  void TerminateSimulation();
};
#endif