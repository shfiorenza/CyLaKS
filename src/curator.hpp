#ifndef _CYLAKS_CURATOR_HPP_
#define _CYLAKS_CURATOR_HPP_
#include "definitions.hpp"
#include "filament_manager.hpp"
#include "file_utilities.hpp"
#include "protein_manager.hpp"
#include "rng_manager.hpp"
#include "system_parameters.hpp"

class Curator {
private:
  FileManager files_;
  ProteinManager proteins_;
  FilamentManager filaments_;

  SysTimepoint start_time_;

public:
  size_t verbosity_{0};
  bool sim_running_{true};
  bool sim_equilibrated_{true};

  char *sim_name_{nullptr};
  char *test_mode_{nullptr};

  size_t n_sim_species_{0};
  size_t n_sim_objs_{0};

  double t_tot_{0.0};
  double t_equil_{0.0};
  double t_snapshot_{0.0};

  size_t i_step_{0};
  size_t i_datapoint_{0};

  size_t n_steps_pre_equil_{0};
  size_t n_steps_equil_{0};
  size_t n_steps_tot_{0};
  size_t n_steps_snapshot_{0};

  SysParameters params_;
  RandomNumberManager gsl_;

private:
  void CheckArgs(char *agrv[]);
  void GenerateLogFile();
  void GenerateDataFiles();
  void ParseParameters();
  void SetLocalParameters();
  void InitializeSimObjects();
  void CheckEquilibration();
  void PrintProgress();
  void OutputData();
  void OutputSimDuration();
  void CloseDataFiles();

public:
  Curator(char *argv[]) {
    CheckArgs(argv);
    GenerateLogFile();
    ParseParameters();
    SetLocalParameters();
    InitializeSimObjects();
    GenerateDataFiles();
  }
  void InitializeSimulation(char *argv[], system_properties *properties,
                            system_parameters *parameters);

  void ErrorExit(const char *function_name);
  template <typename... Args> void Log(const char *msg, const Args... args) {
    // This is technically a horrendous vulnerability, but we don't care about
    // 'hackers' in our sim; also should never be linked to input
    int chars_printed{printf(msg, args..., "MISSING STRING")};
    int chars_written{fprintf(log_file_, msg, args..., "MISSING STRING")};
    if (chars_printed < 0 or chars_written < 0) {
      printf("Fatal error in Curator::Log()\n *** EXITING ***\n");
      fprintf(log_file_, "Fatal error in Curator::Log()\n *** EXITING ***\n");
      exit(1);
    }
  }
  // Overloaded for verbosity-based alerts
  template <typename... Args>
  void Log(int alert_lvl, const char *msg, const Args... args) {
    if (verbosity_ < alert_lvl) {
      return;
    }
    Log(msg, args...);
  }
  void EvolveSimulation();
  void PrintMicrotubules();
  void PrintMicrotubules(double pause_duration);
  void PauseSim(double duration);
  void StartDataCollection();
  void TerminateSimulation();
  void CleanUp();
};
#endif