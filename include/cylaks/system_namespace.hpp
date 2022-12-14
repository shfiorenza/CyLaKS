#ifndef _CYLAKS_SYSTEM_NAMESPACE_HPP_
#define _CYLAKS_SYSTEM_NAMESPACE_HPP_
#include <filesystem>
#include <string>

// Sys namespace: Used to easily share variables and functions across CyLaKS
// classes
namespace Sys {

inline std::string sim_name_;  // Name of simulation
inline std::string test_mode_; // Name of test mode, if any
inline std::string yaml_file_; // Name of parameter file

inline int teth_x_min_; // temporary trash variable
inline int teth_x_max_; // temporary trash variable

inline FILE *log_file_;      // Pointer to log fle
inline size_t verbosity_{0}; // How much info is output to log; higher is more

inline bool running_{true};       // If main kmc-bd loop is running
inline bool equilibrating_{true}; // If proteins are still equilibrating

inline size_t n_objects_{0}; // No. of total objects across all species
inline size_t n_species_{0}; // No. of different protein species

inline size_t n_steps_pre_equil_{0}; // kMC-BD steps before proteins are checked
inline size_t n_steps_equil_{0};     // kMC-BD steps to complete protein equil.
inline size_t n_steps_run_{0};       // kMC-BD steps to run post-equil.

inline size_t i_step_{0};      // Current kMC-BD step
inline size_t i_datapoint_{0}; // Current datapoint index
inline size_t n_datapoints_max_{0};

inline size_t n_runs_recorded_{0}; // No. of motor runs post-equilibration
inline bool early_exit_triggered_{false}; // Exit after next data output?

// Auxiliary variables related to test_modes
inline size_t ablation_step_{0}; // For 'filament_ablation'

inline int n_xlinks_{-1};       // For 'filament_separation' and 'forced_slide'
inline int binding_active_{-1}; // bool+1; -1 for null; 0 for false; 1 for true
inline bool constant_velocity_{true};
inline double slide_velocity_{-1.0}; // For 'forced_slide'
inline int i_pause_{-1};             // when to pause force clamp
inline int i_resume_{-1};            // when to resume force clamp
inline bool rescale_times_{false};

inline double p_mutant_{-1.0};         // For 'hetero_tubulin'
inline double binding_affinity_{-1.0}; // For 'hetero_tubulin'

// Simple look-up tables for Boltzmann factors used throughout simulation
inline std::vector<double> weight_neighb_bind_;   // index: [n_neighbs]
inline std::vector<double> weight_neighb_unbind_; // index: [n_neighbs]
inline size_t lattice_cutoff_; // Range of long-range Gaussian in n_sites
inline std::vector<double> weight_lattice_bind_;       // index: [delta]
inline std::vector<double> weight_lattice_unbind_;     // index: [delta]
inline std::vector<double> weight_lattice_bind_max_;   // index: [n_neighbs]
inline std::vector<double> weight_lattice_unbind_max_; // index: [n_neighbs]

// Log function: prints to both terminal and sim_name.log file
template <typename... Args>
inline void Log(const char *msg, const Args... args) {
  // Tag each log pintout in terminal w/ simulation name
  printf("[%s] ", sim_name_.c_str());
  // This is technically a horrendous vulnerability, but we don't really care
  // about 'hackers' in our sim; also Log() is never explicitly linked to input
  int chars_printed{printf(msg, args..., "MISSING STRING")};
  int chars_written{fprintf(log_file_, msg, args..., "MISSING STRING")};
  if (chars_printed < 0 or chars_written < 0) {
    printf("Fatal error in Curator::Log()\n *** EXITING ***\n");
    fprintf(log_file_, "Fatal error in Curator::Log()\n *** EXITING ***\n");
    exit(1);
  }
}
// Log function but w/ verbosity: only log msg if set verbosity is high enough
template <typename... Args>
inline void Log(size_t tier, const char *msg, const Args... args) {
  if (verbosity_ < tier) {
    return;
  }
  Log(msg, args...);
}
// If an error occurs, report function name and stats needed for analysis
inline void ErrorExit(const char *fn_name) {
  Log("Fatal error in simulation.\n");
  Log("   Function name: %s\n", fn_name);
  Log("   I_STEP = %zu\n", i_step_);
  Log("   N_STEPS = %zu\n", i_step_ - n_steps_equil_);
  Log("   N_DATAPOINTS = %zu\n", i_datapoint_);
  exit(1);
}
// After sufficient number of kinesin runs, report stats needed for analysis
inline void EarlyExit() {
  Log("Run terminated after sufficient data collection (%i kinesin runs)\n",
      n_runs_recorded_);
  Log("   N_STEPS = %zu\n", i_step_ - n_steps_equil_);
  Log("   N_DATAPOINTS = %zu\n", i_datapoint_);
  running_ = false;
}
template <typename DATA_T1, typename DATA_T2>
inline void OverrideParam(std::string name, DATA_T1 *param, DATA_T2 value) {
  if (*param != value) {
    Log("    %s = %s\n", name.c_str(), std::to_string(value).c_str());
    *param = value;
  }
}
/* Useful data structures */
struct DataFile {
  std::string name_{"example"};
  std::string filename_{"simName_example.file"};
  FILE *fileptr_;
  DataFile() {}
  DataFile(std::string name) : name_{name} {
    filename_ = Sys::sim_name_ + "_" + name_ + ".file";
    fileptr_ = fopen(filename_.c_str(), "w");
    if (fileptr_ == nullptr) {
      printf("Error; cannot open '%s'\n", filename_.c_str());
      exit(1);
    }
  }
  template <typename DATA_T> void Write(DATA_T *array, size_t count) {
    size_t n_chars_written{fwrite(array, sizeof(DATA_T), count, fileptr_)};
    if (n_chars_written < count) {
      printf("Error writing to '%s'\n", filename_.c_str());
      exit(1);
    }
  }
};
struct ProbEntry {
  std::string event_name_;
  double val_;
  std::vector<std::vector<std::vector<double>>> vals_;
  ProbEntry() {}
  ProbEntry(Str name, double val) : event_name_{name}, val_{val} {}
  ProbEntry(Str name, std::vector<std::vector<std::vector<double>>> vals)
      : event_name_{name}, vals_{vals} {}
  double GetVal() { return val_; }
  double GetVal(size_t i) { return vals_[0][0][i]; }
  double GetVal(size_t i, size_t j) { return vals_[0][i][j]; }
  double GetVal(size_t i, size_t j, size_t k) { return vals_[i][j][k]; }
};
struct BoltzmannFactor {
  std::string effect_name_;
  int i_start_{0};
  size_t size_{0};
  Vec<double> bind_;
  Vec<double> unbind_;
  BoltzmannFactor() {}
  BoltzmannFactor(std::string name, size_t sz) : BoltzmannFactor(name, sz, 0) {}
  BoltzmannFactor(std::string name, size_t size, int i_start)
      : effect_name_{name}, size_{size}, i_start_{i_start} {
    bind_.resize(size);
    unbind_.resize(size);
  }
};
}; // namespace Sys
#endif