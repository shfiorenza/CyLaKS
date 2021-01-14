#ifndef _CYLAKS_SYSTEM_NAMESPACE_HPP_
#define _CYLAKS_SYSTEM_NAMESPACE_HPP_
#include <cstring>
#include <filesystem>

namespace Sys {

inline std::string sim_name_;
inline std::string test_mode_;
inline std::string yaml_file_;

inline size_t ablation_step_{0};

inline FILE *log_file_;
inline size_t verbosity_{0};

inline bool running_{true};
inline bool equilibrating_{true};

inline size_t n_unique_objects_{0};
inline size_t n_unique_species_{0};

inline size_t n_steps_pre_equil_{0};
inline size_t n_steps_equil_{0};
inline size_t n_steps_run_{0};

inline size_t i_step_{0};
inline size_t i_datapoint_{0};

inline std::vector<double> weight_neighb_bind_;   // [n_neighbs]
inline std::vector<double> weight_neighb_unbind_; // [n_neighbs]

inline size_t lattice_cutoff_;
inline std::vector<double> weight_lattice_bind_;       // [delta]
inline std::vector<double> weight_lattice_unbind_;     // [delta]
inline std::vector<double> weight_lattice_bind_max_;   // [n_neighbs]
inline std::vector<double> weight_lattice_unbind_max_; // [n_neighbs]

template <typename... Args>
inline void Log(const char *msg, const Args... args) {
  // Tag each log pintout in terminal w/ simulation name
  printf("[%s] ", sim_name_.c_str());
  // This is technically a horrendous vulnerability, but we don't care about
  // 'hackers' in our sim; also Log() is never explicitly linked to input
  int chars_printed{printf(msg, args..., "MISSING STRING")};
  int chars_written{fprintf(log_file_, msg, args..., "MISSING STRING")};
  if (chars_printed < 0 or chars_written < 0) {
    printf("Fatal error in Curator::Log()\n *** EXITING ***\n");
    fprintf(log_file_, "Fatal error in Curator::Log()\n *** EXITING ***\n");
    exit(1);
  }
}
template <typename... Args>
inline void Log(size_t tier, const char *msg, const Args... args) {
  if (verbosity_ < tier) {
    return;
  }
  Log(msg, args...);
}
inline void ErrorExit(const char *fn_name) {
  Log("Fatal error in simulation.\n");
  Log("   Function name: %s\n", fn_name);
  Log("   I_STEP = %zu\n", i_step_);
  Log("   N_STEPS = %zu\n", i_step_ - n_steps_equil_);
  Log("   N_DATAPOINTS = %zu\n", i_datapoint_);
  exit(1);
}
inline void EarlyExit() {
  Log("Run terminated after sufficient data collection\n");
  Log("   N_STEPS = %zu\n", i_step_ - n_steps_equil_);
  Log("   N_DATAPOINTS = %zu\n", i_datapoint_);
  running_ = false;
}
}; // namespace Sys
#endif