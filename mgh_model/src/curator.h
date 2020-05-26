#ifndef _CURATOR_
#define _CURATOR_
#include <chrono>
#include <string.h>
#include <yaml-cpp/yaml.h>
struct system_parameters;
struct system_properties;

using sys_clock = std::chrono::steady_clock;
using sys_timepoint = sys_clock::time_point;
class Curator {
private:
  FILE *log_file_{nullptr};
  unsigned long data_threshold_{0};
  unsigned long n_steps_recorded_{0};
  unsigned long n_steps_per_output_{0};
  unsigned long equil_milestone_{0};
  unsigned long data_milestone_{0};

  system_parameters *parameters_;
  system_properties *properties_;

public:
  int verbosity_{0};
  double t_motors_[4];
  double t_xlinks_[4];
  double t_MTs_[4];
  char *param_file_{nullptr};
  char *sim_name_{nullptr};
  char *test_mode_{nullptr};

  struct timespec pause_dur_;

  sys_timepoint start_;
  sys_timepoint finish_;

private:
  FILE *OpenFile(const char *file_name, const char *type);
  bool FileExists(std::string file_name);
  void CheckArgs(char *agrv[]);
  void GenerateLogFile();
  void GenerateDataFiles();
  void ParseParameters();
  void SetLocalParameters();
  void InitializeSimObjects();
  void OutputData();
  void OutputSimDuration();
  void CloseDataFiles();

public:
  Curator();
  void InitializeSimulation(char *argv[], system_properties *properties,
                            system_parameters *parameters);
  void ErrorExit(const char *function_name);
  template <typename... Args> void Log(const char *msg, const Args... args) {
    // This is technically a horrendous vulnerability, but we don't care about
    // 'hackers' in our sim; also should never be linked to input
    int chars_printed{printf(msg, args..., "MISSING STRING")};
    int chars_written{fprintf(log_file_, msg, args..., "MISSING STRING")};
    if (chars_printed < 0 or chars_written < 0) {
      printf("Fatal error in Curator::Log()\n");
      printf(" *** EXITING ***\n");
      fprintf(log_file_, "Fatal error in Curator::Log()\n");
      fprintf(log_file_, " *** EXITING ***\n");
      exit(1);
    }
  }
  template <typename... Args>
  // Overloaded for verbosity-based alerts
  void Log(int alert_level, const char *msg, const Args... args) {
    if (verbosity_ < alert_level) {
      return;
    }
    Log(msg, args...);
  }
  void UpdateTimestep(unsigned long i_step);
  void PrintMicrotubules();
  void PrintMicrotubules(double pause_duration);
  void PauseSim(double duration);
  void CleanUp();
};
#endif