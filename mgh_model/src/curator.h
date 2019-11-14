#ifndef _CURATOR_
#define _CURATOR_
#include <chrono>
#include <yaml-cpp/yaml.h>
struct system_parameters;
struct system_properties;

using sys_clock = std::chrono::steady_clock;
using sys_timepoint = sys_clock::time_point;
// FIXME; depricate this (used in kinesin & mt MGMT)
using t_unit = std::chrono::duration<double, std::nano>;

class Curator {
private:
  int data_threshold_{0};
  int n_steps_recorded_{0};
  int n_steps_per_output_{0};
  int equil_milestone_{0};
  int data_milestone_{0};

public:
  double t_motors_[4];
  double t_xlinks_[8];
  double t_MTs_[4];

  struct timespec pause_dur_;

  sys_timepoint start_;
  sys_timepoint finish_;

  system_parameters *parameters_;
  system_properties *properties_;

private:
  FILE *OpenFile(const char *file_name, const char *type);
  bool FileExists(std::string file_name);
  void CheckArgs(char *exe_name, int argc);
  void GenerateLogFile(char *sim_name);
  void GenerateDataFiles(char *sim_name);
  void ParseParameters(system_parameters *params, char *param_file);
  void SetLocalParameters();
  void InitializeSimObjects();
  void OutputData();
  void OutputSimDuration();
  void CloseDataFiles();

public:
  Curator();
  void InitializeSimulation(char *exe_name, char *param_file, char *sim_name,
                            int argc, system_properties *properties,
                            system_parameters *parameters);
  void ErrorExit(const char *function_name);
  template <typename... Args> void Log(const char *msg, const Args... args);
  void UpdateTimestep(int i_step);
  void PrintMicrotubules();
  void PrintMicrotubules(double pause_duration);
  void PauseSim(double duration);
  void CleanUp();
};
#endif