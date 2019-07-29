#ifndef _CURATOR_H
#define _CURATOR_H
#include <chrono>
#include <yaml-cpp/yaml.h>
struct system_parameters;
struct system_properties;

using sys_clock = std::chrono::steady_clock;
using sys_time = sys_clock::time_point;
using t_unit = std::chrono::duration<double, std::nano>;
const double n_per_sec_ = 1e9; 

class Curator{
	private:

	public:
		int data_threshold_;
		int range_of_data_;
		int n_pickup_;
		int equil_milestone_;
		int data_milestone_;

		double sim_duration_; 	// For all, indices correspond to:
		double t_motors_[4]; 		//  - 0: net time
		double t_xlinks_[5]; 	//  - 1: time spent updating
		double t_MTs_[4];			//  - 3: time spent executing events

		struct timespec pause_dur_;

		sys_time start_;
		sys_time finish_;

		system_parameters *parameters_;
		system_properties *properties_;

	private:
		void OpenLog(char *sim_name);
		void CheckArguments(char *exe_name, int argc);
		void ParseParameters(system_parameters *params, char *param_file); 
		void SetParameters();
		void SetExperimentalStage();
		void GenerateDataFiles(char *sim_name);
		bool FileExists(std::string file_name);
		FILE* OpenFile(const char *file_name, const char *type);  

		void OutputData();
		void OutputSimDuration();
		void CloseDataFiles();

	public:
		Curator();
		void InitializeSimulation(char *exe_name, char *param_file, 
				char *sim_name, int argc, system_properties *properties, 
				system_parameters *parameters);

		void Log(const char *msg);
		template <typename T> void Log(const char *msg, T arg);
		template <typename T1, typename T2> void Log(const char *msg,
				T1 arg1, T2 arg2);
		void UpdateTimestep(int i_step);
		void PrintMicrotubules();
		void PrintMicrotubules(double pause_duration);
		void PauseSim(double duration); 

		void CleanUp();
};
#endif
