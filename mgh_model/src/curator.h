#ifndef _CURATOR_H
#define _CURATOR_H
#include <chrono>
#include <yaml-cpp/yaml.h>
struct system_parameters;
struct system_properties;

using sys_clock = std::chrono::high_resolution_clock;
using sys_time = sys_clock::time_point;
using t_microsec = std::chrono::microseconds;

class Curator{
	private:

	public:
		int data_threshold_;
		int range_of_data_;
		int n_pickup_;
		int equil_milestone_;
		int data_milestone_;

		double sim_duration_; 			// For all, indices correspond to:
		double t_motors_[4]; 			//  - 0: net time
		double t_xlinks_kmc_[4]; 		//  - 1: time spent updating
		double t_xlinks_dif_[4]; 		//  - 2: time spent correcting stats
		double t_MTs_[4];				//  - 3: time spent executing events

//		FILE *stream_;

		struct timespec pause_dur_;

		sys_time start_;
		sys_time finish_;

		system_parameters *parameters_;
		system_properties *properties_;

	private:
	
	public:
		Curator();
		void ParseParameters(system_parameters *parameters, 
							 char *param_file); 
		void InitializeSimulation(system_properties *properties);
		void SetParameters();
		void SetExperimentalStage();

		void CheckArguments(char *sim_name, int argc);

		void GenerateDataFiles(char *sim_name);
		FILE* OpenFile(const char *file_name, const char *type);  
		bool FileExists(std::string file_name);

		void PrintMicrotubules();
		void PrintMicrotubules(double pause_duration);
		void OutputData();
		void UpdateTimestep(int i_step);
		void PauseSim(double duration); 
		void OutputSimDuration();
		void CloseDataFiles();
};
#endif
