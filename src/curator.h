#ifndef _CURATOR_H
#define _CURATOR_H
#include <yaml-cpp/yaml.h>
#include <thread>
struct system_parameters;
struct system_properties;

class Curator{
	private:

	public:
		int data_threshold_;
		int range_of_data_;
		int n_pickup_;
		int equil_milestone_;
		int data_milestone_;
		double start_;
		double finish_;
		double sim_duration_;
		FILE *stream_;
		struct timespec pause_dur_;

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

		void CheckMPI(int threads_provided);
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
