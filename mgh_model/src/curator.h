#ifndef _CURATOR_H
#define _CURATOR_H
#include <iostream>
#include <sys/stat.h>
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
		double sim_duration_;
		clock_t start_, finish_;
		FILE *stream_;
		struct timespec pause_dur_;

		system_parameters *parameters_;
		system_properties *properties_;
	private:
	
	public:
		Curator();
		void InitializeSimulation(system_parameters *parameters, 
								  system_properties *properties);
		void SetParameters();
		void SetExperimentalStage();
		void OutputSimDetails();

		void OpenFiles(char *sim_name);

		void PrintMicrotubules();
		void PrintMicrotubules(double pause_duration);
		void OutputData();
		void UpdateTimestep(int i_step);
		void PauseSim(double duration); 
		void OutputSimDuration();
		void CleanUp();

		bool FileExists(std::string file_name);
};
#endif
