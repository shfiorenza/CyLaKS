#ifndef _CURATOR_H
#define _CURATOR_H

#include <iostream>
#ifndef _PARAMETERS_H
typedef struct system_parameters system_parameters;
#endif
#ifndef _SYSTEM_PROPERTIES_H
typedef struct system_properties system_properties;
#endif

class Curator{
	private:

	public:
		system_parameters *parameters_;
		system_properties *properties_;
	private:
	
	public:
		Curator();
		void Initialize(system_parameters *parameters, 
						system_properties *properties);

		void OutputSimDetails();
		void PrintMicrotubules();
		void OutputData(FILE *occupancy_file, FILE *ID_file);
};
#endif
