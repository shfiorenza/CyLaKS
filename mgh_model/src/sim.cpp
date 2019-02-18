#include "master_header.h"

int main(int argc, char *argv[]){

	system_parameters parameters;
	system_properties properties;

	// Check that input has the correct number of arguments
	properties.wallace.CheckArguments(argv[0], argc);
	// Parse parameters from yaml file into simulation's parameter struct
	properties.wallace.ParseParameters(&parameters, argv[1]);
	// Initialize sim objects (MTs, kinesin, MAPs, etc.)
	properties.wallace.InitializeSimulation(&properties);
	// Generate data files
	properties.wallace.GenerateDataFiles(argv[2]);

	// Main KMC loop
	for(int i_step = 0; i_step < parameters.n_steps; i_step++){
		properties.wallace.UpdateTimestep(i_step);
		properties.kinesin4.RunKMC();
		properties.prc1.RunKMC();
		properties.prc1.RunDiffusion();
		properties.microtubules.RunDiffusion();
	}

	// Cleanup stuff
	properties.wallace.OutputSimDuration();
	properties.wallace.CloseDataFiles();
	properties.gsl.CleanUp();

	return 0;
}
