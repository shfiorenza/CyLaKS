/* Main simulation file */
#include "master_header.h"

int main(int argc, char *argv[]){

	char param_file[160];
   
	system_parameters parameters;
	system_properties properties;

	// Check that the user has input the correct number of arguments 
	if(argc != 3){
		fprintf(stderr, "\nWrong number of command-line arguments in main\n");
		fprintf(stderr, "Usage: %s parameters.yaml sim_name\n\n", argv[0]);
		exit(1);
	}
	// Copy inputted parameter file name
	strcpy(param_file, argv[1]);
	// Parse through input parameter file and copy values to sim's internal parameter structure
	parse_parameters(param_file, &parameters);

	// Use our experimental curator, Wallace, to open appropriate files and initialize classes/etc. used in the simulation
	properties.wallace.InitializeSimulation(&parameters, &properties);
	properties.wallace.OpenFiles(argv[2]);

	// Run kinetic Monte Carlo loop n_steps times 
	for(int i_step = 0; i_step < parameters.n_steps; i_step++){
		// Wallace keeps track of outputting data, etc
		properties.wallace.UpdateTimestep(i_step);
		// Explicit KMC actions (binding, stepping, etc)
		properties.kinesin4.RunKMC();
		properties.prc1.RunKMC();
		// Diffusion
		properties.kinesin4.RunDiffusion();
		properties.prc1.RunDiffusion();
		// MTs go last because they sum up all the forces and stuff
//		properties.microtubules.RunDiffusion();
		// Some good ole-fashioned ASCII printout
		if(parameters.mt_printout == true)
			if(i_step % 1000 == 0)
				properties.wallace.PrintMicrotubules(0.5);
	}

	properties.wallace.OutputSimDuration();
	properties.wallace.CleanUp();
	return 0;
}
