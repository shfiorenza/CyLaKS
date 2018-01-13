/* Main simulation file */
#include "master_header.h"

int main(int argc, char *argv[]){

	char param_file[160], occupancy_file[160], motor_ID_file[160], xlink_ID_file[160], MT_coord_file[160];
	system_parameters parameters;
	system_properties properties;

	// Check that the user has input the correct number of arguments 
	if(argc != 3){
		fprintf(stderr, "\nWrong number of command-line arguments in main\n");
		fprintf(stderr, "Usage: %s parameters.yaml sim_name\n\n", argv[0]);
		exit(1);
	}
	strcpy(param_file, argv[1]);
	// Generate names of output files based on the input simulation name
	sprintf(occupancy_file, "%s_occupancy.file", argv[2]);
	sprintf(motor_ID_file, "%s_motorID.file", argv[2]);
	sprintf(xlink_ID_file, "%s_xlinkID.file", argv[2]);	
	sprintf(MT_coord_file, "%s_MTcoord.file", argv[2]);
	// Parse through input parameter file and copy values to sim's internal parameter structure
	parse_parameters(param_file, &parameters);
	// Open occupancy file, which stores the species ID of each occupant (or -1 for none) for all MT sites during data collection stage (DCS)
	properties.occupancy_file_ = gfopen(occupancy_file, "w");
	// Open motor ID file, which stores the unique ID of all bound motors (unbound are not tracked) and their respective site indices during DCS
	properties.motor_ID_file_ = gfopen(motor_ID_file, "w");
	// Open xlink ID file, which does the same as the motor ID file but for xlinks
	properties.xlink_ID_file_ = gfopen(xlink_ID_file, "w");
	// Open MT coord file, which stores the coordinates of the left-most edge of each microtubule during DCS
	properties.MT_coord_file_ = gfopen(MT_coord_file, "w");

	// Initialize the experimental curator, Wallace; he does things such as output data, print ASCII models of microtubules, etc.
	properties.wallace.Initialize(&parameters, &properties);
	// Initialize the general science library (gsl) class; just an easy way of sampling distributions and referencing the RNG
	properties.gsl.Initialize(parameters.seed);
	// Initialize microtubules, kinesin4, and prc1 to set the experimental stage 
	properties.microtubules.Initialize(&parameters, &properties);
	properties.kinesin4.Initialize(&parameters, &properties); 
	properties.prc1.Initialize(&parameters, &properties);

	// Temporary way of starting MTs with an offset (in sites)
	properties.microtubules.mt_list_[1].coord_ = 5;

	// Run kinetic Monte Carlo loop n_steps times 
	for(int i_step = 0; i_step < parameters.n_steps; i_step++){
		properties.wallace.UpdateTimestep(i_step);
//		properties.kinesin4.RunKMC();
		properties.prc1.RunKMC();
//		properties.microtubules.RunDiffusion();
//		properties.kinesin4.RunDiffusion();
		properties.prc1.RunDiffusion();
		properties.wallace.PrintMicrotubules(0.1);
	}
	properties.wallace.OutputSimDuration();
	properties.wallace.CleanUp();

	return 0;
}
