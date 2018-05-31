/* Main simulation file */
#include "master_header.h"

int main(int argc, char *argv[]){

	char param_file[160], occupancy_file[160], 
		 motor_ID_file[160], xlink_ID_file[160], 
		 tether_coord_file[160], mt_coord_file[160], 
		 motor_extension_file[160], xlink_extension_file[160],
		 motor_force_file[160], xlink_force_file[160], total_force_file[160];

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
	sprintf(tether_coord_file, "%s_tether_coord.file", argv[2]);
	sprintf(mt_coord_file, "%s_mt_coord.file", argv[2]);
	sprintf(motor_extension_file, "%s_motor_extension.file", argv[2]);
	sprintf(xlink_extension_file, "%s_xlink_extension.file", argv[2]);
	sprintf(motor_force_file, "%s_motor_force.file", argv[2]);
	sprintf(xlink_force_file, "%s_xlink_force.file", argv[2]);
	sprintf(total_force_file, "%s_total_force.file", argv[2]);
	// Parse through input parameter file and copy values to sim's internal parameter structure
	parse_parameters(param_file, &parameters);
	// Open occupancy file, which stores the species ID of each occupant (or -1 for none) for all MT sites during data collection stage (DCS)
	properties.occupancy_file_ = gfopen(occupancy_file, "w");
	// Open motor ID file, which stores the unique ID of all bound motors (unbound are not tracked) and their respective site indices during DCS
	properties.motor_ID_file_ = gfopen(motor_ID_file, "w");
	// Open xlink ID file, which does the same as the motor ID file but for xlinks
	properties.xlink_ID_file_ = gfopen(xlink_ID_file, "w");
	// Open tether coord file, which stores the coordinates of the anchor points of tethered motors
	properties.tether_coord_file_ = gfopen(tether_coord_file, "w");
	// Open mt coord file, which stores the coordinates of the left-most edge of each microtubule during DCS
	properties.mt_coord_file_ = gfopen(mt_coord_file, "w");
	// Open motor extension file, which stores the number of motors with a certain tether extension for all possible extensions
	properties.motor_extension_file_ = gfopen(motor_extension_file, "w");
	// Open xlink extension file, which stores the number of stage-2 xlinks at a certain extension for all possible extensions
	properties.xlink_extension_file_ = gfopen(xlink_extension_file, "w");
	// Open motor force file, which stores the sum of forces coming from motor tether extensions
	properties.motor_force_file_ = gfopen(motor_force_file, "w");
	// Open xlink force file, which stores the sum of forces coming from xlink extensions
	properties.xlink_force_file_ = gfopen(xlink_force_file, "w");
	// Open total force file, which stores the sum of ALL forces coming from xlink and motor tether extensions
	properties.total_force_file_ = gfopen(total_force_file, "w");

	// Initialize the experimental curator, Wallace; he does things such as output data, print ASCII models of microtubules, etc.
	properties.wallace.Initialize(&parameters, &properties);
	// Initialize the general science library (gsl) class; just an easy way of sampling distributions and referencing the RNG
	properties.gsl.Initialize(parameters.seed);
	// Initialize microtubules, kinesin4, and prc1 to set the experimental stage 
	properties.microtubules.Initialize(&parameters, &properties);
	properties.kinesin4.Initialize(&parameters, &properties); 
	properties.prc1.Initialize(&parameters, &properties);

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
		properties.microtubules.RunDiffusion();  // XXX fix to let both MTs move
		// Some good ole-fashioned ASCII printout
		if(parameters.mt_printout == true)
			if(i_step % 1000 == 0)
				properties.wallace.PrintMicrotubules(0.5);
	}
	properties.wallace.OutputSimDuration();
	properties.wallace.CleanUp();

	return 0;
}
