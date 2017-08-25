#include "master_header.h"

int main(int argc, char *argv[]){

	long seed;
	int n_steps, data_threshold, n_datapoints, range_of_data, n_pickup, equil_milestone, data_milestone;
	double duration;
	char param_file[160];
	FILE *occupancy_file, *ID_file, *stream; 
	clock_t start, finish;

	system_parameters parameters;
	system_properties properties;

	start = clock();
	// Get command-line input
	if(argc != 4){
		fprintf(stderr, "Wrong number of command-line arguments in main\n");
		fprintf(stderr, "Usage: %s parameters.yaml main_output.file ID_output.file\n", argv[0]);
		exit(1);
	}
	strcpy(param_file, argv[1]);
	occupancy_file = gfopen(argv[2], "w");
	ID_file = gfopen(argv[3], "w");

	// Parse parameters and set local variables
	parse_parameters(param_file, &parameters);
	seed = parameters.seed;											// Seed for the RNG
	n_steps = parameters.n_steps;									// Total number of steps in one run
	data_threshold = parameters.data_threshold;						// Step at which data collection starts 
	n_datapoints = parameters.n_datapoints;							// Number of data points per MT to be written to file during sim

	range_of_data = n_steps - data_threshold;
	n_pickup = range_of_data/n_datapoints;		  					// Number of timesteps between data collection events 	
	equil_milestone = data_threshold/10;
	data_milestone = range_of_data/10;

	// Initialize wallace so he print out some details about the sim 
	properties.wallace.Initialize(&parameters, &properties);
	properties.wallace.OutputSimDetails();
	// Initialize "experimental stage" for simulation
	properties.gsl.Initialize(seed);
	properties.microtubules.Initialize(&parameters, &properties);
	properties.kinesin4.Initialize(&parameters, &properties); 

	// Main kinetic Monte Carlo (KMC) simulation
	for(int i_step = 0; i_step < n_steps; i_step++){

//		printf("VVVVV STEP %i VVVVV\n", i_step);
		properties.current_step_ = i_step;
		properties.kinesin4.RunKMC();
//		printf("%g_%g\n", properties.p_bind_cum_, properties.p_unbind_cum_);
//		printf("%i_%i\n", properties.n_binds_, properties.n_unbinds_);
//		properties.wallace.PrintMicrotubules();

		//TODO: let wallace do the below
		// Give updates on equilibrium process (every 10%)
		if(i_step < data_threshold && i_step%equil_milestone == 0){
			printf("Equilibration is %i percent complete (step number %i)\n", (int)(i_step/equil_milestone)*10, i_step);
		}
		// Start data collection at appropriate step threshold
		else if(i_step >= data_threshold){
			int delta = i_step - data_threshold;
			// Collect data every n_pickup timesteps
			if(delta%n_pickup == 0){ 
				properties.wallace.OutputData(occupancy_file, ID_file);
			}
			if(delta%data_milestone == 0){ 
				printf("Data collection is %i percent complete (step number %i)\n", (int)(delta/data_milestone)*10, i_step);
				fflush(stdout);
			}
			else if(delta == range_of_data - 1){
				printf("Done!\n");  
			}
		}
	}
	fclose(occupancy_file);
	fclose(ID_file);

	double p_bind_avg = properties.p_bind_cum_/properties.n_binds_/parameters.delta_t;
	double p_unbind_avg = properties.p_unbind_cum_/properties.n_unbinds_/parameters.delta_t;
	double ratio = properties.n_binds_/(double)properties.n_unbinds_;
	printf("Average binding probability: %g\n", p_bind_avg);
	printf("Average unbinding probability: %g\n", p_unbind_avg);
	printf("%g bind event(s) [%i] for every unbind event [%i].\n", ratio, properties.n_binds_, properties.n_unbinds_);

	// TODO: let wallace track sim duation
	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC;

	stream = fopen("time_sim.dat","w");
	fprintf(stream, "Time to execute main code: %f seconds\n", duration);
	fclose(stream);

	return 0;
}
