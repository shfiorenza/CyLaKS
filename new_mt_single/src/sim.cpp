/*Overlaping microtubules.*/

#include "master_header.h"

int main(int argc, char *argv[]){

	long seed;
	int length_of_microtubule, n_microtubules, n_steps, data_threshold, range_of_data, pickup_time, equil_milestone, c_motor;
	double delta_t, duration, k_on, k_off, switch_rate, motor_speed, motor_speed_eff, p_bind, p_unbind, p_switch, p_move;

	char param_file[160];
	system_parameters parameters;

	clock_t start, finish;
	FILE *f_config, *output_file, *stream; 

	// Initialize and set generator type for RNG
	const gsl_rng_type *generator_type; 
	generator_type = gsl_rng_mt19937;   
	gsl_rng *rng; 

	// Get command-line input
	if(argc != 3){
		fprintf(stderr, "Wrong number of command-line arguments in main\n");
		fprintf(stderr, "Usage: %s parameters.yaml output_file.file\n", argv[0]);
		exit(1);
	}
	strcpy(param_file, argv[1]);
	output_file = gfopen(argv[2], "w");

	start = clock();

	// Parse through parameter file 
	parse_parameters(param_file, &parameters);

	// Initialize local variables 
	seed = parameters.seed;											// Seed for the RNG

	n_steps = parameters.n_steps;									// Total number of steps in one run
	data_threshold = parameters.data_threshold;						// Step at which data collection starts 
	delta_t = parameters.delta_t;									// Duration of one timestep in seconds

	n_microtubules = parameters.n_microtubules;						// Number of microtubules to be used in simulation
	length_of_microtubule = parameters.length_of_microtubule;		// Length of each microtubule
	
	k_on = parameters.k_on;											// Binding rate in units of 1/(nanomolar*second)
	c_motor = parameters.c_motor;									// Motor concentration in nanomolar
	k_off = parameters.k_off;										// Unbinding rate in units of 1/second
	switch_rate = parameters.switch_rate;							// Switch rate in units of 1/second
	motor_speed = parameters.motor_speed;							// Motor speed in units of micrometers/second
	motor_speed_eff = 125*motor_speed;								// Convert units to sites/second (each site is ~8nm long)

	p_bind = k_on*c_motor*delta_t;									// Probability to bind per unit timestep 
	p_unbind = k_off*delta_t;										// Probability to unbind per unit timestep
	p_switch = switch_rate*delta_t;									// Probability to switch per unit timestep
	p_move = motor_speed_eff*delta_t;								// Probability to move one site per unit timestep

	range_of_data = n_steps - data_threshold;
	pickup_time = range_of_data/100000;								// We want ~10,000 timesteps to average from
	equil_milestone = data_threshold/20;							// Runs are long so we'd like some update on progress during equilibration

	// Allocate memory for and seed the RNG	
	rng = gsl_rng_alloc(generator_type);
	gsl_rng_set(rng, seed);

	// Output parameters in "sim units"  
	printf("Motor velocity: %f sites per timestep\n", p_move);
	printf("Motor switching frequency: %f per timestep\n", p_switch);
	printf("Motor binding frequency: %f per timestep\n", p_bind);
	printf("Motor unbinding frequency: %f per timestep\n", p_unbind);
	printf("Total simulation duration: %f seconds\n", delta_t*n_steps);
	printf("Timestep duration: %f seconds\n\n", delta_t);
	fflush(stdout);

	// Array of microtubules that acts as experimental stage
	microtubule mt_array[n_microtubules];	
	// Populate the array and initialize microtubules
	for(int i_mt = 0; i_mt < n_microtubules; i_mt++){
		if(i_mt % 2 == 0){
			mt_array[i_mt].polarity = 0;
			mt_array[i_mt].plus_end = length_of_microtubule - 1;
			mt_array[i_mt].minus_end = 0;
			mt_array[i_mt].delta_x = 1;
			mt_array[i_mt].mt_index_adj = i_mt + 1;		// only allows switching between PAIRS; fix
		}
		else if(i_mt % 2 == 1){
			mt_array[i_mt].polarity = 1;
			mt_array[i_mt].plus_end = 0;
			mt_array[i_mt].minus_end = length_of_microtubule - 1;
			mt_array[i_mt].delta_x = -1;
			mt_array[i_mt].mt_index_adj = i_mt - 1;		// only allows switching between PAIRS; fix
		}
		mt_array[i_mt].length = length_of_microtubule;
		mt_array[i_mt].coord = 0;
		mt_array[i_mt].n_bound = 0;
		mt_array[i_mt].track.resize(length_of_microtubule, nullptr);
	}

	// List of motors that acts as a resevoir throughout simulation
	std::vector<motor> motor_list;
	// Populate motor_list with new motors, each having its own unique ID
	for(int i_motor = 0; i_motor < n_microtubules*length_of_microtubule; i_motor++){
		motor new_motor;
		new_motor.ID = i_motor; 
		new_motor.global_coord = i_motor;	// Initially 'ID' and 'global_coord' correspond, but global_coord will be shuffled around; ID is static
		motor_list.push_back(new_motor);
	}

	// List of pointers that stores the memory addresses of bound motors
	std::vector<motor*> bound_list;

	// List of integers that stores raw coordinates of unoccupied sites
	std::vector<int> unbound_list;
	// Populate unbound_list (exclude boundary sites)
	for(int i_mt = 0; i_mt < n_microtubules; i_mt++){
		for(int i_site = 1; i_site < length_of_microtubule - 1; i_site++){      // Sites in range [0, length_of_mt) correspond to first MT; sites in
			unbound_list.push_back((i_mt*length_of_microtubule) + i_site);   // range [length_of_mt, 2*length_of_mt) correspond to second MT, etc.
		}
	}

	// Main kMC loop 
	for(int i_step = 0; i_step < n_steps; i_step++){

		double n_avg;
		int n_unoccupied, n_to_bind, n_bound, n_to_unbind, n_switchable, n_to_switch, n_movable, n_to_move, n_events, i_mt_adj;

		// Find total number of unoccupied sites on all microtubules (excluding boundary sites)
		n_unoccupied = unbound_list.size();
		// Calculate average number of binding events that occur per timestep (over many timesteps)
		n_avg = p_bind*n_unoccupied;
		// Sample poisson distribution to 'predict' exact number of binding events in this particular timestep
		n_to_bind = gsl_ran_poisson(rng, n_avg);

		// Find total number of bound motors on all microtubules
		n_bound = bound_list.size();
		// Sample binomial distribution to 'predict' exact  number of unbinding events in this particular timestep
		n_to_unbind = gsl_ran_binomial(rng, p_unbind, n_bound);

		// Run through bound list to determine how many motors are capable of stepping and switching for this timestep
		n_switchable = 0;
		n_movable = 0;
		for(int entry = 0; entry < bound_list.size(); entry++){
			// Pick motor from our bound list; get its site coord and MT index
			int site_coord = bound_list[entry]->site_coord;
			int mt_index = bound_list[entry]->mt_index;

			// Get relevant parameters from microtubule 
			int delta_x = mt_array[mt_index].delta_x;
			int mt_index_adj = mt_array[mt_index].mt_index_adj;
			int plus_end = mt_array[mt_index].plus_end;
			int minus_end = mt_array[mt_index].minus_end;

			// Make sure everything is fine and dandy in our bound list
			if(mt_array[mt_index].track[site_coord] == NULL){
				printf("Error with bound_list in main loop!!\n");
				exit(1);
			}

			// Check if motor at minus end is capable of stepping; switching only occurs in bulk so don't check that
			if(site_coord == minus_end){
				if(mt_array[mt_index].track[minus_end + delta_x] == NULL){
					n_movable++;
				}
			}		
			// Otherwise, only look at sites in the bulk (not minus_end or plus_end)
			else if(site_coord != plus_end){
				// If site coordinate on adjacent MT isn't occupied, count motor as switchable
				if(mt_array[mt_index_adj].track[site_coord] == NULL){
					n_switchable++;
				}
				// If adjacent site coordinate on same MT isn't occupied, count motor as movable 
				if(mt_array[mt_index].track[site_coord + delta_x] == NULL){
					n_movable++;
				}
			}
		}
		// Sample binomial distribution to calculate expected number of switch and move events in this particular timestep
		n_to_switch = gsl_ran_binomial(rng, p_switch, n_switchable);
		n_to_move = gsl_ran_binomial(rng, p_move, n_movable);
		n_events = n_to_bind + n_to_unbind + n_to_switch + n_to_move;

		// kMC_list stores the number of expected kMC events in an array for shuffling; type of event is encoded in integer value
		int kmc_list[n_events];
		for(int i_list = 0; i_list < n_to_bind; i_list++){
			kmc_list[i_list] = 0;		// "0" means the kMC event corresponds to binding
		}
		for(int i_list = n_to_bind; i_list < (n_to_bind + n_to_unbind); i_list++){
			kmc_list[i_list] = 1;		// "1" means the kMC event corresponds to unbinding
		}
		for(int i_list = (n_to_bind + n_to_unbind); i_list < (n_events - n_to_move); i_list++){
			kmc_list[i_list] = 2;		// "2" means the kMC event corresponds to switching
		}
		for(int i_list = (n_events - n_to_move); i_list < n_events; i_list++){
			kmc_list[i_list] = 3;		// "3" means the kMC event corresponds to moving
		}
		// Only attempt to shuffle kmc_array[i_mt] if it exists
		if(n_events > 0){
			gsl_ran_shuffle(rng, kmc_list, n_events, sizeof(int));
		}

		// Iterates through kmc_list and executes the expected number of events 
		for(int list_entry = 0; list_entry < n_events; list_entry++){

			// Pick kmc_index from next entry in kmc_list; still random because kmc_list was shuffled
			int kmc_index = kmc_list[list_entry];

			switch(kmc_index){
				case 0:{
						   motors_bind(&parameters, mt_array, motor_list, bound_list, unbound_list, rng);
						   break;
				}
				case 1:{
						   motors_unbind(&parameters, mt_array, motor_list, bound_list, unbound_list, rng);
						   break;
				}
				case 2:{
						   motors_switch(&parameters, mt_array, motor_list, bound_list, unbound_list, rng);
						   break;
				}
				case 3:{
						   motors_move(&parameters, mt_array, motor_list, bound_list, unbound_list, rng);
						   break;
				}
			}
		}
		motors_boundaries(&parameters, mt_array, motor_list, bound_list, unbound_list, rng);
//		print_microtubules(&parameters, mt_array);
//		printf("%i__%i\n", n_movable, n_to_move);
//		printf("\n");

		// Gives updates on equilibration process (every 5%)
		if(i_step < data_threshold){
			if(i_step%equil_milestone == 0){
				printf("Equilibration is %i percent complete. (step number %i)\n", (int)(i_step/equil_milestone)*5, i_step);
			}
		}
		// Starts data collection at threshold
		else if(i_step >= data_threshold){
			int delta = i_step - data_threshold;
			if(delta%pickup_time == 0){
				printf("Step %d\n", (i_step));
//				print_microtubules(&parameters, mt_array);
				output_data(&parameters, mt_array, output_file);	
				fflush(stdout);
			}
		}
	}
	fclose(output_file);

	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC;

	stream = fopen("time_sim.dat","w");
	fprintf(stream, "Time to execute main code: %f seconds\n", duration);
	fclose(stream);

	return 0;
}
