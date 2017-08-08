/*Main simulation*/

#include "master_header.h"

int main(int argc, char *argv[]){

	long seed;
	int length_of_microtubule, n_microtubules, n_steps, data_threshold, range_of_data, pickup_time, n_datapoints, data_milestone, equil_milestone;
	double delta_t, duration, p_mutant, mutant_affinity, k_on, c_motor, k_off, motor_speed, motor_speed_eff, p_bind_n, p_unbind_n, p_move_nn, 
		   p_move_nm, p_bind_m, p_unbind_m, p_move_mm, p_move_mn;

	char param_file[160];
	system_parameters parameters;

	clock_t start, finish;
	FILE *output_file, *ID_file, *stream; 

	// Initialize and set generator type for RNG
	const gsl_rng_type *generator_type; 
	generator_type = gsl_rng_mt19937;   
	gsl_rng *rng; 

	// Get command-line input
	if(argc != 4){
		fprintf(stderr, "Wrong number of command-line arguments in main\n");
		fprintf(stderr, "Usage: %s parameters.yaml main_output.file ID_output.file\n", argv[0]);
		exit(1);
	}
	strcpy(param_file, argv[1]);
	output_file = gfopen(argv[2], "w");
	ID_file = gfopen(argv[3], "w");

	start = clock();

	// Parse through parameter file 
	parse_parameters(param_file, &parameters);

	// Initialize local variables 
	seed = parameters.seed;											// Seed for the RNG

	n_steps = parameters.n_steps;									// Total number of steps in one run
	data_threshold = parameters.data_threshold;						// Step at which data collection starts 
	n_datapoints = parameters.n_datapoints;							// Number of data points per MT to be written to file during sim
	delta_t = parameters.delta_t;									// Duration of one timestep in seconds	

	n_microtubules = parameters.n_microtubules;						// Number of microtubules to be used in simulation
	length_of_microtubule = parameters.length_of_microtubule;		// Length of each microtubule
	
	k_on = parameters.k_on;											// Binding rate in units of 1/(nanomolar*second)
	c_motor = parameters.c_motor;									// Motor concentration in nanomolar
	k_off = parameters.k_off;										// Unbinding rate in units of 1/second
	motor_speed = parameters.motor_speed;							// Motor speed in units of micrometers/second
	motor_speed_eff = motor_speed;									// Convert units to sites/second (each site is ~8nm long)
	p_mutant = parameters.p_mutant;									// Probability that any given tubulin dimer is a mutant
	mutant_affinity = parameters.mutant_affinity;					// Mutant tubulin binding affinity relative to normal tubulin

	p_bind_n = k_on*c_motor*delta_t;								// Probability to bind to normal tubulin per unit timestep
	p_unbind_n = k_off*delta_t;										// Probability to unbind per unit timestep for motors on normal tubulin
	p_move_nn = motor_speed_eff*delta_t;							// Probability to step one site per unit timestep for normal->normal tubulin moves
	p_move_nm = p_move_nn*mutant_affinity;							// Probability to step one site per unit timestep for normal->mutant tubulin moves

	p_bind_m = p_bind_n*mutant_affinity;							// Probability to bind to mutant tubulin per unit timestep
	p_unbind_m = p_unbind_n/mutant_affinity;						// Probability to unbind per unit timestep for motors on mutant tubulin
	p_move_mn = p_move_nn*(mutant_affinity/(mutant_affinity + 0.05));	// Some guess
	p_move_mm = p_move_mn*mutant_affinity;							// Probability to step one site per unit timestep for mutant->mutant tubulin moves

	range_of_data = n_steps - data_threshold;
	pickup_time = range_of_data/n_datapoints;				
	equil_milestone = data_threshold/10;							// Runs are long so we'd like some update on progress during equilibration
	data_milestone = range_of_data/10;								// Ditto for data collection

	// Allocate memory for and seed the RNG	
	rng = gsl_rng_alloc(generator_type);
	gsl_rng_set(rng, seed);

	// Output parameters in "sim units"  
	printf("Simulation elapses %i seconds overall; each timestep is %g seconds.\n", (int)(delta_t*n_steps), delta_t);
	printf("Microtubules have %i tubulin sites each; ~%i percent of sites are mutant.\n", length_of_microtubule, (int)(100*p_mutant));
	printf("Mutant tubulin is assumed to have %g the binding affinity of normal tubulin.\n\n", mutant_affinity);
	printf("Motor binding frequency: %f per timestep for normal tubulin\n", p_bind_n);
	printf("			 %f per timestep for mutant tubulin\n\n", p_bind_m);
	printf("Motor unbinding frequency: %f per timestep on normal tubulin\n", p_unbind_n);
	printf("			   %f per timestep on mutant tubulin\n\n", p_unbind_m);
	printf("Motor velocity: %f sites per timestep for normal->normal steps\n", p_move_nn);
	printf("		%f sites per timestep for normal->mutant steps\n", p_move_nm);
	printf("		%f sites per timestep for mutant->normal steps\n", p_move_mn);
	printf("		%f sites per timestep for mutant->mutant steps\n\n", p_move_mm);
	fflush(stdout);

	// Array of microtubules that acts as the experimental stage
	microtubule mt_array[n_microtubules];	
	// Initialize individual parameters based on polarity
	for(int i_mt = 0; i_mt < n_microtubules; i_mt++){
		if(i_mt % 2 == 0){
			mt_array[i_mt].polarity = 0;
			mt_array[i_mt].plus_end = length_of_microtubule - 1;
			mt_array[i_mt].minus_end = 0;
			mt_array[i_mt].delta_x = 1;
		}
		else if(i_mt % 2 == 1){
			mt_array[i_mt].polarity = 1;
			mt_array[i_mt].plus_end = 0;
			mt_array[i_mt].minus_end = length_of_microtubule - 1;
			mt_array[i_mt].delta_x = -1;
		}
		mt_array[i_mt].length = length_of_microtubule;
		// Populate lattice with tubulin sites; also store tubulin
		// site type in array that will later be used to write to file
		int tubulin_distribution[length_of_microtubule];
		int *array_ptr = tubulin_distribution;
		for(int i_site = 0; i_site < length_of_microtubule; i_site++){
			tubulin new_site;
			new_site.mt_index = i_mt;
			new_site.site_coord = i_site;
			// Ensure boundary sites are normal tubulin 
			if(i_site == 0 || i_site == (length_of_microtubule - 1)){
				new_site.mutant = false;
				tubulin_distribution[i_site] = 0;			// "0" corresponds to normal tubulin
			}
			// In bulk, insert mutant tubulin sites with probability p_mutant
			else{
				double ran = gsl_rng_uniform(rng);
				if(ran < p_mutant){
					new_site.mutant = true;
					tubulin_distribution[i_site] = 1;		// "1" corresponds to mutant tubulin
				}
				else{
					new_site.mutant = false;
					tubulin_distribution[i_site] = 0;		// "0" corresponds to normal tubulin
				}
			}
			mt_array[i_mt].lattice.push_back(new_site);	
		}
		// Generate a tubulin_distribution for each MT in the simulation
		char file_name[160];
		sprintf(file_name, "tubulin_distribution_%i.file", i_mt);
		stream = fopen(file_name, "w");
		fwrite(array_ptr, sizeof(int), length_of_microtubule, stream);
		fclose(stream);
	}
	// List of motors that acts as a resevoir throughout simulation
	std::vector<motor> motor_list;
	// Populate motor_list with new motors, each having its own unique ID
	for(int i_motor = 0; i_motor < n_microtubules*length_of_microtubule; i_motor++){
		motor new_motor;
		new_motor.ID = i_motor; 
		motor_list.push_back(new_motor);
	}
	// List that stores pointers to unbound tubulin sites
	std::vector<tubulin*> normal_unbound_list;
	std::vector<tubulin*> mutant_unbound_list;
	// Populate unbound_list with pointers to all tubulin sites (excluding boundaries) to set initial conditions 
	for(int i_mt = 0; i_mt < n_microtubules; i_mt++){
		for(int i_site = 1; i_site < length_of_microtubule - 1; i_site++){
			tubulin *site_address = &mt_array[i_mt].lattice[i_site];	
			if(mt_array[i_mt].lattice[i_site].mutant == true){
				mutant_unbound_list.push_back(site_address);
			}
			else if(mt_array[i_mt].lattice[i_site].mutant == false){
				normal_unbound_list.push_back(site_address);
			}	
			else{
				printf("what.\n");
			}
		}
	}
	// List that stores pointers to bound motors
	std::vector<motor*> normal_bound_list;
	std::vector<motor*> mutant_bound_list;

	// Begin main simulation loop
	for(int i_step = 0; i_step < n_steps; i_step++){

		int n_unoccupied, n_to_bind, n_bound, n_to_unbind, n_movable_n, n_to_move_n, n_movable_m, n_to_move_m, n_events,
			m_unoccupied, m_to_bind, m_bound, m_to_unbind, m_movable_m, m_to_move_m, m_movable_n, m_to_move_n, m_events; 
		double n_avg, m_avg;

		// Find total number of unoccupied sites on all microtubules (excludes boundary sites by design)
		n_unoccupied = normal_unbound_list.size();
		m_unoccupied = mutant_unbound_list.size();
		// Calculate average number of binding events that occur per timestep (over many timesteps)
		n_avg = p_bind_n*n_unoccupied;
		m_avg = p_bind_m*m_unoccupied;
		// Sample poisson distribution and determine the exact number of binding events to execute during this timestep
		n_to_bind = gsl_ran_poisson(rng, n_avg);
		m_to_bind = gsl_ran_poisson(rng, m_avg);

		// Find total number of bound motors on all microtubules
		n_bound = normal_bound_list.size();
		m_bound = mutant_bound_list.size();
		// Check boundary sites for motors and remove their contribution from n_bound if they exist
		for(int i_mt = 0; i_mt < n_microtubules; i_mt++){
			int plus_end = mt_array[i_mt].plus_end;
			if(mt_array[i_mt].lattice[plus_end].occupant != NULL){
				n_bound--;		// Boundary tubulin will always be normal by design
			}
			int minus_end = mt_array[i_mt].minus_end;
			if(mt_array[i_mt].lattice[minus_end].occupant != NULL){
				n_bound--;		// Boundary tubulin will always be normal by design
			}
		}
		// Sample binomial distribution and determine exact number of unbinding events to execute during this timestep
		n_to_unbind = gsl_ran_binomial(rng, p_unbind_n, n_bound);
		m_to_unbind = gsl_ran_binomial(rng, p_unbind_m, m_bound);

		// Iterate through normal_bound_list and determine how many motors are capable of stepping at the start of this timestep
		n_movable_n = 0;
		n_movable_m = 0;
		for(int normal_index = 0; normal_index < normal_bound_list.size(); normal_index++){
			// Pick motor from normal_bound_list; get its site coord and MT index
			int site_coord = normal_bound_list[normal_index]->site_coord;
			int mt_index = normal_bound_list[normal_index]->mt_index;
			// Get relevant parameters from microtubule 
			int delta_x = mt_array[mt_index].delta_x;
			int plus_end = mt_array[mt_index].plus_end;
			int minus_end = mt_array[mt_index].minus_end;
			// Make sure everything is fine and dandy in normal_bound_list
			if(mt_array[mt_index].lattice[site_coord].occupant == NULL){
				printf("Error with normal_bound_list in main loop!!\n");
				exit(1);
			}
			if(mt_array[mt_index].lattice[site_coord].mutant == true){
				printf("Error with normal_bound_list (type two) in main loop!!\n");
				exit(1);
			}
			// Exclude the plus end from statistics (end-pausing prevents stepping at plus_end)
			if(site_coord != plus_end){
				// Check if adjacent tubulin is unoccupied before proceeding
				if(mt_array[mt_index].lattice[site_coord + delta_x].occupant == NULL){
					// If adjacent tubulin is normal, count this site as n->n movable
					if(mt_array[mt_index].lattice[site_coord + delta_x].mutant == false){
						n_movable_n++;
					}
					// Else if adjacent tubulin is mutated, count this site as n->m movable
					else if(mt_array[mt_index].lattice[site_coord + delta_x].mutant == true){
						n_movable_m++;
					}
					else{
						printf("what. TWO!\n");
					}
				}
			}		
		}
		// Iterate through mutant_bound_list and determine how many motors are capable of stepping at the start of this timestep
		m_movable_n = 0;
		m_movable_m = 0;
		for(int index = 0; index < mutant_bound_list.size(); index++){
			// Pick motor from normal_bound_list; get its site coord and MT index
			int site_coord = mutant_bound_list[index]->site_coord;
			int mt_index = mutant_bound_list[index]->mt_index;
			// Get relevant parameters from microtubule 
			int delta_x = mt_array[mt_index].delta_x;
			int plus_end = mt_array[mt_index].plus_end;
			int minus_end = mt_array[mt_index].minus_end;
			// Make sure everything is fine and dandy in normal_bound_list
			if(mt_array[mt_index].lattice[site_coord].occupant == NULL){
				printf("Error with mutant_bound_list in main loop!!\n");
				exit(1);
			}
			if(mt_array[mt_index].lattice[site_coord].mutant == false){
				printf("Error with mutant_bound_list (type two) in main loop!!\n");
				exit(1);
			}
			// Exclude the plus end from statistics (end-pausing prevents stepping at plus_end)
			if(site_coord != plus_end){
				// Check if adjacent tubulin is unoccupied before proceeding
				if(mt_array[mt_index].lattice[site_coord + delta_x].occupant == NULL){
					// If adjacent tubulin is normal, count this site as n->n movable
					if(mt_array[mt_index].lattice[site_coord + delta_x].mutant == false){
						m_movable_n++;
					}
					// Else if adjacent tubulin is mutated, count this site as n->m movable
					else if(mt_array[mt_index].lattice[site_coord + delta_x].mutant == true){
						m_movable_m++;
					}
					else{
						printf("what. THREEEE!\n");
					}
				}
			}		
		}
		// Sample binomial distribution and determine exact number of move events to execute in this timestep
		n_to_move_n = gsl_ran_binomial(rng, p_move_nn, n_movable_n);
		n_to_move_m = gsl_ran_binomial(rng, p_move_nm, n_movable_m);
		m_to_move_n = gsl_ran_binomial(rng, p_move_mn, m_movable_n);
		m_to_move_m = gsl_ran_binomial(rng, p_move_mm, m_movable_m);

		// Calculate total number of kMC events
		n_events = n_to_bind + n_to_unbind + n_to_move_m + n_to_move_n;
		m_events = m_to_bind + m_to_unbind + m_to_move_n + m_to_move_m;		


		// Store the number of expected kMC events in a 1-D array for shuffling; type of event is encoded in integer value
		int kmc_list[n_events + m_events];
		// Even kMC events correspond to normal tubulin
		for(int i_list = 0; i_list < n_to_bind; i_list++){
			kmc_list[i_list] = 0;		// "0" corresponds to the kMC event of binding to normal tubulin
		}
		for(int i_list = n_to_bind; i_list < (n_to_bind + n_to_unbind); i_list++){
			kmc_list[i_list] = 2;		// "2" corresponds to the kMC event of unbinding from normal tubulin
		}
		for(int i_list = (n_to_bind + n_to_unbind); i_list < (n_events - n_to_move_n); i_list++){
			kmc_list[i_list] = 4;		// "4" corresponds to the kMC event of stepping from normal tubulin to mutant tubulin
		}
		for(int i_list = (n_events - n_to_move_n); i_list < n_events; i_list++){
			kmc_list[i_list] = 6;		// "6" corresponds to the kMC event of stepping from normal tubulin to normal tubulin
		}
		// Odd kMC events correspond to mutant tubulin
		for(int i_list = n_events; i_list < (n_events + m_to_bind); i_list++){
			kmc_list[i_list] = 1;		// "1" corresponds to the kMC event of binding to mutant tubulin
		}
		for(int i_list = (n_events + m_to_bind); i_list < (n_events + m_to_bind + m_to_unbind); i_list++){
			kmc_list[i_list] = 3;		// "3" corresponds to the kMC event of unbinding from mutant tubulin
		}
		for(int i_list = (n_events + m_to_bind + m_to_unbind); i_list < (n_events + m_events - m_to_move_m); i_list++){
			kmc_list[i_list] = 5;		// "5" corresponds to the kMC event of stepping from mutant tubulin to normal tubulin
		}
		for(int i_list = (n_events + m_events - m_to_move_m); i_list < (n_events + m_events); i_list++){
			kmc_list[i_list] = 7;		// "7" corresponds to the kMC event of stepping from mutant tubulin to mutant tubulin
		}

/*		printf("STATISTICS FOR STEP %i:\n", i_step);
		printf("  Normal: %i events; %i binds, %i unbinds, %i (out of %i) n->m steps, %i (out of %i) n->n steps || %lu bound, %lu unbound\n", 
				n_events, n_to_bind, n_to_unbind, n_to_move_m, n_movable_m, n_to_move_n, n_movable_n, normal_bound_list.size(), normal_unbound_list.size());
		printf("    Normal portion of kmc_list: ");
		for(int i_list = 0; i_list < n_events; i_list++){
			printf("%i", kmc_list[i_list]);
		}
		printf("\n");
		printf("  Mutant: %i events; %i binds, %i unbinds, %i (out of %i) m->n steps, %i (out of %i) m->m steps || %lu bound, %lu unbound\n", 
				m_events, m_to_bind, m_to_unbind, m_to_move_n, m_movable_n, m_to_move_m, m_movable_m, mutant_bound_list.size(), mutant_unbound_list.size());
		printf("    Mutant portion of kmc_list: ");
		for (int i_list = n_events; i_list < n_events + m_events; i_list++){
			printf("%i", kmc_list[i_list]);
		}
		printf("\n");
		printf("Initial state:\n");
		print_microtubules(&parameters, mt_array);
*/

		// If one or more kMC events occur, execute them accordingly
		if((n_events + m_events) > 0){
			// Shuffle kmc_list to randomize event order
			gsl_ran_shuffle(rng, kmc_list, (n_events + m_events), sizeof(int));
			// Iterate through kmc_list and perform appropriate algorithms  
			for(int list_entry = 0; list_entry < (n_events + m_events); list_entry++){

				int kmc_index = kmc_list[list_entry];

				switch(kmc_index){
					case 0:{
							motors_bind(&parameters, mt_array, motor_list, normal_bound_list, normal_unbound_list, rng);
/*							printf("After normal binding event:\n");
							print_microtubules(&parameters, mt_array);
*/	
							break;
					}
					case 1:{
							motors_bind(&parameters, mt_array, motor_list, mutant_bound_list, mutant_unbound_list, rng);
/*							printf("After mutant binding event:\n");
							print_microtubules(&parameters, mt_array);
*/	
							break;
					}
					case 2:{
							motors_unbind(&parameters, mt_array, motor_list, normal_bound_list, normal_unbound_list, rng);
/*							printf("After normal unbinding event:\n");
							print_microtubules(&parameters, mt_array);
*/	
							break;
					}
					case 3:{
							motors_unbind(&parameters, mt_array, motor_list, mutant_bound_list, mutant_unbound_list, rng);
/*							printf("After mutant unbinding event:\n");
							print_microtubules(&parameters, mt_array);
*/
							break;
					}
					case 4:{
							motors_move(&parameters, mt_array, motor_list, normal_bound_list, normal_unbound_list, 
										 mutant_bound_list, mutant_unbound_list, true, rng);
/*							printf("After n->m stepping event:\n");
							print_microtubules(&parameters, mt_array);
*/	
							break;
					}
					case 5:{
							motors_move(&parameters, mt_array, motor_list, mutant_bound_list, mutant_unbound_list, 
										 normal_bound_list, normal_unbound_list, false, rng);
/*							printf("After m->n stepping event:\n");
							print_microtubules(&parameters, mt_array);
*/							
							break;
					}
					case 6:{
							motors_move(&parameters, mt_array, motor_list, normal_bound_list, normal_unbound_list, false, rng);
/*							printf("After n->n stepping event:\n");
							print_microtubules(&parameters, mt_array);
*/
							break;
					}
					case 7:{
							motors_move(&parameters, mt_array, motor_list, mutant_bound_list, mutant_unbound_list, true, rng);
/*							printf("After m->m stepping event:\n");
							print_microtubules(&parameters, mt_array);
*/
							break;
					}
				}
				motors_boundaries(&parameters, mt_array, motor_list, normal_bound_list, rng, (n_events + m_events));
			}
		}
		// If no kMC events occur, simply update boundary conditions
		else if((n_events + m_events) == 0){
			motors_boundaries(&parameters, mt_array, motor_list, normal_bound_list, rng, 1);
		}

		// Give updates on equilibration process (every 10%)
		if(i_step < data_threshold){
			if(i_step%equil_milestone == 0){
				printf("Equilibration is %i percent complete (step number %i)\n", (int)(i_step/equil_milestone)*10, i_step);
				fflush(stdout);
			}
		}
		// Start data collection at appropriate step threshold
		else if(i_step >= data_threshold){
			int delta = i_step - data_threshold;
			// Collect data at integer multiples of pickup_time
			if(delta%pickup_time == 0){
				output_data(&parameters, mt_array, output_file, ID_file);	
			}
			// Give update on data collection process (every 10%)
			if(delta%data_milestone == 0){
				printf("Data collection is %i percent complete (step number %i)\n", (int)(delta/data_milestone)*10, i_step);
				fflush(stdout);
			}
		}
		else if(i_step == n_steps){
			printf("Done!\n");	
		}
	}
	fclose(output_file);
	fclose(ID_file);

	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC;

	stream = fopen("time_sim.dat","w");
	fprintf(stream, "Time to execute main code: %f seconds\n", duration);
	fclose(stream);

	return 0;
}
