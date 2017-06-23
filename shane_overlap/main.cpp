/*Overlaping microtubules.*/

#include "master_header.h"

int main(int argc, char *argv[])
{
	system_parameters parameters;
	microtubule *mt_array, *mt_array_initial;
	environment *timestep; 

	char param_file[160];
	int length_of_microtubule, n_runs, n_steps, pickup_time, n_protofilaments;
	double random;
	long seed;
	double v_motor_g, f_turning, c_kon, koff, duration, c_eff, random_check;
	FILE *f_config, *output_file, *stream, *detail_file; 
	clock_t start, finish;

	// Initialize timespec pointer for nanosleep function used below (also why time.h is needed)
	struct timespec ts; 
	ts.tv_sec = 0;								//how many seconds to sleep
	ts.tv_nsec = 60000000;						//how many nanonseconds to sleep (added to above for total sleep time)

	// Temporary local initalization of variables 
	n_steps = 100;
	pickup_time = 1;
	//	seed = 2338238232322;
	seed = 2382323223142322;

	//seeds the random number generator used in <algorithm> for shuffling
	std::srand(seed);

	// Get command-line input.
	if (argc != 4) {
		fprintf(stderr, "Wrong number of command-line arguments in main\n");
		fprintf(stderr, "Usage: %s variables.param init.config output_file\n", argv[0]);
		exit(1);
	}

	start = clock();
	/*
	   strcpy(param_file, argv[1]);

	// Read in input parameters. 
	parse_parameters(param_file, &parameters);

	// Initialize local variables. 
	n_steps = parameters.n_steps;	
	pickup_time = parameters.pickup_time;								//time interval between data collection 
	seed = parameters.seed;

	// Unit conversion, 1um=1000nm, 1 tubulin length ~ 8 nm. The standard units are showed in parameters.h file. 
	parameters.v_motor_g = parameters.v_motor_g*1000.0/8.0/2000.0;		// change the unit from um/s to tubulin/0.005s
	parameters.f_turning = parameters.f_turning/60.0/2000.0;			// change the unit from 1/min to 1/0.005s
	c_eff = parameters.c_motor;                                     	// effective c
	parameters.koff = parameters.koff/60.0/2000.0;	 					// change the unit to 1/0.005s

	printf("Motor moving velocity : %f tubulin/0.005s\nMotor turning frequency : %f 1/0.005s\n", parameters.v_motor_g, parameters.f_turning);
	printf("Dissociation frequency : %f 1/0.005s\n", parameters.koff);
	printf("Overall frequency %f\n", parameters.v_motor_g+parameters.f_turning+parameters.koff);
	fflush(stdout);
	 */
	// Read initial state 
	f_config = gfopen(argv[2], "rb");
	fprintf(stdout, "\nReading from config file: %s\n", argv[2]);
	fread(&n_runs, sizeof(int), 1, f_config);
	fread(&n_protofilaments, sizeof(int), 1, f_config);
	fread(&length_of_microtubule, sizeof(int), 1, f_config);
	fprintf(stdout, "   n_runs = %d\n", n_runs);
	fprintf(stdout, "   n_protofilaments = %d\n", n_protofilaments);
	fprintf(stdout, "   length = %d\n\n", length_of_microtubule);

	// Keeps length_of_microtubule blow 2000...relic from Hui-Shun's code; may remove later 
	if ((length_of_microtubule >= N_SITES_MAX)){ 
		fprintf(stderr, "The length is too long to do the simulation, change N_SITES_MAX in code.\n");
		exit(1);
	}

	// By allocating dynamic memory to pointers, we can effectively use them as arrays	

	// Allocate dynamic memory for the MT array at all timesteps 
	timestep = (environment*) malloc((n_steps/pickup_time)*sizeof(environment));
	for(int i = 0; i < (n_steps/pickup_time); i++){
		timestep[i].mt_array = (microtubule*) malloc(n_protofilaments*sizeof(microtubule));
		for(int j = 0; j < n_protofilaments; j++){
			timestep[i].mt_array[j].track = (site*) malloc(length_of_microtubule*sizeof(site));
			for(int k = 0; k < length_of_microtubule; k++){
				timestep[i].mt_array[j].track[k].occupancy = (int*) malloc(sizeof(int));			
				timestep[i].mt_array[j].track[k].coord = (int*) malloc(sizeof(int));			
			}
			timestep[i].mt_array[j].polarity = (int*) malloc(sizeof(int));
			timestep[i].mt_array[j].n_bound = (int*) malloc(sizeof(int));
		}
	}

	// Allocates dynamic memory for initial MT pointer, then populates it from build file
	mt_array = (microtubule*) malloc(n_protofilaments*sizeof(microtubule));
	mt_array_initial = (microtubule*) malloc(n_protofilaments*sizeof(microtubule));	
	for (int j = 0; j < n_protofilaments; j++){	
		mt_array[j].track = (site*) malloc(length_of_microtubule*sizeof(site));
		mt_array_initial[j].track = (site*) malloc(length_of_microtubule*sizeof(site));
		for(int k = 0; k < length_of_microtubule; k++){
			mt_array[j].track[k].occupancy = (int*) malloc(sizeof(int));
			mt_array_initial[j].track[k].occupancy = (int*) malloc(sizeof(int));
			fread(mt_array_initial[j].track[k].occupancy, sizeof(int), 1, f_config);

			mt_array[j].track[k].coord = (int*) malloc(sizeof(int));
			mt_array_initial[j].track[k].coord = (int*) malloc(sizeof(int));
			fread(mt_array_initial[j].track[k].coord, sizeof(int), 1, f_config);
		}
		mt_array[j].polarity = (int*) malloc(sizeof(int));
		mt_array_initial[j].polarity = (int*) malloc(sizeof(int));
		fread(mt_array_initial[j].polarity, sizeof(int), 1, f_config);

		mt_array[j].coord = (int*) malloc(sizeof(int));
		mt_array_initial[j].coord = (int*) malloc(sizeof(int));
		fread(mt_array_initial[j].coord, sizeof(int), 1, f_config);

		mt_array[j].n_bound = (int*) malloc(sizeof(int));
		mt_array_initial[j].n_bound = (int*) malloc(sizeof(int));
		fread(mt_array_initial[j].n_bound, sizeof(int), 1, f_config);
	}
	fclose(f_config);


	detail_file = gfopen(argv[3], "w");
	for (int run = 0; run < 1; run++){

		// Use vectors rather than linked lists for their built-in functions and multi-dimensionality
		// Vectors will store the ABSOLUTE COORD (MT array index + coord of MT edge) of relevant items 
		std::vector< std::vector<int> > bound_list;
		std::vector< std::vector<int> > unbound_list;

		// Populates vectors from initial conditions 
		for (int j = 0; j < n_protofilaments; j++){

			// Resets the MT array to what was read from the config file 
			// MAKE CLEAR_ARRAY FUNCTION TO ENABLE MULTIPLE RUNS
			memcpy(&mt_array[j], &mt_array_initial[j], sizeof(microtubule));

			std::vector<int> temp_bound, temp_unbound;
			for(int k = 0; k < length_of_microtubule; k++){
				if(mt_array[j].track[k].occupancy[0] == 0){
					temp_unbound.push_back(k + mt_array[j].coord[0]);
				}
				else if(mt_array[j].track[k].occupancy[0] == 2){
					temp_bound.push_back(k + mt_array[j].coord[0]);
				}
			}
			unbound_list.push_back(temp_unbound);
			bound_list.push_back(temp_bound);

			temp_unbound.clear();
			temp_bound.clear();
		}

		// kMC loop runs for n_steps and stores mt_array into timestep[i] for integer multiples of pickup_time
		for (int i_step = 0; i_step < n_steps; i_step++){
			if (i_step%pickup_time == 0){
				printf("\nRun %d, step %d:\n", (run+1), (i_step+1));
				print_microtubules(mt_array);
				for (int j = 0; j < n_protofilaments; j++){
					timestep[(int)(i_step/pickup_time)].mt_array[j] = mt_array[j];
				}
				nanosleep(&ts, NULL);
			}

			motors_bind(mt_array, bound_list, unbound_list);	
			printf("motors_bind:\n");
			print_microtubules(mt_array);
			printf("\n");

			motors_unbind(mt_array, bound_list, unbound_list);
			printf("motors_unbind:\n");
			print_microtubules(mt_array);
			printf("\n");

			motors_move(mt_array, bound_list, unbound_list);
			printf("motors_move:\n");
			print_microtubules(mt_array);
			printf("\n");

			motors_switch(mt_array, bound_list, unbound_list);
			printf("motors_switch: \n");
			print_microtubules(mt_array);
			printf("\n");
			fflush(stdout);

		}
		printf("\n");
	}
	//output(microtubule_final, (int)(n_steps/pickup_time), n_protofilaments, output_file, detail_file);
	fclose(detail_file);
	fflush(stdout);

	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC;

	stream = fopen("time_main.dat","w");
	fprintf(stream, "Time to execute main code: %f seconds\n", duration);
	fclose(stream);

	return 0;
}

void print_microtubules(microtubule *mt_array){

    int n_protofilaments = 2;
    int length_of_microtubule = 120;

			for (int j = 0; j < n_protofilaments; j++){
				for (int k = 0; k < length_of_microtubule; k++){
					switch (mt_array[j].track[k].occupancy[0]){
						case 0:
							printf("=");
							break;
						case 1:
							printf("x");
							break;
						case 2: 
							printf("M");
							break;
						case 8: 
							printf("?");
							break;
					}
				}
				printf(" %i", mt_array[j].polarity[0]);
				printf("\n");
				fflush(stdout);
			}

return;

}
