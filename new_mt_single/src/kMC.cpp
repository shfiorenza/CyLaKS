#include "master_header.h"

void motors_bind(system_parameters *parameters, microtubule *mt_array, std::vector<motor> &motor_list, std::vector<motor*> &bound_list, std::vector<int> &unbound_list, gsl_rng *rng){

	int length_of_microtubule = parameters->length_of_microtubule;
	int n_microtubules = parameters->n_microtubules;

	// Check to make sure at least one unbound site exists
	if(unbound_list.empty() != true){
	
		int unbound_entry, global_coord, site_coord, mt_index;
		int n_attempts = 0;
		bool failure = false;

		// Ensure (within reason) that boundary sites aren't included
		do{	if(n_attempts > unbound_list.size()){
				failure = true;
				break;
			}
			// Randomly pick an unbound list entry and get the corresponding global coordinate
			unbound_entry = gsl_rng_uniform_int(rng, unbound_list.size());
			global_coord = unbound_list[unbound_entry];
			// Calculate the site coordinate and mt index associated with this global coordinate
			site_coord = (global_coord % length_of_microtubule);
			mt_index = (int)(global_coord / length_of_microtubule);
			n_attempts++;
		}while(site_coord == 0 || site_coord == (length_of_microtubule -1));

		// Check to make sure an eligible site index was found
		if(failure != true){
			// Check to make sure the site is actually unoccupied 
			if(mt_array[mt_index].track[site_coord] != NULL){
				printf("Error in motor binding code at site %i on microtubule %i.\n", site_coord, mt_index);
				exit(1); 
			}
			
			// Find a motor to bind
			int motor_entry = 0;
			while(motor_list[motor_entry].bound == true){
				motor_entry++;
			}
			// Update motor details
			motor_list[motor_entry].bound = true;
			motor_list[motor_entry].mt_index = mt_index;
			motor_list[motor_entry].site_coord = site_coord;
			motor_list[motor_entry].global_coord = global_coord;
			motor_list[motor_entry].motor_entry = motor_entry;
			// Update microtubule details
			mt_array[mt_index].track[site_coord] = &motor_list[motor_entry];
			mt_array[mt_index].n_bound++;
			// Remove global coordinate of this site from unbound_list
			unbound_list.erase(unbound_list.begin() + unbound_entry);
			// Get memory location of this motor and add it to bound_list
			motor* motor_address = &motor_list[motor_entry];
			bound_list.push_back(motor_address);
		}
		// Announce if an eligible site index was NOT found
		else{
			printf("gnarly bro; failed to bind @ %i_%i!\n", mt_index, site_coord);
		}
	}
}

void motors_unbind(system_parameters *parameters, microtubule *mt_array, std::vector<motor> &motor_list, std::vector<motor*> &bound_list, std::vector<int> &unbound_list, gsl_rng *rng){

	int length_of_microtubule = parameters->length_of_microtubule;
	int n_microtubules = parameters->n_microtubules;

	// Check to make sure at least one bound motor exists
	if(bound_list.empty() != true){
	
		int bound_entry, site_coord, mt_index;
		int n_attempts = 0;
		bool failure = false;

		// Ensure (within reason) that boundary sites aren't included
		do{	if(n_attempts > bound_list.size()){
				failure = true;
				break;
			}
			// Randomly pick a bound list entry and get the site coordinate/mt index of the motor it points to
			bound_entry = gsl_rng_uniform_int(rng, bound_list.size());
			site_coord = bound_list[bound_entry]->site_coord;
			mt_index = bound_list[bound_entry]->mt_index;
			n_attempts++;

		}while(site_coord == 0 || site_coord == (length_of_microtubule -1));

		// Check to make sure an eligible bound motor was found
		if(failure != true){
			// Check to make sure the site actually has a motor on it
			if(mt_array[mt_index].track[site_coord] == NULL){
				printf("Error in motor unbinding code at site %i on microtubule %i.\n", site_coord, mt_index);
				exit(1); 
			}
 
			// Get motor list entry that points to this motor
			int motor_entry = bound_list[bound_entry]->motor_entry;			
			// Update motor details
			motor_list[motor_entry].bound = false;
			// Update microtubule details
			mt_array[mt_index].track[site_coord] = NULL;
			mt_array[mt_index].n_bound--;
			// Get global coordinate of the site this motor was attached to
			int global_coord = motor_list[motor_entry].global_coord;
			// Remove the entry that points to this motor from bound_list
			bound_list.erase(bound_list.begin() + bound_entry);
			// Add global coordinate of this site to unbound_list
			unbound_list.push_back(global_coord);
		}
		// Announce if an eligible site index was NOT found
		else{
			printf("gnarly bro; failed to unbind @ %i_%i!\n", mt_index, site_coord);
		}
	}
}

void motors_switch(system_parameters *parameters, microtubule *mt_array, std::vector<motor> &motor_list, std::vector<motor*> &bound_list, std::vector<int> &unbound_list, gsl_rng *rng){
	
	int n_microtubules = parameters->n_microtubules;
	int length_of_microtubule = parameters->length_of_microtubule;

	if(bound_list.empty() != true){

		// Randomly pick a list entry and get the corresponding bound motor coordinates
		int bound_entry, site_coord, mt_index, mt_index_adj;

		int n_attempts = 0;
		bool failure = false;
		// Ensures (within reason) that we obtain the index of a switchable motor that isn't on the boundary
		do{	if(n_attempts > 2*bound_list.size()){
				failure = true;
				break;
			}
			bound_entry = gsl_rng_uniform_int(rng, bound_list.size());
			site_coord = bound_list[bound_entry]->site_coord;
			mt_index = bound_list[bound_entry]->mt_index;
			mt_index_adj = mt_array[mt_index].mt_index_adj;
			n_attempts++;

		}while(mt_array[mt_index_adj].track[site_coord] != NULL || site_coord == 0 || (site_coord == length_of_microtubule - 1));
		
		// Only attempt to switch if a suitable motor was found
		if(failure != true){
			// Verifies there is a motor capable of switching
			if(mt_array[mt_index].track[site_coord] == NULL){
				printf("Error in motor switching code at site %i on microtubule %i.\n", site_coord, mt_index);
				exit(1);
			}
			// Verifies that neighbor site on adjacent MT is empty
			if(mt_array[mt_index_adj].track[site_coord] != NULL){
				printf("Error in motor switching code (type two) at site %i on microtubule %i.\n", site_coord, mt_index);
				exit(1);
			}

			// Find motor list entry that corresponds to this motor
			int motor_entry = bound_list[bound_entry]->motor_entry;
			// Get old global coordinate of motor
			int old_global_coord = motor_list[motor_entry].global_coord;
			// Calculate new global coordinate of motor
			int new_global_coord = mt_index_adj*length_of_microtubule + site_coord;
			// Update motor
			motor_list[motor_entry].mt_index = mt_index_adj;
			motor_list[motor_entry].global_coord = new_global_coord;
			// Update microtubules
			mt_array[mt_index].track[site_coord] = NULL;
			mt_array[mt_index].n_bound--;
			mt_array[mt_index_adj].track[site_coord] = &motor_list[motor_entry];
			mt_array[mt_index_adj].n_bound++;
			// Remove new global coordinate from unbound_list
			auto unbound_entry = std::find(unbound_list.begin(), unbound_list.end(), new_global_coord);
			int entry_index = std::distance(unbound_list.begin(), unbound_entry);
			unbound_list.erase(unbound_list.begin() + entry_index);
			// Add old global coordinate to unbound_list
			unbound_list.push_back(old_global_coord);
		}
		// Announce if an eligible site index was NOT found
		else{
			printf("gnarly bro; failed to switch @ %i_%i!\n", mt_index, site_coord);
		}
	}
}

void motors_move(system_parameters *parameters, microtubule *mt_array, std::vector<motor> &motor_list, std::vector<motor*> &bound_list, std::vector<int> &unbound_list, gsl_rng *rng){

	int n_microtubules = parameters->n_microtubules;
	int length_of_microtubule = parameters->length_of_microtubule;

	if(bound_list.empty() != true){

		int bound_entry, site_coord, mt_index, plus_end, minus_end, delta_x;
		
		int n_attempts = 0;
		bool failure = false;
		// Ensure (within reason) we find a motor that isn't blocked
		do{ if(n_attempts > 3*bound_list.size()){
				failure = true;
				break;
			}
			bound_entry = gsl_rng_uniform_int(rng, bound_list.size());
			site_coord = bound_list[bound_entry]->site_coord;
			mt_index = bound_list[bound_entry]->mt_index;
			plus_end = mt_array[mt_index].plus_end;
			minus_end = mt_array[mt_index].minus_end;
			delta_x = mt_array[mt_index].delta_x;
			n_attempts++;
			// Ensure (within reason) that motors don't step off the plus_end (i.e. enforce end-pausing)
			while(site_coord == plus_end){
				if(n_attempts > 3*bound_list.size()){
					failure = true;
					break;
				}
				bound_entry = gsl_rng_uniform_int(rng, bound_list.size());
				site_coord = bound_list[bound_entry]->site_coord;
				mt_index = bound_list[bound_entry]->mt_index;
				plus_end = mt_array[mt_index].plus_end;
				minus_end = mt_array[mt_index].minus_end;
				delta_x = mt_array[mt_index].delta_x;
				n_attempts++;
			}
		}while(mt_array[mt_index].track[site_coord + delta_x] != NULL);
		
		if(failure != true){
			// Verify motiry isn't blocked before attempting to step
			if(mt_array[mt_index].track[site_coord+delta_x] != NULL){
				printf("Error in motor moving code at site %i on microtubule %i.\n", site_coord, mt_index);
				exit(1);
			}

			// Find entry in motor list that corresponds to this motor
			int motor_entry = mt_array[mt_index].track[site_coord]->motor_entry;
			// Get old global coordinate of this motor
			int old_global_coord = motor_list[motor_entry].global_coord;
			// Calculate new global coordinate of this motor
			int new_global_coord = mt_index*length_of_microtubule + site_coord + delta_x;
			// Update motor details
			motor_list[motor_entry].global_coord = new_global_coord;
			motor_list[motor_entry].site_coord = site_coord + delta_x;
			// Update microtubule details
			mt_array[mt_index].track[site_coord] = NULL;
			mt_array[mt_index].track[site_coord + delta_x] = &motor_list[motor_entry];
			// Check if new site corresponds to the plus end, in which case its global coord is already absent from unbound_list 
			if(site_coord + delta_x != plus_end){
				// Find entry in unbound_list that corresponds to new global_coord, then remove it
				auto unbound_entry = std::find(unbound_list.begin(), unbound_list.end(), new_global_coord);
				int entry_index = std::distance(unbound_list.begin(), unbound_entry);
				unbound_list.erase(unbound_list.begin() + entry_index);
			}
			// Check if old site corresponds to the minus end, in which case its global coord should be suppressed from unbound_list
			if(site_coord != minus_end){
				// Add old global coordinate to unbound_list			
				unbound_list.push_back(old_global_coord);
			}

		}
		// Announce if an eligible site index was NOT found
		else if(failure == true){
			printf("gnarly bro; failed to move @ %i_%i\n", mt_index, site_coord);
		}
	}
}
void motors_boundaries(system_parameters *parameters, microtubule *mt_array, std::vector<motor> &motor_list, std::vector<motor*> &bound_list, std::vector<int> &unbound_list, gsl_rng *rng){

	int n_microtubules, length_of_microtubule, plus_end, minus_end, delta_x;
	double alpha, beta, p_move, p_switch;

	n_microtubules = parameters->n_microtubules;
	length_of_microtubule = parameters->length_of_microtubule;
	alpha = parameters->alpha;
	beta = parameters->beta;
	p_move = 125*parameters->motor_speed*parameters->delta_t;		
	p_switch = parameters->switch_rate*parameters->delta_t;
	// For BC to work, alpha and beta must be scaled so that motor speed is unity 
	double alpha_eff = alpha*p_move;
	double beta_eff = beta*p_move;		

	for(int i_mt = 0; i_mt < n_microtubules; i_mt++){
		// Get appropriate plus_ and minus_end site coordinates		
		plus_end = mt_array[i_mt].plus_end;
		minus_end = mt_array[i_mt].minus_end;
		delta_x = mt_array[i_mt].delta_x;

		// Generate two random numbers in the range [0.0, 1.0)
		double ran1 = gsl_rng_uniform(rng);
		double ran2 = gsl_rng_uniform(rng);		

		// If minus end is unoccupied, insert a motor onto the microtubule track with probability alpha_eff
		if(mt_array[i_mt].track[minus_end] == NULL && ran1 < alpha_eff){
			// Find motor to bind
			int motor_entry = 0;
			while(motor_list[motor_entry].bound == true){
				motor_entry++;
			}
			// Get global_coord of this site
			int global_coord = i_mt*length_of_microtubule + minus_end;
			// Update motor details
			motor_list[motor_entry].bound = true;
			motor_list[motor_entry].mt_index = i_mt;
			motor_list[motor_entry].site_coord = minus_end;
			motor_list[motor_entry].global_coord = global_coord;
			motor_list[motor_entry].motor_entry = motor_entry;
			// Update microtubule
			mt_array[i_mt].track[minus_end] = &motor_list[motor_entry];
			mt_array[i_mt].n_bound++;
			// Get memory location of this motor and add it to bound_list
			motor* motor_address = &motor_list[motor_entry];
			bound_list.push_back(motor_address);
		}
		// If a motor is at the plus end, remove it from the microtubule track with probability beta_eff 
		if(mt_array[i_mt].track[plus_end] != NULL && ran2 < beta_eff){
			// Find motor list entry that corresponds to this motor
			int motor_entry = mt_array[i_mt].track[plus_end]->motor_entry;
			// Update motor
			motor_list[motor_entry].bound = false; 
			// Update microtubule
			mt_array[i_mt].track[plus_end] = NULL;	
			mt_array[i_mt].n_bound--;
			// Find bound list entry that points to the memory location of this motor, then remove it
			motor* motor_address = &motor_list[motor_entry];
			auto bound_entry = std::find(bound_list.begin(), bound_list.end(), motor_address);		// Generate iterator pointing to entry
			int entry_index = std::distance(bound_list.begin(), bound_entry);						// Find distance from start of vector to iterator
			bound_list.erase(bound_list.begin() + entry_index);
		}
	}
}


