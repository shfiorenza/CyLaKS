#include "master_header.h"

void motors_bind(system_parameters *parameters, microtubule *mt_array, std::vector< std::vector<int> > &bound_list, std::vector< std::vector<int> > &unbound_list, int mt_index, gsl_rng *rng){

	int length_of_microtubule = parameters->length_of_microtubule;

	int list_size = unbound_list[mt_index].size();
	if(list_size > 0){

		int list_entry = gsl_rng_uniform_int(rng, list_size);
		int index = unbound_list[mt_index][list_entry];
		
		int n_attempts = 0;
		bool failure = false;
		// Ensures (within reason) that we don't touch boundaries
		while(index == 0 || index == (length_of_microtubule -1)){
			if(n_attempts > 2*list_size){
				failure = true;
				break;
			}
			list_entry = gsl_rng_uniform_int(rng, list_size);
			index = unbound_list[mt_index][list_entry];
			n_attempts++;
		}

		// Only attempt to bind if an eligible index was found
		if(failure != true){
			// Checks to make sure site is actually unoccupied 
			if(mt_array[mt_index].track[index].occupancy[0] != 0){
				printf("Error in motor binding code: %i_%i.\n", mt_index, index);
				exit(1); 
			}
			// Update lists and microtubule details
			mt_array[mt_index].track[index].occupancy[0] = 2;
			mt_array[mt_index].n_bound[0]++;
			// Need to randomize where we insert (index) in bound_list to avoid correlations
			int other_list_size = bound_list[mt_index].size();
			int pushback_offset = 0;
			if(other_list_size > 0){	
				pushback_offset = gsl_rng_uniform_int(rng, other_list_size);
			}
			bound_list[mt_index].insert(bound_list[mt_index].begin() + pushback_offset, index);
			unbound_list[mt_index].erase(std::remove(unbound_list[mt_index].begin(), unbound_list[mt_index].end(), index), unbound_list[mt_index].end());
		}
	}
}

void motors_unbind(system_parameters *parameters, microtubule *mt_array, std::vector< std::vector<int> > &bound_list, std::vector< std::vector<int> > &unbound_list, int mt_index, gsl_rng *rng){

	int length_of_microtubule = parameters->length_of_microtubule;

	int list_size = bound_list[mt_index].size();
	if(list_size > 0){

		int list_entry = gsl_rng_uniform_int(rng, list_size);
		int index = bound_list[mt_index][list_entry];

		int n_attempts = 0;
		bool failure = false;
		// Ensures (within reason) that we don't touch boundaries
		while(index == 0 || index == (length_of_microtubule - 1)){
			if(n_attempts > 2*list_size){
				failure = true;
				break;
			}
			list_entry = gsl_rng_uniform_int(rng, list_size);
			index = bound_list[mt_index][list_entry];
			n_attempts++;
		}

		// Only attempt to bind if an eligible index was found
		if(failure != true){
			// Checks to make sure the site actually has a motor on it
			if(mt_array[mt_index].track[index].occupancy[0] != 2){
				printf("Error in motor unbinding code: %i_%i.\n", mt_index, index);
				exit(1); 
			}
			// Update lists and microtubule details
			mt_array[mt_index].track[index].occupancy[0] = 0;
			mt_array[mt_index].n_bound[0]--;
			// Need to randomize where we insert (index) in unbound_list to avoid correlations
			int other_list_size = unbound_list[mt_index].size();
			int pushback_offset = 0;
			if(other_list_size > 0){
				pushback_offset = gsl_rng_uniform_int(rng, other_list_size);
			}
			unbound_list[mt_index].insert(unbound_list[mt_index].begin() + pushback_offset, index);
			bound_list[mt_index].erase(std::remove(bound_list[mt_index].begin(), bound_list[mt_index].end(), index), bound_list[mt_index].end());
		}
	}
}

void motors_switch(system_parameters *parameters, microtubule *mt_array, std::vector< std::vector<int> > &bound_list, std::vector< std::vector<int> > &unbound_list, int mt_index, gsl_rng *rng){

	int n_microtubules, length_of_microtubule, mt_index_adj;
	n_microtubules = parameters->n_microtubules;
	length_of_microtubule = parameters->length_of_microtubule;
	// Assigns the "relative" index (above or below, essentially)  of a MT's partner based on its polarity 	
	if(mt_array[mt_index].polarity[0] == 0){
		mt_index_adj = mt_index + 1;
	}	
	else if(mt_array[mt_index].polarity[0] == 1){
		mt_index_adj = mt_index - 1;
	}

	int list_size = bound_list[mt_index].size();
	if(list_size > 0){

		int list_entry = gsl_rng_uniform_int(rng, list_size);
		int index = bound_list[mt_index][list_entry];

		int n_attempts = 0;
		bool failure = false;
		// Ensures (within reason) that we obtain the index of a switchable motor that isn't on the boundary
		while(mt_array[mt_index_adj].track[index].occupancy[0] != 0 ||index == 0 || index == (length_of_microtubule - 1)){
			if(n_attempts > 2*list_size){
				failure = true;
				break;
			}
			list_entry = gsl_rng_uniform_int(rng, list_size);
			index = bound_list[mt_index][list_entry];
			n_attempts++;
		}

		// Only attempt to switch if an eligible index was found
		if(failure != true){
			// Verifies that there is a motor capable of switching
			if(mt_array[mt_index].track[index].occupancy[0] != 2){	
				printf("Error in motor switching code: %i_%i.\n", mt_index, index);
				exit(1);
			}
			// Verifies that the adjacent site on the neighboring MT is empty
			if(mt_array[mt_index_adj].track[index].occupancy[0] != 0){
				printf("Error in motor switching code; type two: %i_%i.\n", mt_index, index);
				exit(1);
			}

			// Update lists and microtubule details
			mt_array[mt_index].track[index].occupancy[0] = 0;
			mt_array[mt_index].n_bound[0]--;
			mt_array[mt_index_adj].track[index].occupancy[0] = 2;
			mt_array[mt_index_adj].n_bound[0]++;
			// Randomize where we insert (index) in the bound_list of the adjacent MT
			int other_list_size1 = bound_list[mt_index_adj].size();
			int pushback_offset1 = 0;
			if(other_list_size1 > 0){
				pushback_offset1 = gsl_rng_uniform_int(rng, other_list_size1);
			}
			bound_list[mt_index_adj].insert(bound_list[mt_index_adj].begin() + pushback_offset1, index);
			// Erase (index) from unbound list of the adjacent MT
			unbound_list[mt_index_adj].erase(std::remove(unbound_list[mt_index_adj].begin(), unbound_list[mt_index_adj].end(), 
														                             index), unbound_list[mt_index_adj].end());
			// Randomize where we insert (index) in the unbound_list of the 'active' MT
			int other_list_size2 = unbound_list[mt_index].size();
			int pushback_offset2 = 0;
			if(other_list_size2 > 0){
				pushback_offset2 = gsl_rng_uniform_int(rng, other_list_size2);
			}
			unbound_list[mt_index].insert(unbound_list[mt_index].begin() + pushback_offset2, index);
			// Erase (index) from bound_list of the 'active' MT
			bound_list[mt_index].erase(std::remove(bound_list[mt_index].begin(), bound_list[mt_index].end(), index), bound_list[mt_index].end());
		}
	}
}

void motors_move(system_parameters *parameters, microtubule *mt_array, std::vector< std::vector<int> > &bound_list, std::vector< std::vector<int> > &unbound_list, int mt_index, gsl_rng *rng){

	int length_of_microtubule, delta_x;
	length_of_microtubule = parameters->length_of_microtubule;
	// Assigns the "direction" of motor movement based on the MT's polarity 	
	if(mt_array[mt_index].polarity[0] == 0){
		delta_x = 1;								// Plus end is on right; moves towards increasing coord
	}	
	else if(mt_array[mt_index].polarity[0] == 1){
		delta_x = -1;								// Plus end is on left; moves towards decreasing coord
	}
	
	int list_size = bound_list[mt_index].size();
	if(list_size > 0){

		int list_entry;
		int index;
		int n_attempts = 0;
		bool failure = false;
		// Ensures (within reason) that we obtain the index of a motor that isn't blocked
		do{ if(n_attempts > 2*list_size){
				failure = true;
				break;
			}
			list_entry = gsl_rng_uniform_int(rng, list_size);
			index = bound_list[mt_index][list_entry];
			n_attempts++;
			// If motors try to step into the 'void,' i.e. out of the overlap, simply treat it as being blocked
			while(((index + delta_x) < 0) || ((index + delta_x) > (length_of_microtubule -1))){
				if(n_attempts > 2*list_size){
					failure = true;
					// Sets index to middle of MT to avoid seg_fault; doesn't accually matter since move will not execute
					index = (int)(length_of_microtubule/2);
					break;
				}
				list_entry = gsl_rng_uniform_int(rng, list_size);
				index = bound_list[mt_index][list_entry];
				n_attempts++;
			}
		}while(mt_array[mt_index].track[(index+delta_x)].occupancy[0] !=0);
		// Only attempt to move if an eligible index was found
		if(failure != true){
			// Verify that the motor isn't blocked before attempting to execute move
			if(mt_array[mt_index].track[(index+delta_x)].occupancy[0] != 0){
				printf("Error in motor moving code: %i_%i.\n", mt_index, index);
				exit(1);
			}
			mt_array[mt_index].track[index].occupancy[0] = 0;
			mt_array[mt_index].track[(index + delta_x)].occupancy[0] = 2;
			// Randomize where we place (index + delta_x) in bound_list
			int pushback_offset1 = gsl_rng_uniform_int(rng, list_size);
			bound_list[mt_index].insert(bound_list[mt_index].begin() + pushback_offset1, (index + delta_x));
			// Removes (index) from bound list
			bound_list[mt_index].erase(std::remove(bound_list[mt_index].begin(), bound_list[mt_index].end(), index), bound_list[mt_index].end());
			// Randomize where we place (index) in unbound_list
			int other_list_size = unbound_list[mt_index].size();
			int pushback_offset2 = 0;
			if(other_list_size > 0){
				pushback_offset2 = gsl_rng_uniform_int(rng, other_list_size);
			}
			unbound_list[mt_index].insert(unbound_list[mt_index].begin() + pushback_offset2, index);
			// Removes (index + delta_x) from unbound_list
			unbound_list[mt_index].erase(std::remove(unbound_list[mt_index].begin(), unbound_list[mt_index].end(), (index + delta_x)), unbound_list[mt_index].end());
		}
		else if(failure == true){
			printf("gnarly bro\n");
		}
	}
}

void motors_boundaries(system_parameters *parameters, microtubule *mt_array, std::vector< std::vector<int> > &bound_list, std::vector< std::vector<int> > &unbound_list, gsl_rng *rng){

	int n_microtubules, length_of_microtubule, plus_end, minus_end;
	double delta_t, alpha, beta, p_plus, p_minus;

	n_microtubules = parameters->n_microtubules;
	length_of_microtubule = parameters->length_of_microtubule;
	delta_t = parameters->delta_t;
	alpha = parameters->alpha;
	beta = parameters->beta;

	p_plus = (1 - beta);		// Probability of finding a motor at the plus end after one timestep 
	p_minus = alpha;			// Probability of finding a motor at the minus end per time step

	for(int j_mt = 0; j_mt < n_microtubules; j_mt++){
		// Sets appropriate polarity; see parameter.h file for description
		if(mt_array[j_mt].polarity[0] == 0){
			plus_end = (length_of_microtubule - 1);
			minus_end = 0;
		}
		else if(mt_array[j_mt].polarity[0] == 1){
			plus_end = 0;
			minus_end = (length_of_microtubule - 1);
		}
		// Generates two random numbers in the range [0.0, 1.0)
		double random1 = gsl_rng_uniform(rng);
		double random2 = gsl_rng_uniform(rng);		

		// Enforces boundary condition for plus-end
		if(random1 >= p_plus && mt_array[j_mt].track[plus_end].occupancy[0] == 2){
			mt_array[j_mt].track[plus_end].occupancy[0] = 0;	
			mt_array[j_mt].n_bound[0]--;
			// Removes (plus_end) from bound list
			bound_list[j_mt].erase(std::remove(bound_list[j_mt].begin(), bound_list[j_mt].end(), plus_end), bound_list[j_mt].end());
			// Randomize where we place (plus_end) in unbound_list
			int list_size = unbound_list[j_mt].size();
			int pushback_offset = 0;
			if(list_size > 0){
				pushback_offset = gsl_random_uniform_int(rng, list_size);
			}
			unbound_list[j_mt].insert(unbound_list[j_mt].begin() + pushback_offset, plus_end);
		}
		else if(random1 < p_plus && mt_array[j_mt].track[plus_end].occupancy[0] == 0){
			mt_array[j_mt].track[plus_end].occupancy[0] = 2;
			mt_array[j_mt].n_bound[0]++;			
			// Removes (plus_end) from unbound list
			unbound_list[j_mt].erase(std::remove(unbound_list[j_mt].begin(), unbound_list[j_mt].end(), plus_end), unbound_list[j_mt].end());
			// Randomize where we place (plus_end) in bound_list
			int list_size = bound_list[j_mt].size();
			int pushback_offset = 0;
			if(list_size > 0){
				pushback_offset = gsl_random_uniform_int(rng, list_size);
			}
			bound_list[j_mt].insert(bound_list[j_mt].begin() + pushback_offset, plus_end);
		}
		// Enforces boundary condition for minus-end
		if(random2 >= p_minus && mt_array[j_mt].track[minus_end].occupancy[0] == 2){
			mt_array[j_mt].track[minus_end].occupancy[0] = 0;	
			mt_array[j_mt].n_bound[0]--;
			// Removes (minus_end) from bound list
			bound_list[j_mt].erase(std::remove(bound_list[j_mt].begin(), bound_list[j_mt].end(), minus_end), bound_list[j_mt].end());
			// Randomize where we place (minus_end) in unbound_list
			int list_size = unbound_list[j_mt].size();
			int pushback_offset = 0;
			if(list_size > 0){
				pushback_offset = gsl_random_uniform_int(rng, list_size);
			}
			unbound_list[j_mt].insert(unbound_list[j_mt].begin() + pushback_offset, minus_end);
		}
		else if(random2 < p_minus && mt_array[j_mt].track[minus_end].occupancy[0] == 0){
			mt_array[j_mt].track[minus_end].occupancy[0] = 2;
			mt_array[j_mt].n_bound[0]++;			
			// Removes (minus_end) from unbound list; adds it to bound list
			unbound_list[j_mt].erase(std::remove(unbound_list[j_mt].begin(), unbound_list[j_mt].end(), minus_end), unbound_list[j_mt].end());
			// Randomize where we place (minus) end in bound_list
			int list_size = bound_list[j_mt].size();
			int pushback_offset = 0;
			if(list_size > 0){
				pushback_offset = gsl_random_uniform_int(rng, list_size);
			}
			bound_list[j_mt].push_back(minus_end);
		}
	}
}
