#include "master_header.h"

void motors_bind(microtubule *mt_array, std::vector< std::vector<int> > &bound_list, std::vector< std::vector<int> > &unbound_list){

	int n_steps, pickup_time, n_protofilaments, length_of_microtubule, delta, n_bound, n_to_bind, n;
	double n_avg, k0;
	long seed;	

	// Decloration/initialization of pointers needed for GSL functions
	const gsl_rng_type *T;
	gsl_rng *r;	
	gsl_rng_env_setup();
	T = gsl_rng_default;

	// Temporary local assignment of parameters
	n_steps = 100;
	pickup_time = 1;
	n_protofilaments = 2;
	length_of_microtubule = 120;
	k0 = 0.25;

	for(int j = 0; j < n_protofilaments; j++){

		// This goes here so that the random numbers of each MT are not coupled
		r = gsl_rng_alloc (T);
		// Average number of attachments in one time step (measured over many time steps)
		n_avg = k0*(length_of_microtubule - mt_array[j].n_bound[0]);
		// Sample from Poisson distribution to estimate how many will bind in this particular time step
		// FIGURE OUT HOW TO SEED THIS RNG
		n_to_bind = gsl_ran_poisson(r, n_avg);

		std::cout << " n_bound: " << mt_array[j].n_bound[0];
		std::cout << ", n_avg: " << n_avg << ", n_to_bind: " << n_to_bind;
		std::cout << ", (un)bound list size: " << bound_list[j].size() << " (" << unbound_list[j].size() << ")" << std::endl;
/*
		std::sort(unbound_list[j].begin(), unbound_list[j].end());
		for(int k = 0; k < 10; k++){
			std::cout << unbound_list[j][k] << " ";
		}
		std::cout << std::endl;
*/
		// temp error
		if(unbound_list[j].empty() == true){
			std::cout << "WHOOPDIE DOO TWO" << std::endl;
			exit(1);
		}
		// Shuffles the list of unbound sites
		std::random_shuffle(unbound_list[j].begin(), unbound_list[j].end());

		// Caps number of motors to unbind at unbound_list.size()
		if(n_to_bind > unbound_list[j].size()){
			n = unbound_list[j].size();
		}
		else{
			n = n_to_bind;
		}
		for(int k = 0; k < n; k++){

			int	index = unbound_list[j].back();
			// temp error
			if(index > (length_of_microtubule - 1)){
				std::cout << "WHOOPDIE DOO" << std::endl;
				exit(1);
			}
			// Checks to make sure site is indeed unoccupied
			if(mt_array[j].track[index].occupancy[0] != 0){
				std::cout << "Error in motor binding code: " << j << "_" << index <<  std::endl;
				exit(1); 
			}
			mt_array[j].track[index].occupancy[0] = 2;
			mt_array[j].n_bound[0]++;
			bound_list[j].push_back(index);
			unbound_list[j].pop_back();
		}
	}
	return;
}
void motors_unbind(microtubule *mt_array, std::vector< std::vector<int> > &bound_list, std::vector< std::vector<int> > &unbound_list){

	int n_steps, pickup_time, n_protofilaments, length_of_microtubule, delta, n_bound, n_to_unbind, n;
	double k0, p;
	long seed;	

	// Decloration/initialization of pointers needed for GSL functions
	const gsl_rng_type *T;
	gsl_rng *r;	
	gsl_rng_env_setup();
	T = gsl_rng_default;

	// Temporary local assignment of parameters
	n_steps = 100;
	pickup_time = 1;
	n_protofilaments = 2;
	length_of_microtubule = 120;
	k0 = 0.5;
	p = 0.169;

	for(int j = 0; j < n_protofilaments; j++){

		// This goes here so that the random numbers of each MT are not coupled
		r = gsl_rng_alloc (T);
		// Probability to unbind per unit time step
		n_bound = mt_array[j].n_bound[0];
		// Sample binomial distribution for expected number of unbinds out of total occupied sites
		// FIGURE OUT HOW TO SEED THIS RNG
		n_to_unbind = gsl_ran_binomial(r, p, n_bound);		

		// temp error
		if(bound_list[j].empty() == true){
			std::cout << "WHOOPDIE DOO TRES" << std::endl;
			exit(1);
		}
		// Shuffles the list of bound sites
		std::random_shuffle(bound_list[j].begin(), bound_list[j].end());

		// Caps number of motors to unbinde at bound_list.size()
		if(n_to_unbind > bound_list[j].size()){
			n = bound_list[j].size();
		}
		else {
			n = n_to_unbind;
		}

		for(int k = 0; k < n; k++){

			int	index = bound_list[j].back();
			// temp error
			if(index > (length_of_microtubule - 1)){
				std::cout << "WHOOPDIE DOO" << std::endl;
				exit(1);
			}
			// Checks to make sure the site actually has a motor on it
			if(mt_array[j].track[index].occupancy[0] != 2){
				std::cout << "Error in motor unbinding code: " << j << "_" << index << std::endl;
				exit(1); 
			}
			mt_array[j].track[index].occupancy[0] = 0;
			mt_array[j].n_bound[0]--;
			unbound_list[j].push_back(index);
			bound_list[j].pop_back();
		}
	}
	return;
}
void motors_switch(microtubule *mt_array, std::vector< std::vector<int> > &bound_list, std::vector< std::vector<int> > &unbound_list){

	int n_steps, pickup_time, n_protofilaments, length_of_microtubule, delta, n_bound, n_to_switch, n_switched, n_bound_eff,  n, t;
	double k0, switch_freq;
	long seed;	

	std::vector<int> switchable_list, just_switched, switchable_list_temp;

	// Decloration/initialization of pointers needed for GSL functions
	const gsl_rng_type *T;
	gsl_rng *r;	
	gsl_rng_env_setup();
	T = gsl_rng_default;

	// Temporary local assignment of parameters
	n_steps = 100;
	pickup_time = 1;
	n_protofilaments = 2;
	length_of_microtubule = 120;
	k0 = 0.5;
	switch_freq = 0.44;

	for(int j = 0; j < n_protofilaments; j++){

		// This goes here so that the random numbers of each MT are not coupled
		r = gsl_rng_alloc (T);

		n_bound = mt_array[j].n_bound[0];
		n_switched = just_switched.size();
		n_bound_eff = n_bound - n_switched;

		// Assigns the "relative" index (above or below, essentially)  of a MT's partner based on its polarity 	
		if(mt_array[j].polarity[0] == 0){
			t = j + 1;
		}	
		else if(mt_array[j].polarity[0] == 1){
			t = j - 1;
		}

		//Stores the coordinates of motors that switched prior and then clears the switchable list
		// Resets the switchable list and populates it based on neighboring MT	
		for(int k = 0; k < n_bound; k++){
			int index = bound_list[j][k];
			// temp error
			if(mt_array[j].track[index].occupancy[0] !=2 ){
				std::cout << "WHAT DID YA DOOO" << std::endl;
				exit(1);
			}
			if(std::find(bound_list[t].begin(), bound_list[t].end(), index) == bound_list[t].end()){
				switchable_list_temp.push_back(index);
			}
		}
		// Sorts the temp_switch list and just_switch lists (which
		// is necessary to use the set_difference function), then
		// uses it to subtract the elements of just_switched from 
		// switchable_list_temp to form the final switchable_list
		switchable_list.clear();
		std::sort(switchable_list_temp.begin(), switchable_list_temp.end());
		std::sort(just_switched.begin(), just_switched.end());
		std::set_difference(switchable_list_temp.begin(), switchable_list_temp.end(), just_switched.begin(), just_switched.end(), std::inserter(switchable_list, switchable_list.begin()));
		switchable_list_temp.clear();
		just_switched.clear();
/*
		// Output of each list for troubleshooting purposes
		for(int k = 0; k < switchable_list.size(); k++){
			std::cout << switchable_list[k] << " ";
		}
		std::cout << "(" << switchable_list.size() << ")" << std::endl;
*/
		// Sample binomial distribution for expected number of unbinds out of total occupied sites
		// FIGURE OUT HOW TO SEED THIS RNG
		n_to_switch = gsl_ran_binomial(r, switch_freq, n_bound_eff);
		// Shuffles the list (ALREADY SEEDED)
		std::random_shuffle(switchable_list.begin(), switchable_list.end());

		// Caps number of motors to switch at switchable_list.size()
		if (n_to_switch > switchable_list.size()){
			n = switchable_list.size();
		}
		else{
			n = n_to_switch;
		}
		for(int k = 0; k < n; k++){

			int index = switchable_list.back();

			if(mt_array[j].track[index].occupancy[0] != 2){
				std::cout << "Error in motor switching code." << std::endl;
				exit(1);
			}
			if(mt_array[t].track[index].occupancy[0] != 0){
				std::cout << "Error in motor switching code; type two." << std::endl;
			}
			mt_array[j].track[index].occupancy[0] = 0;
			mt_array[j].n_bound[0]--;

			mt_array[t].track[index].occupancy[0] = 2;
			mt_array[t].n_bound[0]++;

			// Adds the switched motor's coords to the appropriate lists
			just_switched.push_back(index);
			bound_list[t].push_back(index);
			unbound_list[j].push_back(index);

			// Removes the switched motor's coords from the respective bound/unbound list
			unbound_list[t].erase(std::remove(unbound_list[t].begin(), unbound_list[t].end(), index), unbound_list[t].end());
			bound_list[j].erase(std::remove(bound_list[j].begin(), bound_list[j].end(), index), bound_list[j].end());
			
			switchable_list.pop_back();
		}	
	}
	return;
}
void motors_move(microtubule *mt_array, std::vector< std::vector<int> > &bound_list, std::vector< std::vector<int> > &unbound_list){

	int n_steps, pickup_time, n_protofilaments, length_of_microtubule, delta, n_bound, n_to_move, n, t;
	double k0, p;
	long seed;

	std::vector<int> move_list;

	// Decloration/initialization of pointers needed for GSL functions
	const gsl_rng_type *T;
	gsl_rng *r;	
	gsl_rng_env_setup();
	T = gsl_rng_default;

	// Temporary local assignment of parameters
	n_steps = 100;
	pickup_time = 1;
	n_protofilaments = 2;
	length_of_microtubule = 120;
	k0 = 0.5;
	p = 0.5;

	for(int j = 0; j < n_protofilaments; j++){

		// This goes here so that the random numbers of each MT are not coupled
		r = gsl_rng_alloc (T);
		n_bound = mt_array[j].n_bound[0];

		if(bound_list[j].empty() == true){
			std::cout << "WHOOPDIE DOO FIEF" << std::endl;
			exit(1);
		}
		// Sample binomial distribution for expected number of unbinds out of total occupied sites
		// FIGURE OUT HOW TO SEED THIS RNG
		n_to_move = gsl_ran_binomial(r, p, n_bound);	

		// Shuffles the list of bound sites
		std::random_shuffle(bound_list[j].begin(), bound_list[j].end());

		// Caps number of motors to move at bound_list.size()
		if(n_to_move > bound_list[j].size()){
			n = bound_list.size();
		}
		else{
			n = n_to_move;
		}
		// Clears the move_list and then repopulates it for the microtubule we're looking at	
		move_list.clear();
		for(int k = 0; k < n; k++){
			int index = bound_list[j][k];
			move_list.push_back(index);
		}
		// Assigns the "direction" of motor movement based on the MT's polarity 	
		if(mt_array[j].polarity[0] == 0){
			t = 1;													// Plus end is on right; moves towards increasing coord
			std::sort(move_list.begin(), move_list.end());			// Want descending order so motors don't block eachother when stepping
		}	
		else if(mt_array[j].polarity[0] == 1){
			t = -1;													// Plus end is on left; moves towards decreasing coord
			std::sort(move_list.rbegin(), move_list.rend());		// Want ascending order for the same reason
		}
/*		for(int k = 0; k < move_list.size(); k++){
			std::cout << move_list[k] << " ";
		}
		std::cout << t << std::endl;
*/
		for(int k = 0; k < move_list.size(); k++){

			int index = move_list.back();

			// Temp cop-out of having to deal with boundary conditions -- FIX THIS
			if((index + t ) <= 1 || (index + t) >= (length_of_microtubule - 1)){
				move_list.pop_back();
				break;
			}

			// Only step if another motor isn't in the way
			if(mt_array[j].track[index + t].occupancy[0] == 0){

				mt_array[j].track[index].occupancy[0] = 0;
				mt_array[j].track[index + t].occupancy[0] = 2;

				// Removes (index) from bound list, replaces it with (index + t)
				bound_list[j].erase(std::remove(bound_list[j].begin(), bound_list[j].end(), index), bound_list[j].end());
				bound_list[j].push_back((index + t));

				// Removes (index + t) from unbound list, replaces it with (index)
				unbound_list[j].erase(std::remove(unbound_list[j].begin(), unbound_list[j].end(), (index + t)), unbound_list[j].end());
				unbound_list[j].push_back(index);
			}
			else if(mt_array[j].track[index + t].occupancy[0] == 2){
				move_list.pop_back();
			}
			else{
				std::cout << "Something V wrong in move_motors, bro" << std::endl;
				exit(1);
			}
		}
	}
	return;
}
