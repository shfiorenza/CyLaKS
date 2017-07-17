#include "master_header.h"

void allocate_memory(system_parameters *parameters, microtubule *mt_array){

	int n_microtubules = parameters->n_microtubules;
	int length_of_microtubule = parameters->length_of_microtubule;

	mt_array = (microtubule*) malloc(n_microtubules*sizeof(microtubule));
	for (int j_mt = 0; j_mt < n_microtubules; j_mt++){	
		// 'Track' of each MT, i.e. what the motors move along (1 site = 1 tubulin)
		mt_array[j_mt].track = (site*) malloc(length_of_microtubule*sizeof(site));
		for(int k_site = 0; k_site < length_of_microtubule; k_site++){
			// Initial occupancy of each site on track (0 means empty and 2 means motor as of now)
			mt_array[j_mt].track[k_site].occupancy = (int*) malloc(sizeof(int));
			// Initial coord (or ID maybe??) distribution of each site on track 			 ***POTENTIALLY UNNECESSARY***
			mt_array[j_mt].track[k_site].coord = (int*) malloc(sizeof(int));
		}
		// Polarity of each MT (0 means plus-end is on RIGHT, 1 means plus-end is on LEFT)
		mt_array[j_mt].polarity = (int*) malloc(sizeof(int));
		// Coordinates of each MT (relative to the left edge, i.e. the 0 index of each array)
		mt_array[j_mt].coord = (int*) malloc(sizeof(int));
		// Number of motors bound to each MT 
		mt_array[j_mt].n_bound = (int*) malloc(sizeof(int));
	}
return;
}
