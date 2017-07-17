#include "master_header.h"

// Writes data to file
void output_data(system_parameters *parameters, microtubule *mt_array, FILE *output_file){

	int n_microtubules = parameters->n_microtubules;
	int length_of_microtubule = parameters->length_of_microtubule;

	for(int j_mt = 0; j_mt < n_microtubules; j_mt++){
		for(int k_site = 0; k_site < length_of_microtubule; k_site++){
			fwrite(mt_array[j_mt].track[k_site].occupancy, sizeof(int), 1, output_file);
		}
	}
	return;
}

// Rudimentary ASCII output function
void print_microtubules(system_parameters *parameters, microtubule *mt_array){

	int n_microtubules = parameters->n_microtubules;
	int length_of_microtubule = parameters->length_of_microtubule;

	for (int j_mt = 0; j_mt < n_microtubules; j_mt++){
		for (int k_site = 0; k_site < length_of_microtubule; k_site++){
			switch (mt_array[j_mt].track[k_site].occupancy[0]){
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
		printf(" %i", mt_array[j_mt].polarity[0]);
		printf("\n");
		fflush(stdout);
	}
	return;
}
