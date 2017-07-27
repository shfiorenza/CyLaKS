#include "master_header.h"

// Writes data to file
void output_data(system_parameters *parameters, microtubule *mt_array, FILE *output_file){

	int n_microtubules = parameters->n_microtubules;
	int length_of_microtubule = parameters->length_of_microtubule;

	for(int j_mt = 0; j_mt < n_microtubules; j_mt++){

		int temp_array[length_of_microtubule];
		int *array_ptr = temp_array;
		for(int i_site = 0; i_site < length_of_microtubule; i_site++){
			if(mt_array[j_mt].track[i_site] == NULL){
				temp_array[i_site] = 0;
			}
			else{
				temp_array[i_site] = 2;
			}
		}
		fwrite(array_ptr, sizeof(int), length_of_microtubule, output_file);
	}
	return;
}

// Rudimentary ASCII output function
void print_microtubules(system_parameters *parameters, microtubule *mt_array){

	int n_microtubules = parameters->n_microtubules;
	int length_of_microtubule = parameters->length_of_microtubule;

	for (int j_mt = 0; j_mt < n_microtubules; j_mt++){
		for (int k_site = 0; k_site < length_of_microtubule; k_site++){
			if(mt_array[j_mt].track[k_site] == NULL)
				printf("=");
			else
				printf("M");
		}
		printf(" %i", mt_array[j_mt].polarity);
		printf("\n");
		fflush(stdout);
	}
	return;
}
