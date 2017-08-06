#include "master_header.h"

// Writes data to file
void output_data(system_parameters *parameters, microtubule *mt_array, FILE *output_file, FILE *ID_file){

	int n_microtubules = parameters->n_microtubules;
	int length_of_microtubule = parameters->length_of_microtubule;

	for(int i_mt = 0; i_mt < n_microtubules; i_mt++){
		// Translate occupancy data into a simple integer array
		int temp_array[length_of_microtubule];
		int *array_ptr = temp_array;
		for(int i_site = 0; i_site < length_of_microtubule; i_site++){
			if(mt_array[i_mt].lattice[i_site].occupant == NULL){
				// Unoccupied sites correspond to 0
				temp_array[i_site] = 0;
			}
			else{
				// Sites occupied with single-headed motors correspond to 2
				temp_array[i_site] = 2;
			}
		}
		// Transcribe ID data into an integer array for easier data writing
		int temp_ID_array[length_of_microtubule];
		int *ID_array_ptr = temp_ID_array;
		for(int i_site = 0; i_site < length_of_microtubule; i_site++){
			// If site is indeed occupied, get the ID of the motor occupying it
			if(mt_array[i_mt].lattice[i_site].occupant != NULL){
				temp_ID_array[i_site] = mt_array[i_mt].lattice[i_site].occupant->ID;
			}
			// Otherwise, set entry to -1 (no real motor will have this ID)
			else{
				temp_ID_array[i_site] = -1;
			}
		}
		// Write two data arays to their appropriate files 
		fwrite(array_ptr, sizeof(int), length_of_microtubule, output_file);
		fwrite(ID_array_ptr, sizeof(int), length_of_microtubule, ID_file);
	}
}

// Rudimentary ASCII output function
void print_microtubules(system_parameters *parameters, microtubule *mt_array){

	int n_microtubules = parameters->n_microtubules;
	int length_of_microtubule = parameters->length_of_microtubule;

	for (int i_mt = 0; i_mt < n_microtubules; i_mt++){
		for (int i_site = 0; i_site < length_of_microtubule; i_site++){
			if(mt_array[i_mt].lattice[i_site].occupant == NULL)
				printf("=");
			else
				printf("M");
		}
		printf(" %i", mt_array[i_mt].polarity);
		printf("\n");
		fflush(stdout);
	}
}
