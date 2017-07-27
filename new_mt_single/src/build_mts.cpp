#include "master_header.h"

int main(int argc, char *argv[])
{
	microtubule *mt_array;
 
	int n_run, n_microtubules, length_of_microtubule;
	FILE *output_file;

	/* Get command-line input. */
	if (argc != 5) {
		printf("Wrong number of command-line arguments in main\n");
		fprintf(stderr, "Usage: %s n_run n_microtubules length_of_microtubule output_file\n", argv[0]);
		exit(1);
	}

	n_run = atoi(argv[1]);
	n_microtubules = atoi(argv[2]);
	length_of_microtubule = atoi(argv[3]);

	// Below is not true; FIX FOR BETTER GENERALITY
	if(n_microtubules%2 != 0){
		printf("Error. Must have even pairs of microtubules for proper overlap\n");
		exit(1);
	}

    // Allocate dynamic memory to microtubule pointer and treat its memory space as an array
	mt_array = (microtubule*) malloc(n_microtubules*sizeof(microtubule));
	for(int i = 0; i < n_microtubules; i++){
		
		// Allocate dynamic memory to the track/polarity of each microtubule 
		mt_array[i].track = (site*) malloc(length_of_microtubule*sizeof(site));
		for(int j = 0; j < length_of_microtubule; j++){
			mt_array[i].track[j].occupancy = (int*) malloc(sizeof(int));
			mt_array[i].track[j].occupancy[0] = 0;
		
			mt_array[i].track[j].coord = (int*) malloc(sizeof(int));
			mt_array[i].track[j].coord[0] = j; 
		}

		mt_array[i].polarity = (int*) malloc(sizeof(int));
		if(i%2 == 0){
			mt_array[i].polarity[0] = 0;
		}
		else{
			mt_array[i].polarity[0] = 1;
		}

		mt_array[i].coord = (int*) malloc(sizeof(int));
		mt_array[i].coord[0] = 0;

		mt_array[i].n_bound = (int*) malloc(sizeof(int));
		mt_array[i].n_bound[0] = 0; 
	}

	output_file = fopen(argv[4], "wb");
	fwrite(&n_run, sizeof(int), 1, output_file);
	fwrite(&n_microtubules, sizeof(int), 1, output_file);
	fwrite(&length_of_microtubule, sizeof(int), 1, output_file);

	for (int i=0;i<n_microtubules;i++){
		for(int j = 0; j < length_of_microtubule; j++){
			fwrite(mt_array[i].track[j].occupancy, sizeof(int), 1 , output_file);
			fwrite(mt_array[i].track[j].coord, sizeof(int), 1, output_file);
		}
		fwrite(mt_array[i].polarity, sizeof(int), 1, output_file);
		fwrite(mt_array[i].coord, sizeof(int), 1, output_file);
		fwrite(mt_array[i].n_bound, sizeof(int), 1, output_file);
	}

	fclose(output_file);

	return 0;
}
