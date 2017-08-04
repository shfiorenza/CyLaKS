#include "hskuan.h"

int main(int argc, char *argv[])
{
  int n_run, n_protofilaments, length_of_microtubule, i, j;
  int **microtubule_right, **microtubule_left;
  FILE *output_file;

  /* Get command-line input. */
  if (argc != 5) {
    fprintf(stderr, "Usage: %s n_run n_protofilaments length_of_microtubule output_file\n", argv[0]);
    error_exit("Wrong number of command-line arguments in main");
  }

  n_run = atoi(argv[1]);
  n_protofilaments = atoi(argv[2]);
  length_of_microtubule = atoi(argv[3]);

  microtubule_right = (int**) allocate_2d_array(n_protofilaments, N_SITES_MAX + 1, sizeof(int));
  microtubule_left = (int**) allocate_2d_array(n_protofilaments, N_SITES_MAX + 1, sizeof(int));

  output_file = fopen(argv[4], "wb");
  fwrite(&n_run, sizeof(int), 1, output_file);
  fwrite(&n_protofilaments, sizeof(int), 1, output_file);
  fwrite(&length_of_microtubule, sizeof(int), 1, output_file);
  for (i=0;i<n_protofilaments;i++){
     for (j=0;j<length_of_microtubule;j++){
        microtubule_right[i][j] = 1;
        microtubule_left[i][j] = 1;	//symmetric in left and right.
     }
     microtubule_right[i][N_SITES_MAX] = length_of_microtubule;
     microtubule_left[i][N_SITES_MAX] = length_of_microtubule;
     fwrite(microtubule_right[i], sizeof(int), N_SITES_MAX + 1, output_file);
     fwrite(microtubule_left[i], sizeof(int), N_SITES_MAX + 1, output_file);
  }
  
  fclose(output_file);

  return 0;
}
