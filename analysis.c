#include "hskuan.h"

int main(int argc, char *argv[])
{
  FILE *file, *output;
  int ***tubulin_array;
  int time = 1e5;
  double **final;
  int i, j, flag;

  /* Get command-line input. */
  if (argc != 3) {
    fprintf(stderr, "Usage: %s file flag\n", argv[0]);
    error_exit("Wrong number of command-line arguments in main");
  }

  flag = atoi(argv[2]);

  tubulin_array = (int***) allocate_3d_array(time, 1, N_SITES_MAX + 1, sizeof(int));
  final = (double**) allocate_2d_array(1, N_SITES_MAX + 1, sizeof(double));

  file = gfopen(argv[1], "rb");
  for (i=0;i<time;i++){
     fread(tubulin_array[i][0], sizeof(int), N_SITES_MAX + 1, file);
  }
  fclose(file);

  if (flag == 1){
     for (i=0;i<time;i++){
        for (j=0;j<N_SITES_MAX + 1;j++){
           if (tubulin_array[i][0][j] == 3){
              final[0][j] += 1;
           }
        }
     }
     for (i=0;i<N_SITES_MAX + 1;i++){
        final[0][j] = final[0][j]/time;
     }

     output = gfopen("data/ana.bin", "wb");
     fwrite(final[0], sizeof(double), N_SITES_MAX + 1, output);
     fclose(output);
  }
  else{
     output = gfopen("data/ana.bin", "wb");
     for (i=0;i<time;i++){
        fwrite(&tubulin_array[i][0][10], sizeof(int), 1, output);
     }
     fclose(output);
  }

  return 0;
}
