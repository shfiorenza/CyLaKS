#include "hskuan.h"

void output(int ****microtubule_final, int time, int n_protofilaments, FILE *output_file, FILE *detail_file)
{
  int i, j;

  for (j=0;j<time;j++){ 
     for (i=0;i<n_protofilaments;i++){
  //      fwrite(&microtubule_final[j][i][0][N_SITES_MAX], sizeof(int), 1, output_file);
    //    fwrite(&microtubule_final[j][i][1][N_SITES_MAX], sizeof(int), 1, output_file);
        //fprintf(output_file, "%d\n", tubulin_array[i][N_SITES_MAX]);
        fwrite(microtubule_final[j][i][0], sizeof(int), N_SITES_MAX, detail_file);
        fwrite(microtubule_final[j][i][1], sizeof(int), N_SITES_MAX, detail_file);
        //printf("%d %d %d\n", N_SITES_MAX+1, N_SITES_MAX+1, 2*(N_SITES_MAX+1));
     }
  }

  return;
}
