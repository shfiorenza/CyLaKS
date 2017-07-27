/* Overlaping microtubules.
  */

#include "hskuan.h"

//extern int temp_end[2];

int main(int argc, char *argv[])
{
   system_parameters parameters;
   char param_file[160];
   int initial_length, n_runs, n_steps, pickup_time, **microtubule_right, **microtubule_left, i, j, i_step, run, **microtubule_right_0,
       **microtubule_left_0, n_protofilaments, ****microtubule_final, *pseudo_end, *temp_end, jj;
   double ***microtubule_right_final, ***microtubule_left_final, random;
   long seed;
   double vg, v_motor_g, f_turning, c_kon, koff, koff_end, duration, c_eff, length_temp, random_check;
   FILE *f_config, *output_file, *stream, *detail_file, *final_config, *final_config2;
   clock_t start, finish;
   
   /* Get command-line input. */
   if (argc != 5) {
      fprintf(stderr, "Usage: %s param_file init.config output_file detail_file\n", argv[0]);
      error_exit("Wrong number of command-line arguments in main");
   }

   start = clock();

   strcpy(param_file, argv[1]);

   /* Read in input parameters. */
   parse_parameters(param_file, &parameters);

   /* Initialize local variables. */
   n_steps = parameters.n_steps;	/* number of simulation runs */
   pickup_time = parameters.pickup_time;/* how frequent to output data */
   seed = parameters.seed;

   /* Unit conversion, 1um=1000nm, 1 tubulin length ~ 8 nm. The standard units are showed in parameters.h file. */
   parameters.vg = parameters.vg*1000.0/8.0/60.0/2000.0;	// change the unit from um/min to tubulin/0.005s
   parameters.shrink = parameters.shrink*1000.0/8.0/60.0/2000.0; // change the unit from um/min to tubulin/0.005s
   parameters.v_motor_g = parameters.v_motor_g*1000.0/8.0/2000.0;	// change the unit from um/s to tubulin/0.005s
   parameters.f_turning = parameters.f_turning/60.0/2000.0;		// change the unit from 1/min to 1/0.005s
   c_eff = parameters.c_motor;                                        // effective c
   //if (parameters.vg != 0){
   //   parameters.vg = 38.53*c_eff/(c_eff+82.61)/8.0/200.0;
   //}
   //c_eff = (sqrt(1.0+8.0*0.4*parameters.c_motor)-1.0)/(4.0*0.4);// assume some dimer forms and only monomers work. The association constant (K) is 0.5. K = c_2/c^2;
   parameters.c_kon = c_eff * parameters.kon/2000.0; // change the unit to /0.005s
   // c_motor is in unit of nM, and kon is in unit of nM-1 s-1
   parameters.koff = parameters.koff/2000.0;	// change the unit to 1/0.005s

   //parameters.shrink = parameters.v_motor_g - parameters.vg;
   length_temp = 0.0;
   if (parameters.vg != 0){
      //parameters.shrink = 0.0286*parameters.c_motor*1000.0/8.0/60.0/2000.0;//parameters.shrink*parameters.c_motor/(parameters.c_motor+84.0);
      //parameters.vg = parameters.vg + 0.0286*parameters.c_motor*1000.0/8.0/60.0/2000.0;
   }

   if (parameters.flux_flag == 2){
      //parameters.vg = (2.4+parameters.c_motor*0.0285)*1000.0/8.0/60.0/2000.0;//parameters.vg + parameters.shrink;	//u-w
      //parameters.shrink = parameters.c_motor*0.0285*1000.0/8.0/60.0/2000.0;//parameters.shrink;	//w
      printf("Calculated flux: %f\n", parameters.shrink*(parameters.vg-parameters.shrink)/parameters.vg/parameters.v_motor_g);
      printf("Calculated end density: %f\n", (parameters.vg-parameters.shrink)/parameters.vg);
      fflush(stdout);
   }

   //parameters.koff_end = parameters.koff_end/60/200;        // change the unit to 1/0.005s
   printf("Growing velocity : %f tubulin/0.005s\nShrinking velocity : %f tubulin/0.005s\nMotor moving velocity : %f tubulin/0.005s\nMotor turning frequency : %f 1/0.005s\n", parameters.vg, parameters.shrink, parameters.v_motor_g, parameters.f_turning);
   printf("kon * c_eff : %f 1/0.005s\nDissociation frequency : %f 1/0.005s\n", parameters.c_kon, parameters.koff);
   printf("rho_0 : %f\nlambda : %f\n", parameters.c_kon/(parameters.c_kon + parameters.koff), parameters.v_motor_g/(parameters.c_kon + parameters.koff));
   printf("Overall frequency %f\n", parameters.v_motor_g+parameters.f_turning+parameters.koff);
   fflush(stdout);

   /* Read initial state */
   f_config = gfopen(argv[2], "rb");
   fprintf(stdout, "reading from config file:\n\n");
   fread(&n_runs, sizeof(int), 1, f_config);
   fread(&n_protofilaments, sizeof(int), 1, f_config);
   fread(&initial_length, sizeof(int), 1, f_config);
   fprintf(stdout, "   n_runs = %d\n", n_runs);
   fprintf(stdout, "   n_protofilaments = %d\n", n_protofilaments);
   fprintf(stdout, "   initial length = %d\n", initial_length);

   microtubule_right = (int**) allocate_2d_array(n_protofilaments, N_SITES_MAX + 1, sizeof(int));
   microtubule_left = (int**) allocate_2d_array(n_protofilaments, N_SITES_MAX + 1, sizeof(int));
   microtubule_right_0 = (int**) allocate_2d_array(n_protofilaments, N_SITES_MAX + 1, sizeof(int));
   microtubule_left_0 = (int**) allocate_2d_array(n_protofilaments, N_SITES_MAX + 1, sizeof(int));
   microtubule_right_final = (double***) allocate_3d_array(n_protofilaments, 2, N_SITES_MAX + 1, sizeof(double));
   microtubule_left_final = (double***) allocate_3d_array(n_protofilaments, 2, N_SITES_MAX + 1, sizeof(double));
   microtubule_final = (int ****) allocate_4d_array((int)(n_steps/pickup_time), n_protofilaments, 2, N_SITES_MAX + 1, sizeof(int));
   pseudo_end = parameters.pseudo_end = (int *) allocate_1d_array(2, sizeof(int));
   temp_end = parameters.temp_end = (int *) allocate_1d_array(4, sizeof(int));
   parameters.n_protofilaments = n_protofilaments;

   stream = fopen ("pseudo_end.dat","r");
   if (stream!=NULL){
      fscanf (stream, "%d %d", &pseudo_end[0], &pseudo_end[1]);
      fclose (stream);
   }
   else{
      pseudo_end[0] = 0;pseudo_end[1] = 0;
   }
   printf("pseudo_end: %d %d\n", pseudo_end[0], pseudo_end[1]);

   for (i=0;i<N_SITES_MAX+1;i++){
      microtubule_right_final[0][0][i] = 0.0;	//right, moving right motor
      microtubule_right_final[0][1][i] = 0.0;     //right, moving left motor
      microtubule_left_final[0][0][i] = 0.0;     //left, moving right motor
      microtubule_left_final[0][1][i] = 0.0;     //left, moving left motor
   }

   for (i=0;i<n_protofilaments;i++){
      fread(microtubule_right_0[i], sizeof(int), N_SITES_MAX+1, f_config);
      fread(microtubule_left_0[i], sizeof(int), N_SITES_MAX+1, f_config);
      if ((microtubule_right_0[i][N_SITES_MAX] >= N_SITES_MAX) || (microtubule_left_0[i][N_SITES_MAX] >= N_SITES_MAX)){
         fprintf(stderr, "The length is too long to do the simulation, change N_SITES_MAX in code.\n");
         exit(1);
      }
   }
   fclose(f_config);

   output_file = gfopen(argv[3], "w");
   detail_file = gfopen(argv[4], "w");
   temp_end[0] = 0;temp_end[1] = 0;temp_end[2] = 0;random_check = 0.0;

   for (run=0;run<n_runs;run++){
      for (i=0;i<n_protofilaments;i++){
         for (j=0;j<N_SITES_MAX+1;j++){
            microtubule_right[i][j] = microtubule_right_0[i][j];
            microtubule_left[i][j] = microtubule_left_0[i][j];
            //if (j<10)
            //   printf("%d %d %d\n", microtubule_right[i][j], microtubule_left[i][j], j);
         }
      }
      for (i_step=0;i_step<n_steps;i_step++){

         //printf("%d\n", i_step);
         //fflush(stdout);

         /* Monte Carlo simulation. Both of them are in overlap_growth.c . */
         for (jj=0;jj<2;jj++){
            random = ran3(&seed);
            if (random < 0.0 || random > 1.0){
               printf("Warning! random number has a problem.\n");
               fflush(stdout);
            }
            random_check += random;

            if (random > 0.5){
               overlap_growth(&parameters, microtubule_right, microtubule_left);
            }
            else{
               motors_move(&parameters, microtubule_right, microtubule_left);
            }
         }
         /* End. */

         //printf("%d\n", i_step);
         //fflush(stdout);

         //printf("%d\n", i_step);
         if (i_step%pickup_time == 0){
            printf("%d", i_step);
            for (i=0;i<n_protofilaments;i++){
               for (j=0;j<(N_SITES_MAX+1);j++){
                  microtubule_final[(int)(i_step/pickup_time)][i][0][j] = microtubule_right[i][j];
                  microtubule_final[(int)(i_step/pickup_time)][i][1][j] = microtubule_left[i][j];
               }
            }
            printf(", length = %d.\n", microtubule_right[0][N_SITES_MAX]+microtubule_left[0][N_SITES_MAX]);
            fflush(stdout);
            //if ((i_step + 1) > (n_steps/2)){
            //   length_temp = length_temp + microtubule_right[0][N_SITES_MAX]+microtubule_left[0][N_SITES_MAX];
            //}
         }
         if ((i_step + 1) == (n_steps/2)){
            temp_end[0] = 0;
            temp_end[1] = 0;
            temp_end[2] = 0;
            temp_end[3] = 0;
         }
         if ((i_step + 1) >= (n_steps/2)){
            length_temp = length_temp + microtubule_right[0][N_SITES_MAX]+microtubule_left[0][N_SITES_MAX];
            if (microtubule_left[0][microtubule_left[0][N_SITES_MAX]-2] == 3 || microtubule_left[0][microtubule_left[0][N_SITES_MAX]-2] == 4){
               //if (microtubule_right[0][microtubule_right[0][N_SITES_MAX]-3] == 2 || microtubule_right[0][microtubule_right[0][N_SITES_MAX]-3] == 4){
                  temp_end[2] = temp_end[2] + 1;
               //}
            }
            if (microtubule_right[0][microtubule_right[0][N_SITES_MAX]-1] == 2 || microtubule_right[0][microtubule_right[0][N_SITES_MAX]-1] == 4){
               temp_end[0] = temp_end[0] + 1;
            }
            if (microtubule_left[0][microtubule_left[0][N_SITES_MAX]-1] == 3 || microtubule_left[0][microtubule_left[0][N_SITES_MAX]-1] == 4){
               if (microtubule_left[0][N_SITES_MAX] == 1){
                  printf("Warning! site-2 doesn't exist properly.\n");
                  fflush(stdout);
               }
            
               if (microtubule_left[0][microtubule_left[0][N_SITES_MAX]-2] == 3 || microtubule_left[0][microtubule_left[0][N_SITES_MAX]-2] == 4){
                  temp_end[1] = temp_end[1] + 1;
               }
            }
            if (microtubule_right[0][microtubule_right[0][N_SITES_MAX]-2] == 2 || microtubule_right[0][microtubule_right[0][N_SITES_MAX]-2] == 4){
               if (microtubule_left[0][N_SITES_MAX] == 1){
                  printf("Warning! site-2 doesn't exist properly.\n");
                  fflush(stdout);
               }

               if (microtubule_right[0][microtubule_right[0][N_SITES_MAX]-3] == 2 || microtubule_right[0][microtubule_right[0][N_SITES_MAX]-3] == 4){
                  temp_end[3] = temp_end[3] + 1;
               }
            }
         }

		// Seems like HSK was attempting to do the averaging in this code at first??

         /* microtubule[i][j] = 1, tubulin only.
         microtubule[i][j] = 1+1, tubulin + motor_right
         microtubule[i][j] = 1+2, tubulin + motor_left
         microtubule[i][j] = 1+1+2, tubulin + motor_right_left *
         for (j=0;j<(N_SITES_MAX+1);j++){
            if (microtubule_right[0][j] == 2 || microtubule_right[0][j] == 4){
               microtubule_right_final[0][0][j] = microtubule_right_final[0][0][j] + (double)(1.0/(double)(n_steps));
               //microtubule_final[(int)(n_steps/pickup_time)][0][0][j] = microtubule_final[(int)(n_steps/pickup_time)][0][0][j]+1;
            }
            else if (microtubule_right[0][j] == 3 || microtubule_right[0][j] == 4){
               microtubule_right_final[0][1][j] = microtubule_right_final[0][1][j] + (double)(1.0/(double)(n_steps));
               //microtubule_final[(int)(n_steps/pickup_time)+1][0][0][j] = microtubule_final[(int)(n_steps/pickup_time)+1][0][0][j]+1;
            }
            if (microtubule_left[0][j] == 2 || microtubule_left[0][j] == 4){
               microtubule_left_final[0][0][j] = microtubule_left_final[0][0][j] + (double)(1.0/(double)(n_steps));
               //microtubule_final[(int)(n_steps/pickup_time)][0][1][j] = microtubule_final[(int)(n_steps/pickup_time)][0][1][j]+1;
            }
            else if (microtubule_left[0][j] == 3 || microtubule_left[0][j] == 4){
               microtubule_left_final[0][1][j] = microtubule_left_final[0][1][j] + (double)(1.0/(double)(n_steps));
               //microtubule_final[(int)(n_steps/pickup_time)+1][0][1][j] = microtubule_final[(int)(n_steps/pickup_time)+1][0][1][j]+1;
            }
         }
         */
      }
   }
   output(microtubule_final, (int)(n_steps/pickup_time), n_protofilaments, output_file, detail_file);
   printf("temp_end: %f %f %f %f\n", (double)((double)(temp_end[0])/(double)(n_steps/2)), (double)((double)(temp_end[1])/(double)(n_steps/2)), 
                                  (double)((double)(temp_end[2])/(double)(n_steps/2)), (double)((double)(temp_end[3])/(double)(n_steps/2)));
   printf("tau_0: %f. tau_0_again: %f. ", (0.00003125*parameters.c_motor-parameters.shrink*(double)((double)(temp_end[1]))/
                                   ((double)(n_steps/2)))/(0.00003125*parameters.c_motor+0.00003125*84.0), 
                                   0.00003125*parameters.c_motor/(0.00003125*parameters.c_motor+0.00003125*84.0+parameters.vg));
   printf("tau_0 tau_1: %f.\n", 0.00003125*parameters.c_motor/(0.00003125*parameters.c_motor+0.00003125*84.0+parameters.vg+parameters.shrink));
   printf("Random number check: %f\n", random_check/2.0/((double)(n_steps)));

   if (parameters.flux_flag == 2){
      //parameters.vg = parameters.vg + parameters.shrink;        //u-w
      //parameters.shrink = parameters.shrink;    //w
      printf("Calculated flux: %f\n", parameters.shrink*(parameters.vg-parameters.shrink)/parameters.vg/parameters.v_motor_g);
      printf("Calculated end density: %f\n", (parameters.vg-parameters.shrink)/parameters.vg);
      fflush(stdout);
   }

   fclose(output_file);
   fclose(detail_file);

   printf("End, final length = %d.\n", microtubule_right[0][N_SITES_MAX]+microtubule_left[0][N_SITES_MAX]);
   length_temp = length_temp/((double)((n_steps/2)));
   fprintf(stdout, "Average length = %f.\n", length_temp);
   fflush(stdout);
/*
   final_config = gfopen("final.config", "wb");
   fwrite(&n_runs, sizeof(int), 1, final_config);
   fwrite(&n_protofilaments, sizeof(int), 1, final_config);
   initial_length = microtubule_right[0][N_SITES_MAX];//+microtubule_left[0][N_SITES_MAX];
   fwrite(&initial_length, sizeof(int), 1, final_config);
   for (i=0;i<n_protofilaments;i++){
      fwrite(microtubule_right[i], sizeof(int), N_SITES_MAX+1, final_config);
      fwrite(microtubule_left[i], sizeof(int), N_SITES_MAX+1, final_config);
   }
   //fflush(final_config);
   fclose(final_config);
*/

/*
   final_config = gfopen("all.config", "wb");
   final_config2 = gfopen("diff.config", "wb");
   double all, diff;
   for (i=0;i<n_protofilaments;i++){
      for (j=0;j<N_SITES_MAX+1;j++){
         all = microtubule_right_final[i][0][j]+microtubule_right_final[i][1][j];
         diff = microtubule_right_final[i][0][j]-microtubule_right_final[i][1][j];
         fwrite(&all, sizeof(double), 1, final_config);
         fwrite(&diff, sizeof(double), 1, final_config2);
      }
      for (j=0;j<N_SITES_MAX+1;j++){
         all = microtubule_left_final[i][0][j]+microtubule_left_final[i][1][j];
         diff = microtubule_left_final[i][0][j]-microtubule_left_final[i][1][j];
         fwrite(&all, sizeof(double), 1, final_config);
         fwrite(&diff, sizeof(double), 1, final_config2);
      }
   }
   //printf("%d\n", microtubule_final[(int)(n_steps/pickup_time)][0][0][0]);
   fclose(final_config);
   fclose(final_config2);
*/   

   finish = clock();

   duration = (double)(start-finish)/CLOCKS_PER_SEC;

   stream = fopen("time_main.dat","w");
   fprintf(stream, "%f seconds\n", duration);
   fclose(stream);
  /* stream = fopen("pseudo_end.dat","w");
   fprintf(stream, "%d %d\n", pseudo_end[0], pseudo_end[1]);
   fclose(stream);
*/
   return 0;
}
