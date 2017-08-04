#include	"hskuan.h"

void overlap_growth(system_parameters * parameters, int **microtubule_right, int **microtubule_left)
{
   int i, k, n_protofilaments, ii, j, site, length_of_microtubule, center, m, flux_flag, gap, *pseudo_end, *temp_end;
   long seed;
   double vg, random, shrink, c_kon, koff, temp, back, back_old, random2, random3, pre_opp, vg_old, shrink_old, random0;
   static int ini_flag = 0;
   static int *microtubule_temp_r, *microtubule_temp_l;

   pseudo_end = parameters->pseudo_end; //pseudo_end[0] left, pseudo_end[1] right. 0 means no motor and 1 means a motor.
   temp_end = parameters->temp_end;
   seed = parameters->seed;
   vg = parameters->vg;
   if (parameters->shrink != 0.0){
      shrink = parameters->shrink;//parameters->vg - parameters->vg*parameters->c_motor/84.0;
   }
   else{
      shrink = parameters->shrink;
   }
   vg = parameters->vg*1.0;//(parameters->c_motor/(84.0+parameters->c_motor));
   //*parameters->c_motor/(parameters->c_motor+84.0);// - (parameters->vg - (parameters->vg*parameters->c_motor/84.0));
   shrink = parameters->shrink;
   c_kon = parameters->c_kon;
   koff = parameters->koff;
   n_protofilaments = parameters->n_protofilaments;
   flux_flag = parameters->flux_flag;
   back = 0.0;//(1.0-parameters->c_motor/(parameters->c_motor+120.0));//82.41;//parameters->right_conc_right;//0.966;//0.0253*10.0;///20.0;
   back_old = back;
   pre_opp = 0.0;//0.008*parameters->c_motor;///(parameters->c_motor+10);//1.0;//parameters->c_motor/(parameters->c_motor+82.41);
   gap = 0;
   vg_old = vg;
   shrink_old = shrink;

   if (ini_flag == 0){
      microtubule_temp_r = (int*) calloc (N_SITES_MAX+1,sizeof(int));
      microtubule_temp_l = (int*) calloc (N_SITES_MAX+1,sizeof(int));
      printf("flux_flag = %d\n", flux_flag);
      printf("vg = %f, shrink = %f, rho_0 = %f\n", vg, shrink, c_kon);
      printf("Old: vg = %f, shrink = %f, back = %f\n", vg_old, shrink_old, back_old);
      printf("back = %f, pre_opp = %f\n", back, pre_opp);
      ini_flag = 1;
   }

   /* microtubule[i][the last one] tells the total length of microtubule, right means the (+) end is at right end, and left means the (+) end
      at the left end. 0 is at the microtubule_right[i][0]. */
   for (ii=0;ii<n_protofilaments;ii++){
      random = ran3(&seed);
      if (random < 0.5){	//This works when n_protofilaments = 2;
         i = 0;
      }
      else{
         i = 1;
      }
      i = ii;

      for (m=0;m<2;m++){	// 2 means right and left.
         random0 = ran3(&seed);
         if (random0 > 0.5){	//right
            /* Right end */
            site = microtubule_right[i][N_SITES_MAX]-1;	//end site
            if (flux_flag != 2){
               if (microtubule_right[i][site-gap] == 2 || microtubule_right[i][site-gap] == 4){
                  if (flux_flag == 1){
                     random = ran3(&seed);
                     if (microtubule_right[i][site] == 2 || microtubule_right[i][site] == 4){
                        //if (random < shrink){
                           //microtubule_right[0][N_SITES_MAX] -= 1;
                           //microtubule_right[1][N_SITES_MAX] -= 1;
                        //}
                        //if (random < (vg)){
                        //   microtubule_right[0][N_SITES_MAX] += 1;
                           //microtubule_right[1][N_SITES_MAX] += 1;
                        //}
                     }
                     else if (random < shrink){
                        microtubule_right[0][N_SITES_MAX] -= 1;
                        //microtubule_right[1][N_SITES_MAX] += 1;
                     }
                  }
                  else{
                     random = ran3(&seed);
                     /*
                     printf("Check: site = %d, microtubule_right[i][site-gap] = %d, gap = %d, shrink = %f, vg = %f\n",
                            site, microtubule_right[i][site-gap], gap, shrink, vg);
                     printf("Check: microtubule_right[i][site-gap-1] = %d, microtubule_right[i][site-gap+1] = %d, n_prota = %d, i = %d, microtubule_right[0][N_SITES_MAX] = %d\n",
                            microtubule_right[i][site-gap-1], microtubule_right[i][site-gap+1], n_protofilaments, i, microtubule_right[0][N_SITES_MAX]);
                     fflush(stdout);
                     */

                     if (random < shrink){
                        microtubule_right[0][N_SITES_MAX] -= 1;

                        /**/
                        random2 = ran3(&seed);
                        if (microtubule_right[i][site-gap+1] == 2 || microtubule_right[i][site-gap+1] == 4){
                           back = 0.0;//back_old;
                        }
                        else{
                           back = 1.0;//back_old;
                        }

                        //back = 1.0-parameters->c_motor/(parameters->c_motor+84.0);

                        if (random2 < back){
                           if ((site-gap) != 0 ){
                              if (microtubule_right[i][site-gap-1] == 1){
                                 microtubule_right[i][site-gap-1] = 2;
                              }
                              else if (microtubule_right[i][site-gap-1] == 3){
                                 microtubule_right[i][site-gap-1] = 4;
                              }
                           }
                           else{
                              if (microtubule_left[i][0] == 1){
                                 microtubule_left[i][0] = 2;
                              }
                              else if (microtubule_left[i][0] == 3){
                                 microtubule_left[i][0] = 4;
                              }
                              if (gap !=0 ){
                                 printf("Warning!!\n");
                              }
                           }
                        }
                        //microtubule_right[i][site-gap] = microtubule_right[i][site];
                        /**/
                        if (microtubule_right[i][site-gap] == 4){
                           microtubule_right[i][site-gap] = 3;
                        }
                        else{
                           microtubule_right[i][site-gap] = 1;
                        }
                        /**/
                        /**/
                     }
                     else if (random < (shrink+vg)){
                        microtubule_right[0][N_SITES_MAX] += 1;
                        microtubule_right[i][site+1] = microtubule_right[i][site];
                        microtubule_right[i][site] = 1;
                        /*
                        random3 = ran3(&seed);
                        if (random3 < pre_opp){
                        //if (pseudo_end[1] == 1){
                           microtubule_right[i][site+1] = 2;
                           //pseudo_end[1] = 0;
                        }
                        else{
                           microtubule_right[i][site+1] = 1;
                        }
                        //microtubule_right[1][N_SITES_MAX] += 1;
                        */
                     }
                  }
               }
               else{
                  /*
                  if (microtubule_right[i][site-gap+1] == 2 || microtubule_right[i][site-gap+1] == 4){
                     shrink = 0.0;
                     vg = 0.0;
                  }
                  else{
                     shrink = shrink_old;
                     vg = vg_old;
                  }
                  */
                  random = ran3(&seed);
                  /**/
                  if (flux_flag == 1){
                     if (microtubule_right[i][site] == 2 || microtubule_right[i][site] == 4){
                        //if (random < shrink){
                        //   microtubule_right[i][N_SITES_MAX] -= 1;
                        //}
                        if (random < (vg)){
                           microtubule_right[i][N_SITES_MAX] += 1;
                           if (microtubule_right[i][site] == 2){
                              microtubule_right[i][site] == 1;
                              microtubule_right[i][site+1] == 2;
                           }
                           else{
                              microtubule_right[i][site] == 3;
                              microtubule_right[i][site+1] == 2;
                           }
                        }
                     }
                  }
                  /**/
                  else{
                     if (random < vg){
                        microtubule_right[0][N_SITES_MAX] += 1;
                        microtubule_right[i][site+1] = microtubule_right[i][site];
                        microtubule_right[i][site] = 1;
                        /*
                        random3 = ran3(&seed);
                        if (random3 < pre_opp){
                        //if (pseudo_end[1] == 1){
                           microtubule_right[i][site+1] = 2;
                           //pseudo_end[1] = 0;
                        }
                        else{
                           microtubule_right[i][site+1] = 1;
                        }
                        //microtubule_right[1][N_SITES_MAX] += 1;
                        */
                     }
                  }

                  if (microtubule_right[i][site] == 0){
                     printf("Fatal error.\n");
                     exit(1);
                  }
               }
            }
            else{
               random = ran3(&seed);
               if (microtubule_right[i][site-gap] == 2 || microtubule_right[i][site-gap] == 4){
                  //printf("Right\n");
                  if (random < shrink){
                     microtubule_right[i][N_SITES_MAX] -= 1;
                     //temp_end[2] = temp_end[2] + 1;
                  }
                  //fflush(stdout);
               }
               else{
                  //printf("Right, %d ", microtubule_right[i][site-gap]);
                  if (random < shrink){
                     microtubule_right[i][N_SITES_MAX] -= 1;
                     //temp_end[2] = temp_end[2] + 1;
                  }
                  else if (random < (shrink+vg)){
                     microtubule_right[i][N_SITES_MAX] += 1;
                  }
                  //printf("%f %d %d\n", random, microtubule_right[i][N_SITES_MAX], microtubule_right[i][site-gap]);
                  //fflush(stdout);
               }
            }

            if (microtubule_right[i][N_SITES_MAX] >= N_SITES_MAX){ //make microtubule[i][N_SITES_MAX-1] always 0.
               fprintf(stderr, "The microtubule length exceeds the boundary, change N_SITES_MAX in code or reset the parameters.\n");
               exit(1);
            }
            else if (microtubule_right[i][N_SITES_MAX]+microtubule_left[i][N_SITES_MAX] <= 0){
               fprintf(stderr, "The microtubule length is 0, please check the parameters.\n");
               exit(1);
               //microtubule_right[i][N_SITES_MAX] = 1;
            }
            else if (microtubule_right[i][N_SITES_MAX] < 3 || microtubule_left[i][N_SITES_MAX] < 3){
               fprintf(stderr, "Warning, either microtubule_right[i][N_SITES_MAX] or microtubule_left[i][N_SITES_MAX] is smaller than 3.\n");
            }
            /* End. */
         }
         else{	//left
            /* Left end */
            site = microtubule_left[i][N_SITES_MAX]-1;       //end site
            if (flux_flag != 2){
               if (microtubule_left[i][site-gap] == 3 || microtubule_left[i][site-gap] == 4){
                  if (flux_flag == 1){
                     random = ran3(&seed);
                     if (microtubule_left[i][site] == 3 || microtubule_left[i][site] == 4){
                        //if (random < shrink){
                        //   microtubule_left[0][N_SITES_MAX] -= 1;
                           //microtubule_left[1][N_SITES_MAX] -= 1;
                        //}
                        //if (random < (vg)){
                           //microtubule_left[0][N_SITES_MAX] += 1;
                           //microtubule_left[1][N_SITES_MAX] += 1;
                        //}
                     }
                     else if (random < shrink){
                        microtubule_left[0][N_SITES_MAX] -= 1;
                        //microtubule_left[1][N_SITES_MAX] += 1;
                     }
                  }
                  else{
                     random = ran3(&seed);
                     /*
                     printf("Check: site = %d, microtubule_left[i][site-gap] = %d, gap = %d, shrink = %f, vg = %f\n",
                            site, microtubule_left[i][site-gap], gap, shrink, vg);
                     printf("Check: microtubule_left[i][site-gap-1] = %d, microtubule_left[i][site-gap+1] = %d, n_prota = %d, i = %d, microtubule_left[0][N_SITES_MAX] = %d\n",
                            microtubule_left[i][site-gap-1], microtubule_left[i][site-gap+1], n_protofilaments, i, microtubule_left[0][N_SITES_MAX]);
                     fflush(stdout);
                     */

                     if (random < shrink){
                        microtubule_left[0][N_SITES_MAX] -= 1;

                        /**/
                        random2 = ran3(&seed);
                        if (microtubule_left[i][site-gap+1] == 3 || microtubule_left[i][site-gap+1] == 4){
                           back = 0.0;//back_old;
                        }
                        else{
                           back = 1.0;//back_old;
                        }

                        //back = 1.0-parameters->c_motor/(parameters->c_motor+84.0);

                        if (random2 < back){
                           if ((site-gap) != 0){
                              if (microtubule_left[i][site-gap-1] == 1){
                                 microtubule_left[i][site-gap-1] = 3;
                              }
                              else if (microtubule_left[i][site-gap-1] == 2){
                                 microtubule_left[i][site-gap-1] = 4;
                              }
                           }
                           else{
                              if (microtubule_right[i][0] == 1){
                                 microtubule_right[i][0] = 3;
                              }
                              else if (microtubule_right[i][0] == 2){
                                 microtubule_right[i][0] = 4;
                              }
                              if (gap != 0){
                                 printf("Warning!\n");
                              }
                           }
                        }
                        //microtubule_left[i][site-gap] = microtubule_left[i][site];
                        /**/
                        if (microtubule_left[i][site-gap] == 4){
                           microtubule_left[i][site-gap] = 2;
                        }
                        else{
                           microtubule_left[i][site-gap] = 1;
                        }
                        /**/
                        /**/
                     }
                     else if (random < (shrink+vg)){
                        microtubule_left[0][N_SITES_MAX] += 1;
                        microtubule_left[i][site+1] = microtubule_left[i][site];
                        microtubule_left[i][site] = 1;
                     }
                  }
               }
               else{
                  random = ran3(&seed);
                  /**/
                  if (flux_flag == 1){
                     if (microtubule_left[i][site] == 3 || microtubule_left[i][site] == 4){
                        //if (random < shrink){
                        //   microtubule_left[i][N_SITES_MAX] -= 1;
                        //}
                        if (random < (vg)){
                           microtubule_left[i][N_SITES_MAX] += 1;
                           if (microtubule_left[i][site] == 3){
                              microtubule_left[i][site] == 1;
                              microtubule_left[i][site+1] == 3;
                           }
                           else{
                              microtubule_left[i][site] == 2;
                              microtubule_left[i][site+1] == 3;
                           }
                        }
                     }
                  }
                  /**/
                  else{
                     if (random < vg){
                        microtubule_left[0][N_SITES_MAX] += 1;
                        microtubule_left[i][site+1] = microtubule_left[i][site];
                        microtubule_left[i][site] = 1;
                     }
                  }

                  if (microtubule_left[i][site] == 0){
                     printf("Fatal error. left\n");
                     exit(1);
                  }
               }
            }
            else{
               random = ran3(&seed);
               if (microtubule_left[i][site-gap] == 3 || microtubule_left[i][site-gap] == 4){
                  //printf("Left\n");
                  if (random < shrink){
                     microtubule_left[i][N_SITES_MAX] -= 1;
                     //temp_end[2] = temp_end[2] + 1;
                  }
                  //fflush(stdout);
               }
               else{
                  //printf("Left, %d ", microtubule_left[i][site-gap]);
                  if (random < shrink){
                     microtubule_left[i][N_SITES_MAX] -= 1;
                     //temp_end[1] = temp_end[1] + 1;
                  }
                  else if (random < (shrink+vg)){
                     microtubule_left[i][N_SITES_MAX] += 1;
                     //temp_end[2] = temp_end[2] + 1;
                  }
                  //printf("%f %d %d\n", random, microtubule_left[i][N_SITES_MAX], microtubule_left[i][site-gap]);
                  //fflush(stdout);
               }
            }

            if (microtubule_left[i][N_SITES_MAX] >= N_SITES_MAX){       //make microtubule[i][N_SITES_MAX-1] always 0.
               fprintf(stderr, "The microtubule length exceeds the boundary, change N_SITES_MAX in code or reset the parameters.\n");
               exit(1);
            }
            else if (microtubule_right[i][N_SITES_MAX]+microtubule_left[i][N_SITES_MAX] <= 0){
               fprintf(stderr, "The microtubule length is 0, please check the parameters.\n");
               //microtubule_left[i][N_SITES_MAX] = 1;
               exit(1);
            }
            else if (microtubule_right[i][N_SITES_MAX] < 3 || microtubule_left[i][N_SITES_MAX] < 3){
               fprintf(stderr, "Warning, either microtubule_right[i][N_SITES_MAX] or microtubule_left[i][N_SITES_MAX] is smaller than 3.\n");
            }
            /* End. */
         }

         /* >0 means the site is occuppied with tubulin, 0 means the site is empty. microtubule[i][N_SITES_MAX] won't be touched. */
         /* 1 means the site is occuppied with the tubulin only, and no other tubulin on it. */
         /**/
         for (k=0;k<n_protofilaments;k++){
         for (j=0;j<microtubule_right[k][N_SITES_MAX];j++){
            if (microtubule_right[k][j] == 0){
               microtubule_right[k][j] = 1;
            }
         }
         for (j=microtubule_right[k][N_SITES_MAX];j<N_SITES_MAX;j++){
            microtubule_right[k][j] = 0;
         }
         for (j=0;j<microtubule_left[k][N_SITES_MAX];j++){
            if (microtubule_left[k][j] == 0){
               microtubule_left[k][j] = 1;
            }
         }
         for (j=microtubule_left[k][N_SITES_MAX];j<N_SITES_MAX;j++){
            microtubule_left[k][j] = 0;
         }
         /**/
         /* Centering */
         if (microtubule_right[k][N_SITES_MAX] < 3 || microtubule_left[k][N_SITES_MAX] < 3){
            for (j=0;j<N_SITES_MAX;j++){
               microtubule_temp_r[j] = 0;
               microtubule_temp_l[j] = 0;
            }
            printf("%d %d %d   ", microtubule_right[k][N_SITES_MAX] + microtubule_left[k][N_SITES_MAX], 
                                  microtubule_right[k][N_SITES_MAX], microtubule_left[k][N_SITES_MAX]);
            fflush(stdout);

            length_of_microtubule = microtubule_right[k][N_SITES_MAX] + microtubule_left[k][N_SITES_MAX];
            center = length_of_microtubule/2;

            //printf("%d %d %d %d\n", microtubule_right[i][N_SITES_MAX], microtubule_left[i][N_SITES_MAX], length_of_microtubule, center);
            //fflush(stdout);

            if (center > microtubule_right[k][N_SITES_MAX]){
               for (j=0;j<center;j++){
                  microtubule_temp_l[j] = microtubule_left[k][microtubule_left[k][N_SITES_MAX]-center+j];
               }
               for (j=0;j<(microtubule_left[k][N_SITES_MAX]-center);j++){
                  microtubule_temp_r[j] = microtubule_left[k][microtubule_left[k][N_SITES_MAX]-center-1-j];
               }
               for (j=0;j<microtubule_right[k][N_SITES_MAX];j++){
                  microtubule_temp_r[j+microtubule_left[k][N_SITES_MAX]-center] = microtubule_right[k][j];
               }
               microtubule_temp_r[N_SITES_MAX] = length_of_microtubule-center;
               microtubule_temp_l[N_SITES_MAX] = center;
            }
            else{
               for (j=0;j<center;j++){
                  microtubule_temp_r[j] = microtubule_right[k][microtubule_right[k][N_SITES_MAX]-center+j];
               }
               for (j=0;j<(microtubule_right[k][N_SITES_MAX]-center);j++){
                  microtubule_temp_l[j] = microtubule_right[k][microtubule_right[k][N_SITES_MAX]-center-1-j];
               }
               for (j=0;j<microtubule_left[k][N_SITES_MAX];j++){
                  microtubule_temp_l[j+microtubule_right[k][N_SITES_MAX]-center] = microtubule_left[k][j];
               }
               microtubule_temp_r[N_SITES_MAX] = center;
               microtubule_temp_l[N_SITES_MAX] = length_of_microtubule-center;
            }
            for (j=0;j<(N_SITES_MAX+1);j++){
               microtubule_right[k][j] = microtubule_temp_r[j];
               microtubule_left[k][j] = microtubule_temp_l[j];
            }
            printf("%d %d\n", length_of_microtubule, microtubule_left[k][N_SITES_MAX]);
            fflush(stdout);
         }
         }
         /* End.*/
      }
   }

   //printf("pseudo_end: %d %d %d %d\n", pseudo_end[0], pseudo_end[1], parameters->pseudo_end[0], parameters->pseudo_end[1]);
   //fflush(stdout);

   return;
}

void motors_move(system_parameters * parameters, int **microtubule_right, int **microtubule_left)
{
   int n_protofilaments, length_of_microtubule, i, j, temp, site, ii, *pseudo_end;
   long seed;
   double random, random2, c_kon, v_motor_g, koff, f_turning, vg, r_c_r, r_c_l, l_c_r, l_c_l, rho0, v_motor_g_boundary, v_motor_g_old, qb, qv, c_kon_end;
   double one_l_c_l, one_r_c_r, shrink, koff_end[2], f_turning_end, random3, c_kon_old, factor, koff_old, factor_off, factor_off2, koff_endd, c_kon_endd;
   static int flag = 0;

   pseudo_end = parameters->pseudo_end; //pseudo_end[0] left, pseudo_end[1] right. 0 means no motor and 1 means a motor.
   //printf("motor move: %d %d %d %d\n", pseudo_end[0], pseudo_end[1], parameters->pseudo_end[0], parameters->pseudo_end[1]);
   //fflush(stdout);
   seed = parameters->seed;
   n_protofilaments = parameters->n_protofilaments;
   c_kon = parameters->c_kon;
   c_kon_old = c_kon;
   v_motor_g = parameters->v_motor_g;
   v_motor_g_boundary = 1.0*v_motor_g;
   //*(1.0+1.0/(parameters->c_motor+82.61)*parameters->right_conc_right);//*parameters->right_conc_right;
   v_motor_g_old = v_motor_g;
   vg = parameters->vg;
   koff = parameters->koff;
   koff_old = koff;
   f_turning = parameters->f_turning;
   rho0 = c_kon/(c_kon+koff);
   shrink = parameters->shrink;
   c_kon_end = c_kon;//1.0*0.00003125*parameters->c_motor;
   //0.0*c_kon;//0.003125*parameters->c_motor;//1000.0*c_kon;//+shrink*vg/(shrink-vg)*(1-parameters->c_motor/(parameters->c_motor+84.0));
   c_kon_endd = 1.0*c_kon;
   //0.0;//c_kon;//parameters->c_motor*vg;//c_kon;//0.001*parameters->c_motor;//c_kon;//0.001*parameters->c_motor;//* parameters->right_conc_right;
   koff_end[0] = koff;//84.0*0.003125;//1000.0*koff;//0.0083;//*(parameters->c_motor/(parameters->c_motor+84.0));
   //parameters->c_motor/1500.0;//koff;//shrink*(parameters->c_motor/(parameters->c_motor+84.0));//koff;//0.08241;//koff;//0.08241;
   koff_end[1] = koff;//84.0*0.003125;//1000.0*koff;//0.0083;//*(parameters->c_motor/(parameters->c_motor+84.0));
   koff_endd = 1.0*koff;
   factor = 1.0;//1.0+100.0/parameters->c_motor;
   factor_off = 1.0;
   factor_off2 = 1.0;//1.0/(1.0*koff)*2.0*(parameters->vg-parameters->shrink)*parameters->vg/(parameters->shrink*(parameters->vg-parameters->shrink));
   f_turning_end = 1.0*f_turning;//0.0;//f_turning;//0.0;//parameters->f_turning;//0.0;
   r_c_r = v_motor_g * parameters->right_conc_right;
   r_c_l = v_motor_g * parameters->right_conc_left;//vg * 80.6 / (parameters->c_motor + 80.6);//v_motor_g * parameters->right_conc_left;
   l_c_r = v_motor_g * parameters->left_conc_right;//vg * 80.6 / (parameters->c_motor + 80.6);//v_motor_g * parameters->left_conc_right;
   l_c_l = v_motor_g * parameters->left_conc_left;
   qb = c_kon*0.0; qv = v_motor_g*0.0; one_l_c_l = 0.0; one_r_c_r = 0.0;

   /*
   if (pseudo_end[0] == 0){
      koff_end[0] = 0;
   }
   else if (pseudo_end[0] == 1){
      koff_end[0] = shrink;
   }
   if (pseudo_end[1] == 0){
      koff_end[1] = 0;
   }
   else if (pseudo_end[1] == 1){
      koff_end[1] = shrink;
   }
   */

   //r_c_r = v_motor_g * rho0;
   //l_c_l = v_motor_g * rho0;
   if (flag ==0){
      printf("r_c_r = %f\n", r_c_r);	//1.0
      printf("r_c_l = %f\n", r_c_l);	//0.0
      printf("l_c_r = %f\n", l_c_r);	//0.0
      printf("l_c_l = %f\n", l_c_l);	//1.0
      printf("c_kon_end = %f, koff_end = %f\n", c_kon_end, koff_end[0]);
      printf("v_motor_g*qv = %f\n", qv);
      printf("v_motor_g_boundary = %f\n", v_motor_g_boundary);
      flag = 1;
   }

   //printf("%f\n", r_c_r);fflush(stdout);

   //printf("bbbbb %d %d\n", n_protofilaments, parameters->n_protofilaments);
   //fflush(stdout);

   /* microtubule[i][j] = 1, tubulin only.
      microtubule[i][j] = 1+1, tubulin + motor_right
      microtubule[i][j] = 1+2, tubulin + motor_left
      microtubule[i][j] = 1+1+2, tubulin + motor_right_left
      For the following for-loop, microtubule[i][length_of_microtubule-1] must greater than 0,
      i.e. it has at least one tubulin; otherwise, the calculations are meaningless. */
   for (ii=0;ii<n_protofilaments;ii++){
      random = ran3(&seed);
      if (random<0.5){	//This works when n_protofilaments = 2.
         i = 0;
      }
      else{
         i = 1;
      }
      i = ii;

      length_of_microtubule = microtubule_right[i][N_SITES_MAX]+microtubule_left[i][N_SITES_MAX];
      //printf("%d %d %d\n", length_of_microtubule, microtubule_right[i][N_SITES_MAX], microtubule_left[i][N_SITES_MAX]);
      //fflush(stdout);
      
      /*
      random = ran3(&seed);
      if (random < parameters->c_motor/(parameters->c_motor+84.0)){
         if (microtubule_left[i][(microtubule_left[i][N_SITES_MAX]-1)] == 2 || microtubule_left[i][(microtubule_left[i][N_SITES_MAX]-1)] == 4){
            microtubule_left[i][(microtubule_left[i][N_SITES_MAX]-1)] = 4;
         }
         else{
            microtubule_left[i][(microtubule_left[i][N_SITES_MAX]-1)] = 3;
         }
      }
      else{
         if (microtubule_left[i][(microtubule_left[i][N_SITES_MAX]-1)] == 2 || microtubule_left[i][(microtubule_left[i][N_SITES_MAX]-1)] == 4){
            microtubule_left[i][(microtubule_left[i][N_SITES_MAX]-1)] = 2;
         }
         else{
            microtubule_left[i][(microtubule_left[i][N_SITES_MAX]-1)] = 1;
         }
      }
      random = ran3(&seed);
      if (random < parameters->c_motor/(parameters->c_motor+84.0)){
         if (microtubule_right[i][(microtubule_right[i][N_SITES_MAX]-1)] == 3 || microtubule_right[i][(microtubule_right[i][N_SITES_MAX]-1)] == 4){
            microtubule_right[i][(microtubule_right[i][N_SITES_MAX]-1)] = 4;
         }
         else{
            microtubule_right[i][(microtubule_right[i][N_SITES_MAX]-1)] = 2;
         }
      }
      else{
         if (microtubule_right[i][(microtubule_right[i][N_SITES_MAX]-1)] == 3 || microtubule_right[i][(microtubule_right[i][N_SITES_MAX]-1)] == 4){
            microtubule_right[i][(microtubule_right[i][N_SITES_MAX]-1)] = 3;
         }
         else{
            microtubule_right[i][(microtubule_right[i][N_SITES_MAX]-1)] = 1;
         }
      }
      */

      /* The end site has a different dissociation constant, so treat it differently. 
         microtubule[i][j]%4 = 3 means the site has tubulin + motor; tubulin_array[i][j]%4 = 0 means the site has nothing. 
         Since microtubule[i][N_SITES_MAX-1] is always 0, I don't have to worry if j+1 = N_SITES_MAX. */
      for (j=0;j<(2*length_of_microtubule);j++){
         random = ran3(&seed);     // The random number is from 0 to 1.
         if (random > 1 || random < 0){
            printf("error. random > 1 || random < 0\n");
         }
         temp = (int)(random*((double)(2*length_of_microtubule)));	// Pick a site, 2 means right or left
         //temp = 24;
         //printf("%d %f %d %d\n", temp, random, (int)(temp/2), (int)((temp-2*microtubule_left[i][N_SITES_MAX])/2));	// Check
         //fflush(stdout);
         
         random2 = ran3(&seed);
         /*
         if (pseudo_end[0] == 0){
            if (random2 < (c_kon/10.0)){
               pseudo_end[0] = 1;
               koff_end[0] = shrink;
            }
         }
         else{
            if (random2 < 84.0/(c_kon/10.0/parameters->c_motor)){
               pseudo_end[0] = 0;
               koff_end[0] = 0;
            }
         }
         random2 = ran3(&seed);
         if (pseudo_end[1] == 0){
            if (random2 < (c_kon/10.0)){
               pseudo_end[1] = 1;
               koff_end[1] = shrink;
            }
         }
         else{
            if (random2 < 84.0/(c_kon/10.0/parameters->c_motor)){
               pseudo_end[1] = 0;
               koff_end[1] = 0;
            }
         }
         */

         if (temp < (2*microtubule_left[i][N_SITES_MAX])){
            site = (int)(temp/2);

            /* Treat boundary independently. */
            if (site == (microtubule_left[i][N_SITES_MAX]-1)){
               random2 = ran3(&seed);
               if (microtubule_left[i][site] == 1){
                  if (temp%2 == 0){	//motors moving toward right
                     /*
                     if (microtubule_left[i][site-1] == 2 || microtubule_left[i][site-1] == 4){
                        c_kon_end = c_kon_old*factor;
                     }
                     else{
                        c_kon_end = c_kon_old;
                     }
                     */

                     if (random2 < l_c_r){
                        microtubule_left[i][site] = 2;
                     }
                     else if (random2 < (l_c_r+c_kon_endd)){
                        microtubule_left[i][site] = 2;
                     }

                     //c_kon_end = c_kon_old;
                  }
                  else{			//motors moving toward left
                     /*
                     if (microtubule_left[i][site-1] == 3 || microtubule_left[i][site-1] == 4){
                        c_kon_end = c_kon_old*factor;
                     }
                     else{
                        c_kon_end = c_kon_old;
                     }
                     */

                     if (random2 < c_kon_end){
                        microtubule_left[i][site] = 3;
                     }

                     //c_kon_end = c_kon_old;
                  }
               }
               else if (microtubule_left[i][site] == 2){
                  if (temp%2 == 0){	//motors moving toward right
                     if (microtubule_left[i][site-1] == 2 || microtubule_left[i][site-1] == 4){
                        koff = factor_off * koff_old;
                     }
                     else{
                        koff = koff_old;
                     }

                     if (random2 < koff_endd){
                        microtubule_left[i][site] = 1;
                     }
                     else if (random2 < (koff_endd+f_turning_end)){             // turn around
                        microtubule_left[i][site] = 3;
                     }
                     else if (microtubule_left[i][site-1] == 1){	// walk toward right
                        if (random2 < (koff_endd+f_turning_end+v_motor_g_boundary)){
                           microtubule_left[i][site] = 1;
                           microtubule_left[i][site-1] = 2;
                        }
                     }
                     else if (microtubule_left[i][site-1] == 3){	// walk toward right, the site is 3 not empty (1)
                        if (random2 < (koff_endd+f_turning_end+v_motor_g_boundary-qv)){		// no motor_right at the targeting site
                           microtubule_left[i][site] = 1;
                           microtubule_left[i][site-1] = 4;
                        }
                     }

                     koff = koff_old;
                  }
                  else{			//motors moving toward left
                     /*
                     if (microtubule_left[i][site-1] == 3 || microtubule_left[i][site-1] == 4){
                        c_kon_end = c_kon_old*factor;
                     }
                     else{
                        c_kon_end = c_kon_old;
                     }
                     */

                     if (random2 < (c_kon_end-qb)){
                        microtubule_left[i][site] = 4;
                     }

                     //c_kon_end = c_kon_old;
                  }
               }
               else if (microtubule_left[i][site] == 3){
                  if (temp%2 == 0){	//motors moving toward right
                     /*
                     if (microtubule_left[i][site-1] == 2 || microtubule_left[i][site-1] == 4){
                        c_kon_end = c_kon_old*factor;
                     }
                     else{
                        c_kon_end = c_kon_old;
                     }
                     */

                     if (random2 < (c_kon_endd-qb)){
                        microtubule_left[i][site] = 4;
                     }
                     else if (random2 < (l_c_r+c_kon_endd-qb)){
                        microtubule_left[i][site] = 4;
                     }

                     //c_kon_end = c_kon_old;
                  }
                  else{

                     if (v_motor_g < l_c_l){
                        one_l_c_l = l_c_l;
                     }
                     else{
                        one_l_c_l = v_motor_g - l_c_l;
                     }

                     if (microtubule_left[i][site-1] == 3 || microtubule_left[i][site-1] == 4){
                        koff = factor_off * koff_old;
                     }
                     else{
                        koff = koff_old;
                     }

                     if (random2 < koff_end[1]){	//(v_motor_g - l_c_l)
                        microtubule_left[i][site] = 1;
                     }
                     else if (random2 < (koff_end[1]+f_turning_end)){         // turn around
                        microtubule_left[i][site] = 2;
                     }
                     else if (random2 < (koff_end[1]+f_turning_end+one_l_c_l)){
                        microtubule_left[i][site] = 1;
                     }
                     else if (random2 < (koff_end[1]+f_turning_end+one_l_c_l+0.0*v_motor_g)){
                        if (pseudo_end[0] == 0){
                           microtubule_left[i][site] = 1;
                           //pseudo_end[0] = 1;
                        }
                     }

                     koff = koff_old;
                  }
               }
               else if (microtubule_left[i][site] == 4){
                  if (temp%2 == 0){     //motors moving toward right
                     if (microtubule_left[i][site-1] == 2 || microtubule_left[i][site-1] == 4){
                        koff = factor_off * koff_old;
                     }
                     else{
                        koff = koff_old;
                     }

                     if (random2 < koff_endd){
                        microtubule_left[i][site] = 3;
                     }
                     //else if (random2 < (koff+qb)){
                     //   microtubule_left[i][site] = 3;
                     //}
                     else if (microtubule_left[i][site-1] == 1){        // walk toward right
                        if (random2 < (koff_endd+v_motor_g_boundary)){
                           microtubule_left[i][site] = 3;
                           microtubule_left[i][site-1] = 2;
                        }
                     }
                     else if (microtubule_left[i][site-1] == 3){        // walk toward right, the site is 3 not empty (1)
                        if (random2 < (koff_endd+v_motor_g_boundary-qv)){           // no motor_right at the targeting site
                           microtubule_left[i][site] = 3;
                           microtubule_left[i][site-1] = 4;
                        }
                     }

                     koff = koff_old;
                  }
                  else{

                     if (microtubule_left[i][site-1] == 3 || microtubule_left[i][site-1] == 4){
                        koff = factor_off * koff_old;
                     }
                     else{
                        koff = koff_old;
                     }

                     if (v_motor_g < l_c_l ){
                        one_l_c_l = l_c_l;
                     }
                     else{
                        one_l_c_l = v_motor_g - l_c_l;
                     }

                     if (random2 < koff_end[1]){
                        microtubule_left[i][site] = 2;
                     }
                     //else if (random2 < (koff+qb)){
                     //   microtubule_left[i][site] = 2;
                     //}
                     else if (random2 < (koff_end[1]+one_l_c_l)){
                        microtubule_left[i][site] = 2;
                     }
                     else if (random2 < (koff_end[1]+one_l_c_l+0.0*v_motor_g)){
                        if (pseudo_end[0] == 0){
                           microtubule_left[i][site] = 2;
                           //pseudo_end[0] = 1;
                        }
                     }

                     koff = koff_old;
                  }
               }
            }
            /* End. */
            else{
               /* Pick microtubule_left[i][(int)(temp/2)] */
               if (microtubule_left[i][site]==1){				// No motor
                  if (temp%2 == 0){						// attach motor_right
                     if (site != 0){
                     //if (site == (microtubule_left[i][N_SITES_MAX]-2)){
                        if (microtubule_left[i][site-1] == 2 || microtubule_left[i][site-1] == 4 ||
                            microtubule_left[i][site+1] == 2 || microtubule_left[i][site+1] == 4){
                           c_kon = c_kon_old*factor;
                        }
                        else{
                           c_kon = c_kon_old;
                        }
                     }
                     else{
                        if (microtubule_right[i][0] == 2 || microtubule_right[i][0] == 4 ||
                           microtubule_left[i][site+1] == 2 || microtubule_left[i][site+1] == 4){
                           c_kon = c_kon_old*factor;
                        }
                        else{
                           c_kon = c_kon_old;
                        }
                     }

                     random = ran3(&seed);
                     if (random < c_kon){
                        microtubule_left[i][site] = 2;
                     }

                     c_kon = c_kon_old;
                  }
                  else {							// attach motor_left
                     if (site != 0){
                        if (microtubule_left[i][site-1] == 3 || microtubule_left[i][site-1] == 4 ||
                            microtubule_left[i][site+1] == 3 || microtubule_left[i][site+1] == 4){
                           c_kon = c_kon_old*factor;
                        }
                        else{
                           c_kon = c_kon_old;
                        }
                     }
                     else{
                        if (microtubule_right[i][0] == 3 || microtubule_right[i][0] == 4 ||
                           microtubule_left[i][site+1] == 3 || microtubule_left[i][site+1] == 4){
                           c_kon = c_kon_old*factor;
                        }
                        else{
                           c_kon = c_kon_old;
                        }
                     }

                     random = ran3(&seed);
                     if (random < c_kon){
                        microtubule_left[i][site] = 3;
                     }

                     c_kon = c_kon_old;
                  }
               }
               else if (microtubule_left[i][site]==2){			// motor_right
                  if (temp%2 == 0){
                     /*
                     if (site != 0){
                        if (microtubule_left[i][site-1] == 2 || microtubule_left[i][site-1] == 4 || 
                            microtubule_left[i][site+1] == 2 || microtubule_left[i][site+1] == 4){
                           koff = factor_off * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }
                     else{
                        if (microtubule_right[i][0] == 2 || microtubule_right[i][0] == 4 ||
                            microtubule_left[i][site+1] == 2 || microtubule_left[i][site+1] == 4){
                           koff = factor_off * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }
                     */
                     if (site == (microtubule_left[i][N_SITES_MAX]-2)){
                        if (microtubule_left[i][site+1] == 2 || microtubule_left[i][site+1] == 4){
                           koff = koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }

                     random = ran3(&seed);
                     if (random < koff){
                        microtubule_left[i][site] = 1;
                     }
                     else if (random < (koff+f_turning)){             // turn around
                        microtubule_left[i][site] = 3;
                     }
                     else if (site != 0){
                        if (microtubule_left[i][site-1] == 1){               // no motor_right at the targeting site
                           if (random < (koff+f_turning+v_motor_g)){
                              microtubule_left[i][site] = 1;
                              microtubule_left[i][site-1] = 2;
                           }
                        }
                        else if (microtubule_left[i][site-1] == 3){               // no motor_right at the targeting site
                           if (random < (koff+f_turning+v_motor_g-qv)){
                              microtubule_left[i][site] = 1;
                              microtubule_left[i][site-1] = 4;
                           }
                        }
                     }
                     else if (site == 0){
                        if (microtubule_right[i][0] == 1){           // no motor_right at the targeting site
                           if (random < (koff+f_turning+v_motor_g)){
                              microtubule_left[i][site] = 1;
                              microtubule_right[i][0] = 2;
                           }
                        }
                        else if (microtubule_right[i][0] == 3){
                           if (random < (koff+f_turning+v_motor_g-qv)){
                              microtubule_left[i][site] = 1;
                              microtubule_right[i][0] = 4;
                           }
                        }
                     }

                     koff = koff_old;
                  }
                  else{							// attach motor_left
                     if (site != 0){
                        if (microtubule_left[i][site-1] == 3 || microtubule_left[i][site-1] == 4 ||
                            microtubule_left[i][site+1] == 3 || microtubule_left[i][site+1] == 4){
                           c_kon = c_kon_old*factor;
                        }
                        else{
                           c_kon = c_kon_old;
                        }
                     }
                     else{
                        if (microtubule_right[i][0] == 3 || microtubule_right[i][0] == 4 ||
                           microtubule_left[i][site+1] == 3 || microtubule_left[i][site+1] == 4){
                           c_kon = c_kon_old*factor;
                        }
                        else{
                           c_kon = c_kon_old;
                        }
                     }

                     random = ran3(&seed);
                     if (random < (c_kon-qb)){
                        microtubule_left[i][site] = 4;
                     }

                     c_kon = c_kon_old;
                  }
               }
               else if (microtubule_left[i][site]==3){			// motor_left
                  if (temp%2 == 0){						// attach motor_right
                     if (site != 0){
                        if (microtubule_left[i][site-1] == 2 || microtubule_left[i][site-1] == 4 ||
                            microtubule_left[i][site+1] == 2 || microtubule_left[i][site+1] == 4){
                           c_kon = c_kon_old*factor;
                        }
                        else{
                           c_kon = c_kon_old;
                        }
                     }
                     else{
                        if (microtubule_right[i][0] == 2 || microtubule_right[i][0] == 4 ||
                           microtubule_left[i][site+1] == 2 || microtubule_left[i][site+1] == 4){
                           c_kon = c_kon_old*factor;
                        }
                        else{
                           c_kon = c_kon_old;
                        }
                     }

                     random = ran3(&seed);
                     if (random < (c_kon-qb)){
                        microtubule_left[i][site] = 4;
                     }

                     c_kon = c_kon_old;
                  }
                  else{
                     /*
                     if (site != 0){
                        if (microtubule_left[i][site-1] == 3 || microtubule_left[i][site-1] == 4 ||
                            microtubule_left[i][site+1] == 3 || microtubule_left[i][site+1] == 4){
                           koff = factor_off * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }
                     else{
                        if (microtubule_right[i][0] == 3 || microtubule_right[i][0] == 4 ||
                            microtubule_left[i][site+1] == 3 || microtubule_left[i][site+1] == 4){
                           koff = factor_off * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }
                     */

                     if (site == (microtubule_left[i][N_SITES_MAX]-2)){
                        if (microtubule_left[i][site+1] == 3 || microtubule_left[i][site+1] == 4){
                           koff = factor_off2 * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }

                     random = ran3(&seed);
                     if (site == (microtubule_left[i][N_SITES_MAX]-2)){
                        v_motor_g = v_motor_g_boundary;
                     }

                     if (random < koff){
                        microtubule_left[i][site] = 1;
                     }
                     else if (random < (koff+f_turning)){         // turn around
                        microtubule_left[i][site] = 2;
                     }
                     else if (microtubule_left[i][site+1] == 1){
                        if (random < (koff+f_turning+v_motor_g)){
                           microtubule_left[i][site] = 1;
                           microtubule_left[i][site+1] = 3;
                        }
                     }
                     else if (microtubule_left[i][site+1] == 2){
                        if (random < (koff+f_turning+v_motor_g-qv)){
                           microtubule_left[i][site] = 1;
                           microtubule_left[i][site+1] = 4;
                        }
                     }

                     if (site == (microtubule_left[i][N_SITES_MAX]-2)){
                        v_motor_g = v_motor_g_old;
                     }

                     koff = koff_old;
                  }
               }
               else if (microtubule_left[i][site]==4){
                  if (temp%2 == 0){
                     /*
                     if (site != 0){
                        if (microtubule_left[i][site-1] == 2 || microtubule_left[i][site-1] == 4 ||
                            microtubule_left[i][site+1] == 2 || microtubule_left[i][site+1] == 4){
                           koff = factor_off * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }
                     else{
                        if (microtubule_right[i][0] == 2 || microtubule_right[i][0] == 4 ||
                            microtubule_left[i][site+1] == 2 || microtubule_left[i][site+1] == 4){
                           koff = factor_off * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }
                     */

                     if (site == (microtubule_left[i][N_SITES_MAX]-2)){
                        if (microtubule_left[i][site+1] == 2 || microtubule_left[i][site+1] == 4){
                           koff = koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }

                     random = ran3(&seed);
                     if (random < koff){
                        microtubule_left[i][site]=3;
                     }
	             //else if (random < (koff+qb)){
                     //   microtubule_left[i][site]=3;
                     //}
                     else if (site != 0){
                        if (microtubule_left[i][site-1] == 1){               // no motor_right at the targeting site
                           if (random < (koff+v_motor_g)){
                              microtubule_left[i][site] = 3;
                              microtubule_left[i][site-1] = 2;
                           }
                        }
                        else if (microtubule_left[i][site-1] == 3){               // no motor_right at the targeting site
                           if (random < (koff+v_motor_g-qv)){
                              microtubule_left[i][site] = 3;
                              microtubule_left[i][site-1] = 4;
                           }
                        }
                     }
                     else if (site == 0){
                        if (microtubule_right[i][0] == 1){           // no motor_right at the targeting site
                           if (random < (koff+v_motor_g)){
                              microtubule_left[i][site] = 3;
                              microtubule_right[i][0] = 2;
                           }
                        }
                        else if (microtubule_right[i][0] == 3){
                           if (random < (koff+v_motor_g-qv)){
                              microtubule_left[i][site] = 3;
                              microtubule_right[i][0] = 4;
                           }
                        }
                     }

                     koff = koff_old;
                  }
                  else{
                     /*
                     if (site != 0){
                        if (microtubule_left[i][site-1] == 3 || microtubule_left[i][site-1] == 4 ||
                            microtubule_left[i][site+1] == 3 || microtubule_left[i][site+1] == 4){
                           koff = factor_off * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }
                     else{
                        if (microtubule_right[i][0] == 3 || microtubule_right[i][0] == 4 ||
                            microtubule_left[i][site+1] == 3 || microtubule_left[i][site+1] == 4){
                           koff = factor_off * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }
                     */

                     if (site == (microtubule_left[i][N_SITES_MAX]-2)){
                        if (microtubule_left[i][site+1] == 3 || microtubule_left[i][site+1] == 4){
                           koff = factor_off2 * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }

                     random = ran3(&seed);

                     if (site == (microtubule_left[i][N_SITES_MAX]-2)){
                        v_motor_g = v_motor_g_boundary;
                     }

                     if (random < koff){
                        microtubule_left[i][site]=2;
                     }
                     //else if (random < (koff+qb)){
                     //   microtubule_left[i][site]=2;
                     //}
                     else if (microtubule_left[i][site+1] == 1){
                        if (random < (koff+v_motor_g)){
                           microtubule_left[i][site] = 2;
                           microtubule_left[i][site+1] = 3;
                        }
                     }
                     else if (microtubule_left[i][site+1] == 2){
                        if (random < (koff+v_motor_g-qv)){
                           microtubule_left[i][site] = 2;
                           microtubule_left[i][site+1] = 4;
                        }
                     }

                     if (site == (microtubule_left[i][N_SITES_MAX]-2)){
                        v_motor_g = v_motor_g_old;
                     }

                     koff = koff_old;
                  }
               }
               /* End. */	
            }
         }
         else{	//temp >= 2*microtubule_left[i][N_SITES_MAX]
            /* Pick microtubule_right[i][(int)((temp-2*microtubule_left[i][N_SITES_MAX])/2)] */
            site = (int)((temp-2*microtubule_left[i][N_SITES_MAX])/2);

            /* Treat boundary independently. */
            if (site == (microtubule_right[i][N_SITES_MAX]-1)){
               random2 = ran3(&seed);
               if (microtubule_right[i][site] == 1){
                  if (temp%2 == 0){	//motors moving toward right
                     /*
                     if (microtubule_right[i][site-1] == 2 || microtubule_right[i][site-1] == 4){
                        c_kon_end = c_kon_old*factor;
                     }
                     else{
                        c_kon_end = c_kon_old;
                     }
                     */

                     if (random2 < c_kon_end){
                        microtubule_right[i][site] = 2;
                     }

                     //c_kon_end = c_kon_old;
                  }
                  else{			//motors moving toward left
                     /*
                     if (microtubule_right[i][site-1] == 3 || microtubule_right[i][site-1] == 4){
                        c_kon_end = c_kon_old*factor;
                     }
                     else{
                        c_kon_end = c_kon_old;
                     }
                     */

                     if (random2 < r_c_l){
                        microtubule_right[i][site] = 3;
                     }
                     else if (random2 < (r_c_l+c_kon_endd)){
                        microtubule_right[i][site] = 3;
                     }

                     //c_kon_end = c_kon_old;
                  }
               }
               else if (microtubule_right[i][site] == 2){
                  if (temp%2 == 0){	//motors moving toward right
                     if (microtubule_right[i][site-1] == 2 || microtubule_right[i][site-1] == 4){
                        koff = factor_off * koff_old;
                     }
                     else{
                        koff = koff_old;
                     }

                     if (v_motor_g < r_c_r ){
                        one_r_c_r = r_c_r;
                     }
                     else{
                        one_r_c_r = v_motor_g - r_c_r;
                     }

                     if (random2 < koff_end[0]){	//(v_motor_g - r_c_r)
                        microtubule_right[i][site] = 1;
                     }
                     else if (random2 < (koff_end[0]+one_r_c_r)){
                        microtubule_right[i][site] = 1;
                     }
                     else if (random2 < (koff_end[0]+f_turning_end+one_r_c_r)){
                        microtubule_right[i][site] = 3;
                     }
                     else if (random2 < (koff_end[0]+f_turning_end+one_r_c_r+0.0*v_motor_g)){
                        if (pseudo_end[1] == 0){
                           microtubule_right[i][site] = 1;
                           //pseudo_end[1] = 1;
                        }
                     }

                     koff = koff_old;
                  }
                  else{
                     /*
                     if (microtubule_right[i][site-1] == 3 || microtubule_right[i][site-1] == 4){
                        c_kon_end = c_kon_old*factor;
                     }
                     else{
                        c_kon_end = c_kon_old;
                     }
                     */

                     if (random2 < (c_kon_endd-qb)){
                        microtubule_right[i][site] = 4;
                     }
                     else if (random2 < (r_c_l+c_kon_endd-qb)){
                        microtubule_right[i][site] = 4;
                     }

                     //c_kon_end = c_kon_old;
                  }
               }
               else if (microtubule_right[i][site] == 3){
                  if (temp%2 == 0){	//motors moving toward right
                     /*
                     if (microtubule_right[i][site-1] == 2 || microtubule_right[i][site-1] == 4){
                        c_kon_end = c_kon_old*factor;
                     }
                     else{
                        c_kon_end = c_kon_old;
                     }
                     */

                     if (random2 < (c_kon_end-qb)){
                        microtubule_right[i][site] = 4;
                     }

                     //c_kon_end = c_kon_old;
                  }
                  else{			//motors moving toward left
                     if (microtubule_right[i][site-1] == 3 || microtubule_right[i][site-1] == 4){
                        koff = factor_off * koff_old;
                     }
                     else{
                        koff = koff_old;
                     }

                     if (random2 < koff_endd){
                        microtubule_right[i][site] = 1;
                     }
                     else if (random2 < (koff_endd+f_turning_end)){
                        microtubule_right[i][site] = 2;
                     }
                     else if (microtubule_right[i][site-1] == 1){
                        if (random2 < (koff_endd+f_turning_end+v_motor_g_boundary)){
                           microtubule_right[i][site] = 1;
                           microtubule_right[i][site-1] = 3;
                        }
                     }
                     else if (microtubule_right[i][site-1] == 2){
                        if (random2 < (koff_endd+f_turning_end+v_motor_g_boundary-qv)){
                           microtubule_right[i][site] = 1;
                           microtubule_right[i][site-1] = 4;
                        }
                     }

                     koff = koff_old;
                  }
               }
               else if (microtubule_right[i][site] == 4){
                  if (temp%2 == 0){     //motors moving toward right
                     if (microtubule_right[i][site-1] == 2 || microtubule_right[i][site-1] == 4){
                        koff = factor_off * koff_old;
                     }
                     else{
                        koff = koff_old;
                     }

                     if (v_motor_g < r_c_r ){
                        one_r_c_r = r_c_r;
                     }
                     else{
                        one_r_c_r = v_motor_g - r_c_r;
                     }

                     if (random2 < koff_end[0]){
                        microtubule_right[i][site] = 3;
                     }
                     //else if (random2 < (koff + qb)){
                     //   microtubule_right[i][site] = 3;
                     //}
                     else if (random2 < (koff_end[0]+one_r_c_r)){
                        microtubule_right[i][site] = 3;
                     }
                     else if (random2 < (koff_end[0]+one_r_c_r+0.0*v_motor_g)){
                        if (pseudo_end[1] == 0){
                           microtubule_right[i][site] = 3;
                           //pseudo_end[1] = 1;
                        }
                     }

                     koff = koff_old;
                  }
                  else{
                     if (microtubule_right[i][site-1] == 3 || microtubule_right[i][site-1] == 4){
                        koff = factor_off * koff_old;
                     }
                     else{
                        koff = koff_old;
                     }

                     if (random2 < koff_endd){
                        microtubule_right[i][site] = 2;
                     }
                     //else if (random2 < (koff+qb)){
                     //   microtubule_right[i][site] = 2;
                     //}
                     else if (microtubule_right[i][site-1] == 1){
                        if (random2 < (koff_endd+v_motor_g_boundary)){
                           microtubule_right[i][site] = 2;
                           microtubule_right[i][site-1] = 3;
                        }
                     }
                     else if (microtubule_right[i][site-1] == 2){
                        if (random2 < (koff_endd+v_motor_g_boundary-qv)){
                           microtubule_right[i][site] = 2;
                           microtubule_right[i][site-1] = 4;
                        }
                     }

                     koff = koff_old;
                  }
               }
            }
            /* End. */
            else{
               if (microtubule_right[i][site]==1){			// No motor
                  if (temp%2 == 0){						// motor_right
                     if (site != 0){
                        if (microtubule_right[i][site-1] == 2 || microtubule_right[i][site-1] == 4 ||
                            microtubule_right[i][site+1] == 2 || microtubule_right[i][site+1] == 4){
                           c_kon = c_kon_old*factor;
                        }
                        else{
                           c_kon = c_kon_old;
                        }
                     }
                     else{
                        if (microtubule_left[i][0] == 2 || microtubule_left[i][0] == 4 ||
                            microtubule_right[i][site+1] == 2 || microtubule_right[i][site+1] == 4){
                           c_kon = c_kon_old*factor;
                        }
                        else{
                           c_kon = c_kon_old;
                        }
                     }

                     random = ran3(&seed);
                     if (random < c_kon){
                        microtubule_right[i][site] = 2;
                     }

                     c_kon = c_kon_old;
                  }
                  else{							// motor_left
                     if (site != 0){
                        if (microtubule_right[i][site-1] == 3 || microtubule_right[i][site-1] == 4 ||
                            microtubule_right[i][site+1] == 3 || microtubule_right[i][site+1] == 4){
                           c_kon = c_kon_old*factor;
                        }
                        else{
                           c_kon = c_kon_old;
                        }
                     }
                     else{
                        if (microtubule_left[i][0] == 3 || microtubule_left[i][0] == 4 ||
                            microtubule_right[i][site+1] == 3 || microtubule_right[i][site+1] == 4){
                           c_kon = c_kon_old*factor;
                        }
                        else{
                           c_kon = c_kon_old;
                        }
                     }

                     random = ran3(&seed);
                     if (random < c_kon){
                        microtubule_right[i][site] = 3;
                     }

                     c_kon = c_kon_old;
                  }
               }
               else if (microtubule_right[i][site]==2){			// motor_right
                  if (temp%2 == 0){                                         // motor_right
                     /*
                     if (site != 0){
                        if (microtubule_right[i][site-1] == 2 || microtubule_right[i][site-1] == 4 ||
                            microtubule_right[i][site+1] == 2 || microtubule_right[i][site+1] == 4){
                           koff = factor_off * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }
                     else{
                        if (microtubule_left[i][0] == 2 || microtubule_left[i][0] == 4 ||
                            microtubule_right[i][site+1] == 2 || microtubule_right[i][site+1] == 4){
                           koff = factor_off * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }
                     */

                     if (site == (microtubule_right[i][N_SITES_MAX]-2)){
                        if (microtubule_right[i][site+1] == 2 || microtubule_right[i][site+1] == 4){
                           koff = factor_off2 * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }

                     random = ran3(&seed);

                     if (site == (microtubule_right[i][N_SITES_MAX]-2)){
                        v_motor_g = v_motor_g_boundary;
                     }

                     if (random < koff){
                        microtubule_right[i][site] = 1;
                     }
                     else if (random < (koff+f_turning)){
                        microtubule_right[i][site] = 3;
                     }
                     else if (microtubule_right[i][site+1] == 1){
                        if (random < (koff+f_turning+v_motor_g)){
                           microtubule_right[i][site] = 1;
                           microtubule_right[i][site+1] = 2;
                        }
                     }
                     else if (microtubule_right[i][site+1] == 3){
                        if (random < (koff+f_turning+v_motor_g-qv)){
                           microtubule_right[i][site] = 1;
                           microtubule_right[i][site+1] = 4;
                        }
                     }

                     if (site == (microtubule_right[i][N_SITES_MAX]-2)){
                        v_motor_g = v_motor_g_old;
                     }

                     koff = koff_old;
                  }
                  else{                                                     // Attach motor_left
                     if (site != 0){
                        if (microtubule_right[i][site-1] == 3 || microtubule_right[i][site-1] == 4 ||
                            microtubule_right[i][site+1] == 3 || microtubule_right[i][site+1] == 4){
                           c_kon = c_kon_old*factor;
                        }
                        else{
                           c_kon = c_kon_old;
                        }
                     }
                     else{
                        if (microtubule_left[i][0] == 3 || microtubule_left[i][0] == 4 ||
                            microtubule_right[i][site+1] == 3 || microtubule_right[i][site+1] == 4){
                           c_kon = c_kon_old*factor;
                        }
                        else{
                           c_kon = c_kon_old;
                        }
                     }

                     random = ran3(&seed);
                     if (random < (c_kon-qb)){
                        microtubule_right[i][site] = 4;
                     }

                     c_kon = c_kon_old;
                  }
               }
               else if (microtubule_right[i][site]==3){			// motor_left
                  if (temp%2 == 0){						// Attach motor_right
                     if (site != 0){
                        if (microtubule_right[i][site-1] == 2 || microtubule_right[i][site-1] == 4 ||
                            microtubule_right[i][site+1] == 2 || microtubule_right[i][site+1] == 4){
                           c_kon = c_kon_old*factor;
                        }
                        else{
                           c_kon = c_kon_old;
                        }
                     }
                     else{
                        if (microtubule_left[i][0] == 2 || microtubule_left[i][0] == 4 ||
                            microtubule_right[i][site+1] == 2 || microtubule_right[i][site+1] == 4){
                           c_kon = c_kon_old*factor;
                        }
                        else{
                           c_kon = c_kon_old;
                        }
                     }

                     random = ran3(&seed);
                     if (random < (c_kon-qb)){
                        microtubule_right[i][site] = 4;
                     }

                     c_kon = c_kon_old;
                  }
                  else{							// motor_left
                     /*
                     if (site != 0){
                        if (microtubule_right[i][site-1] == 3 || microtubule_right[i][site-1] == 4 ||
                            microtubule_right[i][site+1] == 3 || microtubule_right[i][site+1] == 4){
                           koff = factor_off * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }
                     else{
                        if (microtubule_left[i][0] == 3 || microtubule_left[i][0] == 4 ||
                            microtubule_right[i][site+1] == 3 || microtubule_right[i][site+1] == 4){
                           koff = factor_off * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }
                     */

                     if (site == (microtubule_right[i][N_SITES_MAX]-2)){
                        if (microtubule_right[i][site+1] == 3 || microtubule_right[i][site+1] == 4){
                           koff = koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }

                     random = ran3(&seed);
                     if (random < koff){
                        microtubule_right[i][site] = 1;
                     }
                     else if (random < (koff+f_turning)){
                        microtubule_right[i][site] = 2;
                     }
                     else if (site != 0){
                        if (microtubule_right[i][site-1] == 1){
                           if (random < (koff+f_turning+v_motor_g)){
                              microtubule_right[i][site] = 1;
                              microtubule_right[i][site-1] = 3;
                           }
                        }
                        else if (microtubule_right[i][site-1] == 2){
                           if (random < (koff+f_turning+v_motor_g-qv)){
                              microtubule_right[i][site] = 1;
                              microtubule_right[i][site-1] = 4;
                           }
                        }
                     }
                     else if (site == 0){
                        if (microtubule_left[i][0] == 1){
                           if (random < (koff+f_turning+v_motor_g)){
                              microtubule_right[i][site] = 1;
                              microtubule_left[i][0] = 3;
                           }
                        }
                        else if (microtubule_left[i][0] == 2){
                           if (random < (koff+f_turning+v_motor_g-qv)){
                              microtubule_right[i][site] = 1;
                              microtubule_left[i][0] = 4;
                           }
                        }
                     }

                     koff = koff_old;
                  }
               }
               else if (microtubule_right[i][site]==4){
                  if (temp%2 == 0){                                             // Attach motor_right
                     /*
                     if (site != 0){
                        if (microtubule_right[i][site-1] == 2 || microtubule_right[i][site-1] == 4 ||
                            microtubule_right[i][site+1] == 2 || microtubule_right[i][site+1] == 4){
                           koff = factor_off * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }
                     else{
                        if (microtubule_left[i][0] == 2 || microtubule_left[i][0] == 4 ||
                            microtubule_right[i][site+1] == 2 || microtubule_right[i][site+1] == 4){
                           koff = factor_off * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }
                     */

                     if (site == (microtubule_right[i][N_SITES_MAX]-2)){
                        if (microtubule_right[i][site+1] == 2 || microtubule_right[i][site+1] == 4){
                           koff = factor_off2 * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }

                     random = ran3(&seed);

                     if (site == (microtubule_right[i][N_SITES_MAX]-2)){
                        v_motor_g = v_motor_g_boundary;
                     }

                     if (random < koff){
                        microtubule_right[i][site] = 3;
                     }
                     //else if (random < (koff+qb)){
                     //   microtubule_right[i][site] = 3;
                     //}
                     else if (microtubule_right[i][site+1] == 1){
                        if (random < (koff+v_motor_g)){
                           microtubule_right[i][site] = 3;
                           microtubule_right[i][site+1] = 2;
                        }
                     }
                     else if (microtubule_right[i][site+1] == 3){
                        if (random < (koff+v_motor_g-qv)){
                           microtubule_right[i][site] = 3;
                           microtubule_right[i][site+1] = 4;
                        }
                     }

                     if (site == (microtubule_right[i][N_SITES_MAX]-2)){
                        v_motor_g = v_motor_g_old;
                     }

                     koff = koff_old;
                  }
                  else{
                     /*
                     if (site != 0){
                        if (microtubule_right[i][site-1] == 3 || microtubule_right[i][site-1] == 4 ||
                            microtubule_right[i][site+1] == 3 || microtubule_right[i][site+1] == 4){
                           koff = factor_off * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }
                     else{
                        if (microtubule_left[i][0] == 3 || microtubule_left[i][0] == 4 ||
                            microtubule_right[i][site+1] == 3 || microtubule_right[i][site+1] == 4){
                           koff = factor_off * koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }
                     */

                     if (site == (microtubule_right[i][N_SITES_MAX]-2)){
                        if (microtubule_right[i][site+1] == 3 || microtubule_right[i][site+1] == 4){
                           koff = koff_old;
                        }
                        else{
                           koff = koff_old;
                        }
                     }

                     random = ran3(&seed);
                     if (random < koff){
                        microtubule_right[i][site] = 2;
                     }
                     //else if (random < (koff+qb)){
                     //   microtubule_right[i][site] = 2;
                     //}
                     else if (site != 0){
                        if (microtubule_right[i][site-1] == 1){
                           if (random < (koff+v_motor_g)){
                              microtubule_right[i][site] = 2;
                              microtubule_right[i][site-1] = 3;
                           }
                        }
                        else if (microtubule_right[i][site-1] == 2){
                           if (random < (koff+v_motor_g-qv)){
                              microtubule_right[i][site] = 2;
                              microtubule_right[i][site-1] = 4;
                           }
                        }
                     }
                     else if (site == 0){
                        if (microtubule_left[i][0] == 1){
                           if (random < (koff+v_motor_g)){
                              microtubule_right[i][site] = 2;
                              microtubule_left[i][0] = 3;
                           }
                        }
                        else if (microtubule_left[i][0] == 2){
                           if (random < (koff+v_motor_g-qv)){
                              microtubule_right[i][site] = 2;
                              microtubule_left[i][0] = 4;
                           }
                        }
                     }

                     koff = koff_old;
                  }
               }
            }
            /* End. */
         }
      }

      /* Set up boundary condition. 
      * microtubule[i][j] = 1, tubulin only.
      microtubule[i][j] = 1+1, tubulin + motor_right
      microtubule[i][j] = 1+2, tubulin + motor_left
      microtubule[i][j] = 1+1+2, tubulin + motor_right_left *
      random = ran3(&seed);
      if (random < 0.8){	// 0 and 1. set boundary conc = 0 and 0.75.
         microtubule_left[i][microtubule_left[i][N_SITES_MAX]-1] = 3;
      }
      else{
         microtubule_left[i][microtubule_left[i][N_SITES_MAX]-1] = 1;
      }
      random = ran3(&seed);
      if (random < 0.8){      // 0 and 1. set boundary conc = 0 and 0.75.
         microtubule_right[i][microtubule_right[i][N_SITES_MAX]-1] = 2;
      }
      else{
         microtubule_right[i][microtubule_right[i][N_SITES_MAX]-1] = 1;
      }
      */

      //microtubule_left[i][microtubule_left[i][N_SITES_MAX]-1] = 3;   // only have a motor turns left at the left end, and no motor turns right.
      //microtubule_right[i][microtubule_right[i][N_SITES_MAX]-1] = 2;   // only have a motor turns right at the right end, and no motor turns left.
      /* End. */
 
   }

   return;
}
