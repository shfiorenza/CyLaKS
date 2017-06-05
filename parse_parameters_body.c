if (strcmp(param_name, "n_steps") == 0) {
   parameters -> n_steps = atoi(param_value);
   fprintf(stdout, "   n_steps = %d\n", parameters -> n_steps);
}
else if (strcmp(param_name, "seed") == 0) {
   parameters -> seed = atol(param_value);
   fprintf(stdout, "   seed = %ld\n", parameters -> seed);
}
else if (strcmp(param_name, "growing_velocity") == 0){
   parameters -> vg = atof(param_value);
   fprintf(stdout, "   vg = %g\n", parameters -> vg);
}
else if (strcmp(param_name, "shrinking_velocity") == 0){
   parameters -> shrink = atof(param_value);
   fprintf(stdout, "   shrink = %g\n", parameters -> shrink);
}
else if (strcmp(param_name, "motor_velocity") == 0){
   parameters -> v_motor_g = atof(param_value);
   fprintf(stdout, "   v_motor_g = %g\n", parameters -> v_motor_g);
}
else if (strcmp(param_name, "motor_turning_frequency") == 0){
   parameters -> f_turning = atof(param_value);
   fprintf(stdout, "   f_turning = %g\n", parameters -> f_turning);
}
else if (strcmp(param_name, "c_motor") == 0){
   parameters -> c_motor = atof(param_value);
   fprintf(stdout, "   c_motor = %g\n", parameters -> c_motor);
}
else if (strcmp(param_name, "kon") == 0){
   parameters -> kon = atof(param_value);
   fprintf(stdout, "   kon = %g\n", parameters -> kon);
}
else if (strcmp(param_name, "koff") == 0){
   parameters -> koff = atof(param_value);
   fprintf(stdout, "   koff = %g\n", parameters -> koff);
}
else if (strcmp(param_name, "koff_end") == 0){
   parameters -> koff_end = atof(param_value);
   fprintf(stdout, "   koff_end = %g\n", parameters -> koff_end);
}
else if (strcmp(param_name, "pickup_time") == 0){
   parameters -> pickup_time = atoi(param_value);
   fprintf(stdout, "   pickup_time = %d\n", parameters -> pickup_time);
}
else if (strcmp(param_name, "right_conc_right") == 0){
   parameters -> right_conc_right = atof(param_value);
   fprintf(stdout, "   right_conc_right = %g\n", parameters -> right_conc_right);
}
else if (strcmp(param_name, "right_conc_left") == 0){
   parameters -> right_conc_left = atof(param_value);
   fprintf(stdout, "   right_conc_left = %g\n", parameters -> right_conc_left);
}
else if (strcmp(param_name, "left_conc_right") == 0){
   parameters -> left_conc_right = atof(param_value);
   fprintf(stdout, "   left_conc_right = %g\n", parameters -> left_conc_right);
}
else if (strcmp(param_name, "left_conc_left") == 0){
   parameters -> left_conc_left = atof(param_value);
   fprintf(stdout, "   left_conc_left = %g\n", parameters -> left_conc_left);
}
else if (strcmp(param_name, "flux_flag") == 0){
   parameters -> flux_flag = atoi(param_value);
   fprintf(stdout, "   flux_flag = %d\n", parameters -> flux_flag);
}
else {
   fprintf(stderr, "error parsing parameter file %s on line %d:\n", param_file, n_lines);
   fprintf(stderr, "%s", line);
   exit(1);
}
