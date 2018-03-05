if (strcmp(param_name, "seed") == 0) {
	parameters -> seed = atol(param_value);
	fprintf(stdout, "   seed = %ld\n", parameters -> seed);
}
else if (strcmp(param_name, "n_steps") == 0) {
	parameters -> n_steps = atoi(param_value);
	fprintf(stdout, "   n_steps = %i\n", parameters -> n_steps);
}
else if (strcmp(param_name, "n_datapoints") == 0) {
	parameters -> n_datapoints = atoi(param_value);
	fprintf(stdout, "   n_datapoints = %i\n", parameters -> n_datapoints);
}
else if (strcmp(param_name, "data_threshold") == 0) {
	parameters -> data_threshold = atoi(param_value);
	fprintf(stdout, "   data_threshold = %i\n", parameters -> data_threshold);
}
else if (strcmp(param_name, "delta_t") == 0){
	parameters -> delta_t = atof(param_value);
	fprintf(stdout, "   delta_t = %g\n", parameters -> delta_t);
}
else if (strcmp(param_name, "kbT") == 0){
	parameters -> kbT = atof(param_value);
	fprintf(stdout, "   kbT = %g\n", parameters -> kbT);
}
else if (strcmp(param_name, "eta_inverse") == 0){
	parameters -> eta_inverse = atof(param_value);
	fprintf(stdout, "   eta_inverse = %g\n", parameters -> eta_inverse);
}
else if (strcmp(param_name, "n_microtubules") == 0){
	parameters -> n_microtubules = atoi(param_value);
	fprintf(stdout, "   n_microtubules = %i\n", parameters -> n_microtubules);
}
else if (strcmp(param_name, "length_of_microtubule") == 0){
	parameters -> length_of_microtubule = atoi(param_value);
	fprintf(stdout, "   length_of_microtubule = %i\n", parameters -> length_of_microtubule);
}
else if (strcmp(param_name, "mt_radius") == 0){
	parameters -> mt_radius = atof(param_value);
	fprintf(stdout, "   mt_radius = %g\n", parameters -> mt_radius);
}
else if (strcmp(param_name, "mt_height") == 0){
	parameters -> mt_height = atof(param_value);
	fprintf(stdout, "   mt_height = %g\n", parameters -> mt_height);
}
else if (strcmp(param_name, "site_size") == 0){
	parameters -> site_size = atof(param_value);
	fprintf(stdout, "   site_size = %g\n", parameters -> site_size);
}
else if (strcmp(param_name, "k_on_xlink") == 0){
	parameters -> k_on_xlink = atof(param_value);
	fprintf(stdout, "   k_on_xlink = %g\n", parameters -> k_on_xlink);
}
else if (strcmp(param_name, "c_xlink") == 0){
	parameters -> c_xlink = atof(param_value);
	fprintf(stdout, "   c_xlink = %g\n", parameters -> c_xlink);
}
else if (strcmp(param_name, "k_off_xlink_i") == 0){
	parameters -> k_off_xlink_i = atof(param_value);
	fprintf(stdout, "   k_off_xlink_i = %g\n", parameters -> k_off_xlink_i);
}
else if (strcmp(param_name, "k_off_xlink_ii") == 0){
	parameters -> k_off_xlink_ii = atof(param_value);
	fprintf(stdout, "   k_off_xlink_ii = %g\n", parameters -> k_off_xlink_ii);
}
else if (strcmp(param_name, "c_eff_xlink") == 0){
	parameters -> c_eff_xlink = atof(param_value);
	fprintf(stdout, "   c_eff_xlink = %g\n", parameters -> c_eff_xlink);
}
else if (strcmp(param_name, "D_xlink_i") == 0){
	parameters -> D_xlink_i = atof(param_value);
	fprintf(stdout, "   D_xlink_i = %g\n", parameters -> D_xlink_i);
}
else if (strcmp(param_name, "D_xlink_ii") == 0){
	parameters -> D_xlink_ii = atof(param_value);
	fprintf(stdout, "   D_xlink_ii = %g\n", parameters -> D_xlink_ii);
}
else if (strcmp(param_name, "r_0_xlink") == 0){
	parameters -> r_0_xlink = atof(param_value);
	fprintf(stdout, "   r_0_xlink = %g\n", parameters -> r_0_xlink);
}
else if (strcmp(param_name, "k_spring_xlink") == 0){
	parameters -> k_spring_xlink = atof(param_value);
	fprintf(stdout, "   k_spring_xlink = %g\n", parameters -> k_spring_xlink);
}
else if (strcmp(param_name, "k_on_motor") == 0){
	parameters -> k_on_motor = atof(param_value);
	fprintf(stdout, "   k_on_motor = %g\n", parameters -> k_on_motor);
}
else if (strcmp(param_name, "c_motor") == 0){
	parameters -> c_motor = atof(param_value);
	fprintf(stdout, "   c_motor = %g\n", parameters -> c_motor);
}
else if (strcmp(param_name, "c_eff_motor_bind") == 0){
	parameters -> c_eff_motor_bind = atof(param_value);
	fprintf(stdout, "   c_eff_motor_bind = %g\n", 
			parameters -> c_eff_motor_bind);
}
else if (strcmp(param_name, "k_off_motor") == 0){
	parameters -> k_off_motor = atof(param_value);
	fprintf(stdout, "   k_off_motor = %g\n", parameters -> k_off_motor);
}
else if (strcmp(param_name, "k_off_pseudo") == 0){
	parameters -> k_off_pseudo = atof(param_value);
	fprintf(stdout, "   k_off_pseudo = %g\n", parameters -> k_off_pseudo);
}
else if (strcmp(param_name, "k_off_ratio") == 0){
	parameters -> k_off_ratio = atof(param_value);
	fprintf(stdout, "   k_off_ratio = %g\n", parameters -> k_off_ratio);
}
else if (strcmp(param_name, "motor_speed") == 0){
	parameters -> motor_speed = atof(param_value);
	fprintf(stdout, "   motor_speed  = %g\n", parameters -> motor_speed);
}
else if (strcmp(param_name, "failstep_rate") == 0){
	parameters -> failstep_rate = atof(param_value);
	fprintf(stdout, "   failstep_rate = %g\n", parameters -> failstep_rate);
}
else if (strcmp(param_name, "D_motor") == 0){
	parameters -> D_motor = atof(param_value);
	fprintf(stdout, "   D_motor = %g\n", parameters -> D_motor);
}
else if (strcmp(param_name, "k_tether_free") == 0){
	parameters -> k_tether_free = atof(param_value);
	fprintf(stdout, "   k_tether_free = %g\n", parameters -> k_tether_free);
}
else if (strcmp(param_name, "c_eff_motor_teth") == 0){
	parameters -> c_eff_motor_teth = atof(param_value);
	fprintf(stdout, "   c_eff_motor_teth = %g\n", 
			parameters -> c_eff_motor_teth);
}
else if (strcmp(param_name, "k_untether") == 0){
	parameters -> k_untether = atof(param_value);
	fprintf(stdout, "   k_untether = %g\n", parameters->k_untether);
}
else if (strcmp(param_name, "k_untether_free") == 0){
	parameters -> k_untether_free = atof(param_value);
	fprintf(stdout, "   k_untether_free = %g\n", parameters->k_untether_free);
}
else if (strcmp(param_name, "r_0_motor") == 0){
	parameters -> r_0_motor = atof(param_value);
	fprintf(stdout, "   r_0_motor = %g\n", parameters -> r_0_motor);
}
else if (strcmp(param_name, "k_spring_motor") == 0){
	parameters -> k_spring_motor = atof(param_value);
	fprintf(stdout, "   k_spring_motor = %g\n", parameters -> k_spring_motor);
}
else if (strcmp(param_name, "k_slack_motor") == 0){
	parameters -> k_slack_motor = atof(param_value);
	fprintf(stdout, "   k_slack_motor = %g\n", parameters -> k_slack_motor);
}
else if (strcmp(param_name, "stall_force") == 0){
	parameters -> stall_force = atof(param_value);
	fprintf(stdout, "   stall_force = %g\n", parameters -> stall_force);
}
else if (strcmp(param_name, "switch_rate") == 0){
	parameters -> switch_rate = atof(param_value);
	fprintf(stdout, "   switch_rate = %g\n", parameters -> switch_rate);
}
else if (strcmp(param_name, "alpha") == 0){
	parameters -> alpha = atof(param_value);
	fprintf(stdout, "   alpha = %g\n", parameters -> alpha);
}
else if (strcmp(param_name, "beta") == 0){
	parameters -> beta = atof(param_value);
	fprintf(stdout, "   beta = %g\n", parameters -> beta);
}
