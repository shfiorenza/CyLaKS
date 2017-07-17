if (strcmp(param_name, "n_steps") == 0) {
	parameters -> n_steps = atoi(param_value);
	fprintf(stdout, "   n_steps = %d\n", parameters -> n_steps);
}
if (strcmp(param_name, "data_threshold") == 0) {
	parameters -> data_threshold = atoi(param_value);
	fprintf(stdout, "   data_threshold = %d\n", parameters -> data_threshold);
}
if (strcmp(param_name, "seed") == 0) {
	parameters -> seed = atol(param_value);
	fprintf(stdout, "   seed = %ld\n", parameters -> seed);
}
if (strcmp(param_name, "delta_t") == 0){
	parameters -> delta_t = atof(param_value);
	fprintf(stdout, "   delta_t = %g\n", parameters -> delta_t);
}
if (strcmp(param_name, "k_on") == 0){
	parameters -> k_on = atof(param_value);
	fprintf(stdout, "   k_on = %g\n", parameters -> k_on);
}
if (strcmp(param_name, "c_motor") == 0){
	parameters -> c_motor = atoi(param_value);
	fprintf(stdout, "   c_motor = %g\n", parameters -> c_motor);
}
if (strcmp(param_name, "k_off") == 0){
	parameters -> k_off = atof(param_value);
	fprintf(stdout, "   k_off = %g\n", parameters -> k_off);
}
if (strcmp(param_name, "motor_speed") == 0){
	parameters -> motor_speed = atof(param_value);
	fprintf(stdout, "   motor_speed  = %g\n", parameters -> motor_speed);
}
if (strcmp(param_name, "switch_rate") == 0){
	parameters -> switch_rate = atof(param_value);
	fprintf(stdout, "   switch_rate = %g\n", parameters -> switch_rate);
}
if (strcmp(param_name, "alpha") == 0){
	parameters -> alpha = atof(param_value);
	fprintf(stdout, "   alpha = %g\n", parameters -> alpha);
}
if (strcmp(param_name, "beta") == 0){
	parameters -> beta = atof(param_value);
	fprintf(stdout, "   beta = %g\n", parameters -> beta);
}

