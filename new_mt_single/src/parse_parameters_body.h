if (strcmp(param_name, "seed") == 0) {
	parameters -> seed = atol(param_value);
	fprintf(stdout, "   seed = %ld\n\n", parameters -> seed);
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
	fprintf(stdout, "   delta_t = %g\n\n", parameters -> delta_t);
}
else if (strcmp(param_name, "n_microtubules") == 0){
	parameters -> n_microtubules = atof(param_value);
	fprintf(stdout, "   n_microtubules = %i\n", parameters -> n_microtubules);
}
else if (strcmp(param_name, "length_of_microtubule") == 0){
	parameters -> length_of_microtubule = atof(param_value);
	fprintf(stdout, "   length_of_microtubule = %i\n\n", parameters -> length_of_microtubule);
}
else if (strcmp(param_name, "k_on") == 0){
	parameters -> k_on = atof(param_value);
	fprintf(stdout, "   k_on = %g\n\n", parameters -> k_on);
}
else if (strcmp(param_name, "c_motor") == 0){
	parameters -> c_motor = atoi(param_value);
	fprintf(stdout, "   c_motor = %g\n", parameters -> c_motor);
}
else if (strcmp(param_name, "k_off") == 0){
	parameters -> k_off = atof(param_value);
	fprintf(stdout, "   k_off = %g\n", parameters -> k_off);
}
else if (strcmp(param_name, "motor_speed") == 0){
	parameters -> motor_speed = atof(param_value);
	fprintf(stdout, "   motor_speed  = %g\n", parameters -> motor_speed);
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
	fprintf(stdout, "   beta = %g\n\n", parameters -> beta);
}
else if (strcmp(param_name, "p_mutant") == 0){
	parameters -> p_mutant = atof(param_value);
	fprintf(stdout, "   p_mutant = %g\n", parameters -> p_mutant);
}
