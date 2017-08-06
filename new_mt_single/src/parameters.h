#ifndef _PARAMETERS_H
#define _PARAMETERS_H

typedef struct {

	long seed;				/* random number seed */

	int n_steps;			/* number of simulation steps per run */

	int data_threshold;	 	/* step after which data recording begins */
	
	int n_datapoints;		/* number of data points to sample during sim run */

	double delta_t; 		/* how much time passes in one time step (s) */

	int n_microtubules;     /* number of individual microtubules */

	int length_of_microtubule;	/*how many sites are on one MT */

	double k_on;			/* bulk binding rate (1/nM*s) */
	
	double c_motor;         /* bulk motor concentration (nanomolar or nM) */

	double k_off;			/* motor unbinding rate (1/s) */

	double motor_speed;		/* motor velocity (micrometers/s) */

	double switch_rate;		/* motor switching frequency (1/s) */
	
	double alpha;			/* motor flux INTO overlap region */

	double beta;			/* motor flux OUT OF overlap region */

	double p_mutant; 		/* probability of tubulin sites being mutant */

} system_parameters;
#endif
