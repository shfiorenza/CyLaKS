#ifndef _PARAMETERS_H
#define _PARAMETERS_H

struct system_parameters{

	long seed;				/* Random number seed */
	int n_steps;			/* Number of timesteps in simulation */
	int data_threshold;	 	/* Step after which data recording begins */
	int n_datapoints;		/* No. data points to sample during sim run */
	double delta_t; 		/* How much time passes in one time step (s) */
	int n_microtubules;     /* Number of individual microtubules */
	int length_of_microtubule;	/* Length of MT in number of sites */ 
	double k_on;			/* Mulk binding rate (1/nM*s) */
	double c_motor;         /* Mulk motor concentration (nanomolar or nM) */
	double k_off;			/* Motor unbinding rate (1/s) */
	double motor_speed;		/* Motor velocity (sites/s) */		// FIXME
	double switch_rate;		/* Motor switching frequency (1/s) */
	double alpha;			/* Probability of placing motor on minus end*/
	double beta;			/* Probability of removing motor from plus end*/

};
#endif
