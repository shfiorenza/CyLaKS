#ifndef _PARAMETERS_H
#define _PARAMETERS_H
/* The <system_parameters> structure contains thermodynamic, control, and
   configuration parameters for the simulation. */

typedef struct {

	int n_steps;			/* number of simulation steps per run */

	int data_threshold;	 	/* step after which data recording begins */

	long seed;				/* random number seed */

	int n_microtubules;     /* number of individual microtubules */

	int length_of_microtubule;	/*how many sites are on one MT */

	double delta_t; 		/* how much time passes in one time step (s) */

	double k_on;			/* bulk binding rate (1/nM*s) */
	
	double c_motor;         /* bulk motor concentration (nanomolar or nM) */

	double k_off;			/* motor unbinding rate (1/s) */

	double motor_speed;		/* motor velocity (micrometer/s) */

	double switch_rate;		/* motor switching frequency (1/s) */
	
	double alpha;			/* motor flux INTO overlap region */

	double beta;			/* motor flux OUT OF overlap region */


} system_parameters;
#endif
