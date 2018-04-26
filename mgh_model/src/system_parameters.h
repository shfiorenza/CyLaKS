#ifndef _PARAMETERS_H
#define _PARAMETERS_H

struct system_parameters{

	// Simulation properties
	long seed;				/* Random number seed */
	int n_steps;			/* Number of timesteps in simulation */
	int n_datapoints;		/* No. data points to sample during sim run */
	int data_threshold;	 	/* Step after which to begin data recording */
	double delta_t; 		/* How much time passes in one time step (s) */


	// Physical constants
	double kbT;				/* boltzmann constant * temp, in pN * nM */
	double eta_inverse;		/* Inverse viscosity; in [(pN/um^2)*s]^-1 */


	// Microtubule parameters
	int n_microtubules;     /* Number of individual microtubules */
	int length_of_microtubule;		/* Length of MT in number of sites */ 
	int top_mt_start_coord; /* starting coord of top (i = 1) microtubule */
	int bot_mt_start_coord; /* starting coord of bot (i = 0) microtubule */
	double top_mt_imposed_velocity; /* opposes normal sliding; nm/s */
	double mt_radius;		/* For 13 protofilaments; in nm */
	double mt_height;		/* Used in diffusion (dist. from wall); in nm */
	double site_size;		/* Length of tubulin dimer (1 site) in nm */


	// Crosslink parameters
	double k_on_xlink;		/* Xlink bulk binding rate (1/nM*s) */
	double c_xlink;			/* Bulk xlink concentration (nanomolar or nM) */
	double k_off_xlink_i;	/* Unbinding rate for single bound (1 / s) */
	double k_off_xlink_ii;  /* same for double bound (1 / s) */
	double c_eff_xlink;		/* Unitless; used to scale Q (partition func) */

	double D_xlink_i; 	    /* Diffusion const for single-bound; in um^2/s */
	double D_xlink_ii;		/* Diffusion const for double-bound; in um^2/s */

	double r_0_xlink;		/* Rest dist of xlink spring in nm */
	double k_spring_xlink;  /* Spring constant of xlink in pN / nm */


	// Motor parameters
	double k_on_motor;		/* Motor bulk binding rate (1/nM*s) */
	double c_motor;         /* Bulk motor concentration (nanomolar or nM) */
	double c_eff_motor_bind;	/* in nM, used for binding of 2nd head */
	double k_off_motor;		/* Motor unbinding rate (1/s) */
	double k_off_pseudo;	/* Pseudo-bound (1 head) motors unbinding rate */
	double k_off_ratio; 	/* Ratio of stepable to stalled motor unbinds*/
	double motor_speed;		/* Motor velocity (nm/s) */
	double failstep_rate; 	/* Stalled motors can turn pseudo-bound (1/s) */
	double D_motor;			/* Diffusion const for motor; in um^2/s */

	double k_tether_free;   /* Bulk tethering rate from solu. (1/nM*s) */
	double c_eff_motor_teth;	/* in nM, used in Q for motor tethering */
	double k_untether;		/* For tails unbind from xlinks (1/s) */
	double k_untether_free; /* Rate tethered free motors dissoc. (1/s) */

	double r_0_motor;		/* Rest dist of motor tether in nm */
	double k_spring_motor;	/* Spring const of tether in pN / nm */
	double k_slack_motor;	/* EFFECTIVE spring const for neg extension */
	double stall_force;		/* force at which motor stops stepping; in pN */

	double switch_rate;		/* Motor switching frequency (1/s) */
	double alpha;			/* Probability of placing motor on minus end*/
	double beta;			/* Probability of removing motor from plus end*/
};
#endif
