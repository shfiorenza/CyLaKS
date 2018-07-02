#ifndef _KINESIN_PARAMETERS_H
#define _KINESIN_PARAMETERS_H

struct kinesin_parameters{

	double k_on;			// Bulk binding rate; 1/(nM*s)
	double concentration;	// Bulk Concentration; nM
	double conc_eff_bind;	// Effective concentration of second head
							//   when first head is bound; nM
	double k_off_i;			// Unbinding rate for singly-bound; 1/s
	double k_off_ii;		// Unbinding rate for doubly-bound; 1/s
	double velocity;		// Motor stepping velocity; nm/s 
	double diffusion_const;	// Diffusion constant; um^2/s 
	double k_off_ratio;		// Ratio of mobile to stalled motor unbinds; NA
	double failstep_rate;	// Rate at which stalled motors 'failstep' and 
						  	//   go from doubly- to singly-bound 
	double k_tether_free;   // Bulk (not bound) tethering rate; 1/(nM*s)
	double conc_eff_tether;	// Effective concentration of first head
							//   when tethered but not bound; nM
	double k_untether_free; // Untethering rate when free; 1/s
	double k_untether;		// Untethering rate when bound; 1/s
	double r_0;				// Rest length of stalk (or tether); nm
	double k_spring;		// Spring constant of tether; pN/nm
	double k_slack;			// Spring constant of tether when shorter than
							//   rest length--similar to floppy rope; pN/nm
	double stall_force;		// Force at which the motor can no longer step
							//   and has a velocity of 0; pN 
	
/* Below are no longer used in sim -- may remove */
	double alpha;			// Probability of placing motor on minus-end
	double beta;			// Probability of removing motor from plus-end
};
#endif
