#ifndef _KINESIN_PARAMETERS_H
#define _KINESIN_PARAMETERS_H

struct kinesin_parameters{

	double k_on;			// Binding rate of ADP heads to MT; 1/(nM*s)
	double c_bulk;			// Bulk concentration of motors; nM
	double c_eff_bind;		// For 2nd ADP head when 1st is bound; nM
	double k_on_ATP; 		// ATP binding rate to mots; 1/(micromolar*s) 
	double c_ATP;			// Bulk concentration of ATP; micromolar
	double k_phosphorylate; // Rate that ATP->ADPP occurs in motors; 1/s
	double k_off_i;			// Unbinding rate for ADPP-bound heads; 1/s
	double k_off_ii;
	bool endpausing_active;

	bool tethers_active;
	double k_tether_free;   // Bulk (not bound) tethering rate; 1/(nM*s)
	double conc_eff_tether;	// Effective concentration of first head
							//   when tethered but not bound; nM
	double k_untether_free; // Untethering rate when free; 1/s
	double k_untether;		// Untethering rate when bound; 1/s
	double r_0;				// Rest length of stalk (or tether); nm
	double k_spring;		// Spring constant of tether; pN/nm
	double k_slack;			// Spring constant of tether when shorter than
							//  rest length--similar to floppy rope; pN/nm
	double stall_force;		// Force at which the motor can no longer step
							//  and has a velocity of 0; pN 
};
#endif
