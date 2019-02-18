#ifndef _KINESIN_PARAMETERS_H
#define _KINESIN_PARAMETERS_H

struct kinesin_parameters{

	double k_on;			// Binding rate of ADP heads to MT; 1/(nM*s)
	double c_bulk;			// Bulk concentration of motors; nM
	double c_eff_bind;		// For 2nd ADP head when 1st is bound; nM
	double k_on_ATP; 		// ATP binding rate to mots; 1/(micromolar*s) 
	double c_ATP;			// Bulk concentration of ATP; micromolar
	double k_hydrolyze; 	// Rate that ATP->ADPP occurs in motors; 1/s
	double k_off_i;			// Unbinding rate for ADPP-bound heads; 1/s
	double k_off_ii;
	double k_tether;		// Tethering rate; 1/(nM*s)
	double c_eff_tether;	// Effective concentration of free_teth motors
	double k_untether;		// Untethering rate when bound; 1/s
	double r_0;				// Rest length of stalk (or tether); nm
	double k_spring;		// Spring constant of tether; pN/nm
	double k_slack;			// '' but when shorter than rest length
	double stall_force;		// Force at which the motor can't step; pN
	double applied_force;   // Perpetually applied force on motors; pN
	bool tethers_active;
	bool endpausing_active;
};
#endif
