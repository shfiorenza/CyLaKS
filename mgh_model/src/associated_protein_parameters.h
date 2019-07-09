#ifndef _ASSOCIATED_PROTEIN_PARAMETERS_H
#define _ASSOCIATED_PROTEIN_PARAMETERS_H

struct associated_protein_parameters{

	double k_on;				// Bulk binding rate;  1/(nM*s)
	double c_bulk;				// Bulk concentration; nM
	double c_eff_bind;			// Effective conc. of 2nd head binding; nM
	double k_off;				// Unbinding rate;  1/s
	double r_0;					// Rest length of coiled-coil domain; nm
	double k_spring; 			// Spring constant of CC-domain; pN/nm
	double diffusion_coeff;		// Diffusion coefficient; um^2/s
	double interaction_energy;  // Energy between two PRC1 neighbors
};
#endif
