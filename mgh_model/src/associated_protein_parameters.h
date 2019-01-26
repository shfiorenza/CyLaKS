#ifndef _ASSOCIATED_PROTEIN_PARAMETERS_H
#define _ASSOCIATED_PROTEIN_PARAMETERS_H

struct associated_protein_parameters{

	double k_on;				// Bulk binding rate;  1/(nM*s)
	double concentration;		// Bulk concentration; nM
	double conc_eff_bind;		// Effective concentration of second head
	double k_off_i;				// Unbinding rate for singly-bound;  1/s
	double k_off_ii;			// Unbinding rate for doubly-bound;  1/s
								//   when first head is bound; nM
	double diffusion_const_i;	// Diffusion const for singly-bound; um^2/s
	double diffusion_const_ii;	// Diffusion const for doubly-bound; um^2/s
	double r_0;					// Rest length of coiled-coil domain; nm
	double k_spring; 			// Spring constant of CC-domain; pN/nm

};
#endif
