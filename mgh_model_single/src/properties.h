#ifndef _SYSTEM_PROPERTIES
#define _SYSTEM_PROPERTIES
#include "curator.h"
#include "rng_management.h"
#include "microtubule_management.h"
#include "kinesin_management.h"
#include "associated_protein_management.h"

typedef struct system_properties{

	Curator wallace;

	RandomNumberManagement gsl; 

	MicrotubuleManagement microtubules;

	KinesinManagement kinesin4; 

	AssociatedProteinManagement prc1; 
	
	int current_step_;
	
	int n_binds_ = 0;
	double p_bind_cum_= 0;
	int n_unbinds_ = 0;
	double p_unbind_cum_ = 0;

} system_properties;
#endif

