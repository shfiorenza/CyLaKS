#ifndef _SYSTEM_PROPERTIES
#define _SYSTEM_PROPERTIES
#include "curator.h"
#include "rng_management.h"
#include "microtubule_management.h"
#include "kinesin_management.h"
#include "associated_protein_management.h"

struct system_properties{

	Curator wallace;
	RandomNumberManagement gsl; 
	MicrotubuleManagement microtubules;
	KinesinManagement kinesin4; 
	AssociatedProteinManagement prc1; 
	
	int current_step_;
	double t_motors_ = 0;
	double t_xlinks_dif_ = 0;
	double t_xlinks_kmc_ = 0;
	double t_MTs_ = 0;

	FILE *occupancy_file_, 
		 *motor_ID_file_, *xlink_ID_file_, 
		 *tether_coord_file_, *mt_coord_file_, 
		 *motor_extension_file_, *xlink_extension_file_, 
		 *motor_force_file_, *xlink_force_file_, *total_force_file_; 
	
};
#endif

