#ifndef _MICROTUBULE_PARAMETERS_H
#define _MICROTUBULE_PARAMETERS_H
#include <vector>

struct microtubule_parameters{
private:
	using vec_t = std::vector<double>;
public:
	int count;				// Number of MTs in simulation
	int length;				// Length of each MT in sites (tubulin dimers)
	double y_dist;			// Vertical distance between each MT

	double site_size; 		// Length of tubulin dimer; nm
	double radius;			// Outer radius of MT barrel; nm
	double elevation;		// Distance above glass slide; nm
	vec_t start_coord;  	// Coords of each MT at sim start; nm
	vec_t imposed_velocity; // For each MT; nm/s 
	vec_t immobile_until; 	// MTs immobilzed until this time; s

	bool printout_on; 			// Whether or not ASCII printout is enabled 
	bool diffusion_on;			// Whether or not MT diffusion is enabled

};
#endif
