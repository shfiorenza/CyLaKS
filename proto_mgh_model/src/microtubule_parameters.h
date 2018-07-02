#ifndef _MICROTUBULE_PARAMETERS_H
#define _MICROTUBULE_PARAMETERS_H
#include <vector>

struct microtubule_parameters{
	
	int count;				// Number of MTs in simulation
	int length;				// Length of each MT in sites (tubulin dimers)

	double site_size; 		// Length of tubulin dimer; nm
	double radius;			// Outer radius of MT barrel; nm
	double elevation;		// Distance above glass slide; nm
	std::vector<double> start_coord;  // Coords of each MT at sim start; nm
	std::vector<double> imposed_velocity;  // For each MT; nm/s 
	std::vector<double> immobile_until; // MTs immobilzed until this time; s

	bool printout; 			// Whether or not ASCII printout is enabled 


};
#endif
