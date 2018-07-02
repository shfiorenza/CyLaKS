#ifndef _PROTOFILAMENT_H
#define _PROTOFILAMENT_H
#include "tubulin.h"
#include <vector>

struct system_parameters;
struct system_properties; 

class Protofilament{

	private:

	public:
		int index_; 		// Index of PF in MT's barrel array
		int polarity_;		// 0 if plus-end is at start of lattice_; 1 if not
		double z_shift_; 
		
		int n_sites_; 

		Tubulin *plus_end_;
		Tubulin *minus_end_;

		std::vector<Tubulin> lattice_;

		system_parameters *parameters_  = nullptr;
		system_properties *properties_ = nullptr; 
	private:

	public:
		Protofilament(); 
		void Initialize(system_parameters *parameters, 
						system_properties *properties, int i_pf); 

		void SetParameters();
		void GenerateLattice();

};
#endif
