#ifndef _TUBULIN_H
#define _TUBULIN_H

class Microtubule;
class Kinesin;
class AssociatedProtein;
struct system_parameters;
struct system_properties;

class Tubulin{
	private:

	public:
		int index_;		// Index of tubulin site in MT lattice
		int coord_;		// Absolute coord of tubulin site
		int speciesID_ = 0;

		bool occupied_ = false;

		Microtubule *mt_ = nullptr;
		Kinesin *motor_ = nullptr; 
		AssociatedProtein *xlink_ = nullptr;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
	private:

	public:
		Tubulin();
		void Initialize(system_parameters *parameters, 
			system_properties *properties, Microtubule *mt, int i_site);

		// 'springs' refer to prc1 xlinker as well as the kinesin 4
		// tether that attaches to it; the equilibrium position of these 
		// two springs are what we're checking to see if they're on the 
		// same side of the site or not 
		bool SpringEquilOnSameSide();
};
#endif
