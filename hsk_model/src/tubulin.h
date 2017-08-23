#ifndef _TUBULIN_H
#define _TUBULIN_H

#ifndef _PARAMETERS_H
typedef struct system_parameters system_parameters;
#endif
#ifndef _SYSTEM_PROPERTIES_H
typedef struct system_properties system_properties;
#endif

class Microtubule;
class Kinesin;

class Tubulin{
	private:

	public:
		int index_;		// Index of tubulin site in MT lattice
		int coord_;		// Absolute coord of tubulin site
		int speciesID_ = 0;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;

		Microtubule *parent = nullptr;
		Kinesin *occupant = nullptr; 
	private:

	public:
		Tubulin();
		void Initialize(system_parameters *parameters, 
			system_properties *properties, Microtubule *mt, int i_site);
};
#endif
