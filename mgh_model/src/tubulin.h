#ifndef _TUBULIN_H
#define _TUBULIN_H
#include "kinesin.h"

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
		Kinesin::head *motor_head_ = nullptr; 
		AssociatedProtein *xlink_ = nullptr;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
	private:

	public:
		Tubulin();
		void Initialize(system_parameters *parameters, 
			system_properties *properties, Microtubule *mt, int i_site);

		// 'equilibrium' refers to that of the crosslinker itself 
		// and the kinesin 4 tether that attaches to it; this checks 
		// if both are in the same direction w.r.t. this site 
		bool EquilibriumInSameDirection();
};
#endif
