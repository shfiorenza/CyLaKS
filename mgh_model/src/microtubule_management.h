#ifndef _MICROTUBULE_MANAGEMENT_H
#define _MICROTUBULE_MANAGEMENT_H

#include "microtubule.h" 	// Includes <vector> lib as well
#include <iostream>
#ifndef _PARAMETERS_H
typedef struct system_parameters system_parameters;
#endif
#ifndef _SYSTEM_PROPERTIES_H
typedef struct system_properties system_properties;
#endif

class Tubulin;

class MicrotubuleManagement{
	private:

	public:
		int n_unoccupied_;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;

		std::vector<Microtubule> active_mts;
		std::vector<Tubulin*> unoccupied_list;
	private:
	
	public:
		MicrotubuleManagement();

		void Initialize(system_parameters *parameters, 
						system_properties *properties);

		void GenerateActiveMTs();
		void GenerateUnoccupiedList();

		void UnoccupiedCheck(Tubulin *site);
		void UnoccupiedCheck(int i_mt, int i_site);
		void OccupiedCheck(Tubulin *site);
		void OccupiedCheck(int i_mt, int i_site);

		Tubulin* GetUnoccupiedSite();
		void RemoveFromUnoccupiedList(Tubulin* site);
		void AddToUnoccupiedList(Tubulin* site);
};
#endif
