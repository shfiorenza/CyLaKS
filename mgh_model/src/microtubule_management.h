#ifndef _MICROTUBULE_MANAGEMENT_H
#define _MICROTUBULE_MANAGEMENT_H

#include "microtubule.h" 	// Includes <vector> lib as well
//#include <iostream>
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
		int n_sites_tot_ = 0;
		int n_unoccupied_ = 0;
		int n_unoccupied_pairs_ = 0;
		int n_pair_entries_ = 0;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;

		std::vector<Microtubule> mt_list_;
		std::vector<Tubulin*> unoccupied_list_;
		std::vector<Tubulin*> unoccupied_pairs_list_;
	private:
	
	public:
		MicrotubuleManagement();

		void Initialize(system_parameters *parameters, 
						system_properties *properties);

		void SetParameters();		
		void GenerateMicrotubules();

		void UnoccupiedCheck(Tubulin *site);
		void UnoccupiedCheck(int i_mt, int i_site);
		void OccupiedCheck(Tubulin *site);
		void OccupiedCheck(int i_mt, int i_site);

		void UpdateNumUnoccupied();
		void UpdateUnoccupiedList();

		void UpdateNumUnoccupiedPairs();
		void UpdateUnoccupiedPairsList();

		Tubulin* GetUnoccupiedSite();
		Tubulin* GetUnoccupiedPair_1();
		Tubulin* GetUnoccupiedPair_2(Tubulin* first_site);

		void RemoveFromUnoccupiedList(Tubulin* site);
		void AddToUnoccupiedList(Tubulin* site);
};
#endif
