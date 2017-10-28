#ifndef _ASSOCIATED_PROTEIN_H
#define _ASSOCIATED_PROTEIN_H
#include <vector>

#ifndef _PARAMETERS_H
typedef struct system_parameters system_parameters;
#endif
#ifndef _SYSTEM_PROPETIES_H
typedef struct system_properties system_properties;
#endif

class Microtubule;
class Tubulin;
class Kinesin;

class AssociatedProtein{
	private:

	public:
		int ID_;
		int speciesID_ = 1;
		int heads_active_ = 0;
		int n_neighbors_binding_ = 0;
		int n_neighbors_tethering_ = 0;
		
		bool bound_ = false;
		bool tethered_ = false;

		Microtubule *mt_one_ = nullptr;
		Microtubule *mt_two_ = nullptr;
		Tubulin *site_one_ = nullptr;
		Tubulin *site_two_ = nullptr;
		Kinesin *motor_ = nullptr;

		std::vector<Tubulin*> head_neighbor_list_;
		std::vector<Tubulin*> tether_neighbor_list_;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
	private:

	public:
		AssociatedProtein();
		void Initialize(system_parameters *parameters, 
			system_properties *properties, int ID);
		
		void UpdateNeighbors_Head();
		void UpdateNeighbors_Tether();
};
#endif
