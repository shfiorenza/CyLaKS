#ifndef _KINESIN_H
#define _KINESIN_H
#include <vector>

#ifndef _PARAMETERS_H
typedef struct system_parameters system_parameters;
#endif
#ifndef _SYSTEM_PROPERTIES_H
typedef struct system_properties system_properties;
#endif

class Microtubule;
class Tubulin;
class AssociatedProtein;

class Kinesin{
	private:

	public:
		int ID_;						// Unique ID of this kinesin in resevoir
		int speciesID_ = 2;				// Unique ID of this species (kinesin)

		bool bound_ = false; 	
		bool tethered_ = false;

		Microtubule *mt_ = nullptr; 		
		Tubulin *front_site_ = nullptr; 
		Tubulin *rear_site_ = nullptr; 		
		AssociatedProtein *xlink_ = nullptr; 

//		std::vector<AssociatedProtein*> neighbor_xlinks;
//		std::vector<Tubulin*> neighbor_sites;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
	private:

	public:
		Kinesin();
		void Initialize(system_parameters *parameters, 
			system_properties *properties, int ID);

		void UpdateNeighbors();
		void ClearNeighbors();
};
#endif
