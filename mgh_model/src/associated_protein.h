#ifndef _ASSOCIATED_PROTEIN_H
#define _ASSOCIATED_PROTEIN_H
#include <vector>
class Microtubule;
class Tubulin;
class Kinesin;
struct system_properties;
struct system_parameters;

class AssociatedProtein{
	private:

	public:
		int ID_;
		int speciesID_ = 1;
		int heads_active_ = 0;
		int n_neighbor_sites_ = 0;
		int n_neighbor_motors_ = 0;
	
		double kbT_ = 4.11e-21;     //FIXME   in Joules, from bobert
        double k_spring_ = 1.3e-5;  //FIXME   in N/m, from robert
		double r_0_ = 32e-9;			// rest length of prc1 in m
        double site_size_ = 8e-9;   //FIXME   in fuckin m wow 
		double dist_cutoff_ = 12;	//FIXME	  in no. of sites

		bool bound_ = false;
		bool tethered_ = false;

		Tubulin *site_one_ = nullptr;
		Tubulin *site_two_ = nullptr;
		Kinesin *motor_ = nullptr;

		std::vector<Tubulin*> neighbor_sites_;
		std::vector<double> binding_lookup_table_;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
	private:

	public:
		AssociatedProtein();
		void Initialize(system_parameters *parameters, 
			system_properties *properties, int ID);
		void PopulateBindingLookupTable();

		void UpdateNeighborSites();
		void UpdateNeighborMotors();

		Tubulin* GetActiveHeadSite();
		Tubulin* GetNeighborSite();
		double GetBindingWeight(Tubulin *neighbor);
};
#endif
