#ifndef _ASSOCIATED_PROTEIN_H
#define _ASSOCIATED_PROTEIN_H
#include <vector>
class Tubulin;
class Kinesin;
class Microtubule;
struct system_properties;
struct system_parameters;

class AssociatedProtein{
	private:

	public:
		int ID_;
		int speciesID_ = 1;
		int heads_active_ = 0;
		int n_neighbor_sites_ = 0;

		// x_dist_ is used to index the xlink extensions for lookup
		// e.g. x_dist_ = 0 means an extension of 3 nm (35 - 32)
		int x_dist_ = 0; 		// in no. of sites; can only be 0 or pos. 
		int rest_dist_ = 0;		// x_dist at which spring extension is ~0
		int dist_cutoff_ = 6;	// maximum value x_dist_ can have
	
		double kbT_ = 4.114;		// in pN * nm			(sci adv)
        double site_size_ = 8;		// size of tubulin dimer in nm
		double r_0_ = 32;			// spring r_0 in nm  	(radhika)
        double k_spring_ = 0.207; 	// in pN / nm			(sci adv)
		double extension_ = 0;		// in nm; neg. means compression 

		bool tethered_ = false;

		Tubulin *site_one_ = nullptr;
		Tubulin *site_two_ = nullptr;
		Kinesin *motor_ = nullptr;

		std::vector<Tubulin*> neighbor_sites_;
		std::vector<double> binding_weight_lookup_;
		std::vector<double> spring_force_lookup_;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
	private:

	public:
		AssociatedProtein();
		void Initialize(system_parameters *parameters, 
			system_properties *properties, int ID);
		void PopulateBindingLookupTable();

		void UpdateNeighborSites();
		bool NeighborExists(int x_dist);

		void UpdateExtension();
		void ForceUnbind(); 

		double GetAnchorCoordinate();
		double GetBindingWeight(Tubulin *neighbor);
		Tubulin* GetActiveHeadSite();
		Tubulin* GetNeighborSite(int x_dist);
};
#endif
