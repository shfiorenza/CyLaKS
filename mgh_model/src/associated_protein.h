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
		int x_dist_; 			// in no. of sites; can only be 0 or pos. 
		int dist_cutoff_  = 7;		// maximum value x_dist_ can have
		int rest_dist_ = 0;			// x_dist at which spring extension is ~0
	
		double kbT_;			
        double site_size_;	
		double r_0_;		
        double k_spring_;
		double extension_ ;		// current extension of xlink in nm
		double cosine_;			// of angle of xlink w/ respect to MT

		bool tethered_ = false;

		Tubulin *site_one_ = nullptr;
		Tubulin *site_two_ = nullptr;
		Kinesin *motor_ = nullptr;

		std::vector<Tubulin*> neighbor_sites_;
		std::vector<double> binding_weight_lookup_;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
	private:

	public:
		AssociatedProtein();
		void Initialize(system_parameters *parameters, 
			system_properties *properties, int ID);
		void SetParameters();
		void PopulateBindingLookupTable();

		void UpdateNeighborSites();
		bool NeighborExists(int x_dist);

		void UpdateExtension();
		void ForceUnbind(int x_dist_pre); 

		int GetDirectionTowardRest(Tubulin *site);
		int SampleSpringExtension();
		double GetAnchorCoordinate();
		double GetBindingWeight(Tubulin *neighbor);
		double GetExtensionForce(Tubulin *site);
		Tubulin* GetActiveHeadSite();
		Tubulin* GetNeighborSite(int x_dist);
};
#endif
