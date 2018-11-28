#ifndef _KINESIN_H
#define _KINESIN_H
#include <vector>
class Microtubule;
class Tubulin;
class AssociatedProtein;
struct system_properties;
struct system_parameters;

class Kinesin{
	private:

	public:
		int ID_;				// Unique ID of this kinesin in resevoir
		int speciesID_ = 2;		// Unique ID of this species (kinesin)
		int heads_active_ = 0;
		int n_neighbor_sites_ = 0;
		int n_neighbor_xlinks_ = 0;
		
		// x_dist_dub is used to index the tether extension of motors, e.g.
		// an x_dist_dub of 10 means an extension of -80 nm (40 - 120)
		int x_dist_doubled_;		// in no. of sites 
		int dist_cutoff_;			// max value x_dist (not 2x) can be
		int comp_cutoff_;			// min value x_dist (not 2x) can be
		double rest_dist_;		// spring extension is ~0 for this

		double r_0_;			
		double k_spring_;
		double k_slack_;
		double extension_;		// in nm
		double cosine_;			// of motor tether angle w.r.t. horizontal
	
		bool tethered_ = false;

		Microtubule *mt_ = nullptr; 		
		Tubulin *front_site_ = nullptr;
		Tubulin *rear_site_ = nullptr;
		AssociatedProtein *xlink_ = nullptr; 

		// Neighbor xlinks are for when the motor is bound but untethered
		std::vector<AssociatedProtein*> neighbor_xlinks_;
		// Neighbor sites are for when the motor is tethered but unbound
		std::vector<Tubulin*> neighbor_sites_;
		// Indices for lookup tables correspond to distance
		// in (no. of sites)/2, NOT extension of tether 
		// e.g. ...lookup_[1] is an x-dist of 1/2 of a site
		std::vector<double> tethering_weight_lookup_;
		std::vector<double> binding_weight_lookup_;
		std::vector<double> extension_lookup_;
		std::vector<double> cosine_lookup_;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
	private:

	public:
		Kinesin();
		void Initialize(system_parameters *parameters, 
			system_properties *properties, int ID);
		void SetParameters();
		void CalculateCutoffs();
		void InitiateNeighborLists();
		void PopulateTetheringLookupTable();
		void PopulateBindingLookupTable();
		void PopulateExtensionLookups();

		void UpdateNeighborXlinks();
		void UpdateNeighborSites();

		void UpdateExtension();
		void ForceUntether(int x_dub_pre);

		void UntetherSatellite();

		bool AtCutoff();

		int GetDirectionTowardRest();
		double GetRestLengthCoordinate(); 	// coord where ext ~ 0 when bound
		double GetStalkCoordinate(); // tail originates from stalk
		double GetTetheringWeight(AssociatedProtein *xlink);
		double GetBindingWeight(Tubulin *site);
		double GetTetherForce(Tubulin *site);
		Tubulin* GetActiveHeadSite();
		Tubulin* GetSiteCloserToRest();
		Tubulin* GetSiteFartherFromRest();
		Tubulin* GetWeightedNeighborSite();
		AssociatedProtein* GetWeightedNeighborXlink();
};
#endif
