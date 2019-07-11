#pragma once
#include <vector>
#include <string>
class Tubulin;
class Kinesin;
class Microtubule;
struct system_properties;
struct system_parameters;

class AssociatedProtein{
	private:
		// Indices for lookup tables correspond to x_dist
		std::vector<double> cosine_lookup_; 
		std::vector<double> extension_lookup_;
		std::vector<double> binding_weight_lookup_;
		// Indices correspond to x_dist_dub, i.e., 2x teth extension
		std::vector<double> teth_binding_weight_lookup_; 
		// 1st index is current x_dist; 2nd index is proposed x_dist_dub
		std::vector< std::vector<double> > teth_binding_weight_ii_to_; 
		std::vector< std::vector<double> > teth_binding_weight_ii_from_;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;

	public:
		// Neighbor lists; not to be confused with no. of PRC1 neighbs
		std::vector<Tubulin*> neighbor_sites_;
		std::vector<Tubulin*> teth_neighbor_sites_; 
		std::vector<Tubulin*> teth_neighbor_sites_ii_; 
		struct Monomer{
			AssociatedProtein *xlink_;
			Tubulin *site_;
			std::string state_; 
			int GetPRC1NeighbCount(){
				return xlink_->GetPRC1NeighbCount(this);
			};
		};

		// see kinesin header for description of variables
		int ID_;
		int speciesID_ = 1;
		int active_index_;
		int heads_active_ = 0;
		// For neighbor lists; not to be confused with no. of PRC1 neighbs
		int n_neighbor_sites_ = 0;
		int n_teth_neighbor_sites_ = 0;
		int n_teth_neighbor_sites_ii_ = 0;

		// x_dist_ is used to index the xlink extensions for lookup
		// e.g. x_dist_ = 0 means an extension of 3 nm (35 - 32)
		int x_dist_; 			// in no. of sites; can only be 0 or pos. 
		int dist_cutoff_;		// maximum value x_dist_ can have
		int rest_dist_;			// x_dist at which spring extension is ~0
	
		double r_0_;		
        double k_spring_;
		double extension_ ;		// current extension of xlink in nm
		double cosine_;			// of angle of xlink w/ respect to MT

		bool tethered_ = false;
		bool is_outdated_ = false;

		Monomer head_one_ = {this, nullptr, std::string("unbound")}, 
			 	head_two_ = {this, nullptr, std::string("unbound")};

		Kinesin *motor_ = nullptr;

	private:
		void SetParameters();
		void CalculateCutoffs();
		void InitiateNeighborLists();
		void PopulateBindingLookupTable();
		void PopulateTethBindingLookupTable();
		void PopulateTethBindingIILookupTable();
		void PopulateExtensionLookups();

	public:
		AssociatedProtein();
		void Initialize(system_parameters *parameters, 
			system_properties *properties, int ID);

		Monomer* GetActiveHead();
		Tubulin* GetActiveHeadSite();
		double GetAnchorCoordinate();
		int GetPRC1NeighbCount(Monomer* head);

		void UpdateNeighborSites();
		void UpdateTethNeighborSites();
		void UpdateTethNeighborSitesII();
		void UpdateExtension();
		void ForceUnbind(int x_dist_pre); 
		void UntetherSatellite(); 

		int GetDirectionTowardRest(Tubulin *site);
		double GetBindingWeight(Tubulin *neighbor);
		double GetExtensionForce(Tubulin *site);
		double GetTethBindingWeight(Tubulin *neighbor); 
		double GetTethBindingWeightII(Tubulin *neighbor);
		Tubulin* GetSiteCloserToTethRest();
		Tubulin* GetSiteFartherFromTethRest();
		Tubulin* GetWeightedNeighborSite();
		Tubulin* GetWeightedTethNeighborSite();
		Tubulin* GetWeightedTethNeighborSiteII();
};
