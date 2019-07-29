#pragma once
#include <string>
#include <vector>
class Tubulin;
class Kinesin;
class Microtubule;
struct system_properties;
struct system_parameters;

class AssociatedProtein{
	private:
		template<typename DATA_T>
		using Vec = std::vector<DATA_T>;
		// Indices for lookup tables correspond to x_dist
		Vec<double> cosine_lookup_; 
		Vec<double> extension_lookup_;
		// 1st index is n_neighbs (as in PRC1 neighbs); 2nd is x_dist
		Vec<Vec<double>> weights_bind_ii_;
		Vec<Vec<double>> weights_bind_i_teth_;
		Vec<Vec<Vec<double>>> weights_bind_ii_to_teth_;
		Vec<Vec<Vec<double>>> weights_bind_ii_fr_teth_;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;

	public:
		// Neighbor lists; not to be confused with no. of PRC1 neighbs
		Vec<Tubulin*> neighbor_sites_;
		Vec<Tubulin*> teth_neighbor_sites_; 
		Vec<Tubulin*> teth_neighbor_sites_ii_; 
		struct Monomer{
			AssociatedProtein *xlink_;
			Tubulin *site_;
			std::string state_; 
			int GetPRC1NeighbCount(){
				return xlink_->GetPRC1NeighbCount(this);
			};
			int GetDirectionToRest(){
				return xlink_->GetDirectionToRest(this->site_);
			};
			Monomer* GetOtherHead(){
				if(this == &xlink_->head_one_)
					return &xlink_->head_two_;
				else
					return &xlink_->head_one_;
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

		int max_neighbs_ = 0;

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
		void InitializeNeighborLists();
		void InitializeLookupTables(); 

	public:
		AssociatedProtein();
		void Initialize(system_parameters *parameters, 
			system_properties *properties, int ID);

		Monomer* GetActiveHead();
		Tubulin* GetActiveHeadSite();
		double GetAnchorCoordinate();
		int GetPRC1NeighbCount(Monomer* head);

		void UpdateNeighborSites_II();
		void UpdateNeighborSites_I_Teth();
		void UpdateNeighborSites_II_Teth();
		void UpdateExtension();
		void ForceUnbind(int x_dist_pre); 
		void UntetherSatellite(); 

		int GetDirectionToRest(Tubulin *site);
		double GetExtensionForce(Tubulin *site);
		double GetBindingWeight_II(Tubulin *neighbor);
		double GetBindingWeight_I_Teth(Tubulin *neighbor);
		double GetBindingWeight_II_Teth(Tubulin *neighbor);
		Tubulin* GetSiteCloserToTethRest();
		Tubulin* GetSiteFartherFromTethRest();
		Tubulin* GetWeightedSite_Bind_II();
		Tubulin* GetWeightedSite_Bind_I_Teth();
		Tubulin* GetWeightedSite_Bind_II_Teth();
};
