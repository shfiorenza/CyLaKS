#ifndef _KINESIN_H
#define _KINESIN_H
#include <vector>
#include <string>
class Microtubule;
class Tubulin;
class AssociatedProtein;
struct system_properties;
struct system_parameters;

class Kinesin{
	private:
		// Indices for lookup tables correspond to distance
		// in (no. of sites)/2, NOT extension of tether 
		// e.g. ...lookup_[1] is an x-dist of 1/2 of a site
		std::vector<double> cosine_lookup_;
		std::vector<double> extension_lookup_;
		// "creation" refers to creating a state with some extension
		std::vector<double> creation_weight_lookup_;
		std::vector<double> p_step_to_rest_;
		std::vector<double> p_step_fr_rest_;
		// Neighbor sites are for when the motor is tethered but unbound
		std::vector<Tubulin*> neighbor_sites_;
		// Neighbor xlinks are for when the motor is bound but untethered
		std::vector<AssociatedProtein*> neighbor_xlinks_;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
		
	public:
		struct head{
			Kinesin *motor_; 
			Tubulin *site_; 
			bool trailing_; 
			std::string ligand_; 
		};

		int ID_;				// Unique ID of this kinesin in resevoir
		int speciesID_ = 2;		// Unique ID of this species (kinesin)
		int active_index_; 		// index of this motor in active_ list
		int heads_active_ = 0;
		int n_neighbor_sites_ = 0;
		int n_neighbor_xlinks_ = 0;
		
		// x_dist_dub is used to index the tether extension of motors, e.g.
		// an x_dist_dub of 10 means an extension of -80 nm (40 - 120)
		int x_dist_doubled_;		// in no. of sites 
		int dist_cutoff_;			// max value x_dist (not 2x) can be
		int comp_cutoff_;			// min value x_dist (not 2x) can be
		double rest_dist_;			// spring extension is ~0 for this

		double r_0_;			// in nm
		double k_spring_;		// in pN / nm
		double k_slack_;		// in pN / nm
		double extension_;		// in nm
		double cosine_;			// of motor tether angle w.r.t. horizontal

		bool frustrated_ = false;
		bool tethered_ = false;

		head head_one_ = {this, nullptr, false, "ADP"},
			 head_two_ = {this, nullptr, true, "ADP"};

		Microtubule *mt_ = nullptr; 		
		AssociatedProtein *xlink_ = nullptr; 

	private:
		void SetParameters();
		void CalculateCutoffs();
		void InitializeNeighborLists();
		void InitializeExtensionLookup();
		void InitializeCreationWeightLookup();
		void InitializeSteppingProbabilities();

	public:
		Kinesin();
		void Initialize(system_parameters *parameters, 
				system_properties *properties, int ID);
		
		head* GetActiveHead();
		head* GetDockedHead();
		double GetStalkCoordinate(); 	  // tail originates from stalk
		double GetDockedCoordinate();
		void ChangeConformation();

		bool AtCutoff();
		void UpdateNeighborSites();
		void UpdateNeighborXlinks();
		void UpdateExtension();
		void ForceUntether();
		void UntetherSatellite();

		int GetDirectionTowardRest();
		double GetRestLengthCoordinate(); // coord where ext ~ 0 when bound
		double GetTetherForce(Tubulin *site);
		double GetTotalBindingWeight();
		double GetBindingWeight(Tubulin *site);
		double GetTotalTetheringWeight();
		double GetTetheringWeight(AssociatedProtein *xlink);
		Tubulin* GetWeightedNeighborSite();
		AssociatedProtein* GetWeightedNeighborXlink();
};
#endif
