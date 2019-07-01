#ifndef _ASSOCIATED_PROTEIN_MANAGEMENT_H
#define _ASSOCIATED_PROTEIN_MANAGEMENT_H
#include "associated_protein.h"
#include <functional>
#include <string>

template<class T>
using vec = std::vector<T>;

struct system_parameters;
struct system_properties;
class AssociatedProteinManagement{
	private:
		// Structure that holds all pertinent info for a given MC event:
		struct event{
			// Initialization routine
			event(int i, int code, std::string lab, std::string tar, 
					std::function<int(double, int)> p_dist, 
					int *pop, double p): index_(i), kmc_code_(code), 
					label_(lab), target_pop_(tar), prob_dist_(p_dist), 
					pop_ptr_(pop), p_occur_(p) {}
			void SampleStatistics(){
				if(*pop_ptr_ > 0) 
					n_expected_ = prob_dist_(p_occur_, *pop_ptr_);
				else n_expected_ = 0;
			}
		private:
		public:
			std::function<int(double, int)> prob_dist_;
			int index_ = -1; 			// Serialized unique index of event
			int kmc_code_ = -1; 
			std::string label_ = "BRUH";
			std::string target_pop_ = "not me"; 
			double p_occur_ = 0;
			int *pop_ptr_ = nullptr; 	// Pointer to population size variable
			int n_expected_ = 0;
		};
		vec<event> events_;
		vec<vec<int> > events_by_pop_;
		vec<int> kmc_list_;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;

	public:
		int dist_cutoff_;		// see assoc. protein header
		int rest_dist_;			// see assoc. protein header
		int teth_cutoff_; 		// see kinesin header (dist_cutoff_ there)
		int max_neighbs_;
		double interaction_energy_;		// in kBT

		/* Population size trackers */ 
		int n_xlinks_ = 0; 		  // Total number of xlink objects created
		int n_active_ = 0; 		  // No. actively bound; dynamically updated
		int n_free_tethered_ = 0;
		int n_bound_unteth_ = 0;			// needed?
		// First index is number of PRC1 neighbors: [0], [1], or [2]
		// The last entry, [3], includes ALL regardless of n_neighbs
		vec<int> n_bound_i_;
		// Second index is [x] or [x_dub] for base or teth pops.
		vec<vec<int>> n_bound_i_teth_;
		vec<vec<int>> n_sites_ii_;
		// Second index is [x_dub], third & final index is [x]
		vec<vec<vec<int>>> n_sites_ii_teth_same_; 
		vec<vec<vec<int>>> n_sites_ii_teth_oppo_; 

		/* Probabilities of encoded events */
		double p_bind_i_teth_base_;
		double p_bind_ii_base_;
		double p_tether_free_; 
		double p_untether_free_;
		// First index is number of PRC1 neighbors: [0], [1], or [2]
		vec<double> p_bind_i_;			
		vec<double> p_unbind_i_;		
		vec<double> p_diffuse_i_fwd_;
		vec<double> p_diffuse_i_bck_;
		// Second index is [x] or [x_dub] for base or teth pops.
		vec<vec<double>> p_unbind_i_teth_;	
		vec<vec<double>> p_unbind_ii_;
		vec<vec<double>> p_diffuse_i_to_teth_rest_;
		vec<vec<double>> p_diffuse_i_fr_teth_rest_;
		vec<vec<double>> p_diffuse_ii_to_rest_;
		vec<vec<double>> p_diffuse_ii_fr_rest_;
		// Second index is [x_dub]; third & final index is [x]
		vec<vec<vec<double>>> p_unbind_ii_to_teth_;
		vec<vec<vec<double>>> p_unbind_ii_fr_teth_;
		vec<vec<vec<double>>> p_diffuse_ii_to_both_rest_;
		vec<vec<vec<double>>> p_diffuse_ii_fr_both_rest_;
		vec<vec<vec<double>>> p_diffuse_ii_to_self_fr_teth_;
		vec<vec<vec<double>>> p_diffuse_ii_fr_self_to_teth_;

		/* Lists that track different population types */
		vec<AssociatedProtein> xlinks_;				// Actual xlink objects
		vec<AssociatedProtein*> active_;
		vec<AssociatedProtein*> free_tethered_;
		vec<AssociatedProtein*> bound_untethered_;
		// First index is number of PRC1 neighbors: [0], [1], or [2]
		// The last entry, [3], includes ALL regardless of n_neighbs
		// Second index is actual xlink entry
		vec<vec<AssociatedProtein*>> bound_i_;
		// Second index is [x] or [x_dub]; third index is xlink entry
		vec<vec<vec<Tubulin*>>> sites_ii_;
		vec<vec<vec<AssociatedProtein*>>> bound_i_teth_;
		// Second index is [x_dub]; third index is [x]; fourth is entry
		// e.g., [0][16][2][1] -> 2nd xlink w/ x_dub=16, x=2, & 0 neighbs
		vec<vec<vec<vec<Tubulin*>>>> sites_ii_teth_oppo_;
		vec<vec<vec<vec<Tubulin*>>>> sites_ii_teth_same_;	

	private:
		void GenerateXLinks();
		void SetParameters();
		void InitializeLists();
		void InitializeDiffusionEvents();
		void InitializeKMCEvents();

	public:
		AssociatedProteinManagement();
		void Initialize(system_parameters *parameters, 
						system_properties *properties);

		AssociatedProtein* GetFreeXlink();
		AssociatedProtein* GetBoundUntetheredXlink();
		double GetWeight_II(); 
		double GetWeight_I_Teth();
		double GetWeight_II_Teth();

		void Update_All();
		void Update_Bound_I();
		void Update_Bound_I_Teth();
		void Update_Bound_II_Sites();
		void Update_Bound_II_Teth_Sites();
		void Update_Bound_Unteth(); 
		void Update_Free_Teth();

		void GenerateEventList();

		void RunKMC();
		// Diffusion events
		void Diffuse_I(int n_neighbs, int dir);
		void Diffuse_II(int n_neighbs, int x, int dir);
		void Diffuse_I_Teth(int n_neighbs, int x_dub, int dir);
		void Diffuse_II_Teth_Same(int n_neighbs, int x_dub, int x, int dir);
		void Diffuse_II_Teth_Oppo(int n_neighbs, int x_dub, int x, int dir);
		// Kinematic events
		void Bind_I(int n_neighbs);
		void Bind_II();		//XXX add neighb coop in weights
		void Unbind_I(int n_neighbs);
		void Unbind_II(int n_neighbs, int x);
		void Bind_I_Teth(); //XXX add neighb coop in weights
		void Bind_II_Teth(); //XXX add neighb coop in weights
		void Unbind_I_Teth(int n_neighbs, int x_dub); 
		void Unbind_II_To_Teth(int n_neighbs, int x_dub, int x);
		void Unbind_II_From_Teth(int n_neighbs, int x_dub, int x);
		void Tether_Free(); 
		void Untether_Free();
};
#endif
