#ifndef _ASSOCIATED_PROTEIN_MANAGEMENT_H
#define _ASSOCIATED_PROTEIN_MANAGEMENT_H
#include "associated_protein.h"
#include <string>
#include <functional>

struct system_parameters;
struct system_properties;

class AssociatedProteinManagement{
	private:
		// Structure that holds all pertinent info for a given MC event:
		struct event{
			// Initialization routine
			event(int i, int code, std::string l, std::string tar, 
					std::function<int(double, int)> p_dist, 
					int *pop, double p): index_(i), kmc_code_(code), 
					label_(l), target_pop_(tar), prob_dist_(p_dist), 
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
		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;

	public:
		int n_xlinks_ = 0;
		// xlinks actively bound to some MT or motor; dynamically updated
		int n_active_ = 0;
		int n_free_tethered_ = 0;
		int n_bound_i_tot_ = 0; 
		int n_bound_i_tethered_tot_ = 0;		// needed?
		int n_bound_untethered_ = 0;
		// Indexed by number of PRC1 neighbors: [0], [1], or [2]
		std::vector<int> n_bound_i_;
		// Indexed by extension: [x_dub][x] or just [x]
		std::vector<int> n_bound_ii_;
		std::vector<int> n_bound_i_tethered_;
		std::vector< std::vector<int> > n_bound_i_tethered_bindable_;
		std::vector< std::vector<int> > n_bound_ii_tethered_;
		
		// Only one population for singly-bound untethered xlink heads
        int n_sites_i_ = 0;
		// Population for each xlink extension
		std::vector<int> n_sites_ii_;
		// Population for each tether extension
		std::vector<int> n_sites_i_tethered_;
		// Population for each teth AND each xlink extension
		std::vector< std::vector<int> > n_sites_ii_tethered_same_; 
		std::vector< std::vector<int> > n_sites_ii_tethered_oppo_; 

		int teth_cutoff_; 		// see kinesin header (is dist_cutoff_ there)
		int dist_cutoff_;		// see assoc. protein header
		int rest_dist_;			// see assoc. protein header

		// Stand-alone xlink diffusion probabilities
		double p_diffuse_i_fwd_;
		double p_diffuse_i_bck_;
		std::vector<double> p_diffuse_ii_to_rest_;
		std::vector<double> p_diffuse_ii_from_rest_;
		// Tether-dependent xlink diffusion probabilities
		std::vector<double> p_diffuse_i_to_teth_rest_;
		std::vector<double> p_diffuse_i_from_teth_rest_;
		std::vector< std::vector<double> > p_diffuse_ii_to_both_rest_;
		std::vector< std::vector<double> > p_diffuse_ii_from_both_rest_;
		std::vector< std::vector<double> > p_diffuse_ii_to_self_from_teth_;
		std::vector< std::vector<double> > p_diffuse_ii_from_self_to_teth_;

		// Kinematics probabilities
		double p_bind_i_tethered_;
		double p_bind_ii_;
		double p_tether_free_; 
		double p_untether_free_;
		// Indexed by number of PRC1 neighbors: [0], [1], or [2]
		std::vector<double> p_bind_i_;			
		std::vector<double> p_unbind_i_;		
		// Indexed by extension: [x_dub][x] or just [x_dub] / [x]
		std::vector<double> p_unbind_i_tethered_;	
		std::vector<double> p_unbind_ii_;
		std::vector< std::vector<double> > p_unbind_ii_to_teth_;
		std::vector< std::vector<double> > p_unbind_ii_from_teth_;

		// 1-D vectors, index is simply xlink entry
		std::vector<AssociatedProtein> xlinks_;
		std::vector<AssociatedProtein*> active_;
		std::vector<AssociatedProtein*> free_tethered_;
		std::vector<AssociatedProtein*> bound_i_all_; 
		std::vector<AssociatedProtein*> bound_untethered_;
		// 2-D vectors, 1st index is x or x_dub, 2nd is xlink entry
		std::vector< std::vector<AssociatedProtein*> > bound_ii_;
		std::vector< std::vector<AssociatedProtein*> > bound_i_tethered_;
		// 2-D vectors, but 1st index is number of PRC1 neighbors
		std::vector< std::vector<AssociatedProtein*>> bound_i_;
		// 3-D vectors, 1st index is x_dist_dub, 2nd is x_dist, 3rd is entry
		// e.g. vec_[16][2][0] is the 1st xlink w/ x_dist_dub=16, x_dist=2
		std::vector< std::vector< std::vector<AssociatedProtein*> > >
			bound_ii_tethered_;

		// Following vectors refer to sites that xlinks are bound
		// to, but the same indexing convention as above applies
		std::vector<Tubulin*> sites_i_untethered_; 
		std::vector< std::vector<Tubulin*> > sites_ii_untethered_;
		std::vector< std::vector<Tubulin*> > sites_i_tethered_;
		// 'oppo' refers to tether rest and xlink rest being 
		// on opposite sides of the site; likewise for 'same' 
		std::vector< std::vector< std::vector<Tubulin*> > >
			sites_ii_tethered_oppo_;
		std::vector< std::vector< std::vector<Tubulin*> > >
			sites_ii_tethered_same_;	

		std::vector<event> diffu_events_;
		std::vector<std::vector<int> > diffu_by_pop_; 
		std::vector<int> dif_list_;
		std::vector<event> kmc_events_;
		std::vector<std::vector<int> > kmc_by_pop_;
		std::vector<int> kmc_list_;

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

		void UpdateAllLists();
		void UpdateSingleBoundList();
		void UpdateBoundITethered();
		void UpdateDoubleBoundList();
		void UpdateBoundIITethered();
		void UpdateFreeTetheredList();
		void UpdateBoundUntethered();
		
		void UpdateAllSiteLists();
		void UpdateSingleUntetheredSites();
		void UpdateDoubleUntetheredSites();
		void UpdateSingleTetheredSites();
		void UpdateDoubleTetheredSites();

		AssociatedProtein* GetFreeXlink();
		AssociatedProtein* GetBoundUntetheredXlink();

		void GenerateDiffusionList();

		void RunDiffusion();
		void Diffuse_I_Forward(int n_neighbs);
		void Diffuse_I_Backward(int n_neighbs);
		void Diffuse_II_ToRest(int n_neighbs, int x);
		void Diffuse_II_FromRest(int n_neighbs, int x);
		void Diffuse_I_ToTethRest(int n_neighbs, int x_dub);
		void Diffuse_I_FromTethRest(int n_neighbs, int x_dub);
		void Diffuse_II_ToBothRest(int n_neighbs, int x_dub, int x);
		void Diffuse_II_FromBothRest(int n_neighbs, int x_dub, int x);
		void Diffuse_II_ToSelf_FromTeth(int n_neighbs, int x_dub, int x);
		void Diffuse_II_FromSelf_ToTeth(int n_neighbs, int x_dub, int x);

		void GenerateKMCList();
		double GetWeightBindII(); 
		double GetWeightBindITethered();
		double GetWeightBindIITethered();

		void RunKMC();
		void Bind_I(int n_neighbs);
		void Bind_II();		//XXX add neighb coop in weights
		void Unbind_I(int n_neighbs);
		void Unbind_II(int n_neighbs, int x);
		void Bind_I_Tethered(); //XXX add neighb coop in weights
		void Bind_II_Tethered(); //XXX add neighb coop in weights
		void Unbind_I_Tethered(int n_neighbs, int x_dub); 
		void Unbind_II_To_Teth(int n_neighbs, int x_dub, int x);
		void Unbind_II_From_Teth(int n_neighbs, int x_dub, int x);
		void Tether_Free(); 
		void Untether_Free();
};
#endif
