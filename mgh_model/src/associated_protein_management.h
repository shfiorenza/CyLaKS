#ifndef _ASSOCIATED_PROTEIN_MANAGEMENT_H
#define _ASSOCIATED_PROTEIN_MANAGEMENT_H
#include "associated_protein.h"
#include <map>
#include <string>
#include <functional>

struct system_parameters;
struct system_properties;

struct pop_t{
	int n_entries_ = -1;
	std::string type_ = std::string("wut");
	int x_dist_ = -1;
	int x_dist_dub_ = -1;
};

class AssociatedProteinManagement{
	private:
		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;

	public:
		int n_xlinks_ = 0;
		// xlinks actively bound to some MT or motor; dynamically updated
		int n_active_ = 0;

		int n_bound_i_ = 0;
		int n_free_tethered_ = 0;
		int n_bound_untethered_ = 0;
		int n_bound_i_tethered_tot_ = 0;
		std::vector<int> n_bound_i_tethered_;
		// Each possible spring extension has its own double bound list
		std::vector<int> n_double_bound_;
		std::vector< std::vector<int> > n_bound_ii_tethered_;
		
		// Only one population for singly-bound untethered xlink heads
        int n_sites_i_untethered_ = 0;
		// Population for each xlink extension
		std::vector<int> n_sites_ii_untethered_;
		// Population for each tether extension
		std::vector<int> n_sites_i_tethered_;
		// Population for each teth AND each xlink extension
		std::vector< std::vector<int> > n_sites_ii_tethered_;
		std::vector< std::vector<int> > n_sites_ii_tethered_same_; 
		std::vector< std::vector<int> > n_sites_ii_tethered_oppo_; 
		// Only need tether extension when at 'self-rest', i.e. x_dist = 0
		std::vector<int> n_sites_ii_tethered_self_rest_;

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

		// characteristic time to diffuse one site (~8nm) for bound_i/ii
		double tau_i_;	
		double tau_ii_;		

		// Kinematics probabilities
		double p_bind_i_;			
		double p_bind_i_tethered_; 
		double p_bind_ii_; 
		double p_unbind_i_;		
		double p_tether_free_; 
		double p_untether_free_;
		std::vector<double> p_unbind_ii_;
		std::vector<double> p_unbind_i_tethered_;	
		std::vector< std::vector<double> > p_unbind_ii_to_teth_;
		std::vector< std::vector<double> > p_unbind_ii_from_teth_;

		// 1-D vectors, index is simply xlink entry
		std::vector<AssociatedProtein> xlinks_;
		std::vector<AssociatedProtein*> active_;
		std::vector<AssociatedProtein*> bound_i_;
		std::vector<AssociatedProtein*> free_tethered_;
		std::vector<AssociatedProtein*> bound_untethered_;
		// 2-D vectors, 1st index is x_dist, 2nd is xlink entry
		std::vector< std::vector<AssociatedProtein*> > bound_ii_;
		std::vector< std::vector<AssociatedProtein*> > bound_i_tethered_;
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

		// Holds diffusion and KMC events encoded in integer values
		std::vector<int> dif_list_;
		std::vector<int> kmc_list_;
		// Serialized (in regards to self/teth extension) vectors that store
		// number of entries, the population/event label, x, and x_dub
		// (facilitates the parallelization of statistical sampling) 
		std::vector<pop_t> serial_pop_; 
		std::vector<pop_t> serial_dif_;
		std::vector<pop_t> serial_kmc_;
		// Map of population/event label to sampling function;
		// see InitializeFUnctionMap() for function definitions
		std::map<std::string, std::function<int(int)> > dif_sampling_functs_;
		std::map<std::string, std::function<int(int)> > kmc_sampling_functs_;

	private:

	public:
		AssociatedProteinManagement();
		void Initialize(system_parameters *parameters, 
						system_properties *properties);
		void GenerateXLinks();
		void SetParameters();
		void InitiateLists();

		void BoundCheck(AssociatedProtein *xlink);	
		void UntetheredCheck(AssociatedProtein *xlink);
	
		void UpdateSingleBoundList();
		void UpdateBoundITethered();
		void UpdateDoubleBoundList();
		void UpdateBoundIITethered();
		void UpdateFreeTetheredList();
		void UpdateUntethered();
		
		void UpdateSingleUntetheredSites();
		void UpdateDoubleUntetheredSites();
		void UpdateSingleTetheredSites();
		void UpdateDoubleTetheredSites();

		AssociatedProtein* GetFreeXlink();
		AssociatedProtein* GetUntetheredXlink();

		void GenerateDiffusionList();
		int GetNumToStepI_Forward();
		int GetNumToStepI_Backward();
		int GetNumToStepII_ToRest(int x_dist);
		int GetNumToStepII_FromRest(int x_dist);
		// Tether extensions taken into account for the below
		// x_dist_dub refers to tether extension, and is double its
		// real value (can have half-integer values). x_dist is the 
		// actual extension of the crosslink in no of x-sites
		int GetNumToStepI_ToTethRest(int x_dist_dub);
		int GetNumToStepI_FromTethRest(int x_dist_dub);
		int GetNumToStepII_ToBothRest(int x_dist_dub, int x_dist);
		int GetNumToStepII_FromBothRest(int x_dist_dub, int x_dist);
		int GetNumToStepII_ToSelf_FromTeth(int x_dist_dub, int x_dist);
		int GetNumToStepII_FromSelf_ToTeth(int x_dist_dub,int x_dist);

		void RunDiffusion();
		void RunDiffusionI_Forward();
		void RunDiffusionI_Backward();
		void RunDiffusionII_ToRest(int x_dist);
		void RunDiffusionII_FromRest(int x_dist);
		void RunDiffusionI_ToTethRest(int x_dist_dub);
		void RunDiffusionI_FromTethRest(int x_dist_dub);
		void RunDiffusionII_ToBothRest(int x_dist_dub, int x_dist);
		void RunDiffusionII_FromBothRest(int x_dist_dub, int x_dist);
		void RunDiffusionII_ToSelf_FromTeth(int x_dist_dub, int x_dist);
		void RunDiffusionII_FromSelf_ToTeth(int x_dist_dub, int x_dist);

		void GenerateKMCList();
		int GetNumToBind_I();
		int GetNumToBind_I_Tethered();
		int GetNumToBind_II();
		int GetNumToBind_II_Tethered();
		int GetNumToUnbind_I();
		int GetNumToUnbind_I_Tethered(int x_dist_dub);
		int GetNumToUnbind_II(int x_dist);
		int GetNumToUnbind_II_To_Teth(int x_dist_dub, int x_dist);
		int GetNumToUnbind_II_From_Teth(int x_dist_dub, int x_dist);
		int GetNumToTether_Free();
		int GetNumToUntether_Free();

		void RunKMC();
		void RunKMC_Bind_I();
		void RunKMC_Bind_I_Tethered();
		void RunKMC_Bind_II();
		void RunKMC_Bind_II_Tethered();  // XXX 
		void RunKMC_Unbind_I();
		void RunKMC_Unbind_I_Tethered(int x_dist_dub); 
		void RunKMC_Unbind_II(int x_dist);
		void RunKMC_Unbind_II_To_Teth(int x_dist_dub, int x_dist);
		void RunKMC_Unbind_II_From_Teth(int x_dist_dub, int x_dist);
		void RunKMC_Tether_Free();
		void RunKMC_Untether_Free();
};
#endif
