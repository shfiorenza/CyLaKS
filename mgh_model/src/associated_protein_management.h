#ifndef _ASSOCIATED_PROTEIN_MANAGEMENT_H
#define _ASSOCIATED_PROTEIN_MANAGEMENT_H
#include "associated_protein.h"
struct system_parameters;
struct system_properties;

class AssociatedProteinManagement{
	private:

	public:
		int n_xlinks_ = 0;

		int n_single_bound_ = 0;
		std::vector<int> n_bound_i_tethered_;
		int n_bound_i_tethered_tot_ = 0;
		// Each possible spring extension has its own double bound list
		std::vector<int> n_double_bound_;
		std::vector< std::vector<int> > n_bound_ii_tethered_;
		int n_free_tethered_ = 0;
		int n_untethered_ = 0;
		
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
		// Only need different populations for each tether extension
		// when at 'self-rest', i.e. an x_dist of 0
		std::vector<int> n_sites_ii_tethered_self_rest_;

		int dist_cutoff_;		// see assoc. protein header
		int rest_dist_;			// see assoc. protein header

		// Diffusion of xlinks
		double p_diffuse_i_fwd_;
		double p_diffuse_i_bck_;
		std::vector<double> p_diffuse_ii_to_rest_;
		std::vector<double> p_diffuse_ii_from_rest_;
		// Diffusion that is tether-dependent (below)
		std::vector<double> p_diffuse_i_to_teth_rest_;
		std::vector<double> p_diffuse_i_from_teth_rest_;
		std::vector< std::vector<double> > p_diffuse_ii_to_both_rest_;
		std::vector< std::vector<double> > p_diffuse_ii_from_both_rest_;
		std::vector< std::vector<double> > p_diffuse_ii_to_self_from_teth_;
		std::vector< std::vector<double> > p_diffuse_ii_from_self_to_teth_;

		double tau_i_;	
		double tau_ii_;		

		double p_bind_i_;			
		double p_bind_i_tethered_; 
		double p_bind_ii_; 
		double p_unbind_i_;		
		std::vector<double> p_unbind_i_tethered_;	
		std::vector<double> p_unbind_ii_;	// One for each extension
		std::vector< std::vector<double> > p_unbind_ii_to_teth_;
		std::vector< std::vector<double> > p_unbind_ii_from_teth_;
		double p_tether_free_; 
		double p_untether_free_;

		std::vector<AssociatedProtein> xlink_list_;
		std::vector<AssociatedProtein*> single_bound_list_;
		std::vector< std::vector<AssociatedProtein*> > bound_i_tethered_;
		std::vector< std::vector<AssociatedProtein*> > double_bound_list_;
		std::vector< std::vector< std::vector<AssociatedProtein*> > >
			bound_ii_tethered_;
		std::vector<AssociatedProtein*> free_tethered_list_;
		std::vector<AssociatedProtein*> untethered_list_;

		// The following lists all refer to sites of active xlinks
		std::vector<Tubulin*> single_untethered_sites_;
		std::vector< std::vector<Tubulin*> > double_untethered_sites_;
		std::vector< std::vector<Tubulin*> > single_tethered_sites_;
		// 'oppo' refers to tether rest and xlink rest being 
		// on opposite sides of the site; likewise for 'same' 
		std::vector< std::vector< std::vector<Tubulin*> > >
			double_tethered_sites_oppo_;
		std::vector< std::vector< std::vector<Tubulin*> > >
			double_tethered_sites_same_;	

		std::vector<int> diffusion_list_;
		std::vector<int> kmc_list_;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
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
		void UpdateUntetheredList();
		
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
