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
		// Each possible spring extension has its own double bound list
		std::vector<int> n_double_bound_; 
		int n_tethered_ = 0;
		int n_untethered_ = 0;
		
		int dist_cutoff_ = 6;			// for stage2 binding; in sites

		int n_sites_single_tethered_ = 0;
        int n_sites_double_tethered_ = 0;
        int n_sites_single_untethered_ = 0;
        int n_sites_double_untethered_ = 0;

        // all in seconds; from sci adv supp
		// suffix syntax is (binding stage)_(kinesin tether status)
		double tau_single_tethered_ = 0.00192;
		double tau_double_tethered_ = 0.0288;	
		double tau_single_untethered_ = 0.00032;
		double tau_double_untethered_ = 0.00480;

		double c_eff_ = 100;				// unitless 	(educated guess)
		double p_bind_i_;			
		double p_unbind_i_;		
		std::vector<double> p_unbind_ii_;	// One for each extension

		std::vector<AssociatedProtein> xlink_list_;
		std::vector<AssociatedProtein*> single_bound_list_;
		std::vector< std::vector<AssociatedProtein*> > double_bound_list_;
		std::vector<AssociatedProtein*> untethered_list_;

		// The following lists all refer to sites of active xlinks
		std::vector<Tubulin*> single_tethered_sites_;
		std::vector<Tubulin*> double_tethered_sites_;
		std::vector<Tubulin*> single_untethered_sites_;
		std::vector<Tubulin*> double_untethered_sites_;

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
		void UpdateDoubleBoundList();
		void UpdateUntetheredList();
		
		void UpdateSingleTetheredSites();
		void UpdateDoubleTetheredSites();
		void UpdateSingleUntetheredSites();
		void UpdateDoubleUntetheredSites();

		void UpdateExtensions();

		AssociatedProtein* GetUntetheredXlink();

		void RunDiffusion();
		void RunDiffusion_Single_Tethered();
		void RunDiffusion_Double_Tethered();
		void RunDiffusion_Single_Untethered();
		void RunDiffusion_Double_Untethered();

		void GenerateKMCList();
		int GetNumToBind_I();
		int GetNumToBind_II();
		int GetNumToUnbind_I();
		int GetNumToUnbind_II(int x_dist);

		void RunKMC();
		void RunKMC_Bind_I();
		void RunKMC_Bind_II();
		void RunKMC_Unbind_I();
		void RunKMC_Unbind_II(int x_dist);
};
#endif
