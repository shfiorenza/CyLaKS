#ifndef _ASSOCIATED_PROTEIN_MANAGEMENT_H
#define _ASSOCIATED_PROTEIN_MANAGEMENT_H
#include "associated_protein.h"

#ifndef _PARAMETERS_H
typedef struct system_parameters system_parameters;
#endif
#ifndef _SYSTEM_PROPERTIES_H
typedef struct system_properties system_properties;
#endif

class AssociatedProteinManagement{
	private:

	public:
		int n_xlinks_ = 0;
		int n_single_bound_ = 0;
		int n_double_bound_ = 0;
		int n_tethered_ = 0;
		int n_untethered_ = 0;
		int n_sites_occupied_ = 0;

		double p_bind_first_;
		double p_bind_second_;
		double p_unbind_;
		
		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;

		std::vector<AssociatedProtein> xlink_list_;
		std::vector<AssociatedProtein*> tethered_list_;
		std::vector<AssociatedProtein*> untethered_list_;
		std::vector<Tubulin*> occupied_site_list_;
		std::vector<int> kmc_list_;
	private:

	public:
		AssociatedProteinManagement();
		void Initialize(system_parameters *parameters, 
						system_properties *properties);

		void SetParameters();
		void GenerateXLinks();

		void BoundCheck(AssociatedProtein *xlink);	
		void UntetheredCheck(AssociatedProtein *xlink);
	
		void UpdateTetherLists();	
		void UpdateOccupiedSiteList();

		AssociatedProtein* GetUntetheredXlink();

		void RunDiffusion();

		void GenerateKMCList();
		int GetNumToBindFirst();
		int GetNumToBindSecond();
		int GetNumToUnbind();

		void RunKMC();
		void RunKMC_BindFirst();
		void RunKMC_BindSecond();
		void RunKMC_Unbind();
};
#endif
