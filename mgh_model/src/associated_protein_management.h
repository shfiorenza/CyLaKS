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
		int n_double_bound_ = 0;
		int n_tethered_ = 0;
		int n_untethered_ = 0;
		int n_sites_occupied_ = 0;

		/* holy shit motherfucker look at this::
				c = z_m * rho^2 * b^2 * exp(u_o/ kbT)
		   use this jawn for prefactors of partition function */ 

		double kbT_ = 4.11e-21;		//FIXME   in Joules, from bobert
		double k_spring_ = 1.3e-5;	//FIXME   in N/m, from robert
		double site_size_ = 8e-9; 	//FIXME   in fuckin m wow 
		double b_squared = 6.25e-16; 		//FIXME   in m^2, from bober
		double mt_length_squared_;

		double p_bind_i_;			// for first head
		double p_bind_ii_;			// second head  currently not used
		double p_unbind_i_;			// first head (i.e. when single-bound)
		double p_unbind_ii_;		// second head (i.e. when double-bound)

		std::vector<AssociatedProtein> xlink_list_;
		std::vector<AssociatedProtein*> single_bound_list_;
		std::vector<AssociatedProtein*> double_bound_list_;

		std::vector<AssociatedProtein*> tethered_list_;
		std::vector<AssociatedProtein*> untethered_list_;
		std::vector<Tubulin*> occupied_site_list_;
		std::vector<int> kmc_list_;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
	private:

	public:
		AssociatedProteinManagement();
		void Initialize(system_parameters *parameters, 
						system_properties *properties);

		void SetParameters();
		void GenerateXLinks();

		void BoundCheck(AssociatedProtein *xlink);	
		void UntetheredCheck(AssociatedProtein *xlink);
	
		void UpdateSingleBoundList();
		void UpdateDoubleBoundList();
		void UpdateTetheredLists();	
		void UpdateOccupiedSiteList();
		
		void UpdateAllNeighborSiteLists();
		void UpdateAllNeighborMotors();

		AssociatedProtein* GetUntetheredXlink();

		void RunDiffusion();

		void GenerateKMCList();
		int GetNumToBind_I();
		int GetNumToBind_II();
		int GetNumToUnbind_I();
		int GetNumToUnbind_II();

		void RunKMC();
		void RunKMC_Bind_I();
		void RunKMC_Bind_II();
		void RunKMC_Unbind_I();
		void RunKMC_Unbind_II();
};
#endif
