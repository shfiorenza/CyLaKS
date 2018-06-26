#ifndef _KINESIN_MANAGEMENT_H
#define _KINESIN_MANAGEMENT_H
#include "kinesin.h"				// Also includes <vector> lib

struct system_parameters;
struct system_properties;

class KinesinManagement{
	private:

	public:
		int n_motors_ = 0;		

		int n_free_tethered_ = 0;
	 	int	n_pseudo_bound_ = 0;
		int n_eligible_pseudo_ = 0;
		int n_bound_untethered_ = 0;
		int n_stalled_untethered_ = 0;
		int n_stepable_untethered_ = 0; 
		int n_bound_tethered_tot_ = 0; 
		std::vector<int> n_bound_tethered_;
		std::vector<int> n_stepable_to_teth_rest_;
		std::vector<int> n_stepable_from_teth_rest_;
		// switching currently disabled in KMC
		int n_switchable_;		// must be tethered to switch

		int dist_cutoff_;		// see kinesin header
		int comp_cutoff_;		// see kinesin header
		double rest_dist_; 		// see kinesin header

		double p_diffuse_fwd_untethered_; 
		double p_diffuse_bck_untethered_;
		// Diffusing away FROM xlink means the extension increases, 
		// whereas diffusing TO xlink means extension is decreasing
		std::vector<double> p_diffuse_to_tether_rest_;
		std::vector<double> p_diffuse_from_tether_rest_;

		double tau_;		// in seconds, from weird paper
		
		double p_bind_i_free_;
		double p_bind_i_tethered_;
		double p_bind_ii_; 
		double p_unbind_pseudo_;
		double p_unbind_stepable_untethered_;
		double p_unbind_stalled_untethered_;
		double p_unbind_tethered_;
		double p_tether_free_;
		double p_tether_bound_;
		// Indices refer to double the x_distance (in no. of sites)
		// between the stalk of a motor and the anchor of an xlink
		// e.g. p_untether_bound_[23] is for an x_dist of 11.5 sites
		double p_untether_free_;	
		std::vector<double> p_untether_bound_; 	// One for each extension
		double p_step_untethered_;
		double p_failstep_untethered_;
		std::vector<double> p_step_to_teth_rest_;
		std::vector<double> p_step_from_teth_rest_;
		double p_switch_;

		double c_eff_; 		// Used in partition function for tethering
		// boundaries currently disabled in KMC
		double alpha_;		// Prob. to insert onto the minus end
		double beta_;		// Prob. to remove from the plus end
		
		std::vector<Kinesin> motor_list_; 
		std::vector<Kinesin*> free_tethered_list_;
		std::vector<Kinesin*> pseudo_bound_list_;	// Only 1 head bound
		std::vector<Kinesin*> eligible_pseudo_list_;
		std::vector<Kinesin*> bound_untethered_list_;
		std::vector<Kinesin*> stepable_untethered_list_;
		std::vector<Kinesin*> stalled_untethered_list_;
		// For our purposes, 'list' means one-dimensional vector, and
		// the 'table' means two-dimensions: tether extension and index
		std::vector<Kinesin*> bound_tethered_list_;	
		std::vector< std::vector<Kinesin*> > bound_tethered_table_;
		std::vector< std::vector<Kinesin*> > stepable_to_rest_table_;
		std::vector< std::vector<Kinesin*> > stepable_from_rest_table_;
		std::vector<Kinesin*> switchable_list_; 

		std::vector<int> diffusion_list_;
		std::vector<int> kmc_list_;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
	private:

	public:
		KinesinManagement();
		void Initialize(system_parameters *parameters, 
						system_properties *properties);
		void SetParameters();
		void GenerateMotors();
		void InitiateLists();

		void UnboundCheck(Kinesin *motor);
		void BoundCheck(Kinesin *motor);
		bool BoundaryStatus(Kinesin *motor);

		Kinesin* GetFreeMotor();

		void UpdateFreeTetheredList();
		void UpdatePseudoBoundList();
		void UpdateEligiblePseudoList();
		void UpdateBoundUntetheredList();
		void UpdateStepableUntetheredList();
		void UpdateStalledUntetheredList();
		void UpdateBoundTetheredList();
		void UpdateBoundTetheredTable();
		void UpdateStepableTetheredTables();
		void UpdateSwitchableList();
		
		void GenerateDiffusionList();
		int GetNumToStepForward_Unteth();
		int GetNumToStepBackward_Unteth();
		int GetNumToStepTowardRest(int x_dist_doubled);
		int GetNumToStepFromRest(int x_dist_doubled); 

		void RunDiffusion();
		void RunDiffusion_Forward_Untethered();
		void RunDiffusion_Backward_Untethered();
		void RunDiffusion_Toward_Rest(int x_dist_doubled);
		void RunDiffusion_From_Rest(int x_dist_doubled);

		void GenerateKMCList();
		// Roman numerals refer to STAGE of binding, not # of heads to bind
		int GetNumToBind_I_Free();
		int GetNumToBind_I_Tethered();
		int GetNumToBind_II();
		int GetNumToUnbind_Pseudo(); 
		int GetNumToUnbind_Stepable_Untethered();
		int GetNumToUnbind_Stalled_Untethered();
		int GetNumToUnbind_Tethered();
		int GetNumToTether_Free();
		int GetNumToTether_Bound();
		int GetNumToUntether_Bound(int x_dist_doubled);
		int GetNumToUntether_Free();
		int GetNumToStep_Untethered();
		int GetNumToFailstep_Untethered();
		int GetNumToStep_ToTethRest(int x_dist_doubled);
		int GetNumToStep_FromTethRest(int x_dist_doubled);
		int GetNumToSwitch(); 

		void RunKMC();
		void KMC_Bind_I_Free();			// bind from bulk solution
		void KMC_Bind_I_Tethered();		// bind from nearby tether
		void KMC_Bind_II();				// bind 2nd head via 'diffusion'
		void KMC_Unbind_Pseudo();
		void KMC_Unbind_Stepable_Untethered();
		void KMC_Unbind_Stalled_Untethered(); 
		void KMC_Unbind_Tethered();
		void KMC_Tether_Free();
		void KMC_Tether_Bound();
		void KMC_Untether_Free();
		void KMC_Untether_Bound(int x_dist_doubled);
		void KMC_Step_Untethered();
		void KMC_Failstep_Untethered();
		void KMC_Step_ToTethRest(int x_dist_doubled);
		void KMC_Step_FromTethRest(int x_dist_doubled);
		void KMC_Boundaries(int n_events);
		void KMC_Switch();
};
#endif
