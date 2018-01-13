#ifndef _KINESIN_MANAGEMENT_H
#define _KINESIN_MANAGEMENT_H
#include "kinesin.h"				// Also includes <vector> lib

struct system_parameters;
struct system_properties;

class KinesinManagement{
	private:

	public:
		int n_motors_ = 0;		

	 	int	n_single_bound_ = 0;
		int n_double_bound_ = 0;
		// Each possible tether extension has its own bound list
		int n_free_tethered_ = 0;
		std::vector<int> n_bound_tethered_;
		int n_bound_untethered_ = 0; 

		int dist_cutoff_ = 18;

		double tau_ = 0.006; 		// from weird paper on desktop (lol)
		
		double alpha_;			// Prob. to insert onto the minus end
		double beta_;			// Prob. to remove from the plus end
		double c_eff_ = 10; 	// Used in partition function for tethering
		
		double p_bind_i_free_;
		double p_bind_i_tethered_;
		double p_bind_ii_; 
		double p_unbind_;
		double p_unbind_tethered_;
		double p_tether_free_;
		// Indices refer to double the x_distance (in no. of sites)
		// between the stalk of a motor and the anchor of an xlink
		// e.g. p_untether_bound_[23] is for an x_dist of 11.5 sites
		std::vector<double> p_untether_bound_; 		// One for each extension
		double p_untether_free_;	
		std::vector<double> p_step_tethered_;		// One for each extension
		double p_step_untethered_; 
		double p_switch_;

		std::vector<Kinesin> motor_list_; 
		std::vector<Kinesin*> free_tethered_list_;
		std::vector<Kinesin*> pseudo_bound_list_;	// when 1 head is bound
		std::vector<Kinesin*> bound_list_;
		std::vector< std::vector<Kinesin*> > bound_tethered_list_;
		std::vector<Kinesin*> bound_untethered_list_;

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

		void UpdateFreeTetheredList();
		void UpdatePseudoBoundList();
		void UpdateBoundList();
		void UpdateBoundTetheredList();
		void UpdateBoundUntetheredList();

		void UpdateExtensions();
		
		Kinesin* GetFreeMotor();

		void RunDiffusion();
		void RunDiffusion_Bound();

		void GenerateKMCList();
		int GetNumToBind_I_Free();
		int GetNumToBind_I_Tethered();
		int GetNumToBind_II();
		int GetNumToUnbind_Untethered();
		int GetNumToUnbind_Tethered();
		int GetNumToTether_Free();
		int GetNumToTether_Bound();
		int GetNumToSwitch(); 
		int GetNumToUntether_Bound(double x_dist);
		int GetNumToUntether_Free();
		int GetNumToStep_Tethered(double x_dist);
		int GetNumToStep_Untethered();

		void RunKMC();
		void KMC_Bind_I_Free();			// bind from bulk solution
		void KMC_Bind_I_Tethered();		// bind from nearby tether
		void KMC_Bind_II();				// bind 2nd head via 'diffusion'
		void KMC_Unbind_Untethered();
		void KMC_Unbind_Tethered();
		void KMC_Tether_Free();
		void KMC_Tether_Bound();
		void KMC_Switch();
		void KMC_Untether_Bound(double x_dist);
		void KMC_Untether_Free();
		void KMC_Step_Tethered(double x_dist);
		void KMC_Step_Untethered();
		void KMC_Boundaries(int n_events);
};
#endif
