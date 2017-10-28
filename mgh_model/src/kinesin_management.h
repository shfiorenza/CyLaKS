#ifndef _KINESIN_MANAGEMENT_H
#define _KINESIN_MANAGEMENT_H

#include "kinesin.h"	// Also includes <vector> lib
#ifndef _PARAMETERS_H
typedef struct system_parameters system_parameters;
#endif
#ifndef _SYSTEM_PROPERTIES_H
typedef struct system_properties system_properties;
#endif

class KinesinManagement{
	private:

	public:
		int n_motors_ = 0;		
	 	int	n_single_bound_ = 0;
		int n_double_bound_ = 0;
		int n_unbound_but_tethered_ = 0;

		double p_bind_;
		double p_bind_tethered_;
		double p_unbind_;
		double p_tether_unbound_;
		double p_switch_;
		double p_step_;

		double alpha_;
		double beta_;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
		
		std::vector<Kinesin> motor_list_; 
		std::vector<Kinesin*> bound_list_;
		std::vector<Kinesin*> unbound_but_tethered_list_;
		std::vector<int> kmc_list_;
	private:

	public:
		KinesinManagement();
		void Initialize(system_parameters *parameters, 
						system_properties *properties);
		void SetParameters();
		void GenerateMotors();

		void UnboundCheck(Kinesin *motor);
		void BoundCheck(Kinesin *motor);
		bool BoundaryStatus(Kinesin *motor);

		void UpdateLists();
		
		Kinesin* GetUnboundMotor();

		void RunDiffusion();

		void GenerateKMCList();
		int GetNumToBind_I();
		int GetNumToBind_I_Tethered();
		int GetNumToUnbind_I();
		int GetNumToUnbind_II();
		int GetNumToTether_Unbound();
		int GetNumToTether_Bound();
		int GetNumToSwitch(); 
		int GetNumToStep();

		void RunKMC();
		void KMC_Bind_I();
		void KMC_Bind_I_Tethered();
		void KMC_Bind_II();
		void KMC_Unbind_I();
		void KMC_Unbind_II();
		void KMC_Tether_Unbound();
		void KMC_Switch();
		void KMC_Step();
		void KMC_Boundaries(int n_events);
};
#endif
