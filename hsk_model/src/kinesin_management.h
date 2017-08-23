#ifndef _KINESIN_MANAGEMENT_H
#define _KINESIN_MANAGEMENT_H

#include "kinesin.h"
#include <vector>
#ifndef _PARAMETERS_H
typedef struct system_parameters system_parameters;
#endif
#ifndef _SYSTEM_PROPERTIES_H
typedef struct system_properties system_properties;
#endif

class KinesinManagement{
	private:

	public:
		int n_tot_;			// Total number in simulation
		int	n_bound_; 		// Total number of bound kinesin
		int	n_unbound_;		// Total number of unbound kinesin

		double p_bind_;
		double p_unbind_;
		double p_switch_;
		double p_step_;

		double alpha_;
		double beta_;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
		
		std::vector<Kinesin> active_motors;
		std::vector<Kinesin*> bound_list;
		std::vector<Kinesin*> unbound_list;
		std::vector<int> kmc_list;
	private:

	public:
		KinesinManagement();
		void Initialize(system_parameters *parameters, 
						system_properties *properties);

		void SetParameters();
		void GenerateActiveMotors();
		void GenerateUnboundList();

		void UnboundCheck(int ID);
		void BoundCheck(int ID);
		bool BoundaryCheck(Kinesin *motor);

		void GenerateKMCList();
		int GetNumToBind();
		int GetNumToUnbind();
		int GetNumToSwitch(); 
		int GetNumToStep();

		void RunKMC();
		void RunKMC_Bind();
		void RunKMC_Unbind();
		void RunKMC_Switch();
		void RunKMC_Step();
		void RunKMC_Boundaries(int n_events);
};
#endif
