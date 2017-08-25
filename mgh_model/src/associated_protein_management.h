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
		int n_tot_ = 0;
		int n_single_bound_ = 0;
		int n_double_bound_ = 0;

		double p_bind_;
		double p_unbind_;
		
		system_parameters *parameters = nullptr;
		system_properties *properties = nullptr;

		std::vector<AssociatedProtein> xlink_list; 
	private:

	public:
		AssociatedProteinManagement();
		void Initialize(system_parameters *parameters, 
						system_properties *properties);
};
#endif
