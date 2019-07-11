#include "event.h"
#include "associated_protein_management.h" 

template <typename POP_T, typename SITE_T, typename MGMT_T>
std::variant<POP_T, SITE_T> Event<POP_T, SITE_T, MGMT_T>::GetActiveEntry(){

	Entry target;
	if(*n_avail_ > 0) target = target_pool_->at(ran_int_(*n_avail_));
	else target = manager_->CheckScratchFor(target_pop_);
	try{
		if(std::get<POP_T>(target) == nullptr)
			printf("Error in Events:GetActiveEntry BROBRO\n");
	}
	catch(...){
		if(std::get<SITE_T>(target) == nullptr)
			printf("Error TWO in Events:GetActiveEntry BROBRO\n");
	}
	return target;
};

template class 
Event<AssociatedProtein::Monomer*, Tubulin*, AssociatedProteinManagement*>;
