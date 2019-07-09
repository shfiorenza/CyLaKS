#include "event.h"

template <typename POP_T, typename MGMT>
Event<POP_T, MGMT>::Event(MGMT* manager, int ID, int code, std::string name,
	   	std::string target_pop, double p_occur, int* n_avail_ptr, 
		std::vector<POP_T> pop_pool, std::function<int(int)> ran_int,
		std::function<int(double, int)> prob_dist){
//		, std::function<void(POP_T tar)> exe){

	ID_ = ID;
	code_ = code;
	name_ = name;
	target_pop_ = target_pop;
	p_occur_ = p_occur;
	n_avail_ = n_avail_ptr;
	target_pool_ = pop_pool;
	ran_int_ = ran_int;
	prob_dist_ = prob_dist;
//	execute_ = exe;
};

template <typename POP_T, typename MGMT>
void Event<POP_T, MGMT>::SampleStatistics(){

	if(*n_avail_ > 0) n_expected_ = prob_dist_(p_occur_, *n_avail_);
	else n_expected_ = 0;
};

template <typename POP_T, typename MGMT>
void Event<POP_T, MGMT>::Execute(){

	POP_T entry = GetActiveEntry();
//	execute(entry); 
	manager_->KMC_Relay(entry, code_);
};

template <typename POP_T, typename MGMT>
POP_T Event<POP_T, MGMT>::GetActiveEntry(){

	POP_T entry; 
	if(*n_avail_ > 0){
		int i_entry = ran_int_(*n_avail_);
		entry = target_pool_[i_entry];
	}
	else entry = manager_->CheckScratchFor(target_pop_);
	if(entry == nullptr)
		printf("Error in Events:GetActiveEntry BROBRO");
	return entry;
};
