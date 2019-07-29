#pragma once
#include <string>
#include <variant>
#include <iostream>
#include <functional>
struct system_parameters;
struct system_properties;

class AssociatedProteinManagement;

template <typename MGMT_T, typename ENTRY_T>
class Event{
	private:
		MGMT_T manager_;
		std::vector<ENTRY_T>* target_pool_;
		std::function<int(int)> ran_int_;
		std::function<int(double, int)> prob_dist_;

	public:
		int ID_ = -1;
		int code_ = -1;
		std::string name_ = "bruh";
		std::string target_pop_ = "naaw";	
		double p_occur_ = 0;
		int *n_avail_ = nullptr;
		int n_expected_ = 0;

	private:
		ENTRY_T GetActiveEntry(){
			ENTRY_T target;
			if(*n_avail_ > 0) target = target_pool_->at(ran_int_(*n_avail_));
			else target = manager_->CheckScratchFor(target_pop_);
			return target;
		};

	public:
		Event(){};
		Event(MGMT_T manager, int ID, int code, std::string name, 
				std::string target, double p_occur, int* avail, 
				std::vector<ENTRY_T>* pool, std::function<int(int)> ran_int, 
				std::function<int(double, int)> prob_dist){

			manager_ = manager; ID_ = ID; code_ = code; name_ = name;
			target_pop_ = target; p_occur_ = p_occur; n_avail_ = avail;
			target_pool_ = pool; ran_int_ = ran_int; prob_dist_ = prob_dist;
		};
		void SampleStatistics(){
			if(*n_avail_ > 0) n_expected_ = prob_dist_(p_occur_, *n_avail_);
			else n_expected_ = 0;
		};
		void Execute(){
			manager_->Update_List_Relay(name_, target_pop_);
			ENTRY_T target = GetActiveEntry();
			manager_->Execute_Function_Relay(target, code_);
		};
};
