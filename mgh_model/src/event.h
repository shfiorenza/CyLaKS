#pragma once
#include <string>
#include <functional>
struct system_parameters;
struct system_properties;

template <typename POP_T, typename MGMT>
class Event{
	private:
		MGMT* manager_;
		std::vector<POP_T> target_pool_;
		std::function<int(int)> ran_int_;
		std::function<int(double, int)> prob_dist_;
		std::function<void(POP_T target)> execute_; 

	public:
		int ID_ = -1;
		int code_ = -1;
		std::string name_ = "bruh";
		std::string target_pop_ = "naaw";	
		double p_occur_ = 0;
		int *n_avail_ = nullptr;
		int n_expected_ = 0;

	private:
		POP_T GetActiveEntry();

	public:
		Event(MGMT* manager, int ID, int code, std::string name, 
				std::string target_pop, double p_occur, int* pop_ptr, 
				std::vector<POP_T> pop_pool, std::function<int(int)> ran_int, 
				std::function<int(double, int)> prob_dist);
//		, std::function<void(POP_T tar)> exe);
		void SampleStatistics();
		void Execute();
};
