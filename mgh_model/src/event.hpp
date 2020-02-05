#ifndef _EVENT
#define _EVENT
#include <functional>
#include <iostream>
#include <string>
#include <variant>

template <typename MGMT_T, typename ENTRY_T> class Event {
private:
  // Pointer to manager of this event
  MGMT_T manager_;
  // Pointer to list of available specific targets; dynamically updated
  std::vector<ENTRY_T> *target_pool_;
  // Function that updates target_pool_ vector
  std::function<void(void)> update_targets_;
  // Function that actually executes this event
  std::function<void(ENTRY_T)> exe_;
  // Probabilitiy distribution we sample to predict no. of events each timestep
  std::function<int(double, int)> prob_dist_;
  // Generic random integer generator; returns an int from 0 to input N
  std::function<int(int)> ran_int_;

public:
  // Unique ID of this event in event_ list held by MGMT
  int id_ = -1;
  // Expected number of events to occur for current given timestep
  int n_expected_ = 0;
  // Pointer to no. of specific targets this event can act on; dynamic
  int *n_avail_ = nullptr;
  // Probability that this event will occur each timestep
  double p_occur_ = 0.0;
  // Name of this event, e.g., "Bind_II_Teth"
  std::string name_ = "bruh";
  // List of populations that this event targets
  std::vector<std::string> targets_ = {"nope"};
  // Whether or not n_avail_ is redundant in stat correction
  bool is_not_redundant_ = true;

private:
  void UpdateTargetPool() { update_targets_(); }
  ENTRY_T GetActiveEntry() {
    ENTRY_T target;
    if (manager_->verbose_) {
      printf("%i avail for ", *n_avail_);
      std::cout << name_ << std::endl;
    }
    if (*n_avail_ > 0) {
      target = target_pool_->at(ran_int_(*n_avail_));
    } else {
      target = manager_->CheckScratchFor(targets_[0]);
    }
    return target;
  }
  void ExecuteEventOn(ENTRY_T target) { exe_(target); }

public:
  Event(MGMT_T manager, int id, std::string name,
        std::vector<std::string> targets, double p_occur, int *n_avail,
        std::vector<ENTRY_T> *target_pool, std::function<void(void)> update,
        std::function<void(ENTRY_T)> exe_funct,
        std::function<int(double, int)> prob_dist,
        std::function<int(int)> ran_int)
      : manager_{manager}, id_{id}, name_{name}, targets_{targets},
        is_not_redundant_{true}, p_occur_{p_occur}, n_avail_{n_avail},
        target_pool_{target_pool}, update_targets_{update}, exe_{exe_funct},
        prob_dist_{prob_dist}, ran_int_{ran_int} {}
  Event(MGMT_T manager, int id, std::string name,
        std::vector<std::string> targets, bool is_not_redundant, double p_occur,
        int *n_avail, std::vector<ENTRY_T> *target_pool,
        std::function<void(void)> update,
        std::function<void(ENTRY_T)> exe_funct,
        std::function<int(double, int)> prob_dist,
        std::function<int(int)> ran_int)
      : manager_{manager}, id_{id}, name_{name}, targets_{targets},
        is_not_redundant_{is_not_redundant}, p_occur_{p_occur},
        n_avail_{n_avail}, target_pool_{target_pool}, update_targets_{update},
        exe_{exe_funct}, prob_dist_{prob_dist}, ran_int_{ran_int} {}
  int SampleStatistics() {
    if (*n_avail_ > 0) {
      n_expected_ = prob_dist_(p_occur_, *n_avail_);
    } else {
      n_expected_ = 0;
    }
    return n_expected_;
  }
  void Execute() {
    UpdateTargetPool();
    ENTRY_T target = GetActiveEntry();
    ExecuteEventOn(target);
  }
};
#endif