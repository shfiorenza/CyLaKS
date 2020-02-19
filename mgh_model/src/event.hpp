#ifndef _EVENT
#define _EVENT
#include <functional>
#include <string>
#include <variant>

template <typename MGMT_T, typename ENTRY_T> class Event {
private:
  // Pointer to list of available targets to act on; dynamically updated
  std::vector<ENTRY_T> *target_pool_;
  // Function that actually executes this event
  std::function<void(ENTRY_T)> exe_;
  // Probabilitiy distribution we sample to predict no. of events each timestep
  std::function<int(double, int)> prob_dist_;
  // Function that sets n random indices in the range [0, m)
  std::function<void(int *, int, int)> set_ran_indices_;

public:
  // Targets that this event will act on this timestep
  std::vector<ENTRY_T> targets_;
  // Expected number of events to occur for current given timestep
  int n_expected_ = 0;
  // Pointer to no. of specific targets this event can act on; dynamic
  int *n_avail_ = nullptr;
  // Probability that this event will occur each timestep
  double p_occur_ = 0.0;
  // Name of this event, e.g., "Bind_II_Teth"
  std::string name_ = "bruh";

private:
  void SetTargets() {
    if (n_expected_ > targets_.size()) {
      targets_.resize(n_expected_);
    }
    int indices[n_expected_];
    set_ran_indices_(indices, n_expected_, *n_avail_);
    for (int i_entry{0}; i_entry < n_expected_; i_entry++) {
      targets_[i_entry] = target_pool_->at(indices[i_entry]);
    }
  }

public:
  Event(std::string name, double p_occur, int *n_avail,
        std::vector<ENTRY_T> *target_pool,
        std::function<void(ENTRY_T)> exe_funct,
        std::function<int(double, int)> prob_dist,
        std::function<void(int *, int, int)> get_ran_indices)
      : name_{name}, p_occur_{p_occur}, n_avail_{n_avail},
        target_pool_{target_pool}, exe_{exe_funct}, prob_dist_{prob_dist},
        set_ran_indices_{get_ran_indices} {}
  int SampleStatistics() {
    n_expected_ = prob_dist_(p_occur_, *n_avail_);
    SetTargets();
    return n_expected_;
  }
  void RemoveTarget(int index) {
    targets_[index] = targets_[n_expected_ - 1];
    n_expected_--;
  }
  void Execute() {
    exe_(targets_[n_expected_ - 1]);
    n_expected_--;
  }
};
#endif