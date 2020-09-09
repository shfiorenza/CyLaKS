#ifndef _EVENT
#define _EVENT
#include <functional>
#include <string>
#include <variant>

template <typename ENTRY_T> class Event {
private:
  // Whether or not this event uses poisson-based statistics
  bool poisson_based_{false};
  // Poisson based jawns
  double poisson_weight_total_;
  std::vector<double> poisson_weights_;
  std::function<double(ENTRY_T)> get_weight_;
  std::function<double(void)> get_ran_prob_;
  // Pointer to list of available targets to act on; dynamically updated
  std::vector<ENTRY_T> *target_pool_;
  // Function that actually executes this event
  std::function<void(ENTRY_T)> exe_;
  // Probabilitiy distribution we sample to predict no. of events each timestep
  std::function<int(double, int)> prob_dist_;
  // Function that sets n random indices in the range [0, m)
  std::function<void(int *, int, int)> set_ran_indices_;

public:
  unsigned long n_executed_tot_{0};
  unsigned long n_opportunities_tot_{0};
  // Expected number of events to occur for current given timestep
  int n_expected_{0};
  // Pointer to no. of specific targets this event can act on; dynamic
  int *n_avail_{nullptr};
  // Probability that this event will occur each timestep
  double p_occur_{0.0};
  // Name of this event, e.g., "Bind_II_Teth"
  std::string name_{"bruh"};
  // Targets that this event will act on this timestep
  std::vector<ENTRY_T> targets_;

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
  void SampleStatistics_Poisson();
  void SetTargets_Poisson();

public:
  Event(std::string name, double p_occur, int *n_avail,
        std::vector<ENTRY_T> *target_pool,
        std::function<void(ENTRY_T)> exe_funct,
        std::function<int(double, int)> prob_dist,
        std::function<void(int *, int, int)> set_ran_indices)
      : name_{name}, p_occur_{p_occur}, n_avail_{n_avail},
        target_pool_{target_pool}, exe_{exe_funct}, prob_dist_{prob_dist},
        set_ran_indices_{set_ran_indices} {}
  Event(std::string name, double p_occur, int *n_avail,
        std::vector<ENTRY_T> *target_pool,
        std::function<void(ENTRY_T)> exe_funct,
        std::function<double(ENTRY_T)> get_weight,
        std::function<int(double, int)> prob_dist,
        std::function<double(void)> get_ran_prob,
        std::function<void(int *, int, int)> set_ran_indices)
      : name_{name}, p_occur_{p_occur}, n_avail_{n_avail},
        target_pool_{target_pool}, exe_{exe_funct}, prob_dist_{prob_dist},
        get_weight_{get_weight}, get_ran_prob_{get_ran_prob},
        set_ran_indices_{set_ran_indices} {
    poisson_based_ = true;
  }
  int SampleStatistics() {
    n_opportunities_tot_ += *n_avail_;
    // if (name_ == "unbind_ii" and *n_avail_ > 0) {
    //   printf("%i opportunities for unbind_ii\n", *n_avail_);
    // }
    if (poisson_based_) {
      SampleStatistics_Poisson();
    } else {
      n_expected_ = prob_dist_(p_occur_, *n_avail_);
    }
    SetTargets();
    return n_expected_;
  }
  void RemoveTarget(ENTRY_T tar) {
    for (int i_entry{0}; i_entry < n_expected_; i_entry++) {
      if (tar == targets_[i_entry]) {
        targets_[i_entry] = targets_[n_expected_ - 1];
        n_expected_--;
        return;
      }
    }
  }
  void Execute() {
    exe_(targets_[n_expected_ - 1]);
    n_expected_--;
    n_executed_tot_++;
  }
};
#endif