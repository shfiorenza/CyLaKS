#ifndef _EVENT
#define _EVENT
#include <functional>
#include <string>
#include <variant>

template <typename ENTRY_T> class Event {
private:
  // Whether or not this event uses poisson-based statistics
  bool poisson_based_{false};
  // Total poisson weight summed over all entries in target_pool_
  double poisson_weight_total_;
  // Individual poisson weights for each entry in target_pool_
  std::vector<double> poisson_weights_;
  // Returns appropriate poisson weight for given input entry
  std::function<double(ENTRY_T)> get_weight_;
  // Generates a random probability in the range [0, 1)
  std::function<double(void)> get_ran_prob_;
  // Pointer to list of available targets to act on -- dynamically updated
  std::vector<ENTRY_T> *target_pool_;
  // Executes event on given input target
  std::function<void(ENTRY_T)> exe_;
  // Probabilitiy distribution we sample from to predict n_events each timestep
  std::function<int(double, int)> prob_dist_;
  // Sets n random indices in the range [0, m) in given input array
  std::function<void(int *, int, int)> set_ran_indices_;

public:
  // Number of times we actually executed this event over entire simulation
  unsigned long n_executed_tot_{0};
  // Number of opportunities we had to execute this event over entire simulation
  unsigned long n_opportunities_tot_{0};
  // Number of events expected to occur for this timestep
  int n_expected_{0};
  // Pointer to number of targets this event can act on -- dynamically updated
  int *n_avail_{nullptr};
  // Probability that this event will occur each timestep -- DYNAMIC ?!
  double p_occur_{0.0};
  double *p_occur_ptr_{nullptr};
  // Name of this event, e.g., "Bind_II_Teth"
  std::string name_{"bruh"};
  // Targets that this event will act on this timestep
  std::vector<ENTRY_T> targets_;

private:
  // Set n_expected_ targets in random order
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
  // Get total weights and sample from the poisson distribution
  void SamplePoissonStatistics();
  // Choose n_expected_ targets based on poisson weights
  void SetPoissonCandidates();

public:
  // Constructor for binomial-based events
  Event(std::string name, double *p_occur_ptr, int *n_avail,
        std::vector<ENTRY_T> *target_pool,
        std::function<void(ENTRY_T)> exe_funct,
        std::function<int(double, int)> prob_dist,
        std::function<void(int *, int, int)> set_ran_indices)
      : name_{name}, p_occur_ptr_{p_occur_ptr}, n_avail_{n_avail},
        target_pool_{target_pool}, exe_{exe_funct}, prob_dist_{prob_dist},
        set_ran_indices_{set_ran_indices} {
    p_occur_ = *p_occur_ptr_;
  }
  // Constructor for poisson-based events
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
  // Sample appropriate stastistical distribution and return n_expected_
  int SampleStatistics() {
    n_opportunities_tot_ += *n_avail_;
    if (poisson_based_) {
      SamplePoissonStatistics();
    } else {
      p_occur_ = *p_occur_ptr_;
      n_expected_ = prob_dist_(p_occur_, *n_avail_);
    }
    SetTargets();
    return n_expected_;
  }
  // Remove given input target from targets_ array
  void RemoveTarget(ENTRY_T tar) {
    for (int i_entry{0}; i_entry < n_expected_; i_entry++) {
      if (tar == targets_[i_entry]) {
        targets_[i_entry] = targets_[n_expected_ - 1];
        n_expected_--;
        return;
      }
    }
  }
  // Execute event on last entry in targets_ (order is already randomized)
  void Execute() {
    exe_(targets_[n_expected_ - 1]);
    n_expected_--;
    n_executed_tot_++;
  }
};
#endif