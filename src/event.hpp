#ifndef _CYLAKS_EVENT_HPP_
#define _CYLAKS_EVENT_HPP_
#include "system_definitions.hpp"
#include "system_rng.hpp"

class Object;

class Event {
private:
  enum Distribution { Binomial, Poisson };
  size_t n_executed_tot_{0};       // # of times event has been executed
  size_t n_opportunities_tot_{0};  // # of opportunities event had to execute
  Distribution mode_{Binomial};    // Which distribution we sample from
  Vec<Object *> *target_pool_;     // Ptr to list of available targets; dynamic
  Fn<void(Object *)> exe_;         // Function that actually executes this event
  Fn<int(double, int)> prob_dist_; // Sampled to predict n_events each timestep
  struct PoissonToolbox {
    double weight_total_;
    Vec<double> weights_;
    Fn<double(Object *)> get_weight_;
  };
  PoissonToolbox poisson_; // Auxiliary resources for poisson mode
  SysRNG *gsl_;            // Pointer to RNG management structure

public:
  Str name_{"bruh"};         // Name of this event, e.g., "Bind_II_Teth"
  double p_occur_{0.0};      // Probability that event will occur each timestep
  size_t n_expected_{0};     // Expected # of events to occur any given timestep
  size_t *n_avail_{nullptr}; // Ptr to # of targets event can act on; dynamic
  Vec<Object *> targets_;    // Objects this event will act on this timestep

private:
  void SetTargets() {
    if (n_expected_ > targets_.size()) {
      targets_.resize(n_expected_);
    }
    int indices[n_expected_];
    gsl_->SetRanIndices(indices, n_expected_, *n_avail_);
    for (int i_entry{0}; i_entry < n_expected_; i_entry++) {
      targets_[i_entry] = target_pool_->at(indices[i_entry]);
    }
  }
  void SampleStatistics_Poisson();
  void SetTargets_Poisson();

public:
  Event(Str name, double p_occur, size_t *n_avail, Vec<Object *> *target_pool,
        Fn<int(double, int)> prob_dist, Fn<void(Object *)> exe, SysRNG *gsl)
      : gsl_{gsl}, name_{name}, p_occur_{p_occur}, n_avail_{n_avail},
        target_pool_{target_pool}, prob_dist_{prob_dist}, exe_{exe} {}
  Event(Str name, double p_occur, size_t *n_avail, Vec<Object *> *target_pool,
        Fn<int(double, int)> prob_dist, Fn<double(Object *)> weight_fn,
        Fn<void(Object *)> exe, SysRNG *gsl)
      : Event(name, p_occur, n_avail, target_pool, prob_dist, exe, gsl) {
    poisson_.get_weight_ = weight_fn;
    mode_ = Poisson;
  }
  size_t SampleStatistics() {
    n_opportunities_tot_ += *n_avail_;
    if (mode_ == Poisson) {
      SampleStatistics_Poisson();
    } else {
      n_expected_ = prob_dist_(p_occur_, *n_avail_);
    }
    SetTargets();
    return n_expected_;
  }
  void RemoveTarget(Object *tar) {
    for (int i_entry{0}; i_entry < n_expected_; i_entry++) {
      if (tar == targets_[i_entry]) {
        targets_[i_entry] = targets_[--n_expected_];
        return;
      }
    }
  }
  void Execute() {
    exe_(targets_[--n_expected_]);
    n_executed_tot_++;
  }
};
#endif