#include "event.hpp"
#include "associated_protein.h"
#include "kinesin.h"
#include "tubulin.h"

template class Event<
    std::variant<Tubulin *, AssociatedProtein::Monomer *, Kinesin::Monomer *>>;

template <typename ENTRY_T> void Event<ENTRY_T>::SampleStatistics_Poisson() {

  if (*n_avail_ > poisson_weights_.size()) {
    poisson_weights_.resize(*n_avail_);
  }
  poisson_weight_total_ = 0.0;
  for (int i_entry{0}; i_entry < *n_avail_; i_entry++) {
    poisson_weights_[i_entry] = get_weight_(target_pool_->at(i_entry));
    poisson_weight_total_ += poisson_weights_[i_entry];
  }
  n_expected_ = prob_dist_(poisson_weight_total_ * p_occur_, 0);
  if (n_expected_ > 0) {
    SetTargets_Poisson();
  }
}

template <typename ENTRY_T> void Event<ENTRY_T>::SetTargets_Poisson() {

  // Run through potential candidates and remove those with weight = 0.0
  for (int i_entry{0}; i_entry < *n_avail_; i_entry++) {
    if (poisson_weights_[i_entry] == 0.0) {
      int i_last{--*n_avail_};
      if (*n_avail_ == 0) {
        n_expected_ = 0;
        return;
      }
      target_pool_->at(i_entry) = target_pool_->at(i_last);
      poisson_weights_[i_entry] = poisson_weights_[i_last];
      i_entry--;
    }
  }
  // Correct statistics if n_expected > n_avail after above removal
  if (n_expected_ > *n_avail_) {
    n_expected_ = *n_avail_;
  }
  // Scratch array that will hold selected targets before transfer
  ENTRY_T selected_candidates[n_expected_];
  // Select n_expected_ entries at random
  for (int i_set{0}; i_set < n_expected_; i_set++) {
    double p_cum{0.0};
    double ran{get_ran_prob_()};
    for (int i_entry{0}; i_entry < *n_avail_; i_entry++) {
      p_cum += poisson_weights_[i_entry] / poisson_weight_total_;
      if (ran < p_cum) {
        // Add entry to selected candidates
        selected_candidates[i_set] = target_pool_->at(i_entry);
        // Remove selected candidate from target pool so it isn't selected again
        poisson_weight_total_ -= poisson_weights_[i_entry];
        int i_last{--*n_avail_};
        target_pool_->at(i_entry) = target_pool_->at(i_last);
        poisson_weights_[i_entry] = poisson_weights_[i_last];
        break;
      }
    }
  }
  // Set n_avail_ to n_expected_ so that SetTargets() only sees selected
  // candidates; will be updated next timestep regardless
  *n_avail_ = n_expected_;
  // Transfer selected candidates into target pool
  for (int i_entry{0}; i_entry < *n_avail_; i_entry++) {
    target_pool_->at(i_entry) = selected_candidates[i_entry];
  }
}