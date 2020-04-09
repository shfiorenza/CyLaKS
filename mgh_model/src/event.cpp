#include "event.hpp"
#include "associated_protein.h"
#include "kinesin.h"
#include "tubulin.h"

template class Event<
    std::variant<Tubulin *, AssociatedProtein::Monomer *, Kinesin::Monomer *>>;

template <typename ENTRY_T> void Event<ENTRY_T>::SamplePoissonStatistics() {

  // If we have more available targets than the list can hold, resize it
  if (*n_avail_ > poisson_weights_.size()) {
    poisson_weights_.resize(*n_avail_);
  }
  // Set total poisson weight to 0.0 before summing over all entries
  poisson_weight_total_ = 0.0;
  // Run through each entry of target_pool_
  for (int i_entry{0}; i_entry < *n_avail_; i_entry++) {
    // Get individual entry weight and store it in poisson_weights_ array
    poisson_weights_[i_entry] = get_weight_(target_pool_->at(i_entry));
    // Add individual entry weight to poisson_weight_total_
    poisson_weight_total_ += poisson_weights_[i_entry];
  }
  // Calculate average rate based on weights and sample from Poisson dist
  n_expected_ = prob_dist_(poisson_weight_total_ * p_occur_, 0);
  // If we expect at least one event, set target candidates based on weights
  if (n_expected_ > 0) {
    SetPoissonCandidates();
  }
}

template <typename ENTRY_T> void Event<ENTRY_T>::SetPoissonCandidates() {

  // Run through potential candidates and remove those with weight = 0.0
  for (int i_entry{0}; i_entry < *n_avail_; i_entry++) {
    if (poisson_weights_[i_entry] == 0.0) {
      // Swap entry (and respective weight) with last entry in list
      int i_last{--*n_avail_};
      // If this was our last entry, set n_expected_ to 0 and return
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
    // Roll a random probability
    double ran{get_ran_prob_()};
    // Initialize cumulative probability at 0.0
    double p_cum{0.0};
    for (int i_entry{0}; i_entry < *n_avail_; i_entry++) {
      // Probability for this entry is its own weight divided by total weight
      p_cum += poisson_weights_[i_entry] / poisson_weight_total_;
      // Entry for which p_cum > ran will be selected
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