#include "event.hpp"

void Event::SampleStatistics_Poisson() {

  if (*n_avail_ > poisson_.weights_.size()) {
    poisson_.weights_.resize(*n_avail_);
  }
  poisson_.weight_total_ = 0.0;
  for (int i_entry{0}; i_entry < *n_avail_; i_entry++) {
    poisson_.weights_[i_entry] =
        poisson_.get_weight_(target_pool_->at(i_entry));
    poisson_.weight_total_ += poisson_.weights_[i_entry];
  }
  n_expected_ = prob_dist_(poisson_.weight_total_ * p_occur_, 0);
  if (n_expected_ > 0) {
    SetTargets_Poisson();
  }
}

void Event::SetTargets_Poisson() {

  // Run through potential candidates and remove those with weight = 0.0
  for (int i_entry{0}; i_entry < *n_avail_; i_entry++) {
    if (poisson_.weights_[i_entry] == 0.0) {
      int i_last{--*n_avail_};
      if (*n_avail_ == 0) {
        n_expected_ = 0;
        return;
      }
      target_pool_->at(i_entry) = target_pool_->at(i_last);
      poisson_.weights_[i_entry] = poisson_.weights_[i_last];
      i_entry--;
    }
  }
  // Correct statistics if n_expected > n_avail after above removal
  if (n_expected_ > *n_avail_) {
    n_expected_ = *n_avail_;
  }
  // Scratch array that will hold selected targets before transfer
  Object *selected_candidates[n_expected_];
  // Select n_expected_ entries at random
  for (int i_set{0}; i_set < n_expected_; i_set++) {
    double p_cum{0.0};
    double ran{gsl_->GetRanProb()};
    for (int i_entry{0}; i_entry < *n_avail_; i_entry++) {
      p_cum += poisson_.weights_[i_entry] / poisson_.weight_total_;
      if (ran < p_cum) {
        // Add entry to selected candidates
        selected_candidates[i_set] = target_pool_->at(i_entry);
        // Remove selected candidate from target pool so it isn't selected again
        poisson_.weight_total_ -= poisson_.weights_[i_entry];
        int i_last{--*n_avail_};
        target_pool_->at(i_entry) = target_pool_->at(i_last);
        poisson_.weights_[i_entry] = poisson_.weights_[i_last];
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
