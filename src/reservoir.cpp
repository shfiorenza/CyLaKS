#include "reservoir.hpp"
#include "curator.hpp"
#include "motor.hpp"
#include "protein.hpp"

template class Reservoir<Protein>;
template class Reservoir<Motor>;

template <typename ENTRY_T> void Reservoir<ENTRY_T>::SetParameters() {

  if (step_active_ * params_->dt < params_->t_equil + params_->t_run) {
    active_ = true;
  } else {
    return;
  }
  if (!params_->dynamic_equil) {
    equilibrated_ = true;
  }
  size_t window_size{
      (size_t)std::round(params_->dynamic_equil_window / params_->dt)};
  n_bound_.resize(window_size);
}

template <typename ENTRY_T> void Reservoir<ENTRY_T>::CheckEquilibration() {

  if (equilibrated_) {
    return;
  }
  if (wally_->i_step_ < wally_->n_steps_pre_equil_) {
    n_bound_avg_ += double(n_active_entries_) / wally_->n_steps_pre_equil_;
    return;
  }
  if (wally_->i_step_ == wally_->n_steps_pre_equil_) {
    wally_->Log("Run '%s' species %zu ", wally_->sim_name_, species_id_);
    wally_->Log("pre-equililibration ended with n_bound_avg = %.3g\n",
                n_bound_avg_);
    return;
  }
  size_t window_size{n_bound_.size()};
  size_t i_timestep{wally_->i_step_ % window_size};
  if (i_timestep != 0) {
    n_bound_[i_timestep] = n_active_entries_;
    return;
  }
  double window_avg{0.0};
  for (int i_timestep{0}; i_timestep < window_size; i_timestep++) {
    window_avg += double(n_bound_[i_timestep]) / window_size;
  }
  double window_var{0.0};
  for (int i_timestep{0}; i_timestep < window_size; i_timestep++) {
    double diff{window_avg - n_bound_[i_timestep]};
    window_var += diff * diff / (window_size - 1);
  }
  double window_sigma{sqrt(window_var)};
  wally_->Log("Run '%s' species %zu is still ", wally_->sim_name_, species_id_);
  wally_->Log("equilibrating ... (n_bound_avg = %.3g +/- %.1g)\n", window_avg,
              window_sigma);
  double delta{std::fabs(n_bound_avg_ - window_avg)};
  double delta_sigma{sqrt(n_bound_var_ + window_var)};
  if (delta < delta_sigma or delta == n_bound_avg_) {
    equilibrated_ = true;
    wally_->Log("Run '%s' species %zu ", wally_->sim_name_, species_id_);
    wally_->Log("equilibration is complete (delta = %.2g | sigma = %.2g)\n",
                delta, delta_sigma);
  }
  n_bound_avg_ = window_avg;
  n_bound_var_ = window_var;
}

template <typename ENTRY_T> void Reservoir<ENTRY_T>::SortPopulations() {

  if (up_to_date_ or wally_->i_step_ < step_active_) {
    return;
  }
  up_to_date_ = true;
  for (int i_entry{0}; i_entry < n_active_entries_; i_entry++) {
    ENTRY_T *entry{active_entries_[i_entry]};
    for (auto &&pop : sorted_) {
      pop.second.Sort(entry);
    }
  }
}

template <typename ENTRY_T> void Reservoir<ENTRY_T>::Update() {

  CheckEquilibration();
  SortPopulations();
}