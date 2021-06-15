#include "cylaks/reservoir.hpp"
#include "cylaks/motor.hpp"
#include "cylaks/protein.hpp"

template class Reservoir<Motor>;
template class Reservoir<Protein>;

template <typename ENTRY_T>
void Reservoir<ENTRY_T>::GenerateEntries(size_t n_entries) {

  reservoir_.resize(n_entries);
  active_entries_.resize(n_entries);
  for (int i_entry{0}; i_entry < n_entries; i_entry++) {
    reservoir_[i_entry].Initialize(species_id_, Sys::n_objects_++);
  }
  r_min_ = reservoir_[0].spring_.r_min_;
  r_rest_ = reservoir_[0].spring_.r_rest_;
  r_max_ = reservoir_[0].spring_.r_max_;
}

template <typename ENTRY_T> void Reservoir<ENTRY_T>::SetParameters() {

  using namespace Params;
  if (step_active_ * dt < t_equil + t_run) {
    active_ = true;
  }
  if (dynamic_equil_window < 0.0 or !active_) {
    equilibrated_ = true;
  } else {
    size_t window_size{(size_t)std::round(dynamic_equil_window / dt)};
    n_bound_.resize(window_size);
  }
  if (active_ and species_id_ == _id_xlink and Filaments::count > 1) {
    crosslinking_active_ = true;
  }
  if (active_ and species_id_ == _id_motor and Motors::gaussian_range > 0) {
    lattice_coop_active_ = true;
  }
}

template <typename ENTRY_T> void Reservoir<ENTRY_T>::CheckEquilibration() {

  if (equilibrated_) {
    return;
  }
  if (Sys::i_step_ < Sys::n_steps_pre_equil_) {
    n_bound_avg_ += double(n_active_entries_) / Sys::n_steps_pre_equil_;
    return;
  }
  using namespace Sys;
  if (Sys::i_step_ == Sys::n_steps_pre_equil_) {
    Log("Species %zu pre-equililibration ended with n_bound_avg = %.3g\n",
        species_id_, n_bound_avg_);
    return;
  }
  size_t window_size{n_bound_.size()};
  size_t i_timestep{Sys::i_step_ % window_size};
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
  Log("Species %zu is still equilibrating ... (n_bound_avg = %.3g +/- %.2g)\n",
      species_id_, window_avg, window_sigma);
  double delta{std::fabs(n_bound_avg_ - window_avg)};
  double delta_sigma{sqrt(n_bound_var_ + window_var)};
  if (delta < delta_sigma or delta == n_bound_avg_) {
    equilibrated_ = true;
    Log("Species %zu equilibration is complete (delta = % .2g + / - % .2g)\n",
        species_id_, delta, delta_sigma);
  }
  n_bound_avg_ = window_avg;
  n_bound_var_ = window_var;
}

template <typename ENTRY_T> void Reservoir<ENTRY_T>::SortPopulations() {

  if (up_to_date_) {
    return;
  }
  up_to_date_ = true;
  for (auto &&pop : sorted_) {
    pop.second.ZeroOut();
  }
  for (int i_entry{0}; i_entry < n_active_entries_; i_entry++) {
    ENTRY_T *entry{active_entries_[i_entry]};
    Sys::Log(1, " entry no %i (ID %i)\n", i_entry, entry->GetID());
    // printf(" %i HEADS\n", entry->GetNumHeadsActive());
    // entry->UpdateExtension();
    for (auto &&pop : sorted_) {
      Sys::Log(1, "  sorting into %s\n", pop.second.name_.c_str());
      pop.second.Sort(entry);
    }
  }
}
