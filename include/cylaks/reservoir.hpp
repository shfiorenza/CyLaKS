#ifndef _CYLAKS_RESERVOIR_HPP_
#define _CYLAKS_RESERVOIR_HPP_
#include "population.hpp"
#include "system_namespace.hpp"
#include "system_parameters.hpp"
#include "system_rng.hpp"

class Object;

// Reservoir: Holds and organizes all proteins that will ever exist in sim.
//            Creates N objects, where N is the sum of binding sites in sim.
//            Note: objects are not modeled when fully unbound
template <typename ENTRY_T> struct Reservoir {
protected:
  size_t species_id_;
  Vec<ENTRY_T> reservoir_;

  bool up_to_date_{false};

  double n_bound_avg_{0.0};
  double n_bound_var_{0.0};
  Vec<size_t> n_bound_;

public:
  Vec<ENTRY_T *> active_entries_;

  bool active_{false};
  bool equilibrated_{false};
  size_t step_active_{0};
  size_t n_active_entries_{0};

  bool lattice_coop_active_{false};
  bool tethering_active_{false};
  bool crosslinking_active_{false};

  Map<Str, Sys::ProbEntry> p_event_;
  Map<Str, Sys::BoltzmannFactor> weights_;
  Map<Str, Population<Object>> sorted_;

protected:
  void GenerateEntries(size_t n_entries);
  void SetParameters();
  void CheckEquilibration();
  void SortPopulations();

public:
  Reservoir() {}
  void Initialize(size_t sid, size_t size, size_t i_active) {
    species_id_ = sid;
    step_active_ = i_active;
    GenerateEntries(size);
    SetParameters();
  }
  void AddWeight(Str name, size_t size) {
    weights_.emplace(name, Sys::BoltzmannFactor(name, size));
  }
  void AddProb(Str name, double val) {
    p_event_.emplace(name, Sys::ProbEntry(name, val));
  }
  void AddProb(Str name, double val, Str wt_name, size_t mode) {
    assert(mode == 0 or mode == 1);
    size_t size{weights_.at(wt_name).size_};
    // FIXME verify this padding is correct
    Vec3D<double> vals(1, Vec2D<double>(1, Vec<double>(size, val)));
    // Vec<double> vals(weights_.at(wt_name).size_, val);
    for (int i{weights_.at(wt_name).i_start_}; i < size; i++) {
      if (mode == 0) {
        vals[0][0][i] *= weights_.at(wt_name).bind_[i];
      } else {
        vals[0][0][i] *= weights_.at(wt_name).unbind_[i];
      }
    }
    p_event_.emplace(name, Sys::ProbEntry(name, vals));
  }
  void AddProb(Str name, Vec3D<double> vals) {
    p_event_.emplace(name, Sys::ProbEntry(name, vals));
  }
  void AddPop(Str name, Fn<Vec<Object *>(Object *)> sort) {
    sorted_.emplace(name, Population<Object>(name, sort, reservoir_.size()));
  }
  void AddPop(Str name, Fn<Vec<Object *>(Object *)> sort, Vec<size_t> dimsize,
              Vec<int> i_min, Fn<Vec<int>(Object *)> get_i) {
    Vec<size_t> sz{dimsize[0], dimsize[1], dimsize[2], reservoir_.size()};
    sorted_.emplace(name, Population<Object>(name, sort, sz, i_min, get_i));
  }
  ENTRY_T *GetFreeEntry() {
    size_t i_entry = SysRNG::GetRanInt(reservoir_.size());
    ENTRY_T *entry = &reservoir_[i_entry];
    size_t n_attempts{0};
    while (entry->GetNumHeadsActive() > 0 or entry->IsTethered()) {
      if (n_attempts > reservoir_.size()) {
        return nullptr;
      }
      i_entry++;
      if (i_entry == reservoir_.size()) {
        i_entry = 0;
      }
      entry = &reservoir_[i_entry];
      n_attempts++;
    }
    return entry;
  }
  void AddToActive(ENTRY_T *entry) {
    entry->active_index_ = n_active_entries_;
    active_entries_[n_active_entries_++] = entry;
    FlagForUpdate();
  }
  void RemoveFromActive(ENTRY_T *entry) {
    size_t i_entry{entry->active_index_};
    active_entries_[i_entry] = active_entries_[--n_active_entries_];
    active_entries_[i_entry]->active_index_ = i_entry;
    FlagForUpdate();
  }
  void UpdateLatticeDeformation() {
    if (!lattice_coop_active_) {
      return;
    }
    for (int i_entry{0}; i_entry < n_active_entries_; i_entry++) {
      active_entries_[i_entry]->ApplyLatticeDeformation();
    }
  }
  bool UpdateExtensions() {
    bool force_unbind_occurred{false};
    for (int i_active{0}; i_active < n_active_entries_; i_active++) {
      bool successful{active_entries_[i_active]->UpdateExtension()};
      if (!successful) {
        force_unbind_occurred = true;
      }
    }
    return force_unbind_occurred;
  }
  void FlagForUpdate() { up_to_date_ = false; }
  void PrepForKMC() {
    if (Sys::i_step_ < step_active_) {
      return;
    }
    CheckEquilibration();
    SortPopulations();
  }
};
#endif