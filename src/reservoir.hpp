#ifndef _CYLAKS_RESERVOIR_HPP_
#define _CYLAKS_RESERVOIR_HPP_
#include "population.hpp"
#include "system_namespace.hpp"
#include "system_parameters.hpp"
#include "system_rng.hpp"

template <typename ENTRY_T> class Reservoir {
private:
  size_t species_id_;
  Vec<ENTRY_T> reservoir_;

  bool up_to_date_{false};

  double n_bound_avg_{0.0};
  double n_bound_var_{0.0};
  Vec<size_t> n_bound_;

  size_t n_active_entries_{0};
  Vec<ENTRY_T *> active_entries_;

  UMap<Str, double> param_vals_;

  SysRNG *gsl_{nullptr};

protected:
  struct ProbEntry {
    Str event_name_;
    double val_;
    Vec3D<double> vals_;
    ProbEntry() {}
    double GetVal() { return val_; }
    double GetVal(size_t i) { return vals_[0][0][i]; }
    double GetVal(size_t i, size_t j) { return vals_[0][i][j]; }
    double GetVal(size_t i, size_t j, size_t k) { return vals_[i][j][k]; }
    void InitVals(size_t i) { InitVals(1, 1, i); }
    void InitVals(size_t i, size_t j) { InitVals(1, i, j); }
    void InitVals(size_t i, size_t j, size_t k) {
      vals_.resize(i);
      for (int i_entry{0}; i_entry < i; i_entry++) {
        vals_[i_entry].resize(j);
        for (int j_entry{0}; j_entry < j; j_entry++) {
          vals_[i_entry][j_entry].resize(k);
        }
      }
    }
  };
  struct BoltzmannFactor {
    Str effect_name_;
    size_t i_start_{0};
    Vec<double> bind_;
    Vec<double> unbind_;
    BoltzmannFactor() {}
    BoltzmannFactor(Str name, size_t sz) : BoltzmannFactor(name, sz, 0) {}
    BoltzmannFactor(Str name, size_t sz, size_t i)
        : effect_name_{name}, i_start_{i} {
      bind_.resize(sz);
      unbind_.resize(sz);
    }
  };

public:
  bool equilibrated_{false};

  bool active_{false};
  size_t step_active_{0};

  bool tethering_active_{false};
  bool crosslinking_active_{false};

  Map<Str, Population<ENTRY_T>> sorted_;
  Map<Str, ProbEntry> p_event_;
  UMap<Str, BoltzmannFactor> weights_;

private:
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
  ENTRY_T *GetFreeEntry() {
    size_t i_entry = gsl_->GetRanInt(reservoir_.size());
    ENTRY_T *entry = &reservoir_[i_entry];
    while (entry->n_heads_active_ > 0 or entry->tethered_) {
      i_entry++;
      if (i_entry == reservoir_.size()) {
        i_entry = 0;
      }
      entry = &reservoir_[i_entry];
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
  void AddWeight(Str name, size_t size) {
    weights_.emplace(name, BoltzmannFactor(name, size));
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