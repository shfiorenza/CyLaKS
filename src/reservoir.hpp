#ifndef _CYLAKS_RESERVOIR_HPP_
#define _CYLAKS_RESERVOIR_HPP_
#include "population.hpp"
#include "system_namespace.hpp"
#include "system_parameters.hpp"
#include "system_rng.hpp"

class Object;

template <typename ENTRY_T> struct Reservoir {
private:
  size_t species_id_;
  Vec<ENTRY_T> reservoir_;

  bool up_to_date_{false};

  double n_bound_avg_{0.0};
  double n_bound_var_{0.0};
  Vec<size_t> n_bound_;

  Vec<ENTRY_T *> active_entries_;

protected:
  struct ProbEntry {
    Str event_name_;
    double val_;
    Vec3D<double> vals_;
    ProbEntry() {}
    ProbEntry(Str name, double val) : event_name_{name}, val_{val} {}
    ProbEntry(Str name, Vec3D<double> vals) : event_name_{name}, vals_{vals} {}
    double GetVal() { return val_; }
    double GetVal(size_t i) { return vals_[0][0][i]; }
    double GetVal(size_t i, size_t j) { return vals_[0][i][j]; }
    double GetVal(size_t i, size_t j, size_t k) { return vals_[i][j][k]; }
  };
  struct BoltzmannFactor {
    Str effect_name_;
    int i_start_{0};
    size_t size_{0};
    Vec<double> bind_;
    Vec<double> unbind_;
    BoltzmannFactor() {}
    BoltzmannFactor(Str name, size_t sz) : BoltzmannFactor(name, sz, 0) {}
    BoltzmannFactor(Str name, size_t size, int i_start)
        : effect_name_{name}, size_{size}, i_start_{i_start} {
      bind_.resize(size);
      unbind_.resize(size);
    }
  };

public:
  size_t n_active_entries_{0};
  bool equilibrated_{false};

  bool active_{false};
  size_t step_active_{0};

  bool tethering_active_{false};
  bool crosslinking_active_{false};

  int comp_cutoff_{0};
  int ext_cutoff_{0};
  Vec<double> possible_extensions_;
  Vec<double> possible_cosines_;

  Map<Str, ProbEntry> p_event_;
  Map<Str, BoltzmannFactor> weights_;
  Map<Str, Population<Object>> sorted_;

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
  void AddWeight(Str name, size_t size) {
    weights_.emplace(name, BoltzmannFactor(name, size));
  }
  void AddProb(Str name, double val) {
    p_event_.emplace(name, ProbEntry(name, val));
  }
  void AddProb(Str name, double val, Str wt_name, size_t mode) {
    assert(mode == 0 or mode == 1);
    size_t size{weights_.at(wt_name).size_};
    Vec3D<double> vals(1, Vec2D<double>(1, Vec<double>(size, val)));
    // Vec<double> vals(weights_.at(wt_name).size_, val);
    for (int i{weights_.at(wt_name).i_start_}; i < size; i++) {
      if (mode == 0) {
        vals[0][0][i] *= weights_.at(wt_name).bind_[i];
      } else {
        vals[0][0][i] *= weights_.at(wt_name).unbind_[i];
      }
    }
    p_event_.emplace(name, ProbEntry(name, vals));
  }
  void AddProb(Str name, Vec3D<double> vals) {
    p_event_.emplace(name, ProbEntry(name, vals));
  }
  void AddPop(Str name, Fn<Object *(Object *)> sort) {
    sorted_.emplace(name, Population<Object>(name, sort, reservoir_.size()));
  }
  void AddPop(Str name, Fn<Object *(Object *)> sort, Vec<size_t> dimsize,
              Vec<int> i_min, Fn<Vec<int>(Object *)> get_i) {
    Vec<size_t> sz{dimsize[0], dimsize[1], dimsize[2], reservoir_.size()};
    sorted_.emplace(name, Population<Object>(name, sort, sz, i_min, get_i));
  }
  ENTRY_T *GetFreeEntry() {
    size_t i_entry = SysRNG::GetRanInt(reservoir_.size());
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