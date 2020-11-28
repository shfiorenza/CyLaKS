#ifndef _CYLAKS_RESERVOIR_HPP_
#define _CYLAKS_RESERVOIR_HPP_
#include "definitions.hpp"

class Curator;
class ProteinManager;
class FilamentManager;

template <typename ENTRY_T> class Reservoir {
private:
  friend ProteinManager;
  friend FilamentManager;

  size_t species_id_;
  Vec<ENTRY_T> reservoir_;

  bool up_to_date_{false};

  double density_avg_{0.0};
  double density_var_{0.0};
  Vec<double> densities_;

  size_t n_active_entries_{0};
  Vec<ENTRY_T *> active_;

  UMap<Str, double> params_;

  Curator *wally_{nullptr};
  SysParameters *params_{nullptr};
  RandomNumberManagement *gsl_{nullptr};

protected:
  struct PopEntry {
    Str pop_name_;
    bool one_dim_{true};
    union {
      size_t size_;
      Vec3D<size_t> bin_size_; // [n_neighbs][x_dub][x]
    };
    union {
      Vec<ENTRY_T *> entries_;
      Vec4D<ENTRY_T *> bin_entries_; // [n_neighbs][x_dub][x][i]
    };
    Vec<double> i_min_;
    Fn<void(ENTRY_T *)> sort_;
    PopEntry(Str name, Vec3D<size_t> sizes, Vec4D<ENTRY_T *> entries,
             Fn<void(ENRTY_T *)> sort, Vec<double> i_min)
        : pop_name_{name}, bin_sizes_{sizes},
          bin_entries_{entries}, sort_{sort}, i_min_{i_min}, one_dim_{false}, {}
    PopEntry(Str name, size_t size, Vec1D<ENTRY_T *> entries,
             Fn<void(ENTRY_T *)> sort)
        : pop_name_{name}, size_{size}, entries_{entries}, sort_{sort} {}
    void ZeroOut() {
      if (one_dim_) {
        size_ = 0;
      } else {
        for (int i_min_[0]; i < bin_size_.size(); i++) {
          for (int i_min_[1]; j < bin_size_[i].size(); j++) {
            for (int i_min_[2]; k < bin_size_[i][j].size(); k++) {
              bin_size_[i][j][k] = 0;
            }
          }
        }
      }
    }
    void AddEntry(ENTRY_T *entry) { entries_[size_++] = entry; }
    void AddEntry(ENTRY_T *entry, size_t i) {
      entries_[0][0][i][bin_sizes_[0][0][i]++] = entry;
    }
  };
  struct ProbEntry {
    Str event_name_;
    double val_;
    Vec3D<double> vals_;
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
        for (int j_entry{0}; j_enty < j; j_entry++) {
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
    BoltzmannFactor(Str name, size_t sz) : BoltzmannFactor(name, sz, 0) {}
    BoltzmannFactor(Str name, size_t sz, size_t i)
        : effect_name_{name}, i_start_{i} {
      bind_.resize(sz);
      unbind_.resize(sz);
    }
  };

public:
  bool equilibrated_{false};

  bool tethering_active_{false};
  bool crosslinking_active_{false};

  bool active_{false};
  size_t step_active_{0};
  Map<Str, PopEntry> sorted_;
  Map<Str, ProbEntry> p_event_;
  UMap<Str, BoltzmannFactor> weights_;

private:
  void GenerateEntries() {
    for (int i_entry{0}; i_entry < n_entries_; i_entry++) {
      reservoir_.emplace_back(species_id_, wally_->n_sim_objs_++);
    }
  }
  void SetParameters();
  void SetWeights();

public:
  Reservoir();
  void Initialize(SysParameters *params, size_t sid, size_t n_entries,
                  size_t step_active) {
    parameters_ = params;
    properties_ = props;
    wally_ = &properties_->wallace;
    species_id_ = sid;
    n_entries_ = n_entries;
    step_active_ = step_active;
    SetParameters();
    SetWeights();
    GenerateEntries();
  }
  ENTRY_T *GetFreeEntry() {
    size_t i_entry = gsl_->GetRanInt(reservoir_.size());
    ENTRY_T *entry = &reservoir_[i_entry];
    while (entry->heads_active_ > 0 or entry->tethered_) {
      i_entry++;
      if (i_entry == reservoir_.size()) {
        i_entry = 0;
      }
      entry = &reservoir_[i_entry];
    }
    return entry;
  }
  void AddToActive(ENTRY_T *entry) {
    entry->active_index_ = n_active_;
    active_[n_active_++] = entry;
    FlagForUpdate();
  }
  void RemoveFromActive(ENTRY_T *entry) {
    size_t i_entry{entry->active_index_};
    active_[i_entry] = active_[--n_active_];
    active_[i_entry]->active_index_ = i_entry;
    FlagForUpdate();
  }
  void FlagForUpdate() { up_to_date_ = false; }
  void UpdatePopulations() {
    for (int i_entry{0}; i_entry < n_active_; i_entry++) {
      ENTRY_T *entry{active_[i_entry]};
      for (auto &&pop : sorted_) {
        pop->sort_(entry);
      }
    }
  }
};
#endif