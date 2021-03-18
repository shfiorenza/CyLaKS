#ifndef _CYLAKS_POPULATION_HPP_
#define _CYLAKS_POPULATION_HPP_
#include "definitions.hpp"
#include "system_namespace.hpp"

template <typename ENTRY_T> struct Population {
private:
  bool one_d_{true};
  // 1-d stuff
  Fn<Vec<ENTRY_T *>(ENTRY_T *)> get_members_;
  // multi-dim stuff
  Vec<int> min_indices_;
  Fn<Vec<int>(ENTRY_T *)> get_bin_indices_;

  void AddEntry(ENTRY_T *entry) { entries_[size_++] = entry; }
  void AddEntry(ENTRY_T *entry, Vec<int> indices) {
    int k{indices[0]};
    int j{indices.size() > 1 ? indices[1] : 0};
    int i{indices.size() > 2 ? indices[2] : 0};
    Sys::Log(2, "entry added w/ ijk = %i%i%i\n", i, j, k);
    bin_entries_[i][j][k][bin_size_[i][j][k]++] = entry;
    if (entry->GetNumHeadsActive() == 2) {
      bin_entries_[i][j][k][bin_size_[i][j][k]++] = entry->GetOtherHead();
    }
    Sys::Log(2, "bin size = %i\n", bin_size_[i][j][k]);
  }

public:
  Str name_;
  size_t size_{0};
  Vec<ENTRY_T *> entries_;
  Vec3D<size_t> bin_size_;       // [n_neighbs][x_dub][x]
  Vec4D<ENTRY_T *> bin_entries_; // [n_neighbs][x_dub][x][i]
  Population() {}
  Population(Str name, Fn<Vec<ENTRY_T *>(ENTRY_T *)> getmems, size_t size_ceil)
      : name_{name}, get_members_{getmems} {
    entries_.resize(size_ceil);
  }
  Population(Str name, Fn<Vec<ENTRY_T *>(ENTRY_T *)> getmems,
             Vec<size_t> size_ceil, Vec<int> i_min,
             Fn<Vec<int>(ENTRY_T *)> getindices)
      : name_{name}, get_members_{getmems}, min_indices_{i_min},
        get_bin_indices_{getindices}, one_d_{false} {
    assert(size_ceil.size() == 4);
    assert(min_indices_.size() == 3);
    bin_size_.resize(size_ceil[0]);
    bin_entries_.resize(size_ceil[0]);
    for (int i{0}; i < bin_entries_.size(); i++) {
      bin_size_[i].resize(size_ceil[1]);
      bin_entries_[i].resize(size_ceil[1]);
      for (int j{0}; j < bin_entries_[i].size(); j++) {
        bin_size_[i][j].resize(size_ceil[2]);
        bin_entries_[i][j].resize(size_ceil[2]);
        for (int k{0}; k < bin_entries_[i][j].size(); k++) {
          bin_size_[i][j][k] = 0;
          bin_entries_[i][j][k].resize(size_ceil[3]);
        }
      }
    }
  }
  void ZeroOut() {
    if (one_d_) {
      size_ = 0;
    } else {
      for (int i{min_indices_[0]}; i < bin_size_.size(); i++) {
        for (int j{min_indices_[1]}; j < bin_size_[i].size(); j++) {
          for (int k{min_indices_[2]}; k < bin_size_[i][j].size(); k++) {
            bin_size_[i][j][k] = 0;
          }
        }
      }
    }
  }
  void Sort(ENTRY_T *entry) {
    Vec<ENTRY_T *> members{get_members_(entry)};
    // printf("%i MEMBERS\n", members.size());
    for (auto const &member : members) {
      if (one_d_) {
        AddEntry(member);
      } else {
        AddEntry(member, get_bin_indices_(member));
      }
    }
  }
};
#endif