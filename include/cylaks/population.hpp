#ifndef _CYLAKS_POPULATION_HPP_
#define _CYLAKS_POPULATION_HPP_
#include "system_definitions.hpp"
#include "system_namespace.hpp"

// Population: Used to divide active proteins into different bins
//             Stores pointers to proteins that fit given critera
//             E.g., motors with 1 neighbor that can unbind
template <typename ENTRY_T> struct Population {
protected:
  // Whether this pop. is 1-D (a serial list) or multi-dim (a sorted array)
  bool one_d_{true};
  // Sorting function that returns members of this pop. when given a protein
  // Will return an empty vector if no protein components are valid members
  Fn<Vec<ENTRY_T *>(ENTRY_T *)> get_members_;

  // (M-D use) Index to start at in each dimension, often but not always 0
  Vec<int> min_indices_;
  // (M-D use) Returns apt. bin indices <i,j,k> when given a member of the pop.
  Fn<Vec<int>(ENTRY_T *)> get_bin_indices_;

  // (1-D use) Add member to population and increase population size
  void AddEntry(ENTRY_T *entry) { entries_[size_++] = entry; }
  // (M-D use) Add member to apt pop. bin at [i][j][k] and increase bin size
  // Indices are input via vector in format <k>, <j, k>, or <i, j, k>
  void AddEntry(ENTRY_T *entry, Vec<int> indices) {
    int k{indices[0]};
    int j{indices.size() > 1 ? indices[1] : 0};
    int i{indices.size() > 2 ? indices[2] : 0};
    Sys::Log(3, "entry added w/ ijk = %i%i%i\n", i, j, k);
    bin_entries_[i][j][k][bin_size_[i][j][k]++] = entry;
    if (entry->GetNumHeadsActive() == 2) {
      bin_entries_[i][j][k][bin_size_[i][j][k]++] = entry->GetOtherHead();
    }
    Sys::Log(3, "bin size = %i\n", bin_size_[i][j][k]);
  }

public:
  Str name_;       // name of this population, e.g. "bound_i" for singly bound
  size_t size_{0}; // current size of population; 0 when none are active
  Vec<ENTRY_T *> entries_;       // (1-D use) list of active pop entries
  Vec3D<size_t> bin_size_;       // (M-D use) pop size; up to 3 dims.
  Vec4D<ENTRY_T *> bin_entries_; // (M-D use) pop entries; 1 list for each dim.
  Population() {}
  // Constructor for simple 1-D serial population list
  Population(Str name, Fn<Vec<ENTRY_T *>(ENTRY_T *)> getmems, size_t size_ceil)
      : name_{name}, get_members_{getmems} {
    entries_.resize(size_ceil);
  }
  // Constructor for multi-dim. (up to 3-D) sorted population array
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
  // Reset all population sizes to zero
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
  // Sort potential member into either 1-D list or M-D array
  void Sort(ENTRY_T *entry) {
    Vec<ENTRY_T *> members{get_members_(entry)};
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