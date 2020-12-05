#ifndef _CYLAKS_POPULATION_HPP_
#define _CYLAKS_POPULATION_HPP_
#include "system_definitions.hpp"

template <typename ENTRY_T> struct Population {
  Str name_;
  bool one_d_{true};
  Vec<int> min_indices_;
  size_t size_;
  Vec3D<size_t> bin_size_; // [n_neighbs][x_dub][x]
  Vec<ENTRY_T *> entries_;
  Vec4D<ENTRY_T *> bin_entries_; // [n_neighbs][x_dub][x][i]
  Fn<void(ENTRY_T *)> sort_;
  Population() {}
  Population(Str name, Vec3D<size_t> sizes, Vec4D<ENTRY_T *> entries,
             Fn<void(ENTRY_T *)> sort, Vec<double> i_min)
      : name_{name}, bin_size_{sizes}, bin_entries_{entries}, sort_{sort},
        min_indices_{i_min}, one_d_{false} {}
  Population(Str name, size_t size, Vec<ENTRY_T *> entries,
             Fn<void(ENTRY_T *)> sort)
      : name_{name}, size_{size}, entries_{entries}, sort_{sort} {}
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
  void AddEntry(ENTRY_T *entry) { entries_[size_++] = entry; }
  void AddEntry(ENTRY_T *entry, size_t i) {
    bin_entries_[0][0][i][bin_size_[0][0][i]++] = entry;
  }
  void Sort(ENTRY_T *entry) { sort_(entry); }
};
#endif