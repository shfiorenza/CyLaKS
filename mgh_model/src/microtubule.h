#ifndef _MICROTUBULE_H
#define _MICROTUBULE_H
#include "tubulin.h"
#include <vector>
class Curator;
struct system_parameters;
struct system_properties;

class Microtubule {
private:
  system_parameters *parameters_{nullptr};
  system_properties *properties_{nullptr};
  Curator *wally_{nullptr};

public:
  int index_;     // Index of MT in mt_list
  int polarity_;  // 0 for plus end on right; 1 for left
  int plus_end_;  // Index of plus_end in lattice
  int minus_end_; // Index of minus_end in lattice
  int delta_x_;   // Direction motors step (+/- 1)

  int mt_index_adj_; // Index of adjacent microtubule

  int n_sites_; // Number of tubulin sites on MT
  int coord_;   // Absolute coordinate of left-most edge of MT

  int steps_per_iteration_{0};
  double gamma_; // in units of (pN*s) / nm

  Microtubule *neighbor_ = nullptr; // FIXME for 1+ neighbors

  std::vector<Tubulin> lattice_; // All tubulin sites

private:
  void SetParameters();
  void GenerateLattice();

public:
  Microtubule();
  void Initialize(system_parameters *parameters, system_properties *properties,
                  int i_mt);

  void UpdateExtensions();

  double GetNetForce();
  double GetNetForce_Motors();
  double GetNetForce_Xlinks();
};
#endif
