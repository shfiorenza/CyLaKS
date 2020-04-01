#ifndef _MICROTUBULE_H
#define _MICROTUBULE_H
#include "tubulin.h"
#include <vector>
class Curator;
struct system_parameters;
struct system_properties;

class Microtubule {
private:
  double applied_force_{0.0};
  Curator *wally_{nullptr};
  system_parameters *parameters_{nullptr};
  system_properties *properties_{nullptr};

public:
  int index_;                  // Index of MT in mt_list
  int plus_end_;               // Index of plus_end in lattice
  int minus_end_;              // Index of minus_end in lattice
  int delta_x_;                // Direction motors step (+/- 1)
  int n_sites_;                // Number of tubulin sites on MT
  double coord_;               // Absolute coordinate of left-most edge of MT
  double immobile_until_{0.0}; // Time at which MT can move freely

  // Drag coefficient
  double gamma_{0.0}; // in units of (pN*s) / nm
  // Sigma of Gaussian noise
  double sigma_{0.0};

  Microtubule *neighbor_ = nullptr; // FIXME for 1+ neighbors

  std::vector<Tubulin> lattice_; // All tubulin sites

private:
  void SetParameters();
  void GenerateLattice();

public:
  Microtubule();
  void Initialize(system_parameters *parameters, system_properties *properties,
                  int i_mt);

  double GetNetForce();
  double GetNetForce_Motors();
  double GetNetForce_Xlinks();
};
#endif
