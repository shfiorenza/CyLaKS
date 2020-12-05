#ifndef _CYLAKS_RIGID_ROD_H_
#define _CYLAKS_RIGID_ROD_H_
#include "object.hpp"

class RigidRod : public Object {
protected:
  double length_{0.0}; // In microns
  Vec<double> force_;  // In pN; zero'd out every timestep
  Vec<double> orientation_;

public:
  RigidRod(size_t sid, size_t id, double length)
      : length_{length}, Object(sid, id) {
    force_.resize(Sys::_n_dims_max);
    orientation_.resize(Sys::_n_dims_max);
  }
  void UpdateForces();
};
#endif