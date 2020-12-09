#ifndef _CYLAKS_RIGID_ROD_HPP_
#define _CYLAKS_RIGID_ROD_HPP_
#include "object.hpp"

class RigidRod : public Object {
protected:
  double length_{0.0}; // In microns
  Vec<double> force_;  // In pN; zero'd out every timestep
  Vec<double> orientation_;

public:
  RigidRod(size_t sid, size_t id, double length)
      : Object(sid, id), length_{length} {
    force_.resize(_n_dims_max);
    orientation_.resize(_n_dims_max);
  }
  void UpdateForces();
};
#endif