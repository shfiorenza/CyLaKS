#ifndef _CYLAKS_RIGID_ROD_HPP_
#define _CYLAKS_RIGID_ROD_HPP_
#include "object.hpp"

class RigidRod : public Object {

protected:
  double length_{0.0}; // In nm
  Vec<double> gamma_;  // Indicies: 0->par, 1->perp, 2->rot
  Vec<double> sigma_;  // Indicies: 0->par, 1->perp, 2->rot

public:
  Vec<double> force_;       // In pN; zero'd out every timestep
  Vec<double> orientation_; // Unit vector

private:
  void SetParameters() {
    gamma_.resize(3);
    sigma_.resize(3);
    force_.resize(_n_dims_max);
    orientation_.resize(_n_dims_max);
  }

public:
  RigidRod() {}
  void Initialize(size_t sid, size_t id) {
    Object::Initialize(sid, id);
    SetParameters();
  }
  void UpdateForces();
};
#endif