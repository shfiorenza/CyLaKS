#ifndef _CYLAKS_PROTOFILAMENT_HPP_
#define _CYLAKS_PROTOFILAMENT_HPP_
#include "binding_site.hpp"
#include "rigid_rod.hpp"
#include "system_namespace.hpp"
#include "system_parameters.hpp"
#include "system_rng.hpp"

class Protofilament : public RigidRod {
private:
  size_t polarity_{0};
  double center_index_{0.0};

  SysRNG *gsl_{nullptr};

public:
  size_t index_{0};
  size_t immobile_until_{0};

  int dx_{0}; // Towards plus end
  Vec<BindingSite> sites_;

  BindingSite *plus_end_{nullptr};
  BindingSite *minus_end_{nullptr};
  Protofilament *neighbor_{nullptr};

private:
  void SetParameters();
  void GenerateSites();

  double GetNetForce();

public:
  Protofilament() {}
  void Initialize(size_t sid, size_t id, size_t index) {
    RigidRod::Initialize(sid, id);
    index_ = index;
    SetParameters();
    GenerateSites();
    UpdateSitePositions();
    Sys::Log("plus-end: (%g, %g)\n", plus_end_->pos_[0], plus_end_->pos_[1]);
    Sys::Log("minus_end: (%g, %g)\n", minus_end_->pos_[0], minus_end_->pos_[1]);
  }
  void UpdateSitePositions();
  void UpdatePosition(double dt) {
    // FIXME only for x-dim as of now
    // double velocity{GetNetForce() / gamma_par_};
    // double noise{gsl_->GetGaussianNoise(sigma_par_)};
    // pos_[0] += velocity * dt + noise;
  }
};
#endif