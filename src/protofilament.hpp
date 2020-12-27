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
  double dt_eff_{0.0};
  double center_index_{0.0};

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

  void UpdateRodPosition();
  void UpdateSitePositions();

public:
  Protofilament() {}
  void Initialize(size_t sid, size_t id, size_t index) {
    RigidRod::Initialize(sid, id);
    index_ = index;
    SetParameters();
    GenerateSites();
    UpdateSitePositions();
  }
  void AddForce(BindingSite *location, Vec<double> f_applied) {
    if (Sys::i_step_ < immobile_until_) {
      return;
    }
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      // printf("added %g pN at site %i\n", f_applied[i_dim], location->index_);
      force_[i_dim] += f_applied[i_dim];
    }
    if (Params::Filaments::rotation_enabled) {
      // Construct vector pointing from rod COM to site
      Vec<double> r(_n_dims_max, 0.0);
      for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
        r[i_dim] = location->pos_[i_dim] - pos_[i_dim];
      }
      torque_ += Cross(r, f_applied);
    }
  }
  BindingSite *GetNeighb(BindingSite *site, int delta);
  void UpdatePosition() {
    if (Sys::i_step_ < immobile_until_) {
      return;
    }
    UpdateRodPosition();
    UpdateSitePositions();
  }
};
#endif