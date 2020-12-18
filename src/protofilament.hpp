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