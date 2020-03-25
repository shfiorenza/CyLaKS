#ifndef _TUBULIN
#define _TUBULIN
#include "associated_protein.h"
#include "kinesin.h"

class Microtubule;
struct system_parameters;
struct system_properties;

class Tubulin {
private:
  system_parameters *parameters_{nullptr};
  system_properties *properties_{nullptr};

public:
  int index_; // Index of tubulin site in MT lattice
  int species_id_ = 0;

  bool occupied_{false};

  double weight_bind_{0.0};
  double weight_unbind_{0.0};

  Microtubule *mt_{nullptr};
  Kinesin::Monomer *motor_head_{nullptr};
  AssociatedProtein::Monomer *xlink_head_{nullptr};

private:
public:
  Tubulin();
  void Initialize(system_parameters *parameters, system_properties *properties,
                  Microtubule *mt, int i_site);

  // 'equilibrium' refers to that of the crosslinker itself
  // and the kinesin 4 tether that attaches to it; this checks
  // if both are in the same direction w.r.t. this site
  bool EquilibriumInSameDirection();

  int GetPRC1NeighborCount();
  int GetKif4ANeighborCount();

  void UpdateWeights_Kinesin();
};
#endif