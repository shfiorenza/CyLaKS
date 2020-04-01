#ifndef _TUBULIN_
#define _TUBULIN_
#include "associated_protein.h"
#include "kinesin.h"

class Microtubule;
struct system_parameters;
struct system_properties;

class Tubulin {
private:
  double site_size_{0.0};
  system_parameters *parameters_{nullptr};
  system_properties *properties_{nullptr};

public:
  int species_id_{0};
  int index_{-1}; // Index of tubulin site in MT lattice
  Microtubule *mt_{nullptr};

  bool occupied_{false};
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

  double GetCoord();

  int GetPRC1NeighborCount();
  int GetKif4ANeighborCount();
};
#endif