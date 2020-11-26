#include "protofilament.hpp"
#include "master_header.hpp"

void Protofilament::GenerateLattice() {

  for (int i_site{0}; i_site < length_; i_site++) {
    lattice_.emplace_back(_id_tubulin, wally_->n_sim_objs_++, _r_tubulin, this);
  }
}
