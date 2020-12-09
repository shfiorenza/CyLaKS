#include "protofilament.hpp"
#include "curator.hpp"

void Protofilament::GenerateSites() {

  sites_.emplace_back(_id_site, Sys::n_unique_objects_++, _r_site, this);
  printf("lmfao %zu & %g\n", _id_site, _r_site);
}