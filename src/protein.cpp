#include "protein.hpp"
#include "binding_site.hpp"

void Protein::InitializeNeighborList() {}

bool Protein::HasSatellite() {}

void Protein::UntetherSatellite() {}

void Protein::UpdateNeighborList() {}

Object *Protein::GetWeightedNeighbor() {}

BindingHead *Protein::GetActiveHead() {}

bool Protein::Bind(BindingSite *site, BindingHead *head) {

  if (site->occupant_ != nullptr) {
    return false;
  }
  site->occupant_ = head;
  head->site_ = site;
  n_heads_active_++;
}

bool Protein::Unbind(BindingHead *head) {

  BindingSite *site{head->site_};
  site->occupant_ = nullptr;
  head->site_ = nullptr;
  n_heads_active_--;
}

bool Protein::Tether() {}

bool Protein::Untether() {}