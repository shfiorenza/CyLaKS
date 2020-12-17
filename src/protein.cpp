#include "protein.hpp"
#include "binding_site.hpp"
#include "system_namespace.hpp"

void Protein::InitializeNeighborList() {}

bool Protein::HasSatellite() {}

void Protein::UntetherSatellite() {}

void Protein::UpdateNeighborList() {}

Object *Protein::GetWeightedNeighbor() {}

BindingHead *Protein::GetActiveHead() {

  assert(n_heads_active_ == 1);
  if (head_one_.site_ != nullptr) {
    return &head_one_;
  } else if (head_two_.site_ != nullptr) {
    return &head_two_;
  } else {
    Sys::ErrorExit("AssociatedProtein::GetActiveHead");
  }
}

bool Protein::Bind(BindingSite *site, BindingHead *head) {

  if (site->occupant_ != nullptr) {
    return false;
  }
  site->occupant_ = head;
  head->site_ = site;
  n_heads_active_++;
  return true;
}

bool Protein::Unbind(BindingHead *head) {

  BindingSite *site{head->site_};
  site->occupant_ = nullptr;
  head->site_ = nullptr;
  n_heads_active_--;
  return true;
}

bool Protein::Tether() {}

bool Protein::Untether() {}