#include "binding_head.hpp"
#include "binding_site.hpp"
#include "protein.hpp"

int BindingHead::GetDirectionTowardRest() {
  return parent_->GetDirectionTowardRest(this);
}

int BindingHead::GetNumNeighborsOccupied() {
  return site_->GetNumNeighborsOccupied();
}

int BindingHead::GetNumHeadsActive() { return parent_->n_heads_active_; }

void BindingHead::AddForce(Vec<double> f) { site_->AddForce(f); }

void BindingHead::UntetherSatellite() { parent_->UntetherSatellite(); }

bool BindingHead::Unbind() { return parent_->Unbind(this); }

double BindingHead::GetWeight_Diffuse(int dir) {
  return parent_->GetWeight_Diffuse(this, dir);
}

bool BindingHead::Diffuse(int dir) { return parent_->Diffuse(this, dir); }
