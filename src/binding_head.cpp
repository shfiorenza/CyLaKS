#include "cylaks/binding_head.hpp"
#include "cylaks/binding_site.hpp"
#include "cylaks/protein.hpp"
#include "cylaks/protofilament.hpp"

int BindingHead::GetDirectionTowardRest() {
  return parent_->GetDirectionTowardRest(this);
}

int BindingHead::GetNumNeighborsOccupied() {
  return site_->GetNumNeighborsOccupied();
}

int BindingHead::GetNumNeighborsOccupied_Side() {
  return site_->GetNumNeighborsOccupied_Side();
}

int BindingHead::GetNumHeadsActive() { return parent_->n_heads_active_; }

Vec<double> BindingHead::GetBoundObjectOrientation() {
  return site_->filament_->GetPolarOrientation();
}

void BindingHead::AddForce(Vec<double> f) { site_->AddForce(f); }
void BindingHead::AddTorque(double tq) { site_->AddTorque(tq); }

bool BindingHead::UntetherSatellite() { return parent_->UntetherSatellite(); }

double BindingHead::GetWeight_Unbind_II() {
  return parent_->GetWeight_Unbind_II(this);
}

bool BindingHead::Unbind() { return parent_->Unbind(this); }

double BindingHead::GetWeight_Diffuse(int dir) {
  return parent_->GetWeight_Diffuse(this, dir);
}

bool BindingHead::Diffuse(int dir) { return parent_->Diffuse(this, dir); }

bool BindingHead::Diffuse_Side(int dir) {
  return parent_->Diffuse_Side(this, dir);
}