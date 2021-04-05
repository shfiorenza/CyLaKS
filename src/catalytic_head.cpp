#include "cylaks/catalytic_head.hpp"
#include "cylaks/binding_site.hpp"
#include "cylaks/motor.hpp"
#include "cylaks/protofilament.hpp"

int CatalyticHead::GetNumNeighborsOccupied() {
  return site_->GetNumNeighborsOccupied();
}

int CatalyticHead::GetNumHeadsActive() { return parent_->n_heads_active_; }

double CatalyticHead::GetWeight_Unbind_II() {
  return parent_->GetWeight_Unbind_II(this);
}

bool CatalyticHead::Unbind() { return parent_->Unbind(this); }

double CatalyticHead::GetWeight_Diffuse(int dir) {
  return parent_->GetWeight_Diffuse(this, dir);
}

bool CatalyticHead::Diffuse(int dir) { return parent_->Diffuse(this, dir); }
