#include "cylaks/catalytic_head.hpp"
#include "cylaks/binding_site.hpp"
#include "cylaks/motor.hpp"
#include "cylaks/protofilament.hpp"

#include "cylaks/protein.hpp"

void CatalyticHead::Initialize(size_t sid, size_t id, double radius,
                               Motor *parent_ptr,
                               CatalyticHead *other_head_ptr) {
  BindingHead::Initialize(sid, id, radius);
  parent_ = parent_ptr;
  BindingHead::parent_ = dynamic_cast<Protein *>(parent_);
  other_head_ = other_head_ptr;
  BindingHead::other_head_ = dynamic_cast<BindingHead *>(other_head_);
}

int CatalyticHead::GetNumNeighborsOccupied_Motor() {
  return site_->GetNumNeighborsOccupied_Motor();
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
