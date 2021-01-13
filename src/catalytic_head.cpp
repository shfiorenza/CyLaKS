#include "catalytic_head.hpp"
#include "binding_site.hpp"
#include "motor.hpp"
#include "protofilament.hpp"

double CatalyticHead::GetWeight_Unbind_II() {
  return parent_->GetWeight_Unbind_II(this);
}

bool CatalyticHead::Unbind() { return parent_->Unbind(this); }
