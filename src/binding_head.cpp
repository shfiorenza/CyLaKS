#include "binding_head.hpp"
#include "binding_site.hpp"
#include "protein.hpp"

int BindingHead::GetDirectionTowardRest() {

  size_t n_heads_active{parent_->n_heads_active_};
  //   printf("  %zu\n", n_heads_active);
  if (n_heads_active == 1) {
    return 1;
  } else if (n_heads_active == 2) {
  }
}

int BindingHead::GetNeighborCount() { return site_->GetNeighborCount(); }

int BindingHead::GetNumHeadsActive() { return parent_->n_heads_active_; }

void BindingHead::UntetherSatellite() { parent_->UntetherSatellite(); }

bool BindingHead::Unbind() { return parent_->Unbind(this); }

bool BindingHead::Diffuse(int dir) { return parent_->Diffuse(this, dir); }
