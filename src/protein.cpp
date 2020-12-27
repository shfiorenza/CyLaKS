#include "protein.hpp"
#include "binding_head.hpp"
#include "protofilament.hpp"
#include "system_namespace.hpp"

void Protein::InitializeNeighborList() {}

bool Protein::HasSatellite() {}

void Protein::UntetherSatellite() {}

// void Protein::UpdateNeighborList() {}

// Object *Protein::GetWeightedNeighbor() {}

BindingHead *Protein::GetActiveHead() {

  assert(n_heads_active_ == 1);
  // printf("get_active: %i heads active\n", n_heads_active_);
  if (head_one_.site_ != nullptr) {
    return &head_one_;
  } else if (head_two_.site_ != nullptr) {
    return &head_two_;
  } else {
    Sys::ErrorExit("AssociatedProtein::GetActiveHead");
  }
}

void Protein::UpdateNeighbors_Bind_II() {

  n_neighbors_bind_ii_ = 0;
  BindingSite *site{GetActiveHead()->site_};
  Protofilament *neighb_fil{site->filament_->neighbor_};
  for (int delta{-dist_cutoff_}; delta <= dist_cutoff_; delta++) {
    BindingSite *neighb = neighb_fil->GetNeighb(site, delta);
    if (neighb == nullptr) {
      continue;
    }
    if (neighb->occupant_ == nullptr) {
      neighbors_bind_ii_[n_neighbors_bind_ii_++] = neighb;
    }
  }
}

double Protein::GetWeight_Bind_II(BindingSite *neighb) {

  double lambda{0.5};
  BindingSite *site{GetActiveHead()->site_};
  double r_x{neighb->pos_[0] - site->pos_[0]};
  double r_y{neighb->pos_[1] - site->pos_[1]};
  double r{sqrt(Square(r_x) + Square(r_y))};
  double dr{r - Params::Xlinks::r_0};
  double dE{0.5 * Params::Xlinks::k_spring * Square(dr)};
  return exp(-(1.0 - lambda) * dE / Params::kbT);
}

double Protein::GetTotalWeight_Bind_II() {

  double tot_weight{0.0};
  UpdateNeighbors_Bind_II();
  for (int i_neighb{0}; i_neighb < n_neighbors_bind_ii_; i_neighb++) {
    tot_weight += GetWeight_Bind_II(neighbors_bind_ii_[i_neighb]);
  }
  return tot_weight;
}

BindingSite *Protein::GetNeighbor_Bind_II() {

  double weight_tot{GetTotalWeight_Bind_II()};
  double ran{SysRNG::GetRanProb()};
  double p_cum{0.0};
  for (int i_neighb{0}; i_neighb < n_neighbors_bind_ii_; i_neighb++) {
    BindingSite *neighb{neighbors_bind_ii_[i_neighb]};
    p_cum += GetWeight_Bind_II(neighb) / weight_tot;
    if (ran < p_cum) {
      return neighb;
    }
  }
  return nullptr;
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

double Protein::GetWeight_Unbind_II(BindingHead *head) {

  if (n_heads_active_ != 2) {
    Sys::ErrorExit("wut\n");
  }

  double lambda_neighb{1.0};
  double lambda_spring{0.5};
  BindingSite *site{head->site_};
  BindingSite *other_site{head->GetOtherHead()->site_};
  double r_x{site->pos_[0] - other_site->pos_[0]};
  double r_y{site->pos_[1] - other_site->pos_[1]};
  double r{sqrt(Square(r_x) + Square(r_y))};
  // printf("r = %g\n", r);
  double dr{r - Params::Xlinks::r_0};
  double E_spring{0.5 * Params::Xlinks::k_spring * Square(dr)};
  // printf("E = %g\n", E_spring);
  double weight_spring{exp(lambda_spring * fabs(E_spring) / Params::kbT)};
  double E_neighb{site->GetNumNeighborsOccupied() * -1.0 *
                  Params::Xlinks::neighb_neighb_energy};
  double weight_neighb{exp(lambda_neighb * E_neighb)};
  // printf("wt is %g * %g\n", weight_spring, weight_neighb);
  return weight_spring * weight_neighb;
}

double Protein::GetWeight_Diffuse(BindingHead *head, int dir) {

  double lambda_neighb{1.0};
  double lambda_spring{0.5};
  if (n_heads_active_ != 2) {
    Sys::ErrorExit("Protein::GetWeight_diffuse");
  }
  BindingSite *old_site{head->site_};
  int dx{dir * head->GetDirectionTowardRest()};
  // For xlinks exactly at rest,
  if (dx == 0) {
    // Diffuse from rest in a random direction
    if (dir == -1) {
      if (SysRNG::GetRanProb() < 0.5) {
        dx = 1;
      } else {
        dx = -1;
      }
      // Impossible to diffuse toward rest
    } else {
      return 0.0;
    }
  }
  BindingSite *new_site{head->site_->GetNeighbor(dx)};
  if (new_site == nullptr) {
    return 0.0;
  }
  if (new_site->occupant_ != nullptr) {
    return 0.0;
  }
  BindingSite *other_site{head->GetOtherHead()->site_};
  if (old_site->filament_ == other_site->filament_) {
    Sys::ErrorExit("Protein::GetWeight_diffuse [2]");
  }
  double r_x{old_site->pos_[0] - other_site->pos_[0]};
  double r_y{old_site->pos_[1] - other_site->pos_[1]};
  double r{sqrt(Square(r_x) + Square(r_y))};
  double dr{r - Params::Xlinks::r_0};
  // FIXME incorporate k_slack
  double E_old{0.5 * Params::Xlinks::k_spring * Square(dr)};
  double r_x_new{new_site->pos_[0] - other_site->pos_[0]};
  double r_y_new{new_site->pos_[1] - other_site->pos_[1]};
  double r_new{sqrt(Square(r_x_new) + Square(r_y_new))};
  if (r_new > 45.0) {
    return 0.0;
  }
  double dr_new{r_new - Params::Xlinks::r_0};
  double E_new{0.5 * Params::Xlinks::k_spring * Square(dr_new)};
  double dE_spring{E_new - E_old};
  double weight_spring{0.0};
  // Diffusing towards rest is considered an unbinding-type event in
  // regards to Boltzmann factors, since both events let the spring relax
  if (dE_spring < 0.0) {
    weight_spring = exp(lambda_spring * fabs(dE_spring) / Params::kbT);
  }
  // Diffusing away from rest is considered a binding-type event in
  // regards to Boltzmann factors, since both events stretch the spring out
  else {
    weight_spring = exp(-(1.0 - lambda_spring) * fabs(dE_spring) / Params::kbT);
  }
  double E_neighb{old_site->GetNumNeighborsOccupied() * -1.0 *
                  Params::Xlinks::neighb_neighb_energy};
  double weight_neighb{exp(lambda_neighb * E_neighb)};
  // printf("%g & %g\n", weight_spring, weight_neighb);
  return weight_spring * weight_neighb;
}

bool Protein::Diffuse(BindingHead *head, int dir) {

  int dx{dir * head->GetDirectionTowardRest()};
  // FIXME ran num MUST be synchronized with one in GetWeight() above
  // For xlinks exactly at rest,
  if (dx == 0) {
    // Diffuse from rest in a random direction
    if (dir == -1) {
      if (SysRNG::GetRanProb() < 0.5) {
        dx = 1;
      } else {
        dx = -1;
      }
      // Impossible to diffuse toward rest
    } else {
      return false;
    }
  }
  // printf("dx: %i\n", dx);
  BindingSite *old_site = head->site_;
  // printf("no\n");
  int i_new{old_site->index_ + dx};
  // printf("i_old: %i | i_new: %i\n", old_site->index_, i_new);
  if (i_new < 0 or i_new > old_site->filament_->sites_.size() - 1) {
    return false;
  }
  // printf("chaching\n");
  BindingSite *new_site{&old_site->filament_->sites_[i_new]};
  if (new_site->occupant_ != nullptr) {
    return false;
  }
  old_site->occupant_ = nullptr;
  new_site->occupant_ = head;
  head->site_ = new_site;
  // printf("frfr\n\n");
  return true;
}

bool Protein::Tether() {}

bool Protein::Untether() {}