#include "cylaks/protein.hpp"
#include "cylaks/motor.hpp"
#include "cylaks/protofilament.hpp"

bool Protein::HasSatellite() {

  if (teth_partner_ == nullptr) {
    return false;
  }
  if (teth_partner_->GetNumHeadsActive() == 0) {
    return true;
  }
  return false;
}

bool Protein::UntetherSatellite() {

  if (IsTethered()) {
    if (teth_partner_->GetNumHeadsActive() == 0) {
      teth_partner_->teth_partner_ = nullptr;
      teth_partner_->Object::teth_partner_ = nullptr;
      teth_partner_->head_one_.teth_partner_ = nullptr;
      teth_partner_->head_two_.teth_partner_ = nullptr;
      teth_partner_ = nullptr;
      Object::teth_partner_ = nullptr;
      head_one_.teth_partner_ = nullptr;
      head_two_.teth_partner_ = nullptr;
      return true;
    } else {
      Sys::ErrorExit("Protein::UntetherSatellite\n");
    }
  }
  return false;
}

double Protein::GetSoloWeight_Bind_I_Teth(BindingSite *target) {

  double r_x{teth_partner_->GetAnchorCoordinate(0) - target->pos_[0]};
  double x{std::fabs(r_x) / Params::Filaments::site_size};
  if (x < Sys::teth_x_min_ or x > Sys::teth_x_max_) {
    return 0.0;
  }
  double r_y{teth_partner_->GetAnchorCoordinate(1) - target->pos_[1]};
  double r{sqrt(Square(r_x) + Square(r_y))};
  double dr{Params::Motors::r_0 - r};
  double k{dr > 0.0 ? Params::Motors::k_spring : Params::Motors::k_slack};
  double weight_teth{
      exp(-(1.0 - _lambda_spring) * 0.5 * k * Square(dr) / Params::kbT)};
  if (weight_teth > _max_weight) {
    printf("uhhhh\n");
    return 0.0;
  }
  // double weight_teth{tether_.GetWeight_Bind(r)};
  double weight_site{target->GetWeight_Bind()};
  return weight_teth * weight_site;
}

void Protein::UpdateNeighbors_Bind_I_Teth() {

  if (n_heads_active_ != 0) {
    Sys::ErrorExit("Protein::UpdateNeighbors_Bind_I_Teth() [1]");
  }
  if (!IsTethered()) {
    Sys::ErrorExit("Protein::UpdateNeighbors_Bind_I_Teth() [2]");
  }
  if (teth_partner_->GetNumHeadsActive() == 0) {
    Sys::ErrorExit("Protein::UpdateNeighbors_Bind_I_Teth() [3]");
  }
  n_neighbors_bind_i_teth_ = 0;
  // ! FIXME add n_mt = 2 case
  BindingSite *anchor{teth_partner_->GetHeadOne()->site_};
  if (anchor == nullptr) {
    anchor = teth_partner_->GetHeadTwo()->site_;
  }
  for (int dx{Sys::teth_x_min_}; dx <= Sys::teth_x_max_; dx++) {
    for (int dir{-1}; dir <= 1; dir += 2) {
      int i_neighb{(int)anchor->index_ + dir * dx};
      // printf("%i + (%i)(%i) = %i\n", (int)anchor->index_, dir, dx, i_neighb);
      if (i_neighb < 0 or i_neighb >= anchor->filament_->sites_.size()) {
        continue;
      }
      BindingSite *neighb{&anchor->filament_->sites_[i_neighb]};
      if (!neighb->IsOccupied()) {
        // printf("%zu\n", neighbors_bind_i_teth_.size());
        // printf("%i\n", n_neighbors_bind_i_teth_);
        neighbors_bind_i_teth_[n_neighbors_bind_i_teth_++] = neighb;
      }
    }
  }
}

BindingSite *Protein::GetNeighbor_Bind_I_Teth() {

  double weight_tot{GetWeight_Bind_I_Teth()};
  double ran{SysRNG::GetRanProb()};
  double p_cum{0.0};
  Sys::Log(2, "%i NEIGHBS\n", n_neighbors_bind_i_teth_);
  Sys::Log(2, "ran = %g\n", ran);
  for (int i_neighb{0}; i_neighb < n_neighbors_bind_i_teth_; i_neighb++) {
    BindingSite *neighb{neighbors_bind_i_teth_[i_neighb]};
    p_cum += GetSoloWeight_Bind_I_Teth(neighb) / weight_tot;
    Sys::Log(2, "p_cum = %g\n", p_cum);
    if (ran < p_cum) {
      return neighb;
    }
  }
  return nullptr;
}

double Protein::GetWeight_Bind_I_Teth() {

  double tot_weight{0.0};
  UpdateNeighbors_Bind_I_Teth();
  for (int i_neighb{0}; i_neighb < n_neighbors_bind_i_teth_; i_neighb++) {
    tot_weight += GetSoloWeight_Bind_I_Teth(neighbors_bind_i_teth_[i_neighb]);
  }
  // if (tot_weight != 0.0) {
  //   printf("%g\n", tot_weight);
  // }
  return tot_weight;
}

bool Protein::UpdateExtension() {
  if (n_heads_active_ != 2) {
    return true;
  }
  // Update head positions
  for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
    head_one_.pos_[i_dim] = head_one_.site_->pos_[i_dim];
    head_two_.pos_[i_dim] = head_two_.site_->pos_[i_dim];
  }
  // printf("r1 = (%g, %g)\n", head_one_.pos_[0], head_one_.pos_[1]);
  // printf("r2 = (%g, %g)\n", head_two_.pos_[0], head_two_.pos_[1]);
  // Update spring position
  bool spring_attached{spring_.UpdatePosition()};
  if (!spring_attached) {
    return false;
  }
  // bool head_one_attached{pivot_one_.UpdatePosition()};
  // if (!head_one_attached) {
  //   return false;
  // }
  // bool head_two_attached{pivot_two_.UpdatePosition()};
  // if (!head_two_attached) {
  //   return false;
  // }
  // If spring is still attached after updates, apply forces & torques
  spring_.ApplyForces();
  // pivot_one_.ApplyForces();
  // pivot_two_.ApplyForces();
  // printf("**\n\n");
  return true;
}

int Protein::GetDirectionTowardRest(BindingHead *head) {

  if (n_heads_active_ == 1) {
    return 1;
  } else if (n_heads_active_ == 2) {
    double x_static{head->GetOtherHead()->site_->pos_[0]};
    double x_fwd{head->site_->pos_[0] + Params::Filaments::site_size};
    double x_bck{head->site_->pos_[0] - Params::Filaments::site_size};
    double r_x_fwd{x_fwd - x_static};
    double r_x_bck{x_bck - x_static};
    double r_x{head->site_->pos_[0] - x_static};
    double r_y{head->site_->pos_[1] - head->GetOtherHead()->site_->pos_[1]};
    double r{sqrt(Square(r_x) + Square(r_y))};
    double r_fwd{sqrt(Square(r_x_fwd) + Square(r_y))};
    double r_bck{sqrt(Square(r_x_bck) + Square(r_y))};
    double dr{r - spring_.r_rest_};
    double dr_fwd{r_fwd - spring_.r_rest_};
    double dr_bck{r_bck - spring_.r_rest_};
    if (Square(dr_fwd) < Square(dr)) {
      return 1;
    }
    if (Square(dr_bck) < Square(dr)) {
      return -1;
    }
    return 0;
    /*
    // printf("%g > %g?\n", head->site_->pos_[0], GetAnchorCoordinate(0));
    if (head->site_->pos_[0] > GetAnchorCoordinate(0)) {
      return -1;
    } else if (head->site_->pos_[0] < GetAnchorCoordinate(0)) {
      return 1;
    } else {
      return 0;
    }
    */
  } else {
    Sys::ErrorExit("Protein::GetDirToRest()\n");
  }
  return 0;
}

double Protein::GetAnchorCoordinate(int i_dim) {

  if (n_heads_active_ == 0) {
    Sys::ErrorExit("Protein::GetAnchorCoord()");
  }
  if (n_heads_active_ == 1) {
    return GetActiveHead()->site_->pos_[i_dim];
  }
  // printf("%i (%i)\n", n_heads_active_, species_id_);
  return (head_one_.site_->pos_[i_dim] + head_two_.site_->pos_[i_dim]) / 2;
}

void Protein::UpdateNeighbors_Bind_II() {

  n_neighbors_bind_ii_ = 0;
  BindingSite *site{GetActiveHead()->site_};
  Protofilament *neighb_fil{site->filament_->neighbor_};
  double r_y{site->filament_->pos_[1] - neighb_fil->pos_[1]};
  double r_x_max{sqrt(Square(spring_.r_max_) - Square(r_y))};
  int delta_max{(int)std::ceil(r_x_max / Params::Filaments::site_size)};
  for (int delta{-delta_max}; delta <= delta_max; delta++) {
    BindingSite *neighb{neighb_fil->GetNeighb(site, delta)};
    if (neighb == nullptr) {
      continue;
    }
    if (neighb->occupant_ == nullptr) {
      neighbors_bind_ii_[n_neighbors_bind_ii_++] = neighb;
    }
  }
}

double Protein::GetSoloWeight_Bind_II(BindingSite *neighb) {

  BindingSite *site{GetActiveHead()->site_};
  double r_x{neighb->pos_[0] - site->pos_[0]};
  double r_y{neighb->pos_[1] - site->pos_[1]};
  double r{sqrt(Square(r_x) + Square(r_y))};
  if (r < spring_.r_min_ or r > spring_.r_max_) {
    return 0.0;
  }
  double weight_spring{spring_.GetWeight_Bind(r)};
  double weight_site{neighb->GetWeight_Bind()};
  return weight_spring * weight_site;
}

BindingSite *Protein::GetNeighbor_Bind_II() {

  double weight_tot{GetWeight_Bind_II()};
  double ran{SysRNG::GetRanProb()};
  double p_cum{0.0};
  Sys::Log(2, "%i NEIGHBS\n", n_neighbors_bind_ii_);
  Sys::Log(2, "ran = %g\n", ran);
  for (int i_neighb{0}; i_neighb < n_neighbors_bind_ii_; i_neighb++) {
    BindingSite *neighb{neighbors_bind_ii_[i_neighb]};
    p_cum += GetSoloWeight_Bind_II(neighb) / weight_tot;
    Sys::Log(2, "p_cum = %g\n", p_cum);
    if (ran < p_cum) {
      Sys::Log(2, "*** chose neighb %i ***\n\n", neighb->index_);
      return neighb;
    }
  }
  return nullptr;
}

double Protein::GetWeight_Diffuse(BindingHead *head, int dir) {

  if (n_heads_active_ != 2) {
    Sys::ErrorExit("Protein::GetWeight_diffuse");
  }
  int dx{dir * head->GetDirectionTowardRest()};
  // For xlinks exactly at rest,
  if (dx == 0) {
    // Diffuse from rest in a random direction
    if (dir == -1) {
      ran_ = SysRNG::GetRanProb();
      if (ran_ < 0.5) {
        dx = 1;
      } else {
        dx = -1;
      }
      // Impossible to diffuse toward rest
    } else {
      return 0.0;
    }
  }
  BindingSite *new_loc{head->site_->GetNeighbor(dx)};
  // ! FIXME temporary hacky solution for forced_slide test mode
  if (new_loc == nullptr) {
    if (Sys::test_mode_.empty()) {
      return 0.0;
    }
    // Rather than use ghost sites, just use other crosslinker head to get
    // spring weight (should be the exact same dE)
    BindingHead *other_head{head->GetOtherHead()};
    new_loc = other_head->site_->GetNeighbor(-1 * dx);
    // If other head is also trying to diffuse off the end, we're S.O.L.
    if (new_loc == nullptr) {
      return 0.0;
    }
    BindingSite *old_loc{other_head->site_};
    BindingSite *static_loc{head->site_};
    if (old_loc->filament_ == static_loc->filament_) {
      Sys::ErrorExit("Protein::GetWeight_diffuse [3]");
    }
    spring_.UpdatePosition();
    double weight_spring{spring_.GetWeight_Shift(static_loc, old_loc, new_loc)};
    double weight_neighb{head->site_->GetWeight_Unbind()};
    return weight_spring * weight_neighb;
  }
  if (new_loc->occupant_ != nullptr) {
    return 0.0;
  }
  BindingSite *old_loc{head->site_};
  if (old_loc->GetNumNeighborsOccupied() == _n_neighbs_max) {
    return 0.0;
  }
  BindingSite *static_loc{head->GetOtherHead()->site_};
  if (old_loc->filament_ == static_loc->filament_) {
    Sys::ErrorExit("Protein::GetWeight_diffuse [2]");
  }
  spring_.UpdatePosition();
  double weight_spring{spring_.GetWeight_Shift(static_loc, old_loc, new_loc)};
  double weight_neighb{head->site_->GetWeight_Unbind()};
  // printf("WT[%i] = %g\n", dx, weight_spring * weight_neighb);
  return weight_spring * weight_neighb;
}

double Protein::GetWeight_Bind_II() {

  // UpdateExtension();
  double tot_weight{0.0};
  UpdateNeighbors_Bind_II();
  for (int i_neighb{0}; i_neighb < n_neighbors_bind_ii_; i_neighb++) {
    tot_weight += GetSoloWeight_Bind_II(neighbors_bind_ii_[i_neighb]);
  }
  return tot_weight;
}

double Protein::GetWeight_Unbind_II(BindingHead *head) {

  if (n_heads_active_ != 2) {
    Sys::ErrorExit("Protein::GetWeight_Unbind_II()\n");
  }
  double weight_spring{spring_.GetWeight_Unbind()};
  double weight_site{head->site_->GetWeight_Unbind()};
  return weight_spring * weight_site;
}

bool Protein::Diffuse(BindingHead *head, int dir) {

  int dx{dir * head->GetDirectionTowardRest()};
  // For xlinks exactly at rest,
  if (dx == 0) {
    // Diffuse from rest in a random direction
    if (dir == -1) {
      if (ran_ < 0.5) {
        dx = 1;
      } else {
        dx = -1;
      }
      // Impossible to diffuse toward rest
    } else {
      return false;
    }
  }
  BindingSite *old_site = head->site_;
  int i_new{(int)old_site->index_ + dx};
  if (i_new < 0 or i_new > old_site->filament_->sites_.size() - 1) {
    return false;
  }
  BindingSite *new_site{&old_site->filament_->sites_[i_new]};
  if (new_site->occupant_ != nullptr) {
    return false;
  }
  old_site->occupant_ = nullptr;
  new_site->occupant_ = head;
  head->site_ = new_site;
  return true;
}

bool Protein::Diffuse_Side(BindingHead *head, int dir) {

  if (dir != 1 and dir != -1) {
    Sys::ErrorExit("Whoops in Diffuse_Side [Protein]");
  }
  BindingSite *old_site{head->site_};
  BindingSite *new_site{old_site->GetNeighbor_Side(dir)};
  if (new_site == nullptr) {
    return false;
  }
  if (new_site->occupant_ != nullptr) {
    return false;
  }
  old_site->occupant_ = nullptr;
  new_site->occupant_ = head;
  head->site_ = new_site;
  return true;
}

bool Protein::Bind(BindingSite *site, BindingHead *head) {

  if (site == nullptr) {
    return false;
  }
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

bool Protein::Tether(Protein *target) {

  if (target->IsTethered()) {
    printf("target alrdy tethered\n");
    return false;
    // Sys::ErrorExit("Protein::Tether() [1]");
  }
  if (IsTethered()) {
    printf("i am alrdy tethered\n");
    return false;
    // Sys::ErrorExit("Protein::Tether() [2]");
  }
  teth_partner_ = target;
  Object::teth_partner_ = dynamic_cast<Object *>(target);
  head_one_.teth_partner_ = dynamic_cast<Object *>(target);
  head_two_.teth_partner_ = dynamic_cast<Object *>(target);
  target->teth_partner_ = this;
  target->Object::teth_partner_ = dynamic_cast<Object *>(this);
  target->head_one_.teth_partner_ = dynamic_cast<Object *>(this);
  target->head_two_.teth_partner_ = dynamic_cast<Object *>(this);
  return true;
}

bool Protein::Untether() {

  if (IsTethered() and !HasSatellite()) {
    teth_partner_->teth_partner_ = nullptr;
    teth_partner_->Object::teth_partner_ = nullptr;
    teth_partner_->head_one_.teth_partner_ = nullptr;
    teth_partner_->head_two_.teth_partner_ = nullptr;
    teth_partner_ = nullptr;
    Object::teth_partner_ = nullptr;
    head_one_.teth_partner_ = nullptr;
    head_two_.teth_partner_ = nullptr;
    return true;
  } else {
    printf("untether wot m8\n");
    // Sys::ErrorExit("untether wot m8\n");
  }
  return false;
}