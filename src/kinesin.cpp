#include "kinesin.h"
#include "master_header.h"

int Kinesin::Monomer::GetKif4ANeighborCount() {

  if (site_ == nullptr or motor_->properties_->kinesin4.max_neighbs_ == 0) {
    return 0;
  }
  int n_neighbs{site_->GetKif4ANeighborCount()};
  if (motor_->heads_active_ == 2) {
    n_neighbs--;
  }
  return n_neighbs;
}

Kinesin::Kinesin() {}

void Kinesin::Initialize(system_parameters *parameters,
                         system_properties *properties, int id) {

  id_ = id;
  wally_ = &properties->wallace;
  parameters_ = parameters;
  properties_ = properties;
  SetParameters();
  InitializeNeighborLists();
}

void Kinesin::SetParameters() {

  teth_cutoff_ = properties_->kinesin4.teth_cutoff_;
  comp_cutoff_ = properties_->kinesin4.comp_cutoff_;
}

void Kinesin::InitializeNeighborLists() {

  int n_mts = parameters_->microtubules.count;
  // Serialize this bitch so we just roll one random number
  neighbors_bind_i_teth_.resize(n_mts * (2 * teth_cutoff_ + 1));
  neighbors_tether_bound_.resize(n_mts * (2 * teth_cutoff_ + 1));
}
void Kinesin::ChangeConformation() {

  if (heads_active_ == 0) {
    wally_->ErrorExit("Kinesin::ChangeConformation()");
  }
  Tubulin *site{GetActiveHead()->site_};
  if (site->index_ == site->mt_->plus_end_ and
      parameters_->motors.endpausing_active) {
    return;
  }
  head_one_.trailing_ = !head_one_.trailing_;
  head_two_.trailing_ = !head_two_.trailing_;
}

double Kinesin::GetWeight_Bind_II() {

  // Use n_neighbs weight from dock site, use lattice weight from bound
  // site
  Tubulin *dock_site{GetDockSite()};
  if (dock_site->occupied_) {
    return 0.0;
  }
  Monomer *bound_head{GetActiveHead()};
  double corrected_weight{bound_head->site_->weight_bind_};
  int neighbs_dock{dock_site->GetKif4ANeighborCount() - 1};
  int neighbs_bound{bound_head->site_->GetKif4ANeighborCount()};
  corrected_weight /= properties_->kinesin4.weight_neighbs_bind_[neighbs_bound];
  corrected_weight *= properties_->kinesin4.weight_neighbs_bind_[neighbs_dock];
  return corrected_weight;
}

double Kinesin::GetWeight_Unbind_II() {

  Monomer *chosen_head{nullptr};
  bool head_found{false};
  if (head_one_.ligand_ == "ADPP") {
    chosen_head = &head_one_;
    head_found = true;
  } else if (head_two_.ligand_ == "ADPP") {
    chosen_head = &head_two_;
    if (head_found) {
      wally_->ErrorExit("Kinesin::GetWeight_Unbind_II() [1]");
    }
  }
  if (chosen_head == nullptr) {
    wally_->ErrorExit("Kinesin::GetWeight_Unbind_II() [2]");
  }
  double corrected_weight{chosen_head->site_->weight_unbind_};
  if (properties_->kinesin4.no_internal_force_) {
    return corrected_weight;
  }
  corrected_weight *= corrected_weight;
  // Divide out one of the squared neighbor weights
  int n_neighbs{chosen_head->site_->GetKif4ANeighborCount()};
  corrected_weight /= properties_->kinesin4.weight_neighbs_unbind_[n_neighbs];
  return corrected_weight;
}

double Kinesin::GetWeight_BindATP_II() {

  if (properties_->kinesin4.no_internal_force_) {
    return 0.0;
  }
  Monomer *chosen_head{nullptr};
  bool head_found{false};
  if (head_one_.ligand_ == "NULL") {
    chosen_head = &head_one_;
    head_found = true;
  } else if (head_two_.ligand_ == "NULL") {
    chosen_head = &head_two_;
    if (head_found) {
      wally_->ErrorExit("Kinesin::GetWeight_BindATP_II()");
    }
  }
  if (chosen_head == nullptr) {
    wally_->ErrorExit("Kinesin::GetWeight_BindATP_II()");
  }
  double corrected_weight{chosen_head->site_->weight_bind_};
  // Remove contribution from neighb mechanism
  int n_neighbs{chosen_head->site_->GetKif4ANeighborCount()};
  corrected_weight /= properties_->kinesin4.weight_neighbs_bind_[n_neighbs];
  return corrected_weight;
}

Kinesin::Monomer *Kinesin::GetActiveHead() {

  if (heads_active_ != 1) {
    wally_->ErrorExit("Kinesin::GetActiveHead()");
  }
  if (head_one_.site_ != nullptr) {
    return &head_one_;
  } else {
    return &head_two_;
  }
}

Kinesin::Monomer *Kinesin::GetDockedHead() {

  Kinesin::Monomer *active_head{GetActiveHead()};
  if (active_head->ligand_ != "ADPP") {
    wally_->ErrorExit("Kinesin::GetDockedHead() [1]");
  }
  Kinesin::Monomer *docked_head{active_head->GetOtherHead()};
  if (docked_head->ligand_ != "ADP") {
    wally_->ErrorExit("Kinesin::GetDockedHead() [2]");
  }
  return docked_head;
}

Tubulin *Kinesin::GetDockSite() {

  Kinesin::Monomer *active_head{GetActiveHead()};
  Tubulin *site{active_head->site_};
  if (active_head->ligand_ != "ADPP") {
    wally_->ErrorExit("Kinesin::GetDockSite()");
  }
  int dir{-1};
  if (active_head->trailing_) {
    dir = 1;
  }
  int i_dock{site->index_ + dir * site->mt_->delta_x_};
  if (i_dock >= 0 and i_dock <= site->mt_->n_sites_ - 1) {
    return &site->mt_->lattice_[i_dock];
  } else {
    return nullptr;
  }
}

void Kinesin::UpdateExtension() {

  if (heads_active_ == 0 or HasSatellite() or !tethered_) {
    x_dist_doubled_ = -1;
    return;
  }
  double stalk_coord{GetStalkCoordinate()};
  double anchor_coord{xlink_->GetAnchorCoordinate()};
  int x_dub{int(2 * std::fabs(stalk_coord - anchor_coord))};
  if (x_dub > 2 * teth_cutoff_ or x_dub < 2 * comp_cutoff_) {
    ForceUntether();
  } else {
    x_dist_doubled_ = x_dub;
  }
}

void Kinesin::ForceUntether() {

  // Update motor
  tethered_ = false;
  x_dist_doubled_ = 0;
  // Update xlink
  xlink_->motor_ = nullptr;
  xlink_->tethered_ = false;
  xlink_ = nullptr;
  properties_->kinesin4.FlagForUpdate();
  properties_->prc1.FlagForUpdate();
}

void Kinesin::UntetherSatellite() {

  if (!HasSatellite()) {
    return;
  }
  properties_->prc1.RemoveFromActive(xlink_);
  properties_->prc1.FlagForUpdate();
  // Update xlink details
  xlink_->tethered_ = false;
  xlink_->motor_ = nullptr;
  // Update motor details
  tethered_ = false;
  xlink_ = nullptr;
}

bool Kinesin::HasSatellite() {

  if (tethered_) {
    if (xlink_->heads_active_ == 0) {
      return true;
    }
  }
  return false;
}

bool Kinesin::IsStalled() {

  if (tethered_ and heads_active_ == 1) {
    if (!IsSteppingTowardRest() and (x_dist_doubled_ <= 2 * comp_cutoff_ + 1 or
                                     x_dist_doubled_ >= 2 * teth_cutoff_ - 1)) {
      return true;
    }
  }
  return false;
  /*
  Tubulin *site{heads_active_ == 1
                    ? GetActiveHead()->site_
                    : head_one_.trailing_ ? head_two_.site_ : head_one_.site_};
  if (site->index_ == site->mt_->plus_end_) {
    if (parameters_->motors.endpausing_active) {
      return true;
    } else {
      return false;
    }
  } else {
    return site->mt_->lattice_[site->index_ + site->mt_->delta_x_].occupied_;
  }
  */
}

bool Kinesin::IsSteppingTowardRest() {

  if (!tethered_ or HasSatellite()) {
    wally_->ErrorExit("Kinesin::IsSteppingTowardRest()");
  }
  if (GetActiveHead()->site_->mt_->delta_x_ == GetDirectionTowardRest()) {
    return true;
  }
  return false;
}

int Kinesin::GetDirectionTowardRest() {

  if (!tethered_ or HasSatellite()) {
    return 1;
    // wally_->ErrorExit("Kinesin::GetDirectionTowardRest()");
  }
  double stalk_coord{GetStalkCoordinate()};
  double rest_coord{GetRestLengthCoordinate()};
  if (stalk_coord > rest_coord) {
    return -1;
  } else if (stalk_coord < rest_coord) {
    return 1;
  }
  // If stalk is at rest dist, stepping TO xlink is closer
  // to true rest than stepping FROM xlink
  else {
    double anchor_coord{xlink_->GetAnchorCoordinate()};
    if (stalk_coord > anchor_coord) {
      return -1;
    } else {
      return 1;
    }
  }
}

double Kinesin::GetStalkCoordinate() {

  if (heads_active_ == 0) {
    wally_->ErrorExit("Kinesin::GetStalkCoordinate()");
  }
  if (heads_active_ == 1) {
    Monomer *active_head{GetActiveHead()};
    Tubulin *site{active_head->site_};
    double site_coord = site->index_ + site->mt_->coord_;
    if (active_head->trailing_) {
      return site_coord + (double)site->mt_->delta_x_ / 2;
    } else {
      return site_coord - (double)site->mt_->delta_x_ / 2;
    }
  } else {
    int i_one = head_one_.site_->index_;
    int i_two = head_two_.site_->index_;
    double avg_index = (double)(i_one + i_two) / 2;
    return head_one_.site_->mt_->coord_ + avg_index;
  }
}

double Kinesin::GetRestLengthCoordinate() {

  if (!tethered_ or HasSatellite() or heads_active_ == 0) {
    wally_->ErrorExit("Kinesin::GetRestLengthCoordinate() [1]");
  }
  double stalk_coord{GetStalkCoordinate()};
  double anchor_coord{xlink_->GetAnchorCoordinate()};
  if (stalk_coord > anchor_coord) {
    return anchor_coord + properties_->kinesin4.rest_dist_;
  } else if (stalk_coord < anchor_coord) {
    return anchor_coord - properties_->kinesin4.rest_dist_;
  } else {
    wally_->ErrorExit("Kinesin::GetRestLengthCoordinate() [2]");
  }
}

double Kinesin::GetTetherForce(Tubulin *site) {
  return 0.0;
  /*

    if (!tethered_ or HasSatellite() or heads_active_ == 0) {
      wally_->ErrorExit("Kinesin::GetTetherForce()");
    }
    if (!tethered_) {
      return 0.0;
    }
    // Return 0.0 if site == head_two_.site to avoid double-counting
    if (heads_active_ == 2) {
      if (site == head_two_.site_) {
        return 0.0;
      }
    }
    // cosine_ is always positive and doesn't account for relative tether
    position double force_x{-1 * k_spring_ * extension_ * cosine_}; if
    (extension_ < 0) { force_x = -1 * k_slack_ * extension_ * cosine_;
    }
    // Cosine values are defined assuming the objects experiencing the tether
    // force are the RIGHT of the motor, i.e. in the positive direction. Since
    our
    // cosine values do not go negative, we must multiply by -1 if the site is
    // to the LEFT of the motor, i.e., in the negative direction.
    double stalk_coord{GetStalkCoordinate()};
    double anchor_coord{xlink_->GetAnchorCoordinate()};
    double teth_center{(stalk_coord + anchor_coord) / 2};
    double site_coord = site->index_ + site->mt_->coord_;
    if (site_coord < teth_center) {
      force_x *= -1;
    }
    return force_x;
  */
}

void Kinesin::UpdateNeighbors_Bind_I_Teth() {

  if (!tethered_ or heads_active_ != 0) {
    wally_->ErrorExit("Kinesin::UpdateNeighbors_Bind_I_Teth()");
  }
  n_neighbors_bind_i_teth_ = 0;
  double anchor_coord = xlink_->GetAnchorCoordinate();
  // Scan through all potential neighbor sites; add unoccupied to list
  // FIXME this only works for two MTs as of now FIXME
  for (int i_mt = 0; i_mt < parameters_->microtubules.count; i_mt++) {
    Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
    // Scan 1 site past teth_cutoff to account for rounding of stalk_coord
    for (int delta{-teth_cutoff_ + 1}; delta <= teth_cutoff_ + 1; delta++) {
      int i_neighb{(int)(anchor_coord - mt->coord_) + delta};
      // If i_neighb is negative, skip loop iteration
      if (i_neighb < 0) {
        continue;
      }
      // End scan for this MT after end site (mt_length - 1) is checked
      if (i_neighb > mt->n_sites_ - 1) {
        break;
      }
      Tubulin *neighb{&mt->lattice_[i_neighb]};
      if (!neighb->occupied_) {
        int neighb_coord{i_neighb + mt->coord_};
        int x_dub{(int)(2 * fabs(anchor_coord - neighb_coord))};
        if (x_dub >= 2 * comp_cutoff_ and x_dub <= 2 * teth_cutoff_) {
          neighbors_bind_i_teth_[n_neighbors_bind_i_teth_++] = neighb;
        }
      }
    }
  }
}

void Kinesin::UpdateNeighbors_Tether_Bound() {

  if (tethered_ or heads_active_ == 0) {
    wally_->ErrorExit("Kinesin::UpdateNeighbors_Tether_Bound()");
  }
  n_neighbors_tether_bound_ = 0;
  double stalk_coord = GetStalkCoordinate();
  // Scan through all potential neighbor sites; add untethered xlinks
  // FIXME this only works for two MTs as of now FIXME
  for (int i_mt = 0; i_mt < parameters_->microtubules.count; i_mt++) {
    Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
    // Only scan over sites within +/- teth_cutoff_
    for (int delta{-(teth_cutoff_ + 1)}; delta <= (teth_cutoff_ + 1); delta++) {
      int i_neighb{(int)(stalk_coord - mt->coord_) + delta};
      // Start index at first site (0) if site index is <= 0
      if (i_neighb < 0) {
        continue;
      }
      // End scan once last site (mt_length - 1) is reached
      if (i_neighb > mt->n_sites_ - 1) {
        break;
      }
      Tubulin *site{&mt->lattice_[i_neighb]};
      if (site->xlink_head_ == nullptr) {
        continue;
      }
      AssociatedProtein *neighb{site->xlink_head_->xlink_};
      if (neighb->tethered_) {
        continue;
      }
      double anchor_coord{neighb->GetAnchorCoordinate()};
      int x_dub{(int)(2 * fabs(anchor_coord - stalk_coord))};
      if (x_dub >= 2 * comp_cutoff_ and x_dub <= 2 * teth_cutoff_) {
        if (neighb->heads_active_ == 1) {
          neighbors_tether_bound_[n_neighbors_tether_bound_++] = neighb;
        }
        // FIXME this is bad
        else if (neighb->heads_active_ == 2 and i_mt % 2 == 0) {
          neighbors_tether_bound_[n_neighbors_tether_bound_++] = neighb;
        }
      }
    }
  }
}
double Kinesin::GetWeight_Bind_I_Teth(Tubulin *site) {

  // Since only one foot of the motor attaches initially,
  // we only need to consider the coordinate of one site
  // (assuming the diffusion of the other foot avgs out)
  int i_site = site->index_;
  int mt_coord = site->mt_->coord_;
  double site_coord = (double)(mt_coord + i_site);
  // Motor head to bind first is always leading, so stalk coord is 0.5 behind
  double stalk_coord{site_coord - double(site->mt_->delta_x_) / 2};
  double anchor_coord = xlink_->GetAnchorCoordinate();
  double x_dist = fabs(anchor_coord - stalk_coord);
  // Multiply by 2 to guarentee an integer
  int x_dub = x_dist * 2;
  int n_neighbs = site->GetKif4ANeighborCount();
  return (site->weight_bind_ *
          properties_->kinesin4.weight_neighbs_bind_[n_neighbs] *
          properties_->kinesin4.weight_teth_bind_[x_dub]);
}

double Kinesin::GetWeight_Tether_Bound(AssociatedProtein *xlink) {

  double stalk_coord = GetStalkCoordinate();
  double anchor_coord = xlink->GetAnchorCoordinate();
  double x_dist = fabs(anchor_coord - stalk_coord);
  // Multiply by 2 to guarentee an integer
  int x_dub = x_dist * 2;
  // Get tethering weight that corressponds to this adj. x-dist
  return properties_->kinesin4.weight_teth_bind_[x_dub];
}

double Kinesin::GetTotalWeight_Bind_I_Teth() {

  UpdateNeighbors_Bind_I_Teth();
  double tot_weight{0.0};
  for (int i_neighb{0}; i_neighb < n_neighbors_bind_i_teth_; i_neighb++) {
    tot_weight += GetWeight_Bind_I_Teth(neighbors_bind_i_teth_[i_neighb]);
  }
  return tot_weight;
}

double Kinesin::GetTotalWeight_Tether_Bound() {

  double tot_weight{0.0};
  UpdateNeighbors_Tether_Bound();
  for (int i_neighb{0}; i_neighb < n_neighbors_tether_bound_; i_neighb++) {
    tot_weight += GetWeight_Tether_Bound(neighbors_tether_bound_[i_neighb]);
  }
  return tot_weight;
}

Tubulin *Kinesin::GetWeightedSite_Bind_I_Teth() {

  double weight_tot{GetTotalWeight_Bind_I_Teth()};
  double ran{properties_->gsl.GetRanProb()};
  double p_cum{0.0};
  // Scan through neighbor sites; pick one based on normalized weights
  for (int i_neighb{0}; i_neighb < n_neighbors_bind_i_teth_; i_neighb++) {
    Tubulin *neighb{neighbors_bind_i_teth_[i_neighb]};
    p_cum += GetWeight_Bind_I_Teth(neighb) / weight_tot;
    if (ran < p_cum)
      return neighb;
  }
  printf("No neighb for bind_i_teth (MOT)??\n");
  return nullptr;
}

AssociatedProtein *Kinesin::GetWeightedXlink_Tether_Bound() {

  double weight_tot{GetTotalWeight_Tether_Bound()};
  double ran{properties_->gsl.GetRanProb()};
  double p_cum{0.0};
  // Scan through neighbor xlinks; pick one based on normalized weights
  for (int i_neighb{0}; i_neighb < n_neighbors_tether_bound_; i_neighb++) {
    AssociatedProtein *neighb{neighbors_tether_bound_[i_neighb]};
    p_cum += GetWeight_Tether_Bound(neighb) / weight_tot;
    if (ran < p_cum)
      return neighb;
  }
  printf("No neighb for tether_bound??\n");
  return nullptr;
}