#include "associated_protein.h"
#include "master_header.h"

AssociatedProtein::AssociatedProtein() {}

void AssociatedProtein::Initialize(system_parameters *parameters,
                                   system_properties *properties, int id) {

  id_ = id;
  wally_ = &properties->wallace;
  parameters_ = parameters;
  properties_ = properties;
  SetParameters();
  InitializeLookupTables();
  InitializeNeighborLists();
}

void AssociatedProtein::SetParameters() {

  r_0_ = parameters_->xlinks.r_0;
  k_spring_ = parameters_->xlinks.k_spring;
  rest_dist_ = properties_->prc1.rest_dist_;
  dist_cutoff_ = properties_->prc1.dist_cutoff_;
}

void AssociatedProtein::InitializeLookupTables() {

  /* Construct rly basic lookup tables on the fly */
  double r_0{parameters_->xlinks.r_0};
  double r_y{parameters_->microtubules.y_dist};
  double site_size{parameters_->microtubules.site_size};
  cosine_lookup_.resize(dist_cutoff_ + 1);
  extension_lookup_.resize(dist_cutoff_ + 1);
  for (int x_dist{0}; x_dist <= dist_cutoff_; x_dist++) {
    double r_x = site_size * x_dist;
    double r = sqrt(r_x * r_x + r_y * r_y);
    extension_lookup_[x_dist] = r - r_0_;
    cosine_lookup_[x_dist] = r_x / r;
  }
  /* Copy more involved lookup tables from MGMT to ensure consistency */
  int teth_cutoff{properties_->kinesin4.teth_cutoff_};
  int comp_cutoff{properties_->kinesin4.comp_cutoff_};
  int max_neighbs{properties_->prc1.max_neighbs_};
  weight_bind_ii_.resize(max_neighbs + 1);
  weight_bind_i_teth_.resize(max_neighbs + 1);
  weight_bind_ii_to_teth_.resize(max_neighbs + 1);
  weight_bind_ii_fr_teth_.resize(max_neighbs + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs; n_neighbs++) {
    weight_bind_ii_[n_neighbs].resize(dist_cutoff_ + 1);
    for (int x{0}; x <= dist_cutoff_; x++) {
      weight_bind_ii_[n_neighbs][x] =
          properties_->prc1.weight_bind_ii_[n_neighbs][x];
    }
    weight_bind_i_teth_[n_neighbs].resize(2 * teth_cutoff + 1);
    weight_bind_ii_to_teth_[n_neighbs].resize(2 * teth_cutoff + 1);
    weight_bind_ii_fr_teth_[n_neighbs].resize(2 * teth_cutoff + 1);
    for (int x_dub{2 * comp_cutoff}; x_dub <= 2 * teth_cutoff; x_dub++) {
      weight_bind_i_teth_[n_neighbs][x_dub] =
          properties_->prc1.weight_bind_i_teth_[n_neighbs][x_dub];
      weight_bind_ii_to_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      weight_bind_ii_fr_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      for (int x{0}; x <= dist_cutoff_; x++) {
        weight_bind_ii_to_teth_[n_neighbs][x_dub][x] =
            properties_->prc1.weight_bind_ii_to_teth_[n_neighbs][x_dub][x];
        weight_bind_ii_fr_teth_[n_neighbs][x_dub][x] =
            properties_->prc1.weight_bind_ii_fr_teth_[n_neighbs][x_dub][x];
      }
    }
  }
}

void AssociatedProtein::InitializeNeighborLists() {

  int n_mts{parameters_->microtubules.count};
  int teth_cutoff{properties_->kinesin4.teth_cutoff_};
  neighbor_sites_ii_.resize((n_mts - 1) * (2 * dist_cutoff_ + 1));
  neighbor_sites_i_teth_.resize(n_mts * (2 * teth_cutoff + 1));
  neighbor_sites_ii_teth_.resize((n_mts - 1) * (2 * dist_cutoff_ + 1));
}

AssociatedProtein::Monomer *AssociatedProtein::GetActiveHead() {

  if (heads_active_ == 1) {
    if (head_one_.site_ != nullptr) {
      return &head_one_;
    } else if (head_two_.site_ != nullptr) {
      return &head_two_;
    } else {
      wally_->ErrorExit("AssociatedProtein::GetActiveHead [1]");
    }
  } else {
    wally_->ErrorExit("AssociatedProtein::GetActiveHead [2]");
  }
}

double AssociatedProtein::GetAnchorCoordinate() {

  // If single bound, use that head; assume other's diffusion avgs out
  if (heads_active_ == 1) {
    Tubulin *site{GetActiveHead()->site_};
    return (double)(site->mt_->coord_ + site->index_);
  }
  // If double bound, use avg of head's indices
  else if (heads_active_ == 2) {
    int index_one = head_one_.site_->index_;
    int mt_coord_one = head_one_.site_->mt_->coord_;
    double coord_one = (double)(mt_coord_one + index_one);
    int index_two = head_two_.site_->index_;
    int mt_coord_two = head_two_.site_->mt_->coord_;
    double coord_two = (double)(mt_coord_two + index_two);
    return (coord_one + coord_two) / 2;
  } else {
    wally_->ErrorExit("AssociatedProtein::GetAnchorCoordinate()");
  }
}

double AssociatedProtein::GetExtensionForce(Tubulin *site) {

  if (heads_active_ == 2) {
    UpdateExtension();
    // Make sure we didn't force an unbinding event
    if (heads_active_ == 2) {
      double force_mag{extension_ * k_spring_};
      double force{force_mag * cosine_};
      int site_coord{site->index_ + site->mt_->coord_};
      if (site_coord > GetAnchorCoordinate()) {
        force *= -1;
      }
      return force;
    } else {
      return 0;
    }
  } else {
    wally_->ErrorExit("AssociatedProtein::GetExtensionForce()");
  }
}

int AssociatedProtein::GetDirectionTowardRest(Tubulin *site) {

  if (heads_active_ == 0) {
    wally_->ErrorExit("AssociatedProtein::GetDirectionTowardRest()");
  }
  // Xlinks with satellites behave as if untethered (no energy from tether)
  if (!tethered_ or HasSatellite()) {
    // Always return 1 for singly-bound by convention
    if (heads_active_ == 1) {
      return 1;
    } else if (heads_active_ == 2) {
      if (x_dist_ == 0) {
        double ran{properties_->gsl.GetRanProb()};
        if (ran < 0.5) {
          return -1;
        } else {
          return 1;
        }
      } else {
        int site_coord{site->index_ + site->mt_->coord_};
        if (site_coord > GetAnchorCoordinate()) {
          return -1;
        } else {
          return 1;
        }
      }
    }
  }
  // If tethered to a bound motor, use motor's function instead
  else {
    // xlink->GetDirToRest is ill-defined for x = 0, so use motor's
    // function instead (multiply by -1 since xlink is the one moving)
    return -1 * motor_->GetDirectionTowardRest();
  }
}

void AssociatedProtein::UpdateExtension() {

  if (heads_active_ == 2) {
    // Calculate first head's coordinate
    int coord_one{head_one_.site_->mt_->coord_ + head_one_.site_->index_};
    // Calculate second head's coordinate
    int coord_two{head_two_.site_->mt_->coord_ + head_two_.site_->index_};
    // Calculate x distance in # of sites
    int x_dist{abs(coord_one - coord_two)};
    if (x_dist <= dist_cutoff_) {
      x_dist_ = x_dist;
      cosine_ = cosine_lookup_[x_dist];
      extension_ = extension_lookup_[x_dist];
    } else {
      ForceUnbind();
    }
  } else {
    x_dist_ = 0;
    extension_ = 0;
    cosine_ = 0;
  }
  if (tethered_) {
    motor_->UpdateExtension();
  }
}

void AssociatedProtein::ForceUnbind() {

  // Check to see if xlink has a satellite motor (tethered but not bound)
  bool has_satellite{false};
  if (tethered_) {
    if (motor_->heads_active_ == 0) {
      has_satellite = true;
    }
  }
  // Xlinks with satellites behave as if untethered (no energy from tether)
  if (!tethered_ or has_satellite) {
    // If xlink isn't tethered, flip a coin to choose which head to unbind
    double ran{properties_->gsl.GetRanProb()};
    if (ran < 0.5) {
      head_one_.site_->xlink_head_ = nullptr;
      head_one_.site_->occupied_ = false;
      head_one_.site_ = nullptr;
    } else {
      head_two_.site_->xlink_head_ = nullptr;
      head_two_.site_->occupied_ = false;
      head_two_.site_ = nullptr;
    }
  }
  // If xlink is tethered, take energy change of tether into account
  // (Motor must have at least 1 head bound for tether to store energy)
  else if (motor_->heads_active_ > 0) {
    // Unbinding head farthest from teth rest brings teth closer to rest
    Tubulin *site_to{GetSiteFartherFromTethRest()};
    int neighbs_to{site_to->GetPRC1NeighborCount()};
    // Unbinding head closest to teth rest brings teth farther from rest
    Tubulin *site_fr{GetSiteCloserToTethRest()};
    int neighbs_fr{site_fr->GetPRC1NeighborCount()};
    // Calculate relative probabilities of unbinding either head
    int x{x_dist_};
    int x_dub{motor_->x_dist_doubled_};
    double p_to{properties_->prc1.p_unbind_ii_to_teth_[neighbs_to][x_dub][x]};
    double p_fr{properties_->prc1.p_unbind_ii_fr_teth_[neighbs_fr][x_dub][x]};
    double p_tot{p_to + p_fr};
    // Roll a random number to see which should be unbound
    double ran{properties_->gsl.GetRanProb()};
    if (ran < p_to / p_tot) {
      site_to->xlink_head_ = nullptr;
      site_to->occupied_ = false;
      if (site_to == head_one_.site_) {
        head_one_.site_ = nullptr;
      } else {
        head_two_.site_ = nullptr;
      }
    } else {
      site_fr->xlink_head_ = nullptr;
      site_fr->occupied_ = false;
      if (site_fr == head_one_.site_) {
        head_one_.site_ = nullptr;
      } else {
        head_two_.site_ = nullptr;
      }
    }
  }
  heads_active_--;
  x_dist_ = 0;
  extension_ = 0;
  cosine_ = 0;
  if (tethered_) {
    // Update motor extension since anchor coord has changed
    motor_->UpdateExtension();
  }
}

void AssociatedProtein::UntetherSatellite() {

  if (!tethered_) {
    return;
  }
  if (motor_->heads_active_ == 0) {
    // Remove satellite motor from active_ list, replace with last entry
    int i_last{properties_->kinesin4.n_active_ - 1};
    Kinesin *last_entry{properties_->kinesin4.active_[i_last]};
    int i_this{motor_->active_index_};
    properties_->kinesin4.active_[i_this] = last_entry;
    last_entry->active_index_ = i_this;
    properties_->kinesin4.n_active_--;
    // Update motor details
    motor_->tethered_ = false;
    motor_->xlink_ = nullptr;
    // Update xlink details
    tethered_ = false;
    motor_ = nullptr;
  }
}

bool AssociatedProtein::HasSatellite() {

  if (tethered_) {
    if (motor_->heads_active_ == 0) {
      return true;
    }
  }
  return false;
}

Tubulin *AssociatedProtein::GetSiteCloserToTethRest() {

  if (heads_active_ != 2 or !tethered_) {
    wally_->ErrorExit("AssociatedProtein::GetSiteCloserToTethRest()");
  }
  double stalk_coord = motor_->GetStalkCoordinate();
  double site_one_coord =
      head_one_.site_->index_ + head_one_.site_->mt_->coord_;
  int dx_dub_one = 2 * fabs(stalk_coord - site_one_coord);
  double site_two_coord =
      head_two_.site_->index_ + head_two_.site_->mt_->coord_;
  int dx_dub_two = 2 * fabs(stalk_coord - site_two_coord);
  double anchor_coord = GetAnchorCoordinate();
  int x_dub = 2 * fabs(anchor_coord - stalk_coord);
  int rest_dist_dub = 2 * properties_->kinesin4.rest_dist_;
  // twice dx of teth extension if 1 head unbinds
  int dx_teth_dub = x_dist_;
  // twice the distance teth is drom rest
  int dist_from_rest_dub = abs(x_dub - rest_dist_dub);
  // if unbinding causes teth to go from slack to spring or
  // vise versa, always choose the site that makes it go slack
  if (dx_teth_dub > dist_from_rest_dub) {
    if (dx_dub_one < dx_dub_two)
      return head_one_.site_;
    else
      return head_two_.site_;
  } else if (x_dub >= rest_dist_dub) {
    if (dx_dub_one < dx_dub_two)
      return head_one_.site_;
    else
      return head_two_.site_;
  } else {
    if (dx_dub_one > dx_dub_two)
      return head_one_.site_;
    else
      return head_two_.site_;
  }
}

Tubulin *AssociatedProtein::GetSiteFartherFromTethRest() {

  if (heads_active_ == 2 && tethered_ == true) {
    double stalk_coord = motor_->GetStalkCoordinate();
    double site_one_coord =
        head_one_.site_->index_ + head_one_.site_->mt_->coord_;
    int dx_dub_one = 2 * fabs(stalk_coord - site_one_coord);
    double site_two_coord =
        head_two_.site_->index_ + head_two_.site_->mt_->coord_;
    int dx_dub_two = 2 * fabs(stalk_coord - site_two_coord);
    double anchor_coord = GetAnchorCoordinate();
    int x_dub = 2 * fabs(anchor_coord - stalk_coord);
    int rest_dist_dub = 2 * properties_->kinesin4.rest_dist_;
    // twice dx of teth extension if 1 head unbinds
    int dx_teth_dub = x_dist_;
    // twice the distance teth is drom rest
    int dist_from_rest_dub = abs(x_dub - rest_dist_dub);
    // if unbinding causes teth to go from slack to spring or
    // vise versa, always choose the site that makes it taut
    if (dx_teth_dub > dist_from_rest_dub) {
      if (dx_dub_one > dx_dub_two)
        return head_one_.site_;
      else
        return head_two_.site_;
    } else if (x_dub >= rest_dist_dub) {
      if (dx_dub_one > dx_dub_two)
        return head_one_.site_;
      else
        return head_two_.site_;
    } else {
      if (dx_dub_one < dx_dub_two)
        return head_one_.site_;
      else
        return head_two_.site_;
    }
  } else {
    printf("Error in get_site_farther_fr_teth_rest XLINK\n");
    exit(1);
  }
}

void AssociatedProtein::UpdateNeighbors_Bind_II() {

  n_neighbor_sites_ii_ = 0;
  Tubulin *site{GetActiveHead()->site_};
  Microtubule *neighb_mt{site->mt_->neighbor_};
  int site_coord{site->mt_->coord_ + site->index_};
  // Scan through all potential neighbors; only add unoccupied to list
  for (int delta{-dist_cutoff_}; delta <= dist_cutoff_; delta++) {
    int i_neighb = (site_coord - neighb_mt->coord_) + delta;
    // Skip loop iteration if i_neighb is negative
    if (i_neighb < 0) {
      continue;
    }
    // End scan once last site has been checked
    if (i_neighb > neighb_mt->n_sites_ - 1) {
      return;
    }
    Tubulin *neighb = &neighb_mt->lattice_[i_neighb];
    if (!neighb->occupied_) {
      neighbor_sites_ii_[n_neighbor_sites_ii_++] = neighb;
    }
  }
}

void AssociatedProtein::UpdateNeighbors_Bind_I_Teth() {

  if (!tethered_ or heads_active_ != 0) {
    wally_->ErrorExit("AssociatedProtein::UpdateNeighborSites_I_Teth()");
  }
  n_neighbor_sites_i_teth_ = 0;
  int teth_cutoff{properties_->kinesin4.teth_cutoff_};
  int comp_cutoff{properties_->kinesin4.comp_cutoff_};
  double stalk_coord{motor_->GetStalkCoordinate()};
  // Scan through all potential neighbor sites; add unoccupied to list
  for (int i_mt{0}; i_mt < parameters_->microtubules.count; i_mt++) {
    Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
    // Scan 1 site past teth_cutoff to account for rounding of stalk_coord
    for (int delta{-teth_cutoff + 1}; delta <= teth_cutoff + 1; delta++) {
      int i_neighb{(int)(stalk_coord - mt->coord_) + delta};
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
        int x_dub{(int)(2 * fabs(stalk_coord - neighb_coord))};
        if (x_dub >= 2 * comp_cutoff and x_dub <= 2 * teth_cutoff) {
          neighbor_sites_i_teth_[n_neighbor_sites_i_teth_++] = neighb;
        }
      }
    }
  }
}

void AssociatedProtein::UpdateNeighbors_Bind_II_Teth() {

  if (!tethered_ or heads_active_ != 1) {
    wally_->ErrorExit("AssociatedProtein::UpdateNeighborSites_I_Teth()");
  }
  n_neighbor_sites_ii_teth_ = 0;
  int teth_cutoff{properties_->kinesin4.teth_cutoff_};
  int comp_cutoff{properties_->kinesin4.comp_cutoff_};
  double stalk_coord{motor_->GetStalkCoordinate()};
  Tubulin *site{GetActiveHead()->site_};
  Microtubule *neighb_mt{site->mt_->neighbor_};
  int site_coord{site->mt_->coord_ + site->index_};
  // Scan through all potential neighbors; only add unoccupied to list
  for (int delta{-dist_cutoff_}; delta <= dist_cutoff_; delta++) {
    int i_neighb{(site_coord - neighb_mt->coord_) + delta};
    // Start index at first site (0) if i_neighb is 0 or neg
    if (i_neighb < 0) {
      continue;
    }
    // End scan once last site (mt_length - 1) has been checked
    if (i_neighb > neighb_mt->n_sites_ - 1) {
      return;
    }
    Tubulin *neighb{&neighb_mt->lattice_[i_neighb]};
    if (!neighb->occupied_) {
      int neighb_coord{i_neighb + neighb_mt->coord_};
      double new_anchor_coord{(double)(site_coord + neighb_coord) / 2};
      int x_dub{(int)(2 * fabs(new_anchor_coord - stalk_coord))};
      if (x_dub >= 2 * comp_cutoff and x_dub <= 2 * teth_cutoff) {
        neighbor_sites_ii_teth_[n_neighbor_sites_ii_teth_++] = neighb;
      }
    }
  }
}

double AssociatedProtein::GetWeight_Bind_II(Tubulin *neighbor) {

  Tubulin *site{GetActiveHead()->site_};
  int site_coord{site->mt_->coord_ + site->index_};
  int neighbor_coord{neighbor->mt_->coord_ + neighbor->index_};
  int x_dist{abs(site_coord - neighbor_coord)};
  int n_neighbs{neighbor->GetPRC1NeighborCount()};
  return weight_bind_ii_[n_neighbs][x_dist];
}

double AssociatedProtein::GetWeight_Bind_I_Teth(Tubulin *neighbor) {

  double stalk_coord{motor_->GetStalkCoordinate()};
  int site_coord{neighbor->index_ + neighbor->mt_->coord_};
  int x_dub{(int)(2 * fabs(stalk_coord - site_coord))};
  int n_neighbs{neighbor->GetPRC1NeighborCount()};
  return weight_bind_i_teth_[n_neighbs][x_dub];
}

double AssociatedProtein::GetWeight_Bind_II_Teth(Tubulin *neighbor) {

  // Get number of PRC1 neighbors this site has
  int n_neighbs{neighbor->GetPRC1NeighborCount()};
  // Get x_dist of xlink extension if 2nd head were to bind to this neighb
  Tubulin *site{GetActiveHead()->site_};
  int site_coord{site->index_ + site->mt_->coord_};
  int neighb_coord{neighbor->index_ + neighbor->mt_->coord_};
  int x{abs(site_coord - neighb_coord)};
  // Get current x_dist_dub of tether
  double stalk_coord{motor_->GetStalkCoordinate()};
  int x_dub{(int)(2 * fabs(stalk_coord - site_coord))};
  // Get NEW x_dist_dub if 2nd head were to bind to this neighb
  double anchor_coord_post{(double)(neighb_coord + site_coord) / 2};
  int x_dub_post{(int)(2 * fabs(stalk_coord - anchor_coord_post))};
  // If extended, ...
  if (x_dub >= 2 * properties_->kinesin4.rest_dist_) {
    // ... further extension is going away FROM rest
    if (x_dub_post > x_dub) {
      return weight_bind_ii_fr_teth_[n_neighbs][x_dub][x];
    } else {
      return weight_bind_ii_to_teth_[n_neighbs][x_dub][x];
    }
  }
  // Else, if compressed ...
  else {
    // ... further extension is going TO rest
    if (x_dub_post > x_dub) {
      return weight_bind_ii_to_teth_[n_neighbs][x_dub][x];
    } else {
      return weight_bind_ii_fr_teth_[n_neighbs][x_dub][x];
    }
  }
}

double AssociatedProtein::GetTotalWeight_Bind_II() {

  double tot_weight{0.0};
  UpdateNeighbors_Bind_II();
  for (int i_neighb{0}; i_neighb < n_neighbor_sites_ii_; i_neighb++) {
    tot_weight += GetWeight_Bind_II(neighbor_sites_ii_[i_neighb]);
  }
  return tot_weight;
}

double AssociatedProtein::GetTotalWeight_Bind_I_Teth() {

  double tot_weight{0.0};
  UpdateNeighbors_Bind_I_Teth();
  for (int i_neighb{0}; i_neighb < n_neighbor_sites_i_teth_; i_neighb++) {
    tot_weight += GetWeight_Bind_I_Teth(neighbor_sites_i_teth_[i_neighb]);
  }
  return tot_weight;
}

double AssociatedProtein::GetTotalWeight_Bind_II_Teth() {

  double tot_weight{0.0};
  UpdateNeighbors_Bind_II_Teth();
  for (int i_neighb{0}; i_neighb < n_neighbor_sites_ii_teth_; i_neighb++) {
    tot_weight += GetWeight_Bind_II_Teth(neighbor_sites_ii_teth_[i_neighb]);
  }
  return tot_weight;
}

Tubulin *AssociatedProtein::GetWeightedSite_Bind_II() {

  double weight_tot{GetTotalWeight_Bind_II()};
  double ran{properties_->gsl.GetRanProb()};
  double p_cum{0.0};
  for (int i_neighb{0}; i_neighb < n_neighbor_sites_ii_; i_neighb++) {
    Tubulin *neighb{neighbor_sites_ii_[i_neighb]};
    p_cum += GetWeight_Bind_II(neighb) / weight_tot;
    if (ran < p_cum) {
      return neighb;
    }
  }
  printf("No neighb for bind_ii??\n");
  return nullptr;
}

Tubulin *AssociatedProtein::GetWeightedSite_Bind_I_Teth() {

  double weight_tot{GetTotalWeight_Bind_I_Teth()};
  double ran{properties_->gsl.GetRanProb()};
  double p_cum{0.0};
  for (int i_neighb{0}; i_neighb < n_neighbor_sites_i_teth_; i_neighb++) {
    Tubulin *neighb{neighbor_sites_i_teth_[i_neighb]};
    p_cum += GetWeight_Bind_I_Teth(neighb) / weight_tot;
    if (ran < p_cum) {
      return neighb;
    }
  }
  printf("No neighb for bind_i_teth??\n");
  return nullptr;
}

Tubulin *AssociatedProtein::GetWeightedSite_Bind_II_Teth() {

  double weight_tot{GetTotalWeight_Bind_II_Teth()};
  double ran{properties_->gsl.GetRanProb()};
  double p_cum{0.0};
  for (int i_neighb{0}; i_neighb < n_neighbor_sites_ii_teth_; i_neighb++) {
    Tubulin *neighb{neighbor_sites_ii_teth_[i_neighb]};
    p_cum += GetWeight_Bind_II_Teth(neighb) / weight_tot;
    if (ran < p_cum) {
      return neighb;
    }
  }
  printf("No neighb for bind_ii_teth??\n");
  return nullptr;
}