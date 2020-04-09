#include "associated_protein.h"
#include "master_header.h"

int AssociatedProtein::Monomer::GetPRC1NeighborCount() {

  if (site_ == nullptr) {
    return 0;
  } else {
    return site_->GetPRC1NeighborCount();
  }
}

double AssociatedProtein::Monomer::GetCoord() { return site_->GetCoord(); }

AssociatedProtein::AssociatedProtein() {}

void AssociatedProtein::Initialize(system_parameters *parameters,
                                   system_properties *properties, int id) {

  id_ = id;
  wally_ = &properties->wallace;
  parameters_ = parameters;
  properties_ = properties;
  SetParameters();
  InitializeNeighborLists();
}

void AssociatedProtein::SetParameters() {

  // Convert physical quantities from nm to n_sites where applicable
  double site_size{parameters_->microtubules.site_size};
  double kbT{parameters_->kbT};
  // nm * (1 site / X nm) -> n_sites
  r_0_ = parameters_->xlinks.r_0 / site_size;
  // pN / nm * (X nm / 1 site) -> pN / n_sites
  k_spring_ = parameters_->xlinks.k_spring * site_size;
  rest_dist_ = 0;
  dist_cutoff_ = properties_->prc1.dist_cutoff_;
}

void AssociatedProtein::InitializeNeighborLists() {

  int n_mts{parameters_->microtubules.count};
  int teth_cutoff{properties_->kinesin4.teth_cutoff_};
  neighbors_bind_ii_.resize((n_mts - 1) * (2 * dist_cutoff_ + 1));
  neighbors_bind_i_teth_.resize(n_mts * (2 * teth_cutoff + 1));
  neighbors_bind_ii_teth_.resize((n_mts - 1) * (2 * dist_cutoff_ + 1));
}

AssociatedProtein::Monomer *AssociatedProtein::GetActiveHead() {

  if (heads_active_ != 1) {
    wally_->ErrorExit("AssociatedProtein::GetActiveHead [1]");
  }
  if (head_one_.site_ != nullptr) {
    return &head_one_;
  } else if (head_two_.site_ != nullptr) {
    return &head_two_;
  } else {
    wally_->ErrorExit("AssociatedProtein::GetActiveHead [2]");
  }
}

double AssociatedProtein::GetAnchorCoordinate() {

  // If single bound, use that head; assume other's diffusion avgs out
  if (heads_active_ == 1) {
    return GetActiveHead()->GetCoord();
  }
  // If double bound, use avg of head's indices
  else if (heads_active_ == 2) {
    return (head_one_.GetCoord() + head_two_.GetCoord()) / 2;
  } else {
    wally_->ErrorExit("AssociatedProtein::GetAnchorCoordinate()");
  }
}

double AssociatedProtein::GetExtensionForce(Tubulin *site) {

  if (heads_active_ == 0) {
    wally_->ErrorExit("AssociatedProtein::GetExtensionForce()");
  }
  if (heads_active_ == 1) {
    return 0.0;
  }
  double force_x{-1 * k_spring_ * extension_ * cosine_};
  // Cosine values are defined assuming the objects experiencing the spring
  // force are the RIGHT of the xlink head, i.e. in the positive direction.
  // Since our cosine values do not go negative, we must multiply by -1 if the
  // site is to the LEFT of the xlink head, i.e., in the negative direction.
  double site_coord{site->GetCoord()};
  double spring_center{GetAnchorCoordinate()};
  if (site_coord < spring_center) {
    force_x *= -1;
  }
  return force_x;
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
    }
    // Otherwise choose appropriate direction for given site
    else if (heads_active_ == 2) {
      if (x_dist_ == 0) {
        // Need to choose site which has extension of dist_cutoff + 1, not just
        // 1, i.e., need to choose site that will result in larger extension
        double site_offset{site->mt_->coord_ - site->mt_->neighbor_->coord_};
        if (site_offset < 0.0) {
          return -1;
        } else {
          return 1;
        }
      } else {
        double site_coord{site->GetCoord()};
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
    double r_x{fabs(head_one_.GetCoord() - head_two_.GetCoord())};
    double x_dist = std::round(r_x);
    /*
    if (head_one_.site_->index_ == head_two_.site_->index_) {
      printf("i_one: %i | i_two: %i\n", head_one_.site_->index_,
             head_two_.site_->index_);
      printf("coord one: %g | coord two: %g\n", head_one_.GetCoord(),
             head_two_.GetCoord());
      printf("%g -> %g\n", r_x, x_dist);
    }
    */
    if (x_dist > dist_cutoff_) {
      ForceUnbind();
    } else {
      x_dist_ = x_dist;
      // No extension (site_offset): 0
      // Smaller extensions (x_dist - site_offset): 1 to dist_cutoff_
      // Larger extensions (x_dist + site_offset): cutoff + 1 to 2*cutoff
      if (r_x > x_dist and x_dist != 0) {
        x_dist_ += dist_cutoff_;
      }
      cosine_ = properties_->prc1.possible_cosines_[x_dist_];
      extension_ = properties_->prc1.possible_extensions_[x_dist_];
    }
  } else {
    x_dist_ = 0;
    cosine_ = 0.0;
    extension_ = 0.0;
  }
  if (tethered_) {
    motor_->UpdateExtension();
  }
}

void AssociatedProtein::ForceUnbind() {

  // Xlinks with satellites behave as if untethered (no energy from tether)
  if (!tethered_ or HasSatellite()) {
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
    printf("not supported\n");
    wally_->ErrorExit("AP::ForceUnbind()");
    /*
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
    */
  }
  heads_active_--;
  x_dist_ = 0;
  cosine_ = 0;
  extension_ = 0;
  if (tethered_) {
    // Update motor extension since anchor coord has changed
    motor_->UpdateExtension();
  }
  properties_->prc1.FlagForUpdate();
  properties_->kinesin4.FlagForUpdate();
}

void AssociatedProtein::UntetherSatellite() {

  if (!HasSatellite()) {
    return;
  }
  properties_->kinesin4.RemoveFromActive(motor_);
  properties_->kinesin4.FlagForUpdate();
  // Update motor details
  motor_->tethered_ = false;
  motor_->xlink_ = nullptr;
  // Update xlink details
  tethered_ = false;
  motor_ = nullptr;
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
  /*
  double anchor_coord = GetAnchorCoordinate();
  double stalk_coord = motor_->GetStalkCoordinate();
  double x_dub = 2 * fabs(anchor_coord - stalk_coord);
  double x_dub_rest = 2 * properties_->kinesin4.rest_dist_;
  double E_current{};
  if (x_dub > x_dub_rest) {
    E_current = 0.5 * motor_->k_spring_eff_ * (x_dub / 2) * (x_dub / 2);
  } else {
    E_current = 0.5 * motor_->k_slack_eff_ * (x_dub / 2) * (x_dub / 2);
  }
  double x_dub_one{2 * fabs(head_one_.GetCoord() - stalk_coord)};
  double x_dub_two{2 * fabs(head_two_.GetCoord() - stalk_coord)};
  double E_post_one{};
  double E_post_two{};



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
  */
}

Tubulin *AssociatedProtein::GetSiteFartherFromTethRest() {

  /*
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
  */
}

void AssociatedProtein::UpdateNeighbors_Bind_II() {

  n_neighbors_bind_ii_ = 0;
  // printf("dist cutoff is %g\n", dist_cutoff_);
  Monomer *bound_head{GetActiveHead()};
  double bound_coord{bound_head->GetCoord()};
  Microtubule *adjacent_mt{bound_head->site_->mt_->neighbor_};
  // Scan through all potential neighbors; only add unoccupied to list
  for (int delta{-(dist_cutoff_ + 1)}; delta <= dist_cutoff_ + 1; delta++) {
    int i_neighb{int(bound_coord - adjacent_mt->coord_) + delta};
    // printf("i_neighb = %i\n", i_neighb);
    // Skip loop iteration if i_neighb is negative
    if (i_neighb < 0) {
      continue;
    }
    // End scan once last site has been checked
    if (i_neighb > adjacent_mt->n_sites_ - 1) {
      return;
    }
    /*
    printf("delta = %i\n", delta);
    printf("i_neighb = %i\n", i_neighb);
    */
    // printf("occu?\n");
    Tubulin *neighb = &adjacent_mt->lattice_[i_neighb];
    if (!neighb->occupied_) {
      // printf("not occu!\n");
      double r_x{fabs(bound_coord - neighb->GetCoord())};
      double x_dist{std::round(r_x)};
      // printf("x_dist is %g\n", x_dist);
      // printf("r_x = %g\n", r_x);
      if (x_dist <= dist_cutoff_) {
        // printf("n_neighbs = %i\n", n_neighbors_bind_ii_);
        /*
        if (r_x > x_dist and x_dist != 0.0) {
          x_dist += dist_cutoff_;
        }
        printf("added with x_dist = %g\n", x_dist);
        */
        neighbors_bind_ii_[n_neighbors_bind_ii_++] = neighb;
      }
    }
  }
}

void AssociatedProtein::UpdateNeighbors_Bind_I_Teth() {

  /*
  if (!tethered_ or heads_active_ != 0) {
    wally_->ErrorExit("AssociatedProtein::UpdateNeighborSites_I_Teth()");
  }
  n_neighbors_bind_i_teth_ = 0;
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
          neighbors_bind_i_teth_[n_neighbors_bind_i_teth_++] = neighb;
        }
      }
    }
  }
  */
}

void AssociatedProtein::UpdateNeighbors_Bind_II_Teth() {

  /*
  if (!tethered_ or heads_active_ != 1) {
    wally_->ErrorExit("AssociatedProtein::UpdateNeighborSites_I_Teth()");
  }
  n_neighbors_bind_ii_teth_ = 0;
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
        neighbors_bind_ii_teth_[n_neighbors_bind_ii_teth_++] = neighb;
      }
    }
  }
  */
}

double AssociatedProtein::GetWeight_DiffuseToRest(Monomer *head) {

  // Check site we will be diffusing to
  int dx_rest{GetDirectionTowardRest(head->site_)};
  int i_neighb{head->site_->index_ + dx_rest};
  // Check if site exists
  if (i_neighb > head->site_->mt_->n_sites_ - 1 or i_neighb < 0) {
    return 0.0;
  }
  // Check if site is occupied
  if (head->site_->mt_->lattice_[i_neighb].occupied_) {
    return 0.0;
  }
  int x_to{x_dist_ - 1};
  // x_dist == 1 is handled by diffuse_fr for x == 0
  if (x_dist_ == 0) {
    x_to = dist_cutoff_ + 1;
  } else if (x_dist_ == dist_cutoff_ + 1) {
    x_to = 0;
  }
  // We consider diffusing towards rest to be an unbind-like event,
  // since it results in the spring becoming less stretched
  double weight_current{properties_->prc1.weight_spring_unbind_[x_dist_]};
  double weight_post{properties_->prc1.weight_spring_unbind_[x_to]};
  // Want -(U_post - U_current) in exponential, so divide current by post
  double weight_spring{weight_current / weight_post};
  if (properties_->prc1.verbosity_ >= 2) {
    wally_->Log("xlink %i has spring weight %g\n", id_, weight_spring);
  }
  int n_neighbs{head->site_->GetPRC1NeighborCount()};
  return weight_spring * properties_->prc1.weight_neighb_unbind_[n_neighbs];
}

double AssociatedProtein::GetWeight_DiffuseFrRest(Monomer *head) {

  // Check if xlink is at extension cutoff
  if (x_dist_ == dist_cutoff_ or x_dist_ == 2 * dist_cutoff_) {
    return 0.0;
  }
  // Monomer *other_head{head->GetOtherHead()};
  // Check site we will be diffusing to
  int dx_rest{GetDirectionTowardRest(head->site_)};
  int i_neighb{head->site_->index_ - dx_rest};
  // Check if site exists
  if (i_neighb > head->site_->mt_->n_sites_ - 1 or i_neighb < 0) {
    return 0.0;
  }
  // Check if site is occupied
  if (head->site_->mt_->lattice_[i_neighb].occupied_) {
    return 0.0;
  }
  int x_fr{x_dist_ + 1};
  // We consider diffusing from rest to be an bind-like event,
  // since it results in the spring becoming more stretched
  double weight_current{properties_->prc1.weight_spring_bind_[x_dist_]};
  double weight_post{properties_->prc1.weight_spring_bind_[x_fr]};
  // Want U_post - U_current in the exponential; divide accordingly
  double weight_spring{weight_post / weight_current};
  if (properties_->prc1.verbosity_ >= 2) {
    wally_->Log("xlink %i has spring weight %g\n", id_, weight_spring);
  }
  int n_neighbs{head->site_->GetPRC1NeighborCount()};
  return weight_spring * properties_->prc1.weight_neighb_unbind_[n_neighbs];
}

double AssociatedProtein::GetWeight_Unbind_II(Monomer *head) {

  int n_neighbs{head->GetPRC1NeighborCount()};
  return properties_->prc1.weight_spring_unbind_[x_dist_] *
         properties_->prc1.weight_neighb_unbind_[n_neighbs];
}

double AssociatedProtein::GetWeight_Bind_II(Tubulin *neighbor) {

  Tubulin *site{GetActiveHead()->site_};
  if (site->mt_ == neighbor->mt_) {
    wally_->ErrorExit("AP::GetWeight_Bind_II()");
  }
  double r_x{fabs(neighbor->GetCoord() - site->GetCoord())};
  int x_dist{(int)std::round(r_x)};
  // By convention, smaller extensions (w/ x = x_base - site_offset()) range
  // from 0 to dist_cutoff_. Larger extensions from dist_cutoff_ to 2*cutoff_
  if (r_x > x_dist and x_dist != 0) {
    x_dist += dist_cutoff_;
  }
  int n_neighbs{neighbor->GetPRC1NeighborCount()};
  return properties_->prc1.weight_spring_bind_[x_dist] *
         properties_->prc1.weight_neighb_bind_[n_neighbs];
}

double AssociatedProtein::GetWeight_Bind_I_Teth(Tubulin *neighbor) {

  /*
  double stalk_coord{motor_->GetStalkCoordinate()};
  int site_coord{neighbor->index_ + neighbor->mt_->coord_};
  int x_dub{(int)(2 * fabs(stalk_coord - site_coord))};
  int n_neighbs{neighbor->GetPRC1NeighborCount()};
  return weight_bind_i_teth_[n_neighbs][x_dub];
  */
}

double AssociatedProtein::GetWeight_Bind_II_Teth(Tubulin *neighbor) {

  /*
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
  */
}

double AssociatedProtein::GetTotalWeight_Bind_II() {

  double tot_weight{0.0};
  UpdateNeighbors_Bind_II();
  for (int i_neighb{0}; i_neighb < n_neighbors_bind_ii_; i_neighb++) {
    tot_weight += GetWeight_Bind_II(neighbors_bind_ii_[i_neighb]);
  }
  return tot_weight;
}

double AssociatedProtein::GetTotalWeight_Bind_I_Teth() {

  double tot_weight{0.0};
  UpdateNeighbors_Bind_I_Teth();
  for (int i_neighb{0}; i_neighb < n_neighbors_bind_i_teth_; i_neighb++) {
    tot_weight += GetWeight_Bind_I_Teth(neighbors_bind_i_teth_[i_neighb]);
  }
  return tot_weight;
}

double AssociatedProtein::GetTotalWeight_Bind_II_Teth() {

  double tot_weight{0.0};
  UpdateNeighbors_Bind_II_Teth();
  for (int i_neighb{0}; i_neighb < n_neighbors_bind_ii_teth_; i_neighb++) {
    tot_weight += GetWeight_Bind_II_Teth(neighbors_bind_ii_teth_[i_neighb]);
  }
  return tot_weight;
}

Tubulin *AssociatedProtein::GetWeightedSite_Bind_II() {

  // GetTotalWeight also updates relevant neighbor list
  double weight_tot{GetTotalWeight_Bind_II()};
  // printf("Total weight is %g (%i neighbs)\n", weight_tot,
  // n_neighbors_bind_ii_);
  double ran{properties_->gsl.GetRanProb()};
  // printf("Rolled %g\n", ran);
  double p_cum{0.0};
  for (int i_neighb{0}; i_neighb < n_neighbors_bind_ii_; i_neighb++) {
    Tubulin *neighb{neighbors_bind_ii_[i_neighb]};
    p_cum += GetWeight_Bind_II(neighb) / weight_tot;
    // printf("p_cum is %g after entry #%i\n", p_cum, i_neighb);
    if (ran < p_cum) {
      // printf("Chose entry #%i\n", i_neighb);
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
  for (int i_neighb{0}; i_neighb < n_neighbors_bind_i_teth_; i_neighb++) {
    Tubulin *neighb{neighbors_bind_i_teth_[i_neighb]};
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
  for (int i_neighb{0}; i_neighb < n_neighbors_bind_ii_teth_; i_neighb++) {
    Tubulin *neighb{neighbors_bind_ii_teth_[i_neighb]};
    p_cum += GetWeight_Bind_II_Teth(neighb) / weight_tot;
    if (ran < p_cum) {
      return neighb;
    }
  }
  printf("No neighb for bind_ii_teth??\n");
  return nullptr;
}