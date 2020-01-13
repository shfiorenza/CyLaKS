#include "associated_protein.h"
#include "master_header.h"

AssociatedProtein::AssociatedProtein() {}

void AssociatedProtein::Initialize(system_parameters *parameters,
                                   system_properties *properties, int id) {

  id_ = id;
  parameters_ = parameters;
  properties_ = properties;
  SetParameters();
  CalculateCutoffs();
  InitializeNeighborLists();
  InitializeLookupTables();
}

void AssociatedProtein::SetParameters() {

  r_0_ = parameters_->xlinks.r_0;
  k_spring_ = parameters_->xlinks.k_spring;
}

void AssociatedProtein::CalculateCutoffs() {

  int site_size = parameters_->microtubules.site_size;
  double kbT = parameters_->kbT;
  double r_y = parameters_->microtubules.y_dist;
  /* First, calculate rest_dist_ in number of sites */
  int rough_rest_dist = sqrt(r_0_ * r_0_ - r_y * r_y) / site_size;
  double rest_scan[3];
  double scan_force[3];
  for (int i_scan = -1; i_scan <= 1; i_scan++) {
    rest_scan[i_scan + 1] = rough_rest_dist + (i_scan * 0.5);
    double rest_scan_length = rest_scan[i_scan + 1] * site_size;
    double r_scan = sqrt(r_y * r_y + rest_scan_length * rest_scan_length);
    scan_force[i_scan + 1] = (r_scan - r_0_) * k_spring_;
  }
  double min_force = 100;
  for (int i_scan = -1; i_scan <= 1; i_scan++) {
    double force = fabs(scan_force[i_scan + 1]);
    if (force < min_force) {
      min_force = force;
      rest_dist_ = rest_scan[i_scan + 1];
    }
  }
  /* Finally, calculate extension distance cutoff */
  for (int x_dist = (int)rest_dist_; x_dist < 1000; x_dist++) {
    int r_x = x_dist * site_size;
    double r = sqrt(r_y * r_y + r_x * r_x);
    double dr = r - r_0_;
    double U = (k_spring_ / 2) * dr * dr;
    double boltzmann_weight = exp(U / (2 * kbT));
    if (boltzmann_weight > 100) {
      dist_cutoff_ = x_dist;
      break;
    }
  }
}

void AssociatedProtein::InitializeNeighborLists() {

  int n_mts = parameters_->microtubules.count;
  int teth_cutoff = properties_->kinesin4.teth_cutoff_;
  neighbor_sites_.resize((n_mts - 1) * (2 * dist_cutoff_ + 1));
  teth_neighbor_sites_.resize(n_mts * (2 * teth_cutoff + 1));
  teth_neighbor_sites_ii_.resize((n_mts - 1) * (2 * dist_cutoff_ + 1));
}

void AssociatedProtein::InitializeLookupTables() {

  /* Construct rly basic lookup tables first */
  double r_y = parameters_->microtubules.y_dist;
  double site_size = parameters_->microtubules.site_size;
  extension_lookup_.resize(dist_cutoff_ + 1);
  cosine_lookup_.resize(dist_cutoff_ + 1);
  for (int x_dist = 0; x_dist <= dist_cutoff_; x_dist++) {
    double r_x = site_size * x_dist;
    double r = sqrt(r_x * r_x + r_y * r_y);
    extension_lookup_[x_dist] = r - r_0_;
    cosine_lookup_[x_dist] = r_x / r;
  }
  /* Next, construct lookup table for bind_ii weights */
  int max_neighbs = 2;
  double kbT = parameters_->kbT;
  double E_int = -1 * parameters_->xlinks.interaction_energy;
  weights_bind_ii_.resize(max_neighbs + 1);
  for (int n_neighbs(0); n_neighbs <= max_neighbs; n_neighbs++) {
    double E_neighbs = n_neighbs * E_int; // in units of kBT already
    double weight_neighb = exp(-E_neighbs / 2);
    weights_bind_ii_[n_neighbs].resize(dist_cutoff_ + 1);
    for (int x(0); x <= dist_cutoff_; x++) {
      double r_x = (double)x * site_size;
      double r = sqrt(r_y * r_y + r_x * r_x);
      double dr = r - r_0_;
      double U = (k_spring_ / 2) * dr * dr;
      double weight_spring = exp(-U / (2 * kbT));
      double weight_tot = weight_neighb * weight_spring;
      weights_bind_ii_[n_neighbs][x] = weight_tot;
      //			printf("weight is %g for n = %i, x = %i\n",
      // weight_tot, 					n_neighbs, x);
    }
  }
  /* Next, construct lookup table for bind_i_teth weights */
  int teth_cutoff = properties_->kinesin4.teth_cutoff_;
  int comp_cutoff = properties_->kinesin4.comp_cutoff_;
  double k_teth = parameters_->motors.k_spring;
  double k_slack = parameters_->motors.k_slack;
  double r_0_teth = parameters_->motors.r_0;
  double r_y_teth = parameters_->microtubules.y_dist / 2;
  weights_bind_i_teth_.resize(max_neighbs + 1);
  for (int n_neighbs(0); n_neighbs <= max_neighbs; n_neighbs++) {
    double U_neighbs = n_neighbs * E_int;
    double weight_neighbs = exp(-U_neighbs / 2);
    weights_bind_i_teth_[n_neighbs].resize(2 * teth_cutoff + 1);
    for (int x_dub = 0; x_dub <= 2 * teth_cutoff; x_dub++) {
      double r_x = ((double)x_dub) * site_size / 2;
      double r = sqrt(r_y_teth * r_y_teth + r_x * r_x);
      double dr = r - r_0_teth;
      double U_teth(0);
      if (dr < 0)
        U_teth = (k_slack / 2) * dr * dr;
      else
        U_teth = (k_teth / 2) * dr * dr;
      double weight_teth = exp(-U_teth / (2 * kbT));
      if (x_dub < 2 * comp_cutoff)
        weight_teth = 0;
      double weight_tot = weight_neighbs * weight_teth;
      weights_bind_i_teth_[n_neighbs][x_dub] = weight_tot;
      //			printf("weight for 2x = %i is %g\n", x_dub,
      // weight_tot);
    }
  }
  /* Next, construct lookup table for bind_ii_teth weights */
  // XXX input: CURRENT x_dub and PROPOSED x_dist XXX
  double rest_dist_teth = properties_->kinesin4.motors_[0].rest_dist_;
  weights_bind_ii_to_teth_.resize(max_neighbs + 1);
  weights_bind_ii_fr_teth_.resize(max_neighbs + 1);
  for (int n_neighbs(0); n_neighbs <= max_neighbs; n_neighbs++) {
    double U_neighbs = n_neighbs * E_int;
    double weight_neighbs = exp(-U_neighbs / 2);
    weights_bind_ii_to_teth_[n_neighbs].resize(2 * teth_cutoff + 1);
    weights_bind_ii_fr_teth_[n_neighbs].resize(2 * teth_cutoff + 1);
    for (int x_dub(0); x_dub <= 2 * teth_cutoff; x_dub++) {
      // Calc x-distances (in nm) for tether
      double r_x_teth = ((double)x_dub) * site_size / 2;
      // Calc total r values
      double r_teth = sqrt(r_x_teth * r_x_teth + r_y_teth * r_y_teth);
      // Calc tether exts for current dist and stepping to/from rest
      double dr_teth = r_teth - r_0_teth;
      double U_teth;
      if (dr_teth < 0) {
        U_teth = (k_slack / 2) * dr_teth * dr_teth;
      } else {
        U_teth = (k_teth / 2) * dr_teth * dr_teth;
      }
      weights_bind_ii_to_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      weights_bind_ii_fr_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      for (int x(0); x <= dist_cutoff_; x++) {
        double r_x = ((double)x) * site_size;
        double r = sqrt(r_y * r_y + r_x * r_x);
        double dr = r - r_0_;
        double U = (k_spring_ / 2) * dr * dr;
        double weight = exp(-U / (2 * kbT));
        // change in teth extension if 2nd xlink head were to bind
        double dx_teth = ((double)x) / 2;
        double dr_x_teth = dx_teth * site_size;
        double r_x_teth_to, r_x_teth_fr;
        if (dr_teth < 0) {
          r_x_teth_to = r_x_teth + dr_x_teth;
          r_x_teth_fr = r_x_teth - dr_x_teth;
        } else {
          r_x_teth_to = r_x_teth - dr_x_teth;
          r_x_teth_fr = r_x_teth + dr_x_teth;
        }
        double r_teth_to =
            sqrt(r_x_teth_to * r_x_teth_to + r_y_teth * r_y_teth);
        double r_teth_fr =
            sqrt(r_x_teth_fr * r_x_teth_fr + r_y_teth * r_y_teth);
        double dr_teth_to = r_teth_to - r_0_teth;
        double dr_teth_fr = r_teth_fr - r_0_teth;
        double U_teth_to;
        if (dr_teth_to < 0) {
          U_teth_to = (k_slack / 2) * dr_teth_to * dr_teth_to;
        } else {
          U_teth_to = (k_teth / 2) * dr_teth_to * dr_teth_to;
        }
        double U_teth_fr;
        if (dr_teth_fr < 0) {
          U_teth_fr = (k_slack / 2) * dr_teth_fr * dr_teth_fr;
        } else {
          U_teth_fr = (k_teth / 2) * dr_teth_fr * dr_teth_fr;
        }
        double dU_to_teth = U_teth_to - U_teth;
        double dU_fr_teth = U_teth_fr - U_teth;
        double weight_to_teth = exp(-dU_to_teth / (2 * kbT));
        double weight_fr_teth = exp(-dU_fr_teth / (2 * kbT));
        double tot_weight_to(0);
        int x_dub_to = 2 * r_x_teth_to / site_size;
        if (x_dub_to >= 2 * comp_cutoff && x_dub_to <= 2 * teth_cutoff) {
          tot_weight_to = weight_neighbs * weight_to_teth * weight;
        }
        double tot_weight_fr(0);
        int x_dub_fr = 2 * r_x_teth_fr / site_size;
        if (x_dub_fr >= 2 * comp_cutoff && x_dub_fr <= 2 * teth_cutoff) {
          tot_weight_fr = weight_neighbs * weight_fr_teth * weight;
        }
        weights_bind_ii_to_teth_[n_neighbs][x_dub][x] = tot_weight_to;
        weights_bind_ii_fr_teth_[n_neighbs][x_dub][x] = tot_weight_fr;
        /*
        printf("for n=%i, 2x=%i, x=%i, w_to = %g\n",
                        n_neighbs, x_dub, x, weights_bind_ii_to_teth_
                        [n_neighbs][x_dub][x]);
        printf("for n=%i, 2x=%i, x=%i, w_from = %g\n",
                        n_neighbs, x_dub, x, weights_bind_ii_fr_teth_
                        [n_neighbs][x_dub][x]);
        */
      }
    }
  }
}

AssociatedProtein::Monomer *AssociatedProtein::GetActiveHead() {

  if (heads_active_ == 1) {
    if (head_one_.site_ != nullptr)
      return &head_one_;
    else if (head_two_.site_ != nullptr)
      return &head_two_;
    else {
      printf("what in get active head site\n");
      exit(1);
    }
  } else {
    printf("not cool bro...not single bound \n");
    exit(1);
  }
}

Tubulin *AssociatedProtein::GetActiveHeadSite() {

  if (heads_active_ == 1) {
    if (head_one_.site_ != nullptr)
      return head_one_.site_;
    else if (head_two_.site_ != nullptr)
      return head_two_.site_;
    else {
      printf("what in get active head site\n");
      exit(1);
    }
  } else {
    printf("not cool bro...not single bound \n");
    exit(1);
  }
}

double AssociatedProtein::GetAnchorCoordinate() {

  // If single bound, use that head; assume other's diffusion avgs out
  if (heads_active_ == 1) {
    Tubulin *site = GetActiveHeadSite();
    int index = site->index_;
    int mt_coord = site->mt_->coord_;
    double site_coord = (double)(mt_coord + index);
    return site_coord;
  }
  // If double bound, use avg of head's indices
  else if (heads_active_ == 2) {
    int index_one = head_one_.site_->index_;
    int mt_coord_one = head_one_.site_->mt_->coord_;
    double coord_one = (double)(mt_coord_one + index_one);
    int index_two = head_two_.site_->index_;
    int mt_coord_two = head_two_.site_->mt_->coord_;
    double coord_two = (double)(mt_coord_two + index_two);
    double avg_coord = (coord_one + coord_two) / 2;
    return avg_coord;
  } else {
    printf("not NOT cool bro ... cant get anchor index: %i\n", heads_active_);
    exit(1);
  }
}

int AssociatedProtein::GetPRC1NeighbCount(Monomer *head) {

  int n_neighbs = 0;
  Tubulin *site = head->site_;
  if (site != nullptr) {
    int i_plus = site->mt_->plus_end_;
    int i_minus = site->mt_->minus_end_;
    int dx = site->mt_->delta_x_;
    if (site->index_ == i_plus) {
      if (site->mt_->lattice_[site->index_ - dx].xlink_head_ != nullptr)
        n_neighbs++;
    } else if (site->index_ == i_minus) {
      if (site->mt_->lattice_[site->index_ + dx].xlink_head_ != nullptr)
        n_neighbs++;
    } else {
      if (site->mt_->lattice_[site->index_ - dx].xlink_head_ != nullptr)
        n_neighbs++;
      if (site->mt_->lattice_[site->index_ + dx].xlink_head_ != nullptr)
        n_neighbs++;
    }
  }
  return n_neighbs;
}

void AssociatedProtein::UpdateNeighborSites_II() {

  n_neighbor_sites_ = 0;
  int n_mts = parameters_->microtubules.count;
  if (n_mts > 1) {
    Tubulin *site = GetActiveHeadSite();
    int i_site = site->index_;
    Microtubule *mt = site->mt_;
    int mt_length = mt->neighbor_->n_sites_;
    Microtubule *adj_mt = mt->neighbor_;
    int site_coord = mt->coord_ + i_site;
    // Scan through all potential neighbors; only add unoccupied to list
    int i_entry = 0;
    for (int dist = -dist_cutoff_; dist <= dist_cutoff_; dist++) {
      int i_neighbor = (site_coord - adj_mt->coord_) + dist;
      // Start index at first site (0) if i_neighb is 0 or neg
      if (i_neighbor < 0) {
        dist -= (i_neighbor + 1); // - 1;
      }
      // End scan once last site (mt_length - 1) has been checked
      else if (i_neighbor > mt_length - 1) {
        break;
      } else {
        Tubulin *neighbor = &adj_mt->lattice_[i_neighbor];
        if (!neighbor->occupied_) {
          n_neighbor_sites_++;
          neighbor_sites_[i_entry] = neighbor;
          i_entry++;
        }
      }
    }
  }
}

void AssociatedProtein::UpdateNeighborSites_I_Teth() {

  n_teth_neighbor_sites_ = 0;
  if (tethered_ == true && heads_active_ == 0) {
    int n_mts = parameters_->microtubules.count;
    int teth_cutoff = properties_->kinesin4.teth_cutoff_;
    int comp_cutoff = properties_->kinesin4.comp_cutoff_;
    double stalk_coord = motor_->GetStalkCoordinate();
    //		printf("stalk coord is %g\n\n", anchor_coord);
    // Scan through all potential neighbor sites; add unoccupied to list
    int i_entry = 0;
    for (int i_mt = 0; i_mt < n_mts; i_mt++) {
      Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
      int mt_length = mt->n_sites_;
      double mt_coord = mt->coord_;
      int i_stalk = stalk_coord - mt_coord;
      for (int x_dist = -teth_cutoff; x_dist <= teth_cutoff; x_dist++) {
        int i_site = i_stalk + x_dist;
        //				printf("i_site is %i (x_dist %i)\n",
        // i_site, x_dist);
        // Start index at first site (0) if site index is <= 0
        if (i_site < 0) {
          x_dist -= (i_site + 1);
        }
        // End scan at last site (mt_length - 1)
        else if (i_site > mt_length - 1) {
          break;
        } else {
          Tubulin *neighbor = &mt->lattice_[i_site];
          double site_coord = i_site + neighbor->mt_->coord_;
          double x_dist = fabs(stalk_coord - site_coord);
          int x_dist_dub = 2 * x_dist;
          if (x_dist_dub >= 2 * comp_cutoff && x_dist_dub <= 2 * teth_cutoff &&
              neighbor->occupied_ == false) {
            n_teth_neighbor_sites_++;
            teth_neighbor_sites_[i_entry] = neighbor;
            i_entry++;
          }
        }
      }
    }
  } else {
    printf("error in XLINK update teth neighbor sites\n");
    exit(1);
  }
}

void AssociatedProtein::UpdateNeighborSites_II_Teth() {

  n_teth_neighbor_sites_ii_ = 0;
  int n_mts = parameters_->microtubules.count;
  if (n_mts > 1) {
    Tubulin *site = GetActiveHeadSite();
    int i_site = site->index_;
    Microtubule *mt = site->mt_;
    int mt_length = mt->neighbor_->n_sites_;
    Microtubule *adj_mt = mt->neighbor_;
    int site_coord = mt->coord_ + i_site;
    // Scan through all potential neighbors; only add unoccupied to list
    int i_entry = 0;
    for (int dist = -dist_cutoff_; dist <= dist_cutoff_; dist++) {
      int i_neighbor = (site_coord - adj_mt->coord_) + dist;
      // Start index at first site (0) if i_neighb is 0 or neg
      if (i_neighbor < 0) {
        dist -= (i_neighbor + 1);
      }
      // End scan once last bulk site (mt_length - 1) has been checked
      else if (i_neighbor > mt_length - 1) {
        break;
      } else {
        Tubulin *neighbor = &adj_mt->lattice_[i_neighbor];
        if (neighbor->occupied_ == false) {
          n_teth_neighbor_sites_ii_++;
          teth_neighbor_sites_ii_[i_entry] = neighbor;
          i_entry++;
        }
      }
    }
  }
}

void AssociatedProtein::UpdateExtension() {

  if (heads_active_ == 2) {
    int x_dist_pre = x_dist_;
    // Calculate first head's coordinate
    int i_head_one = head_one_.site_->index_;
    int mt_coord_one = head_one_.site_->mt_->coord_;
    int coord_one = mt_coord_one + i_head_one;
    // Calculate second head's coordinate
    int i_head_two = head_two_.site_->index_;
    int mt_coord_two = head_two_.site_->mt_->coord_;
    int coord_two = mt_coord_two + i_head_two;
    // Calculate x_distance in # of sites
    int x_dist = abs(coord_one - coord_two);
    x_dist_ = x_dist;
    if (x_dist <= dist_cutoff_) {
      extension_ = extension_lookup_[x_dist_];
      cosine_ = cosine_lookup_[x_dist_];
    } else {
      ForceUnbind(x_dist_pre);
      //			printf("forced an unbind event >:O\n");
    }
  } else if (heads_active_ == 1) {
    x_dist_ = 0;
    extension_ = 0;
    cosine_ = 0;
  } else {
    printf("some kinda error in assoc. protein update_extension\n");
    exit(1);
  }
}

void AssociatedProtein::ForceUnbind(int x_dist_pre) {

  // different stuff depending on whether or not xlink is tethered
  if (tethered_ == false) {
    double ran = properties_->gsl.GetRanProb();
    if (ran < 0.5) {
      head_one_.site_->xlink_head_ = nullptr;
      head_one_.site_->occupied_ = false;
      head_one_.site_ = nullptr;
    } else {
      head_two_.site_->xlink_head_ = nullptr;
      head_two_.site_->occupied_ = false;
      head_two_.site_ = nullptr;
    }
    heads_active_--;
    x_dist_ = 0;
    extension_ = 0;
    cosine_ = 0;
  } else if (motor_->heads_active_ > 0) {
    // Update motor ext (unbinding 2nd head changes anchor coord)
    Tubulin *site_to = GetSiteFartherFromTethRest();
    int neighbs_to = site_to->GetPRC1NeighborCount();
    Tubulin *site_fr = GetSiteCloserToTethRest();
    int neighbs_fr = site_fr->GetPRC1NeighborCount();

    int x_dub_pre = motor_->x_dist_doubled_;
    double p_unbind_to =
        properties_->prc1
            .p_unbind_ii_to_teth_[neighbs_to][x_dub_pre][x_dist_pre];
    double p_unbind_from =
        properties_->prc1
            .p_unbind_ii_fr_teth_[neighbs_fr][x_dub_pre][x_dist_pre];
    double p_tot = p_unbind_to + p_unbind_from;

    double ran = properties_->gsl.GetRanProb();
    if (ran < p_unbind_to / p_tot) {
      site_to->xlink_head_ = nullptr;
      site_to->occupied_ = false;
      site_to = nullptr;
    } else {
      site_fr->xlink_head_ = nullptr;
      site_fr->occupied_ = false;
      site_fr = nullptr;
    }
    heads_active_--;
    x_dist_ = 0;
    extension_ = 0;
    cosine_ = 0;
    motor_->UpdateExtension();
  }
  // xlinks tethered to free motors diffuse as if untethered
  else {
    double ran = properties_->gsl.GetRanProb();
    if (ran < 0.5) {
      head_one_.site_->xlink_head_ = nullptr;
      head_one_.site_->occupied_ = false;
      head_one_.site_ = nullptr;
    } else {
      head_two_.site_->xlink_head_ = nullptr;
      head_two_.site_->occupied_ = false;
      head_two_.site_ = nullptr;
    }
    heads_active_--;
    x_dist_ = 0;
    extension_ = 0;
    cosine_ = 0;
  }
}

void AssociatedProtein::UntetherSatellite() {

  if (tethered_ == true) {
    // Remove satellite motor from active_ list, replace with last entry
    int i_last = properties_->kinesin4.n_active_ - 1;
    Kinesin *last_entry = properties_->kinesin4.active_[i_last];
    int i_this = motor_->active_index_;
    if (i_this != i_last) {
      properties_->kinesin4.active_[i_this] = last_entry;
      last_entry->active_index_ = i_this;
    }
    properties_->kinesin4.n_active_--;
    // Update motor details
    motor_->tethered_ = false;
    motor_->xlink_ = nullptr;
    // Update xlink details
    tethered_ = false;
    motor_ = nullptr;
  } else {
    printf("Error in xlink UntethSatellite()\n");
    exit(1);
  }
}

int AssociatedProtein::GetDirectionToRest(Tubulin *site) {

  if (heads_active_ == 1)
    return 1;
  else if (heads_active_ == 2) {
    double anchor_coord = GetAnchorCoordinate();
    int site_coord = site->index_ + site->mt_->coord_;
    if (site_coord == anchor_coord) {
      double ran = properties_->gsl.GetRanProb();
      if (ran < 0.5)
        return -1;
      else
        return 1;
    } else if (site_coord > anchor_coord)
      return -1;
    else
      return 1;
  } else
    return 0;
  /*
  else{
          printf("error in get dir. toward rest (xlink)\n");
          exit(1);
  }
  */
}

double AssociatedProtein::GetExtensionForce(Tubulin *site) {

  if (heads_active_ == 2) {
    UpdateExtension();
    // Make sure we didn't force an unbinding event
    if (heads_active_ == 2) {
      double force_mag = extension_ * k_spring_; // in pN
      double site_coord = site->index_ + site->mt_->coord_;
      double anchor_coord = GetAnchorCoordinate();
      double force;
      if (site_coord < anchor_coord)
        force = force_mag * cosine_;
      else
        force = -1 * force_mag * cosine_;
      return force;
    } else
      return 0;
  } else {
    printf("error in get ext force (xlink)\n");
    exit(1);
  }
}

double AssociatedProtein::GetBindingWeight_II(Tubulin *neighbor) {

  Tubulin *site = GetActiveHeadSite();
  Microtubule *mt = site->mt_;
  Microtubule *adj_mt = neighbor->mt_;
  if (adj_mt != mt->neighbor_) {
    printf("adj: %i, neighb: %i\n", adj_mt->index_, mt->neighbor_->index_);
    printf("why the microtubules tho (in assiociated protein GBW)\n");
    exit(1);
  }
  int offset = adj_mt->coord_ - mt->coord_;
  // Calculate distance (in x-dir.) between site and neighb in # of sites
  int i_site = site->index_;
  int i_neighbor = neighbor->index_;
  int x_dist = abs(i_neighbor - i_site + offset);
  // Get number of PRC1 neighbors the neighbor site has
  int n_neighbs = neighbor->GetPRC1NeighborCount();
  //	printf("x: %i, n_neighbs: %i\n", x_dist, n_neighbs);
  // Look up binding weight that corresponds to this x-distance
  double weight = weights_bind_ii_[n_neighbs][x_dist];
  return weight;
}

double AssociatedProtein::GetBindingWeight_I_Teth(Tubulin *neighbor) {

  double stalk_coord = motor_->GetStalkCoordinate();
  double site_coord = neighbor->index_ + neighbor->mt_->coord_;
  double x_dist = fabs(stalk_coord - site_coord);
  int x_dub = 2 * x_dist;
  int n_neighbs = neighbor->GetPRC1NeighborCount();
  double weight = weights_bind_i_teth_[n_neighbs][x_dub];
  return weight;
}

double AssociatedProtein::GetBindingWeight_II_Teth(Tubulin *neighbor) {

  // Get number of PRC1 neighbors this site has
  int n_neighbs = neighbor->GetPRC1NeighborCount();
  // Get x_dist of xlink extension if 2nd head were to bind to this neighb
  Tubulin *site = GetActiveHeadSite();
  double site_coord = site->index_ + site->mt_->coord_;
  double neighb_coord = neighbor->index_ + neighbor->mt_->coord_;
  int x = fabs(site_coord - neighb_coord);
  // Get current x_dist_dub of tether
  double stalk_coord = motor_->GetStalkCoordinate();
  int x_dub = 2 * fabs(stalk_coord - site_coord);
  // Get NEW x_dist_dub if 2nd head were to bind to this neighb
  double anchor_coord_post = (neighb_coord + site_coord) / 2;
  int x_dub_post = 2 * fabs(stalk_coord - anchor_coord_post);
  int x_dub_rest = 2 * properties_->kinesin4.rest_dist_;
  double weight;
  // If extended, ...
  if (x_dub >= x_dub_rest) {
    // ... further extension is going away FROM rest
    if (x_dub_post > x_dub) {
      weight = weights_bind_ii_fr_teth_[n_neighbs][x_dub][x];
    } else {
      weight = weights_bind_ii_to_teth_[n_neighbs][x_dub][x];
    }
  }
  // Else, if compressed ...
  else {
    // ... further extension is going TO rest
    if (x_dub_post > x_dub) {
      weight = weights_bind_ii_to_teth_[n_neighbs][x_dub][x];
    } else {
      weight = weights_bind_ii_fr_teth_[n_neighbs][x_dub][x];
    }
  }
  return weight;
}

Tubulin *AssociatedProtein::GetSiteCloserToTethRest() {

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
  } else {
    printf("Error in get_site_closer_to_teth_rest XLINK\n");
    exit(1);
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

Tubulin *AssociatedProtein::GetWeightedSite_Bind_II() {

  UpdateNeighborSites_II();
  double anch_coord = GetAnchorCoordinate();
  double p_tot = 0;
  for (int i_site = 0; i_site < n_neighbor_sites_; i_site++) {
    Tubulin *site = neighbor_sites_[i_site];
    double site_coord = site->index_ + site->mt_->coord_;
    int x_dist = (int)fabs(anch_coord - site_coord);
    int n_neighbs = site->GetPRC1NeighborCount();
    //		printf("x_dist is %i\n", x_dist);
    p_tot += weights_bind_ii_[n_neighbs][x_dist];
  }
  double ran = properties_->gsl.GetRanProb();
  double p_cum = 0;
  for (int i_site = 0; i_site < n_neighbor_sites_; i_site++) {
    Tubulin *site = neighbor_sites_[i_site];
    double site_coord = site->index_ + site->mt_->coord_;
    int x_dist = (int)fabs(anch_coord - site_coord);
    int n_neighbs = site->GetPRC1NeighborCount();
    p_cum += weights_bind_ii_[n_neighbs][x_dist] / p_tot;
    if (ran < p_cum) {
      return site;
    }
  }
  return nullptr;
}

Tubulin *AssociatedProtein::GetWeightedSite_Bind_I_Teth() {

  UpdateNeighborSites_I_Teth();
  double stalk_coord = motor_->GetStalkCoordinate();
  double p_tot = 0;
  for (int i_site = 0; i_site < n_teth_neighbor_sites_; i_site++) {
    Tubulin *site = teth_neighbor_sites_[i_site];
    double site_coord = site->index_ + site->mt_->coord_;
    double x_dist = fabs(stalk_coord - site_coord);
    int x_dist_dub = 2 * x_dist;
    int n_neighbs = site->GetPRC1NeighborCount();
    p_tot += weights_bind_i_teth_[n_neighbs][x_dist_dub];
  }
  double ran = properties_->gsl.GetRanProb();
  double p_cum = 0;
  for (int i_site = 0; i_site < n_teth_neighbor_sites_; i_site++) {
    Tubulin *site = teth_neighbor_sites_[i_site];
    double site_coord = site->index_ + site->mt_->coord_;
    double x_dist = fabs(stalk_coord - site_coord);
    int x_dist_dub = 2 * x_dist;
    int n_neighbs = site->GetPRC1NeighborCount();
    p_cum += weights_bind_i_teth_[n_neighbs][x_dist_dub] / p_tot;
    if (ran < p_cum) {
      return site;
    }
  }
  return nullptr;
}

Tubulin *AssociatedProtein::GetWeightedSite_Bind_II_Teth() {

  UpdateNeighborSites_II_Teth();
  double stalk_coord = motor_->GetStalkCoordinate();
  double anchor_coord = GetAnchorCoordinate();
  // All neighbor sites have the same initial x_dub
  int x_dub = 2 * fabs(stalk_coord - anchor_coord);
  int x_dub_rest = 2 * properties_->kinesin4.rest_dist_;
  double weight_tot = 0;
  for (int i_site = 0; i_site < n_teth_neighbor_sites_ii_; i_site++) {
    Tubulin *site = teth_neighbor_sites_ii_[i_site];
    double site_coord = site->index_ + site->mt_->coord_;
    int x = (int)fabs(site_coord - anchor_coord);
    double anchor_coord_post = (site_coord + anchor_coord) / 2;
    int x_dub_post = 2 * fabs(stalk_coord - anchor_coord_post);
    int n_neighbs = site->GetPRC1NeighborCount();
    double weight(0);
    // If extended, ...
    if (x_dub >= x_dub_rest) {
      // ... further extension is going away FROM rest
      if (x_dub_post > x_dub) {
        weight = weights_bind_ii_fr_teth_[n_neighbs][x_dub][x];
      } else {
        weight = weights_bind_ii_to_teth_[n_neighbs][x_dub][x];
      }
    }
    // Else, if compressed ...
    else {
      // ... further extension is going TO rest
      if (x_dub_post > x_dub) {
        weight = weights_bind_ii_to_teth_[n_neighbs][x_dub][x];
      } else {
        weight = weights_bind_ii_fr_teth_[n_neighbs][x_dub][x];
      }
    }
    weight_tot += weight;
  }
  double ran = properties_->gsl.GetRanProb();
  double p_cum = 0;
  for (int i_site = 0; i_site < n_teth_neighbor_sites_ii_; i_site++) {
    Tubulin *site = teth_neighbor_sites_ii_[i_site];
    double site_coord = site->index_ + site->mt_->coord_;
    int x = (int)fabs(site_coord - anchor_coord);
    double anchor_coord_post = (site_coord + anchor_coord) / 2;
    int x_dub_post = 2 * fabs(stalk_coord - anchor_coord_post);
    int n_neighbs = site->GetPRC1NeighborCount();
    double weight(0);
    // If extended, ...
    if (x_dub >= x_dub_rest) {
      // ... further extension is going away FROM rest
      if (x_dub_post > x_dub) {
        weight = weights_bind_ii_fr_teth_[n_neighbs][x_dub][x];
      } else {
        weight = weights_bind_ii_to_teth_[n_neighbs][x_dub][x];
      }
    }
    // Else, if compressed ...
    else {
      // ... further extension is going TO rest
      if (x_dub_post > x_dub) {
        weight = weights_bind_ii_to_teth_[n_neighbs][x_dub][x];
      } else {
        weight = weights_bind_ii_fr_teth_[n_neighbs][x_dub][x];
      }
    }
    p_cum += weight / weight_tot;
    if (ran < p_cum) {
      return site;
    }
  }
  return nullptr;
}
