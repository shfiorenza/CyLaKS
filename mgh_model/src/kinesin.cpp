#include "kinesin.h"
#include "master_header.h"

int Kinesin::head::GetAffinity() {

  int n_affs{motor_->properties_->kinesin4.n_affinities_};
  int range{(int)motor_->parameters_->motors.lattice_coop_range};
  int bin_size{0};
  if (n_affs > 1) {
    bin_size = range / (n_affs - 1);
  } else {
    return 0;
  }

  for (int delta{1}; delta <= range; delta++) {
    int i_fwd{site_->index_ + delta};
    if (i_fwd > 0 and i_fwd < site_->mt_->n_sites_) {
      if (site_->mt_->lattice_[i_fwd].motor_head_ != nullptr) {
        if (site_->mt_->lattice_[i_fwd].motor_head_->motor_ != motor_) {
          return ((delta - 1) / bin_size);
        }
      }
    }
    int i_bck{site_->index_ - delta};
    if (i_bck > 0 and i_bck < site_->mt_->n_sites_) {
      if (site_->mt_->lattice_[i_bck].motor_head_ != nullptr) {
        if (site_->mt_->lattice_[i_bck].motor_head_->motor_ != motor_) {
          return ((delta - 1) / bin_size);
        }
      }
    }
  }
  return (n_affs - 1);
}

int Kinesin::head::GetKIF4ANeighbCount() {

  if (site_ == nullptr)
    return 0;
  else {
    int n_neighbs = site_->GetKIF4ANeighborCount();
    if (motor_->heads_active_ == 2)
      n_neighbs--;
    return n_neighbs;
  }
}

Kinesin::Kinesin() {}

void Kinesin::Initialize(system_parameters *parameters,
                         system_properties *properties, int id) {

  id_ = id;
  parameters_ = parameters;
  properties_ = properties;
  SetParameters();
  CalculateCutoffs();
  InitializeNeighborLists();
  InitializeWeightLookup();
  InitializeExtensionLookup();
  InitializeSteppingProbabilities();
}

void Kinesin::SetParameters() {

  r_0_ = parameters_->motors.r_0;
  k_spring_ = parameters_->motors.k_spring;
  k_slack_ = parameters_->motors.k_slack;
}

void Kinesin::CalculateCutoffs() {

  // Only calculate cutoffs if tethering is enabled
  int site_size = parameters_->microtubules.site_size;
  double r_y = parameters_->microtubules.y_dist / 2;
  double kbT = parameters_->kbT;
  // First, calculate rest_dist_ in number of sites
  int roughly_rest = sqrt(r_0_ * r_0_ - r_y * r_y) / site_size;
  double rest_scan[3];
  double scan_force[3];
  for (int i_scan = -1; i_scan <= 1; i_scan++) {
    rest_scan[i_scan + 1] = roughly_rest + ((double)i_scan * 0.5);
    double rest_scan_length = rest_scan[i_scan + 1] * site_size;
    double r_scan = sqrt(r_y * r_y + rest_scan_length * rest_scan_length);
    if (r_scan >= r_0_)
      scan_force[i_scan + 1] = (r_scan - r_0_) * k_spring_;
    else
      scan_force[i_scan + 1] = (r_scan - r_0_) * k_slack_;
  }
  double min_force = 100;
  for (int i_scan = -1; i_scan <= 1; i_scan++) {
    double force = fabs(scan_force[i_scan + 1]);
    if (force < min_force) {
      min_force = force;
      rest_dist_ = rest_scan[i_scan + 1];
    }
  }
  // Next, calculate compression distance cutoff
  comp_cutoff_ = 0;
  for (int x_dub = 2 * rest_dist_; x_dub >= 0; x_dub--) {
    // increment by 0.5x
    int r_x = x_dub * site_size / 2;
    double r = sqrt(r_y * r_y + r_x * r_x);
    double dr = r - r_0_;
    double U;
    if (r < r_0_)
      U = (k_slack_ / 2) * dr * dr;
    else
      U = (k_spring_ / 2) * dr * dr;
    double boltzmann_weight = exp(U / (2 * kbT));
    if (boltzmann_weight > 100) {
      comp_cutoff_ = x_dub / 2;
      break;
    }
  }
  // Finally, calculate extension distance cutoff
  for (int x_dub = 2 * rest_dist_; x_dub < 1000; x_dub++) {
    // increment by 0.5x
    int r_x = x_dub * site_size / 2;
    double r = sqrt(r_y * r_y + r_x * r_x);
    double dr = r - r_0_;
    double U;
    if (r < r_0_)
      U = (k_slack_ / 2) * dr * dr;
    else
      U = (k_spring_ / 2) * dr * dr;
    double boltzmann_weight = exp(U / (2 * kbT));
    if (boltzmann_weight > 100) {
      teth_cutoff_ = x_dub / 2;
      break;
    }
  }
  if (!parameters_->motors.tethers_active) {
    rest_dist_ = 0;
    comp_cutoff_ = 0;
    teth_cutoff_ = 0;
  }
}

void Kinesin::InitializeNeighborLists() {

  int n_mts = parameters_->microtubules.count;
  // Serialize this bitch so we just roll one random number
  neighbor_xlinks_.resize(n_mts * (2 * teth_cutoff_ + 1));
  neighbor_sites_.resize(n_mts * (2 * teth_cutoff_ + 1));
}

void Kinesin::InitializeWeightLookup() {

  double r_y = parameters_->microtubules.y_dist / 2;
  double kbT = parameters_->kbT;
  double site_size = parameters_->microtubules.site_size;
  weight_lookup_.resize(2 * teth_cutoff_ + 1);
  weight_alt_lookup_.resize(2 * teth_cutoff_ + 1);
  for (int x_dub = 0; x_dub <= 2 * teth_cutoff_; x_dub++) {
    double r_x = (double)x_dub * site_size / 2;
    double r = sqrt(r_x * r_x + r_y * r_y);
    double cosine = r_x / r;
    double dr = r - r_0_;
    double k = 0;
    if (dr < 0)
      k = k_slack_;
    else
      k = k_spring_;
    double U_teth = k * dr * dr / 2;
    double weight = exp(-U_teth / (2 * kbT));
    double f_x = fabs(dr * k * cosine);
    double sigma_off = 1.5;
    double weight_alt = exp(-f_x * sigma_off / kbT);
    weight_lookup_[x_dub] = weight;
    weight_alt_lookup_[x_dub] = weight_alt;
  }
}

void Kinesin::InitializeExtensionLookup() {

  double r_y = parameters_->microtubules.y_dist / 2;
  double site_size = parameters_->microtubules.site_size;
  cosine_lookup_.resize(2 * teth_cutoff_ + 1);
  extension_lookup_.resize(2 * teth_cutoff_ + 1);
  for (int x_dub = 0; x_dub <= 2 * teth_cutoff_; x_dub++) {
    double r_x = (double)x_dub * site_size / 2;
    double r = sqrt(r_x * r_x + r_y * r_y);
    cosine_lookup_[x_dub] = r_x / r;
    extension_lookup_[x_dub] = r - r_0_;
  }
}

void Kinesin::InitializeSteppingProbabilities() {

  // Short-range cooperativity
  // ! Assume lambda = 1 & 0 for stepping from & away from neighbors
  double E_int = -1 * parameters_->motors.interaction_energy;
  // If we step AWAY from a single neighbor, overall E increases
  // Due to our assumption about lambda, same weight applies to
  // stepping from a neighbor TO a neighbor
  p_step_fr_neighb_ = exp(E_int);
  // Next, handle applied force
  double r_y = parameters_->microtubules.y_dist / 2;
  double site_size = parameters_->microtubules.site_size;
  double f_stall = parameters_->motors.stall_force;
  if (parameters_->motors.applied_force > 0) {
    double f_app = parameters_->motors.applied_force;
    /*
    double sigma_i = 1.84;
    double p_step = exp(-f_app*sigma_i/kbT);
    applied_p_step_ = p_step;
    */
    double dr = f_app / k_spring_;
    double r = r_0_ + dr;
    double r_x = sqrt(r * r - r_y * r_y);
    double cosine = r_x / r;
    double f_x = f_app * cosine;
    double p = sqrt(1 - (f_x / f_stall));
    applied_p_step_ = p;
    if (id_ == 0)
      printf("p_step scaled from 1 to %g\n", applied_p_step_);
  } else {
    applied_p_step_ = 1.0;
  }
  // Next, handle tethering forces
  p_step_to_rest_.resize(2 * teth_cutoff_ + 1);
  p_step_fr_rest_.resize(2 * teth_cutoff_ + 1);
  for (int x_dub = 0; x_dub <= 2 * teth_cutoff_; x_dub++) {
    double r_x = (double)x_dub * site_size / 2;
    double r = sqrt(r_x * r_x + r_y * r_y);
    double cosine = r_x / r;
    double dr = r - r_0_;
    double k = 0;
    if (dr < 0)
      k = k_slack_;
    else
      k = k_spring_;
    double force_mag = fabs(k * dr);
    double f_x = force_mag * cosine;
    double p = sqrt(1 - (f_x / f_stall));
    if (f_x >= f_stall)
      p = 0;
    if (x_dub < 2 * comp_cutoff_ - 1) {
      p_step_to_rest_[x_dub] = 0;
      p_step_fr_rest_[x_dub] = 0;
    } else if (x_dub <= 2 * comp_cutoff_) {
      p_step_to_rest_[x_dub] = p;
      p_step_fr_rest_[x_dub] = 0;
    } else if (x_dub < 2 * teth_cutoff_ - 1) {
      p_step_to_rest_[x_dub] = p;
      p_step_fr_rest_[x_dub] = p;
    } else {
      p_step_to_rest_[x_dub] = p;
      p_step_fr_rest_[x_dub] = 0;
    }
    //		printf("p_to=%g, p_fr=%g for 2x=%i\n", p_step_to_rest_[x_dub],
    //				p_step_fr_rest_[x_dub], x_dub);
  }
}

Kinesin::head *Kinesin::GetActiveHead() {

  if (heads_active_ == 1) {
    if (head_one_.site_ != nullptr)
      return &head_one_;
    else
      return &head_two_;
  } else {
    printf("error: cannot get active head, my dear brobro\n");
    exit(1);
  }
}

Kinesin::head *Kinesin::GetDockedHead() {

  if (heads_active_ == 1) {
    Kinesin::head *active_head = GetActiveHead();
    if (active_head->ligand_ == "ADPP") {
      Kinesin::head *docked_head(NULL);
      if (active_head == &head_one_)
        docked_head = &head_two_;
      else
        docked_head = &head_one_;
      if (docked_head->ligand_ == "ADP") {
        return docked_head;
      } else {
        printf("error! docked head doesn't have ADP bound??\n");
        exit(1);
      }
    } else {
      printf("error SIR; can't be docked w/o fuggin ADPP\n");
      exit(1);
    }
  } else {
    printf("error; cannot get docked head, madameeeeee\n");
    exit(1);
  }
}

Kinesin::head *Kinesin::StoreDockSite() {

  if (heads_active_ == 1) {
    int i_dock = GetDockedCoordinate() - mt_->coord_;
    Tubulin *dock = &mt_->lattice_[i_dock];
    if (dock->occupied_) {
      return nullptr;
    } else {
      head *docked_head = GetDockedHead();
      docked_head->stored_dock_site_ = dock;
      return docked_head;
    }
  } else {
    printf("error; cannot store dock site\n");
    exit(1);
  }
}

Tubulin *Kinesin::GetDockSite() {

  if (heads_active_ != 1) {
    printf("err 1 in Kinesin::GetDockSite\n");
    printf("*** EXITING ***\n");
    exit(1);
  }
  Kinesin::head *active_head = GetActiveHead();
  if (active_head->ligand_ != "ADPP") {
    printf("err 2 in Kinesin::GetDockSite\n");
    printf("*** EXITING ***\n");
    exit(1);
  }
  int i_dock{-1};
  if (active_head->trailing_) {
    i_dock = active_head->site_->index_ + mt_->delta_x_;
  } else {
    i_dock = active_head->site_->index_ - mt_->delta_x_;
  }
  if (i_dock >= 0 and i_dock < mt_->n_sites_) {
    Tubulin *dock_site = &mt_->lattice_[i_dock];
    if (dock_site->occupied_) {
      return nullptr;
    } else {
      return dock_site;
    }
  } else {
    return nullptr;
  }
}

double Kinesin::GetStalkCoordinate() {

  if (heads_active_ == 1) {
    Tubulin *site = GetActiveHead()->site_;
    double site_coord = site->index_ + mt_->coord_;
    if (GetActiveHead()->trailing_) {
      return site_coord + ((double)mt_->delta_x_) / 2;
    } else {
      return site_coord - ((double)mt_->delta_x_) / 2;
    }
  } else if (heads_active_ == 2) {
    int i_one = head_one_.site_->index_;
    int i_two = head_two_.site_->index_;
    double avg_index = ((double)(i_one + i_two)) / 2;
    double stalk_coord = mt_->coord_ + avg_index;
    return stalk_coord;
  } else {
    printf("error: can NOT get this stalk index bro.\n");
    exit(1);
  }
}

double Kinesin::GetDockedCoordinate() {

  if (heads_active_ == 1) {
    Kinesin::head *active_head = GetActiveHead();
    // finish checking for ADPP on active head here FIXME
    Tubulin *site = active_head->site_;
    double site_coord = site->index_ + mt_->coord_;
    if (active_head->trailing_) {
      return site_coord + mt_->delta_x_;
    } else {
      return site_coord - mt_->delta_x_;
    }
  } else {
    printf("Error in GetDockedCoordinate.\n");
    exit(1);
  }
}

void Kinesin::ChangeConformation() {

  if (heads_active_ == 0) {
    printf("'bruh' moment in kinesin:changeConformation\n");
    exit(1);
  } else if (heads_active_ == 2) {
    // ? Which of these is more physical?
    // Frustrated makes detailed balance easier, but ...
    //// UnbindTrailingHead();
    frustrated_ = true;
    return;
  }
  frustrated_ = false;
  int i_site = GetActiveHead()->site_->index_;
  if (parameters_->motors.endpausing_active and i_site == mt_->plus_end_) {
    return;
  }
  double p_step{1.0};
  p_step *= applied_p_step_;
  if (tethered_) {
    // Motors tethered to satellite xlinks behave as if untethered
    if (xlink_->heads_active_ > 0) {
      int dx_rest = GetDirectionTowardRest();
      if (mt_->delta_x_ == dx_rest) {
        p_step *= p_step_to_rest_[x_dist_doubled_];
      } else {
        p_step *= p_step_fr_rest_[x_dist_doubled_];
      }
    }
  }
  if (i_site != mt_->minus_end_) {
    if (mt_->lattice_[i_site - mt_->delta_x_].motor_head_ != nullptr) {
      p_step *= p_step_fr_neighb_;
    }
  }
  double ran{0.0};
  if (p_step < 1.0) {
    ran = properties_->gsl.GetRanProb();
  }
  if (ran < p_step) {
    head_one_.trailing_ = !head_one_.trailing_;
    head_two_.trailing_ = !head_two_.trailing_;
  }
}

void Kinesin::UnbindTrailingHead() {

  head *trailing_head{nullptr};
  if (head_one_.trailing_) {
    trailing_head = &head_one_;
  } else {
    trailing_head = &head_two_;
  }
  // Update site
  trailing_head->site_->occupied_ = false;
  trailing_head->site_->motor_head_ = nullptr;
  // Update motor
  trailing_head->site_ = nullptr;
  trailing_head->ligand_ = "ADP";
  trailing_head->state_ = "ripped off";
  trailing_head->GetOtherHead()->state_ = "solo stepper";
  heads_active_--;
}

bool Kinesin::IsStalled() {

  if (heads_active_ == 0)
    return false;
  else if (heads_active_ == 1) {
    Kinesin::head *active_head = GetActiveHead();
    int i_site = active_head->site_->index_;
    Microtubule *mt = active_head->site_->mt_;
    int i_plus = mt->plus_end_;
    int dx = mt->delta_x_;
    if (i_site == i_plus)
      return true;
    else
      return mt->lattice_[i_site + dx].occupied_;
  } else {
    Kinesin::head *front_head;
    if (head_one_.trailing_)
      front_head = &head_two_;
    else
      front_head = &head_one_;
    int i_front = front_head->site_->index_;
    Microtubule *mt = front_head->site_->mt_;
    int i_plus = mt->plus_end_;
    int dx = mt->delta_x_;
    if (i_front == i_plus)
      return true;
    else
      return mt->lattice_[i_front + dx].occupied_;
  }
}

bool Kinesin::AtCutoff() {

  int dx = mt_->delta_x_;
  int drest = GetDirectionTowardRest();

  if ((x_dist_doubled_ >= (2 * teth_cutoff_ - 1) && dx == -drest) ||
      (x_dist_doubled_ <= (2 * comp_cutoff_ + 1) && dx == -drest))
    return true;
  else
    return false;
}

void Kinesin::UpdateNeighborSites() {

  n_neighbor_sites_ = 0;
  if (tethered_ == true && heads_active_ == 0) {
    int n_mts = parameters_->microtubules.count;
    double anchor_coord = xlink_->GetAnchorCoordinate();
    // Scan through all potential neighbor sites; add unoccupied to list
    // FIXME this only works for two MTs as of now FIXME
    for (int i_mt = 0; i_mt < n_mts; i_mt++) {
      Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
      int mt_length = mt->n_sites_;
      double mt_coord = mt->coord_;
      int i_anchor = anchor_coord - mt_coord;
      for (int dx = -(teth_cutoff_ + 1); dx <= (teth_cutoff_ + 1); dx++) {
        int i_site = i_anchor + dx;
        // Start index at first site (0) if site index is <= 0
        if (i_site < 0) {
          dx -= (i_site + 1);
        }
        // End scan at last bulk site (mt_length - 1)
        else if (i_site > mt_length - 1) {
          break;
        } else {
          Tubulin *neighbor = &mt->lattice_[i_site];
          double site_coord = i_site + neighbor->mt_->coord_;
          double x = fabs(anchor_coord - site_coord);
          int x_dub = 2 * x;
          if (x_dub >= 2 * comp_cutoff_ && x_dub <= 2 * teth_cutoff_ &&
              neighbor->occupied_ == false) {
            neighbor_sites_[n_neighbor_sites_] = neighbor;
            n_neighbor_sites_++;
          }
        }
      }
    }
  } else {
    printf("error in update neighbor sites\n");
    exit(1);
  }
}

void Kinesin::UpdateNeighborXlinks() {

  n_neighbor_xlinks_ = 0;
  if (tethered_ == false && heads_active_ > 0) {
    int n_mts = parameters_->microtubules.count;
    double stalk_coord = GetStalkCoordinate();
    // Scan through all potential neighbor sites; add untethered xlinks
    // FIXME this only works for two MTs as of now FIXME
    for (int i_mt = 0; i_mt < n_mts; i_mt++) {
      Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
      int mt_length = mt->n_sites_;
      int i_stalk = stalk_coord - mt_->coord_;
      // Only scan over sites within +/- teth_cutoff_
      for (int dx = -(teth_cutoff_ + 1); dx <= (teth_cutoff_ + 1); dx++) {
        int i_site = i_stalk + dx;
        // Start index at first site (0) if site index is <= 0
        if (i_site < 0) {
          dx -= (i_site + 1);
        }
        // End scan once last site (mt_length - 1) is reached
        else if (i_site > mt_length - 1) {
          break;
        } else {
          Tubulin *site = &mt->lattice_[i_site];
          if (site->xlink_head_ != nullptr) {
            AssociatedProtein *xlink = site->xlink_head_->xlink_;
            double anchor_coord = xlink->GetAnchorCoordinate();
            double x = fabs(anchor_coord - stalk_coord);
            int x_dub = 2 * x;
            if (x_dub >= 2 * comp_cutoff_ && x_dub <= 2 * teth_cutoff_ &&
                xlink->tethered_ == false) {
              if (xlink->heads_active_ == 1) {
                neighbor_xlinks_[n_neighbor_xlinks_] = xlink;
                n_neighbor_xlinks_++;
              }
              // FIXME this is bad
              else if (xlink->heads_active_ == 2 && i_mt % 2 == 0) {
                neighbor_xlinks_[n_neighbor_xlinks_] = xlink;
                n_neighbor_xlinks_++;
              }
            }
          }
        }
      }
    }
  } else {
    printf("error in update neighbor xlinks\n");
    exit(1);
  }
}

void Kinesin::UpdateExtension() {

  if (heads_active_ == 0 || tethered_ == false) {
    x_dist_doubled_ = 0;
    extension_ = 0;
  } else if (xlink_->heads_active_ == 0) {
    x_dist_doubled_ = 0;
    extension_ = 0;
  } else {
    double stalk_coord = GetStalkCoordinate();
    double anchor_coord = xlink_->GetAnchorCoordinate();
    double x_dist = fabs(anchor_coord - stalk_coord);
    x_dist_doubled_ = 2 * x_dist;
    if (x_dist_doubled_ > 2 * teth_cutoff_ ||
        x_dist_doubled_ < 2 * comp_cutoff_) {
      ForceUntether();
    } else {
      extension_ = extension_lookup_[x_dist_doubled_];
      cosine_ = cosine_lookup_[x_dist_doubled_];
    }
  }
}

void Kinesin::ForceUntether() {

  // Update motor
  tethered_ = false;
  x_dist_doubled_ = 0;
  extension_ = 0;
  // Update xlink
  xlink_->motor_ = nullptr;
  xlink_->tethered_ = false;
  xlink_ = nullptr;
}

void Kinesin::UntetherSatellite() {

  if (!tethered_) {
    return;
  }
  if (xlink_->heads_active_ == 0) {
    // Remove satellite xlink from active_, replace with last entry
    int i_last = properties_->prc1.n_active_ - 1;
    AssociatedProtein *last_entry = properties_->prc1.active_[i_last];
    int i_this = xlink_->active_index_;
    if (i_this != i_last) {
      properties_->prc1.active_[i_this] = last_entry;
      last_entry->active_index_ = i_this;
    }
    properties_->prc1.n_active_--;
    // Update xlink details
    xlink_->tethered_ = false;
    xlink_->motor_ = nullptr;
    // Update motor details
    tethered_ = false;
    xlink_ = nullptr;
  }
}

int Kinesin::GetDirectionTowardRest() {

  if (tethered_ == true) {
    if (xlink_->heads_active_ > 0) {
      double stalk_coord = GetStalkCoordinate();
      double rest_coord = GetRestLengthCoordinate();
      if (stalk_coord > rest_coord)
        return -1;
      else if (stalk_coord < rest_coord)
        return 1;
      // If stalk is at rest dist, stepping TO xlink is closer
      // to true rest than stepping FROM xlink
      else {
        double anchor_coord = xlink_->GetAnchorCoordinate();
        if (stalk_coord > anchor_coord)
          return -1;
        else
          return 1;
      }
    } else
      return 0;
  } else {
    printf("error in get dir. toward xlink\n");
    exit(1);
  }
}

double Kinesin::GetRestLengthCoordinate() {

  if (tethered_ == true && heads_active_ > 0) {
    if (xlink_->heads_active_ > 0) {
      double stalk_coord = GetStalkCoordinate();
      double anch_coord = xlink_->GetAnchorCoordinate();
      if (stalk_coord > anch_coord) {
        double rest_coord = anch_coord + rest_dist_;
        return rest_coord;
      } else if (stalk_coord < anch_coord) {
        double rest_coord = anch_coord - rest_dist_;
        return rest_coord;
      } else {
        printf("lameness in get rest length coord kinesin\n");
        printf("stalk: %g | anchor: %g\n", stalk_coord, anch_coord);
        exit(1);
      }
    } else {
      printf("error in get_rest_length_coord TWO (kinesin)\n");
      exit(1);
    }
  } else {
    printf("error in get_rest_length_coord (kinesin)\n");
    exit(1);
  }
}

double Kinesin::GetTetherForce(Tubulin *site) {

  if (tethered_ == true && heads_active_ > 0) {
    // Ensure we're not tethered to a satellite xlink
    if (xlink_->heads_active_ > 0) {
      UpdateExtension();
      // Make sure we didn't force an untether event
      if (tethered_ == true) {
        double force_mag;
        if (extension_ < 0)
          force_mag = extension_ * k_slack_;
        else
          force_mag = extension_ * k_spring_;
        double stalk_coord = GetStalkCoordinate();
        double anchor_coord = xlink_->GetAnchorCoordinate();
        double center = (stalk_coord + anchor_coord) / 2;
        double site_coord = site->index_ + site->mt_->coord_;
        double force;
        // Account for Newton's 3rd
        if (site_coord < center)
          force = force_mag * cosine_;
        else
          force = -1 * force_mag * cosine_;
        return force;
      } else
        return 0;
    } else {
      printf("error in get teth force TWO (motor)\n");
      exit(1);
    }
  } else {
    printf("error in get teth force (motor)\n");
    exit(1);
  }
}

double Kinesin::GetTotalBindingWeight() {

  UpdateNeighborSites();
  double tot_wt = 0;
  for (int i_neighb = 0; i_neighb < n_neighbor_sites_; i_neighb++) {
    tot_wt += GetBindingWeight(neighbor_sites_[i_neighb]);
  }
  return tot_wt;
}

double Kinesin::GetBindingWeight(Tubulin *site) {

  // Since only one foot of the motor attaches initially,
  // we only need to consider the coordinate of one site
  // (assuming the diffusion of the other foot avgs out)
  int i_site = site->index_;
  int mt_coord = site->mt_->coord_;
  double site_coord = (double)(mt_coord + i_site);
  double anchor_coord = xlink_->GetAnchorCoordinate();
  double x_dist = fabs(anchor_coord - site_coord);
  // Multiply by 2 to guarentee an integer
  int x_dub = x_dist * 2;
  // Get binding weight that corresponds to this adj. x-dist
  return weight_alt_lookup_[x_dub];
}

double Kinesin::GetTotalTetheringWeight() {

  UpdateNeighborXlinks();
  double tot_wt = 0;
  for (int i_neighb = 0; i_neighb < n_neighbor_xlinks_; i_neighb++) {
    tot_wt += GetTetheringWeight(neighbor_xlinks_[i_neighb]);
  }
  return tot_wt;
}

double Kinesin::GetTetheringWeight(AssociatedProtein *xlink) {

  double stalk_coord = GetStalkCoordinate();
  double anchor_coord = xlink->GetAnchorCoordinate();
  double x_dist = fabs(anchor_coord - stalk_coord);
  // Multiply by 2 to guarentee an integer
  int x_dub = x_dist * 2;
  // Get tethering weight that corressponds to this adj. x-dist
  return weight_lookup_[x_dub];
}

Tubulin *Kinesin::GetWeightedNeighborSite() {

  if (n_neighbor_sites_ == 0) {
    printf("BRUH - no neighbor sites to choose from (MOTOR)\n");
    exit(1);
  }
  // Get total binding probability of all eligible neighbor sites
  double p_tot = 0;
  for (int i_site = 0; i_site < n_neighbor_sites_; i_site++) {
    p_tot += GetBindingWeight(neighbor_sites_[i_site]);
  }
  // Scan through neighbor sites; pick one based on normalized weights
  double ran = properties_->gsl.GetRanProb();
  double p_cum = 0;
  for (int i_site = 0; i_site < n_neighbor_sites_; i_site++) {
    p_cum += GetBindingWeight(neighbor_sites_[i_site]) / p_tot;
    if (ran <= p_cum)
      return neighbor_sites_[i_site];
  }
}

AssociatedProtein *Kinesin::GetWeightedNeighborXlink() {

  if (n_neighbor_xlinks_ == 0) {
    printf("BRUH - no neighbor xlinks to choose from (MOTOR)\n");
    exit(1);
  }
  double p_tot = 0;
  // Get total tethering probability of all eligible neighbor xlinks
  for (int i_xlink = 0; i_xlink < n_neighbor_xlinks_; i_xlink++) {
    p_tot += GetTetheringWeight(neighbor_xlinks_[i_xlink]);
  }
  double ran = properties_->gsl.GetRanProb();
  double p_cum = 0;
  // Scan through neighbor xlinks; pick one based on normalized weights
  for (int i_xlink = 0; i_xlink < n_neighbor_xlinks_; i_xlink++) {
    p_cum += GetTetheringWeight(neighbor_xlinks_[i_xlink]) / p_tot;
    if (ran <= p_cum)
      return neighbor_xlinks_[i_xlink];
  }
}