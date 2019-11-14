#include "associated_protein.h"
#include "master_header.h"

// fuckin come at me bro
void AssociatedProtein::Monomer::UpdateState() {

  std::string root, x, x_dub, neighbs;
  int n_neighbs;
  if (xlink_->heads_active_ == 0) {
    if (xlink_->tethered_) {
      state_ = std::string("free_teth");
    } else {
      state_ = std::string("unbound");
    }
  } else if (xlink_->heads_active_ == 1) {
    if (this == xlink_->GetActiveHead()) {
      neighbs = std::to_string(GetPRC1NeighbCount());
      if (xlink_->tethered_) {
        if (xlink_->motor_->heads_active_ > 0) {
          x_dub = std::to_string(xlink_->motor_->x_dist_doubled_);
          root = std::string("bound_I_") + x_dub;
        } else {
          root = std::string("bound_i");
        }
      } else {
        root = std::string("bound_i");
      }
      state_ = root + std::string("_") + neighbs;
    } else {
      state_ = std::string("docked");
    }
  } else if (xlink_->heads_active_ == 2) {
    neighbs = std::to_string(GetPRC1NeighbCount());
    x = std::to_string(xlink_->x_dist_);
    if (xlink_->tethered_) {
      if (xlink_->motor_->heads_active_ > 0) {
        x_dub = std::to_string(xlink_->motor_->x_dist_doubled_);
        // For x = rest dist, store as "same" regardless
        // Note: this accounted for in CheckScratchFor()
        if (xlink_->x_dist_ == xlink_->rest_dist_) {
          root = "bound_II_same_" + x + "_" + x_dub;
        } else if (site_->EquilibriumInSameDirection()) {
          root = "bound_II_same_" + x + "_" + x_dub;
        } else {
          root = "bound_II_oppo_" + x + "_" + x_dub;
        }
      } else {
        root = std::string("bound_ii_") + x;
      }
    } else {
      root = std::string("bound_ii_") + x;
    }
    state_ = root + std::string("_") + neighbs;
  } else {
    xlink_->wally_->ErrorExit("AssociatedProtein::Monomer::UpdateState()");
  }
}

AssociatedProtein::AssociatedProtein() {}

void AssociatedProtein::Initialize(system_parameters *parameters,
                                   system_properties *properties, int ID) {

  ID_ = ID;
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
  max_neighbs_ = properties_->prc1.max_neighbs_;
}

void AssociatedProtein::InitializeLookupTables() {

  /*
  // Lambda = 1.0 means all of the weight goes into unbinding
  double lambda_neighb{1.0}; // Lambda for neighbor interaction energies
  // Lambda = 0.0 means all of the weight goes into binding
  double lambda_spring{0.0}; // Lambda for spring energies
  // Lambda = 0.5 means the weight is equally split between binding & unbinding
  double lambda_teth{0.5}; // lambda for tether spring energies
  // Calculate neighbor interaction energies
  max_neighbs_ = 2;
  double interaction_energy{-1 * parameters_->xlinks.interaction_energy};
  std::vector<double> neighb_energy(max_neighbs_ + 1, 0.0);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    // Neighbor energies are negative since its an attractive potential
    neighb_energy[n_neighbs] = n_neighbs * interaction_energy; //! in kBT
  }
  if (lambda_neighb != 1.0) {
    printf("Lambda != 1.0 for neighb interactions not implemented yet!\n");
    wally_->ErrorExit("AssociatedProtein::InitializeLookupTables()\n");
  }
  // Calculate spring extension energies
  double k_spring{parameters_->xlinks.k_spring};
  std::vector<double> spring_energy(dist_cutoff_ + 1, 0.0);
  for (int x{0}; x <= dist_cutoff_; x++) {
    double r_x{x * site_size};
    double r{sqrt(r_x * r_x + r_y * r_y)};
    double dr{r - r_0};
    // Spring energies are positive by definition
    spring_energy[x] = 0.5 * k_spring * dr * dr; //! in pN*nm
  }
  */

  /*
  // FIXME: delete below once teth stuff is updated
  // Neighb jawns
  std::vector<double> neighb_weights(max_neighbs_ + 1, 0.0);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    double tot_energy = n_neighbs * interaction_energy;
    // Multiply by to increase binding rate; divide by to reduce unbinding rate
    neighb_weights[n_neighbs] = exp(-tot_energy / 2); // E has units kBT
  }
  */

  double site_size{parameters_->microtubules.site_size};
  double r_y{parameters_->microtubules.y_dist};
  double r_0{parameters_->xlinks.r_0};
  int teth_cutoff{properties_->kinesin4.dist_cutoff_};
  int comp_cutoff{properties_->kinesin4.comp_cutoff_};
  /* Construct rly basic lookup tables first */
  extension_lookup_.resize(dist_cutoff_ + 1);
  cosine_lookup_.resize(dist_cutoff_ + 1);
  for (int x_dist{0}; x_dist <= dist_cutoff_; x_dist++) {
    double r_x = site_size * x_dist;
    double r = sqrt(r_x * r_x + r_y * r_y);
    extension_lookup_[x_dist] = r - r_0_;
    cosine_lookup_[x_dist] = r_x / r;
  }
  weight_bind_ii_.resize(max_neighbs_ + 1);
  weight_bind_i_teth_.resize(max_neighbs_ + 1);
  weight_bind_ii_to_teth_.resize(max_neighbs_ + 1);
  weight_bind_ii_fr_teth_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
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
  /*
  // Next, construct lookup table for bind_ii weights
  double kbT{parameters_->kbT};
  weights_bind_ii_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    double E_i{0.0};
    double E_f{neighb_energy[n_neighbs]};
    double dE = E_f - E_i;
    // dE has units of kbT, so dividing by its numerical value isn't necessary
    double weight_neighb{exp(-(1.0 - lambda_neighb) * dE)};
    weights_bind_ii_[n_neighbs].resize(dist_cutoff_ + 1);
    for (int x{0}; x <= dist_cutoff_; x++) {
      double U_i{0.0};
      double U_f{spring_energy[x]};
      double dU = U_f - U_i;
      double weight_bind{exp(-(1.0 - lambda_spring) * dU / kbT)};
      // printf("weight_bind_[%i][%i] = %g\n", n_neighbs, x, weight_bind);
      weights_bind_ii_[n_neighbs][x] = weight_neighb * weight_bind;
    }
  }
  // Next, construct lookup table for bind_i_teth weights

  double lambda{1.0};
  int teth_cutoff = properties_->kinesin4.dist_cutoff_;
  int comp_cutoff = properties_->kinesin4.comp_cutoff_;
  double k_teth = parameters_->motors.k_spring;
  double k_slack = parameters_->motors.k_slack;
  double r_0_teth = parameters_->motors.r_0;
  double r_y_teth = parameters_->microtubules.y_dist / 2;
  weights_bind_i_teth_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    double weight_neighbs = neighb_weights[n_neighbor_sites_];
    weights_bind_i_teth_[n_neighbs].resize(2 * teth_cutoff + 1);
    for (int x_dub{2 * comp_cutoff}; x_dub <= 2 * teth_cutoff; x_dub++) {
      double r_x = (double)x_dub * site_size / 2;
      double r = sqrt(r_y_teth * r_y_teth + r_x * r_x);
      double dr = r - r_0_teth;
      double U_teth{0};
      if (dr < 0) {
        U_teth = (k_slack / 2) * dr * dr;
      } else {
        U_teth = (k_teth / 2) * dr * dr;
      }
      double weight_teth = exp(-U_teth / (2 * kbT));
      double weight_tot = weight_neighbs * weight_teth;
      weights_bind_i_teth_[n_neighbs][x_dub] = weight_tot;
      //			printf("weight for 2x = %i is %g\n", x_dub,
      // weight_tot);
    }
  }
  // Next, construct lookup table for bind_ii_teth weights
  // XXX input: CURRENT x_dub and PROPOSED x_dist XXX
  double rest_dist_teth = properties_->kinesin4.motors_[0].rest_dist_;
  weights_bind_ii_to_teth_.resize(max_neighbs_ + 1);
  weights_bind_ii_fr_teth_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    double weight_neighbs = neighb_weights[n_neighbs];
    weights_bind_ii_to_teth_[n_neighbs].resize(2 * teth_cutoff + 1);
    weights_bind_ii_fr_teth_[n_neighbs].resize(2 * teth_cutoff + 1);
    for (int x_dub{2 * comp_cutoff}; x_dub <= 2 * teth_cutoff; x_dub++) {
      // Calc x-distances (in nm) for tether
      double r_x_teth = (double)x_dub * site_size / 2;
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
      for (int x{rest_dist_}; x <= dist_cutoff_; x++) {
        double r_x = (double)x * site_size;
        double r = sqrt(r_y * r_y + r_x * r_x);
        double dr = r - r_0_;
        double U = (k_spring_ / 2) * dr * dr;
        double weight = exp(-1 * lambda * U / kbT);
        // change in teth extension if 2nd xlink head were to bind
        double dx_teth = (double)x / 2;
        double dr_x_teth = dx_teth * site_size;
        double r_x_teth_to, r_x_teth_fr;
        if (dr_teth < 0) {
          r_x_teth_to = r_x_teth + dr_x_teth;
          r_x_teth_fr = r_x_teth - dr_x_teth;
        } else {
          r_x_teth_to = r_x_teth - dr_x_teth;
          r_x_teth_fr = r_x_teth + dr_x_teth;
        }
        double r_teth_to{sqrt(r_x_teth_to * r_x_teth_to + r_y_teth * r_y_teth)};
        double r_teth_fr{sqrt(r_x_teth_fr * r_x_teth_fr + r_y_teth * r_y_teth)};
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
        printf("for n=%i, 2x=%i, x=%i, w_to = %g\n",
                        n_neighbs, x_dub, x, weights_bind_ii_to_teth_
                        [n_neighbs][x_dub][x]);
        printf("for n=%i, 2x=%i, x=%i, w_from = %g\n",
                        n_neighbs, x_dub, x, weights_bind_ii_fr_teth_
                        [n_neighbs][x_dub][x]);
      }
    }
  }
  */
}

void AssociatedProtein::InitializeNeighborLists() {

  int n_mts{parameters_->microtubules.count};
  int teth_cutoff{properties_->kinesin4.dist_cutoff_};
  neighbor_sites_.resize((n_mts - 1) * (2 * dist_cutoff_ + 1));
  teth_neighbor_sites_.resize(n_mts * (2 * teth_cutoff + 1));
  teth_neighbor_sites_ii_.resize((n_mts - 1) * (2 * dist_cutoff_ + 1));
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

Tubulin *AssociatedProtein::GetActiveHeadSite() {

  if (heads_active_ == 1) {
    if (head_one_.site_ != nullptr) {
      return head_one_.site_;
    } else if (head_two_.site_ != nullptr) {
      return head_two_.site_;
    } else {
      wally_->ErrorExit("AssociatedProtein::GetActiveHeadSite [1]");
    }
  } else {
    wally_->ErrorExit("AssociatedProtein::GetActiveHeadSite [2]");
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
    wally_->ErrorExit("AssociatedProtein::GetAnchorCoordinate()");
  }
}

int AssociatedProtein::GetPRC1NeighbCount(Monomer *head) {

  int n_neighbs{0};
  Tubulin *site{head->site_};
  if (site != nullptr) {
    int i_site{site->index_};
    int i_last{site->mt_->n_sites_ - 1};
    for (int delta{-1}; delta <= 1; delta += 2) {
      int i_scan = i_site + delta;
      if (i_scan < 0 or i_scan > i_last) {
        continue;
      } else if (site->mt_->lattice_[i_scan].xlink_head_ != nullptr) {
        n_neighbs++;
      }
    }
  }
  return n_neighbs;
}

void AssociatedProtein::UpdateNeighborSites_II() {

  n_neighbor_sites_ = 0;
  Tubulin *site{GetActiveHeadSite()};
  int site_coord{site->mt_->coord_ + site->index_};
  Microtubule *neighb_mt{site->mt_->neighbor_};
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
      neighbor_sites_[n_neighbor_sites_++] = neighb;
    }
  }
}

void AssociatedProtein::UpdateNeighborSites_I_Teth() {

  if (!tethered_ or heads_active_ != 0) {
    wally_->ErrorExit("AssociatedProtein::UpdateNeighborSites_I_Teth()");
  }
  n_teth_neighbor_sites_ = 0;
  int teth_cutoff{properties_->kinesin4.dist_cutoff_};
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
          teth_neighbor_sites_[n_teth_neighbor_sites_++] = neighb;
        }
      }
    }
  }
}

void AssociatedProtein::UpdateNeighborSites_II_Teth() {

  if (!tethered_ or heads_active_ != 1) {
    wally_->ErrorExit("AssociatedProtein::UpdateNeighborSites_I_Teth()");
  }
  n_teth_neighbor_sites_ii_ = 0;
  int teth_cutoff{properties_->kinesin4.dist_cutoff_};
  int comp_cutoff{properties_->kinesin4.comp_cutoff_};
  double stalk_coord{motor_->GetStalkCoordinate()};
  Tubulin *site{GetActiveHeadSite()};
  int site_coord{site->mt_->coord_ + site->index_};
  Microtubule *neighb_mt{site->mt_->neighbor_};
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
        teth_neighbor_sites_ii_[n_teth_neighbor_sites_ii_++] = neighb;
      }
    }
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
      extension_ = extension_lookup_[x_dist];
      cosine_ = cosine_lookup_[x_dist];
    } else {
      ForceUnbind();
    }
  } else {
    x_dist_ = 0;
    extension_ = 0;
    cosine_ = 0;
  }
}

void AssociatedProtein::ForceUnbind() {

  // If xlink isn't tethered, flip a coin to choose which head to unbind
  if (!tethered_) {
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
    heads_active_--;
    x_dist_ = 0;
    extension_ = 0;
    cosine_ = 0;
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
    heads_active_--;
    x_dist_ = 0;
    extension_ = 0;
    cosine_ = 0;
    // Update motor extension since anchor coord has changed
    motor_->UpdateExtension();
  }
  // Xlinks tethered to satellite motors behave as if untethered
  else {
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
    heads_active_--;
    x_dist_ = 0;
    extension_ = 0;
    cosine_ = 0;
  }
}

void AssociatedProtein::UntetherSatellite() {

  if (tethered_) {
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
  } else {
    wally_->ErrorExit("AssociatedProtein::UntetherSatellite()");
  }
}

int AssociatedProtein::GetDirectionToRest(Tubulin *site) {
  if (heads_active_ == 0) {
    return 0;
  } else if (heads_active_ == 1) {
    return 1;
  } else if (heads_active_ == 2) {
    double anchor_coord{GetAnchorCoordinate()};
    int site_coord{site->index_ + site->mt_->coord_};
    if (fabs(anchor_coord - site_coord) < 0.0001) {
      double ran{properties_->gsl.GetRanProb()};
      if (ran < 0.5) {
        return -1;
      } else {
        return 1;
      }
    } else if (site_coord > anchor_coord) {
      return -1;
    } else {
      return 1;
    }
  } else {
    wally_->ErrorExit("AssociatedProtein::GetDirectionToRest()");
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

double AssociatedProtein::GetBindingWeight_II(Tubulin *neighbor) {

  Tubulin *site{GetActiveHeadSite()};
  int site_coord{site->mt_->coord_ + site->index_};
  int neighbor_coord{neighbor->mt_->coord_ + neighbor->index_};
  int x_dist{abs(site_coord - neighbor_coord)};
  int n_neighbs{neighbor->GetPRC1NeighborCount()};
  // printf("x=%i, n_neighbs=%i\n", x_dist, n_neighbs);
  // printf("weight is %g\n", weights_bind_ii_[n_neighbs][x_dist]);
  return weight_bind_ii_[n_neighbs][x_dist];
}

double AssociatedProtein::GetBindingWeight_I_Teth(Tubulin *neighbor) {

  double stalk_coord{motor_->GetStalkCoordinate()};
  int site_coord{neighbor->index_ + neighbor->mt_->coord_};
  int x_dub{(int)(2 * fabs(stalk_coord - site_coord))};
  int n_neighbs{neighbor->GetPRC1NeighborCount()};
  return weight_bind_i_teth_[n_neighbs][x_dub];
}

double AssociatedProtein::GetBindingWeight_II_Teth(Tubulin *neighbor) {

  // Get number of PRC1 neighbors this site has
  int n_neighbs{neighbor->GetPRC1NeighborCount()};
  // Get x_dist of xlink extension if 2nd head were to bind to this neighb
  Tubulin *site{GetActiveHeadSite()};
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

Tubulin *AssociatedProtein::GetWeightedSite_Bind_II() {

  UpdateNeighborSites_II();
  double anch_coord{GetAnchorCoordinate()};
  double p_tot{0.0};
  for (int i_site{0}; i_site < n_neighbor_sites_; i_site++) {
    Tubulin *site{neighbor_sites_[i_site]};
    int site_coord{site->index_ + site->mt_->coord_};
    int x_dist{(int)fabs(anch_coord - site_coord)};
    int n_neighbs{site->GetPRC1NeighborCount()};
    p_tot += weight_bind_ii_[n_neighbs][x_dist];
  }
  double ran{properties_->gsl.GetRanProb()};
  double p_cum{0.0};
  // printf("ran is %g\n", ran);
  for (int i_site{0}; i_site < n_neighbor_sites_; i_site++) {
    Tubulin *site{neighbor_sites_[i_site]};
    int site_coord{site->index_ + site->mt_->coord_};
    int x_dist{(int)fabs(anch_coord - site_coord)};
    int n_neighbs{site->GetPRC1NeighborCount()};
    // printf("p_cum is %g\n", p_cum);
    p_cum += weight_bind_ii_[n_neighbs][x_dist] / p_tot;
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
    p_tot += weight_bind_i_teth_[n_neighbs][x_dist_dub];
  }
  double ran = properties_->gsl.GetRanProb();
  double p_cum = 0;
  for (int i_site = 0; i_site < n_teth_neighbor_sites_; i_site++) {
    Tubulin *site = teth_neighbor_sites_[i_site];
    double site_coord = site->index_ + site->mt_->coord_;
    double x_dist = fabs(stalk_coord - site_coord);
    int x_dist_dub = 2 * x_dist;
    int n_neighbs = site->GetPRC1NeighborCount();
    p_cum += weight_bind_i_teth_[n_neighbs][x_dist_dub] / p_tot;
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
        weight = weight_bind_ii_fr_teth_[n_neighbs][x_dub][x];
      } else {
        weight = weight_bind_ii_to_teth_[n_neighbs][x_dub][x];
      }
    }
    // Else, if compressed ...
    else {
      // ... further extension is going TO rest
      if (x_dub_post > x_dub) {
        weight = weight_bind_ii_to_teth_[n_neighbs][x_dub][x];
      } else {
        weight = weight_bind_ii_fr_teth_[n_neighbs][x_dub][x];
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
        weight = weight_bind_ii_fr_teth_[n_neighbs][x_dub][x];
      } else {
        weight = weight_bind_ii_to_teth_[n_neighbs][x_dub][x];
      }
    }
    // Else, if compressed ...
    else {
      // ... further extension is going TO rest
      if (x_dub_post > x_dub) {
        weight = weight_bind_ii_to_teth_[n_neighbs][x_dub][x];
      } else {
        weight = weight_bind_ii_fr_teth_[n_neighbs][x_dub][x];
      }
    }
    p_cum += weight / weight_tot;
    if (ran < p_cum) {
      return site;
    }
  }
  return nullptr;
}