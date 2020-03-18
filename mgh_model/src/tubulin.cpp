#include "tubulin.h"
#include "master_header.h"

Tubulin::Tubulin() {}

void Tubulin::Initialize(system_parameters *parameters,
                         system_properties *properties, Microtubule *mt,
                         int i_site) {

  parameters_ = parameters;
  properties_ = properties;
  index_ = i_site;
  mt_ = mt;
}

bool Tubulin::EquilibriumInSameDirection() {

  if (occupied_ and xlink_head_ != nullptr) {
    if (xlink_head_->xlink_->tethered_ and
        xlink_head_->xlink_->heads_active_ == 2) {
      Kinesin *motor = xlink_head_->xlink_->motor_;
      double site_coord = index_ + mt_->coord_;
      double xlink_rest = xlink_head_->xlink_->GetAnchorCoordinate();
      double motor_rest = motor->GetRestLengthCoordinate();
      // When tether is extended past rest length, ...
      if (motor->x_dist_doubled_ > 2 * motor->rest_dist_) {
        // ... equils are in same dir. if coords are on same side
        if ((xlink_rest > site_coord && motor_rest > site_coord) ||
            (xlink_rest < site_coord && motor_rest < site_coord)) {
          return true;
        }
        // ... equils are in oppo dir. if coords are on oppo side
        else if ((xlink_rest > site_coord && motor_rest < site_coord) ||
                 (xlink_rest < site_coord && motor_rest > site_coord)) {
          return false;
        } else {
          printf("error in Tubulin::SpringEquilOnSameSide\n");
          exit(1);
        }
      }
      // When tether at or compressed below rest length, ...
      else {
        // ... equils are in same dir. if coords are on oppo side
        if ((xlink_rest > site_coord && motor_rest < site_coord) ||
            (xlink_rest < site_coord && motor_rest > site_coord)) {
          return true;
        }
        // ... equils are in oppo dir. if coords are on same side
        else if ((xlink_rest > site_coord && motor_rest > site_coord) ||
                 (xlink_rest < site_coord && motor_rest < site_coord)) {
          return false;
        } else {
          printf("error in Tubulin::SpringEquilOnSameSide TWO\n");
          exit(1);
        }
      }
    } else {
      printf("error NO. 2 in tubulin: spring equil on same side ??\n");
      exit(1);
    }
  } else {
    printf("error in tubulin:spring equil on same side\n");
    exit(1);
  }
}

/*
void Tubulin::UpdateAffinity() {

  int n_affs{properties_->kinesin4.n_affinities_tot_};
  if (n_affs <= 1) {
    affinity_ = 0;
    return;
  }
  int n_affs_solo{properties_->kinesin4.n_affinities_};
  int range{(int)parameters_->motors.lattice_coop_range};
  int bin_size{range / n_affs_solo};
  if (motor_head_ != nullptr) {
    int n_stacks{motor_head_->GetCometSize() - 1};
    int i_aff{n_affs_solo};
    affinity_ = n_stacks * n_affs_solo + i_aff;
  }
  for (int delta{1}; delta <= range; delta++) {
    int i_fwd{index_ + delta};
    if (i_fwd > 0 and i_fwd < mt_->n_sites_) {
      if (mt_->lattice_[i_fwd].motor_head_ != nullptr) {
        int n_stacks{mt_->lattice_[i_fwd].motor_head_->GetCometSize() - 1};
        int i_aff{n_affs_solo - ((delta - 1) / bin_size)};
        affinity_ = n_stacks * n_affs_solo + i_aff;
        return;
      }
    }
    int i_bck{index_ - delta};
    if (i_bck > 0 and i_bck < mt_->n_sites_) {
      if (mt_->lattice_[i_bck].motor_head_ != nullptr) {
        int n_stacks{mt_->lattice_[i_bck].motor_head_->GetCometSize() - 1};
        int i_aff{n_affs_solo - ((delta - 1) / bin_size)};
        affinity_ = n_stacks * n_affs_solo + i_aff;
        return;
      }
    }
  }
  affinity_ = 0;
}

int Tubulin::GetAffinityExcluding(Kinesin *motor) {

  int n_affs{properties_->kinesin4.n_affinities_tot_};
  if (n_affs <= 1) {
    return 0;
  }
  int n_affs_solo{properties_->kinesin4.n_affinities_};
  int range{(int)parameters_->motors.lattice_coop_range};
  int bin_size{range / n_affs_solo};
  if (motor_head_ != nullptr) {
    if (motor_head_->motor_ != motor) {
      int n_stacks{motor_head_->GetCometSize() - 1};
      int i_aff{n_affs_solo};
      return n_stacks * n_affs_solo + i_aff;
    }
  }
  for (int delta{1}; delta <= range; delta++) {
    int i_fwd{index_ + delta};
    if (i_fwd > 0 and i_fwd < mt_->n_sites_) {
      if (mt_->lattice_[i_fwd].motor_head_ != nullptr) {
        if (mt_->lattice_[i_fwd].motor_head_->motor_ != motor) {
          int n_stacks{mt_->lattice_[i_fwd].motor_head_->GetCometSize() - 1};
          int i_aff{n_affs_solo - ((delta - 1) / bin_size)};
          return n_stacks * n_affs_solo + i_aff;
        }
      }
    }
    int i_bck{index_ - delta};
    if (i_bck > 0 and i_bck < mt_->n_sites_) {
      if (mt_->lattice_[i_bck].motor_head_ != nullptr) {
        if (mt_->lattice_[i_bck].motor_head_->motor_ != motor) {
          int n_stacks{mt_->lattice_[i_bck].motor_head_->GetCometSize() - 1};
          int i_aff{n_affs_solo - ((delta - 1) / bin_size)};
          return n_stacks * n_affs_solo + i_aff;
        }
      }
    }
  }
  return 0;
}
*/

int Tubulin::GetPRC1NeighborCount() {

  if (properties_->prc1.max_neighbs_ == 0) {
    return 0;
  }
  int n_neighbs{0};
  for (int delta{-1}; delta <= 1; delta += 2) {
    int i_scan = index_ + delta;
    if (i_scan < 0 or i_scan > (mt_->n_sites_ - 1)) {
      continue;
    }
    if (mt_->lattice_[i_scan].xlink_head_ != nullptr) {
      n_neighbs++;
    }
  }
  return n_neighbs;
}

int Tubulin::GetKif4ANeighborCount() {

  if (properties_->kinesin4.max_neighbs_ == 0) {
    return 0;
  }
  int n_neighbs{0};
  int mt_end{mt_->n_sites_ - 1};
  for (int delta{-1}; delta <= 1; delta += 2) {
    int i_scan = index_ + delta;
    if (i_scan < 0 or i_scan > mt_end) {
      continue;
    } else if (mt_->lattice_[i_scan].motor_head_ != nullptr) {
      n_neighbs++;
    }
  }
  return n_neighbs;
}

void Tubulin::UpdateWeights_Kinesin() {

  // Get weight from neighbor interactions
  int n_neighbs{GetKif4ANeighborCount()};
  weight_bind_ = properties_->kinesin4.weight_neighbs_bind_[n_neighbs];
  weight_unbind_ = properties_->kinesin4.weight_neighbs_unbind_[n_neighbs];
  // Get weight from lattice deformation
  int cutoff{properties_->kinesin4.lattice_cutoff_};
  double wt_max{weight_bind_ * properties_->kinesin4.weight_lattice_bind_max_};
  for (int delta{0}; delta <= cutoff; delta++) {
    int i_fwd{index_ + delta};
    if (i_fwd >= 0 and i_fwd < mt_->n_sites_) {
      Kinesin::Monomer *head_fwd{mt_->lattice_[i_fwd].motor_head_};
      if (head_fwd != nullptr) {
        if (head_fwd->motor_->heads_active_ == 1) {
          weight_bind_ *= properties_->kinesin4.weight_lattice_bind_[delta];
          weight_unbind_ *= properties_->kinesin4.weight_lattice_unbind_[delta];
        }
        // Forward head is most-recently bound head & thus has a larger effect
        else if (!head_fwd->trailing_) {
          weight_bind_ *= properties_->kinesin4.weight_lattice_bind_[delta];
          weight_unbind_ *= properties_->kinesin4.weight_lattice_unbind_[delta];
        }
        if (weight_bind_ >= wt_max) {
          return;
        }
      }
    }
    if (delta == 0) {
      continue;
    }
    int i_bck{index_ - delta};
    if (i_bck >= 0 and i_bck < mt_->n_sites_) {
      Kinesin::Monomer *head_bck{mt_->lattice_[i_bck].motor_head_};
      if (head_bck != nullptr) {
        if (head_bck->motor_->heads_active_ == 1) {
          weight_bind_ *= properties_->kinesin4.weight_lattice_bind_[delta];
          weight_unbind_ *= properties_->kinesin4.weight_lattice_unbind_[delta];
        } else if (!head_bck->trailing_) {
          weight_bind_ *= properties_->kinesin4.weight_lattice_bind_[delta];
          weight_unbind_ *= properties_->kinesin4.weight_lattice_unbind_[delta];
        }
        if (weight_bind_ >= wt_max) {
          break;
        }
      }
    }
  }
}
