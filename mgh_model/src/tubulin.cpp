#include "tubulin.h"
#include "master_header.h"

Tubulin::Tubulin() {}

void Tubulin::Initialize(system_parameters *parameters,
                         system_properties *properties, Microtubule *mt,
                         int i_site) {

  parameters_ = parameters;
  properties_ = properties;
  index_ = i_site;
  coord_ = i_site;
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

void Tubulin::UpdateAffinity() {

  int n_affs{properties_->kinesin4.n_affinities_};
  int range{(int)parameters_->motors.lattice_coop_range};
  int bin_size{0};
  if (n_affs > 1) {
    bin_size = range / (n_affs - 1);
  } else {
    affinity_ = 0;
    return;
  }
  if (motor_head_ != nullptr) {
    affinity_ = 0;
    return;
  }

  for (int delta{1}; delta <= range; delta++) {
    int i_fwd{index_ + delta};
    if (i_fwd > 0 and i_fwd < mt_->n_sites_) {
      if (mt_->lattice_[i_fwd].motor_head_ != nullptr) {
        affinity_ = (delta - 1) / bin_size;
        return;
      }
    }
    int i_bck{index_ - delta};
    if (i_bck > 0 and i_bck < mt_->n_sites_) {
      if (mt_->lattice_[i_bck].motor_head_ != nullptr) {
        affinity_ = (delta - 1) / bin_size;
        return;
      }
    }
  }
  affinity_ = n_affs - 1;
}

int Tubulin::GetPRC1NeighborCount() {

  int n_neighbs{0};
  int i_plus = mt_->plus_end_;
  int i_minus = mt_->minus_end_;
  int dx = mt_->delta_x_;
  if (index_ == i_plus) {
    if (mt_->lattice_[index_ - dx].xlink_head_ != nullptr)
      n_neighbs++;
  } else if (index_ == i_minus) {
    if (mt_->lattice_[index_ + dx].xlink_head_ != nullptr)
      n_neighbs++;
  } else {
    if (mt_->lattice_[index_ - dx].xlink_head_ != nullptr)
      n_neighbs++;
    if (mt_->lattice_[index_ + dx].xlink_head_ != nullptr)
      n_neighbs++;
  }
  return n_neighbs;
}

int Tubulin::GetKIF4ANeighborCount() {

  int n_neighbs{0};
  int i_plus = mt_->plus_end_;
  int i_minus = mt_->minus_end_;
  int dx = mt_->delta_x_;
  if (index_ == i_plus) {
    if (mt_->lattice_[index_ - dx].motor_head_ != nullptr)
      n_neighbs++;
  } else if (index_ == i_minus) {
    if (mt_->lattice_[index_ + dx].motor_head_ != nullptr)
      n_neighbs++;
  } else {
    if (mt_->lattice_[index_ - dx].motor_head_ != nullptr)
      n_neighbs++;
    if (mt_->lattice_[index_ + dx].motor_head_ != nullptr)
      n_neighbs++;
  }
  return n_neighbs;
}