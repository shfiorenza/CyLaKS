#include "microtubule.h"
#include "master_header.h"

Microtubule::Microtubule() {}

void Microtubule::Initialize(system_parameters *parameters,
                             system_properties *properties, int i_mt) {

  wally_ = &properties->wallace;
  parameters_ = parameters;
  properties_ = properties;
  index_ = i_mt;
  SetParameters();
  GenerateLattice();
}

void Microtubule::SetParameters() {

  applied_force_ = parameters_->microtubules.applied_force;
  n_sites_ = parameters_->microtubules.length[index_];
  coord_ = parameters_->microtubules.start_coord[index_];
  immobile_until_ = parameters_->microtubules.immobile_until[index_];
  if (index_ % 2 == 0) {
    plus_end_ = 0;
    minus_end_ = n_sites_ - 1;
    delta_x_ = -1;
  } else if (index_ % 2 == 1) {
    plus_end_ = n_sites_ - 1;
    minus_end_ = 0;
    delta_x_ = 1;
  }
  double site_size{parameters_->microtubules.site_size};
  double length{n_sites_ * site_size};                   // in nm
  double diameter{2 * parameters_->microtubules.radius}; // in nm
  double eta{parameters_->eta * 1e-06};                  // in (pN*s)/nm^2
  // From The Theory of Polymer Dynamics by Doi & Edwards, pg. 296 eq 8.20
  gamma_ = (2 * M_PI * eta * length) / log(length / diameter);
  double dt_eff{parameters_->delta_t / parameters_->microtubules.n_iterations};
  sigma_ = sqrt(2 * parameters_->kbT * dt_eff / gamma_);     // in nm
  double diffusion_const{parameters_->kbT / gamma_};         // in nm^2/s
  double tau{site_size * site_size / (2 * diffusion_const)}; // in s
  wally_->Log("\nFor microtubule #%i:\n", index_);
  wally_->Log("    Length = %i sites\n", n_sites_);
  wally_->Log("    Gamma = %g (pN*s)/nm\n", gamma_);
  wally_->Log("    D = %g um^2/s\n", diffusion_const * 1e-06);
  wally_->Log("    Tau = %g seconds\n", tau);
}

void Microtubule::GenerateLattice() {

  lattice_.resize(n_sites_);
  for (int i_site = 0; i_site < n_sites_; i_site++) {
    lattice_[i_site].Initialize(parameters_, properties_, this, i_site);
  }
}

double Microtubule::GetNetForce() {

  double forces_summed{applied_force_};
  for (int i_site{0}; i_site < n_sites_; i_site++) {
    Tubulin *site{&lattice_[i_site]};
    // If site is unoccupied, continue on in for loop
    if (!site->occupied_) {
      continue;
    }
    // Check if site is occupied by an xlink head
    if (site->xlink_head_ != nullptr) {
      AssociatedProtein *xlink{site->xlink_head_->xlink_};
      // Add xlink force (returns 0.0 if not doubly-bound)
      forces_summed += xlink->GetExtensionForce(site);
      // If xlink is tethered to non-satellite motor, add tether force too
      if (xlink->tethered_ and !xlink->HasSatellite()) {
        // Only motors tethered to OTHER MTs can exert a force
        if (xlink->motor_->mt_ != site->mt_) {
          forces_summed += xlink->motor_->GetTetherForce(site);
        }
      }
    }
    // Check if site is occupied by a motor head
    else if (site->motor_head_ != nullptr) {
      Kinesin *motor{site->motor_head_->motor_};
      // Only motors tethered to non-satellite xlinks can exert a force
      if (motor->tethered_ and !motor->HasSatellite()) {
        // If tethered xlink is doubly-bound, add reactionary force
        if (motor->xlink_->heads_active_ == 2) {
          // (returns 0.0 if site == head_two_.site to avoid double counting)
          forces_summed += motor->GetTetherForce(site);
        }
        // Otherwise, we need to ensure the xlink isn't on the same MT
        else if (motor->mt_ != motor->xlink_->GetActiveHead()->site_->mt_) {
          // (returns 0.0 if site == head_two_.site to avoid double counting)
          forces_summed += motor->GetTetherForce(site);
        }
      }
    } else {
      wally_->ErrorExit("Microtubule::GetNetForce()");
    }
  }
  return forces_summed;
}

double Microtubule::GetNetForce_Motors() {

  double forces_summed{0.0};
  for (int i_site{0}; i_site < n_sites_; i_site++) {
    Tubulin *site{&lattice_[i_site]};
    // If site is unoccupied, continue on in for loop
    if (!site->occupied_) {
      continue;
    }
    // Check if site is occupied by an xlink head
    if (site->xlink_head_ != nullptr) {
      AssociatedProtein *xlink{site->xlink_head_->xlink_};
      // If xlink is tethered to non-satellite motor, add tether force too
      if (xlink->tethered_ and !xlink->HasSatellite()) {
        // Only motors tethered to OTHER MTs can exert a force
        if (xlink->motor_->mt_ != site->mt_) {
          forces_summed += xlink->motor_->GetTetherForce(site);
        }
      }
    }
    // Check if site is occupied by a motor head
    else if (site->motor_head_ != nullptr) {
      Kinesin *motor{site->motor_head_->motor_};
      // Only motors tethered to non-satellite xlinks can exert a force
      if (motor->tethered_ and !motor->HasSatellite()) {
        // If tethered xlink is doubly-bound, add reactionary force
        if (motor->xlink_->heads_active_ == 2) {
          // (returns 0.0 if site == head_two_.site to avoid double counting)
          forces_summed += motor->GetTetherForce(site);
        }
        // Otherwise, we need to ensure the xlink isn't on the same MT
        else if (motor->mt_ != motor->xlink_->GetActiveHead()->site_->mt_) {
          // (returns 0.0 if site == head_two_.site to avoid double counting)
          forces_summed += motor->GetTetherForce(site);
        }
      }
    } else {
      wally_->ErrorExit("Microtubule::GetNetForce()");
    }
  }
  return forces_summed;
}

double Microtubule::GetNetForce_Xlinks() {

  double forces_summed{0.0};
  for (int i_site{0}; i_site < n_sites_; i_site++) {
    Tubulin *site{&lattice_[i_site]};
    // If site doesn't have an xlink, continue on in for loop
    if (site->xlink_head_ == nullptr) {
      continue;
    }
    AssociatedProtein *xlink{site->xlink_head_->xlink_};
    // Add xlink force (returns 0.0 if not doubly-bound)
    forces_summed += xlink->GetExtensionForce(site);
  }
  return forces_summed;
}
