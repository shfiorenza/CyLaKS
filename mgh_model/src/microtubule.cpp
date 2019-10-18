#include "microtubule.h"
#include "master_header.h"

Microtubule::Microtubule() {}

void Microtubule::Initialize(system_parameters *parameters,
                             system_properties *properties, int i_mt) {

  parameters_ = parameters;
  properties_ = properties;
  index_ = i_mt;

  SetParameters();
  GenerateLattice();
}

void Microtubule::SetParameters() {

  n_sites_ = parameters_->microtubules.length[index_];
  coord_ = parameters_->microtubules.start_coord[index_];
  if (index_ % 2 == 0) {
    polarity_ = 0;
    plus_end_ = 0;
    minus_end_ = n_sites_ - 1;
    delta_x_ = -1;
    mt_index_adj_ = index_ + 1; // FIXME
    neighbor_ = &properties_->microtubules.mt_list_[mt_index_adj_];
  } else if (index_ % 2 == 1) {
    polarity_ = 1;
    plus_end_ = n_sites_ - 1;
    minus_end_ = 0;
    delta_x_ = 1;
    mt_index_adj_ = index_ - 1; // FIXME
    neighbor_ = &properties_->microtubules.mt_list_[mt_index_adj_];
  }
  int mt_length = n_sites_;
  double site_size = parameters_->microtubules.site_size;
  double big_l = mt_length * site_size;
  double radius = parameters_->microtubules.radius;
  double height = parameters_->microtubules.elevation;
  double eta = parameters_->eta;
  // see radhika sliding paper for any of this to make sense
  double numerator = 2 * 3.14159 * big_l * (eta / 1000000);
  ;
  double denom = log(2 * height / radius);
  gamma_ = (numerator / denom);
}

void Microtubule::GenerateLattice() {

  lattice_.resize(n_sites_);
  for (int i_site = 0; i_site < n_sites_; i_site++) {
    lattice_[i_site].Initialize(parameters_, properties_, this, i_site);
  }
}

void Microtubule::UpdateExtensions() {

  // Run through all sites on MT
  for (int i_site = 0; i_site < n_sites_; i_site++) {
    Tubulin *site = &lattice_[i_site];
    // Only update extensions from xlinks (both tether and xlink itself)
    if (site->xlink_head_ != nullptr) {
      AssociatedProtein *xlink = site->xlink_head_->xlink_;
      xlink->UpdateExtension();
      // If xlink is tethered, update tether extension
      if (xlink->tethered_) {
        xlink->motor_->UpdateExtension();
      }
    }
  }
}

double Microtubule::GetNetForce() {

  double forces_summed = 0;
  for (int i_site = 0; i_site < n_sites_; i_site++) {
    Tubulin *site = &lattice_[i_site];
    // Check if site is occupied by xlink head
    if (site->xlink_head_ != nullptr) {
      AssociatedProtein *xlink = site->xlink_head_->xlink_;
      // If doubly-bound, get force from self and potentially teth
      if (xlink->heads_active_ == 2) {
        forces_summed += xlink->GetExtensionForce(site);
        if (xlink->tethered_) {
          Kinesin *motor = xlink->motor_;
          // Only bound motors have valid tether extensions
          if (motor->heads_active_ > 0) {
            // Only motors on OTHER MTs can exert a force
            if (motor->mt_ != site->mt_) {
              forces_summed += motor->GetTetherForce(site);
            }
          }
        }
      }
      // Otherwise if singly-bound, check for tether force
      else if (xlink->tethered_) {
        Kinesin *motor = xlink->motor_;
        // Only bound motors have valid tether extensions
        if (motor->heads_active_ > 0) {
          // Only motors on OTHER MTs can exert a force
          if (motor->mt_ != site->mt_) {
            forces_summed += motor->GetTetherForce(site);
          }
        }
      }
    }
    // Otherwise, check if occupied by motor head
    else if (site->motor_head_ != nullptr) {
      Kinesin *motor = site->motor_head_->motor_;
      // Motors can only exert a force if they are tethered
      if (motor->tethered_) {
        AssociatedProtein *xlink = motor->xlink_;
        // singly-bound xlinks must be on other MTs to exert a force
        if (xlink->heads_active_ == 1) {
          Tubulin *xlink_site = xlink->GetActiveHeadSite();
          if (xlink_site->mt_ != site->mt_) {
            // With 1 head active, no danger of double counting
            if (motor->heads_active_ == 1) {
              forces_summed += motor->GetTetherForce(site);
            }
            // With 2 heads active, only get force from front
            else if (site == motor->head_one_.site_) {
              forces_summed += motor->GetTetherForce(site);
            }
          }
        }
        // If xlink is doubly-bound, get reaction force
        else if (xlink->heads_active_ == 2) {
          // With 1 head active, no danger of double counting
          if (motor->heads_active_ == 1) {
            forces_summed += motor->GetTetherForce(site);
          }
          // With 2 heads active, only get force from front
          else if (site == motor->head_one_.site_) {
            forces_summed += motor->GetTetherForce(site);
          }
        }
      }
    }
  }
  return forces_summed;
}

double Microtubule::GetNetForce_Motors() {

  double forces_summed = 0;
  for (int i_site = 0; i_site < n_sites_; i_site++) {
    Tubulin *site = &lattice_[i_site];
    // Check if site is occupied by xlink head
    if (site->xlink_head_ != nullptr) {
      AssociatedProtein *xlink = site->xlink_head_->xlink_;
      // If doubly-bound, get force from self and potentially teth
      if (xlink->heads_active_ == 2) {
        forces_summed += xlink->GetExtensionForce(site);
        if (xlink->tethered_) {
          Kinesin *motor = xlink->motor_;
          // Only bound motors have valid tether extensions
          if (motor->heads_active_ > 0) {
            // Only motors on OTHER MTs can exert a force
            if (motor->mt_ != site->mt_) {
              forces_summed += motor->GetTetherForce(site);
            }
          }
        }
      }
      // Otherwise if singly-bound, check for tether force
      else if (xlink->tethered_) {
        Kinesin *motor = xlink->motor_;
        // Only bound motors have valid tether extensions
        if (motor->heads_active_ > 0) {
          // Only motors on OTHER MTs can exert a force
          if (motor->mt_ != site->mt_) {
            forces_summed += motor->GetTetherForce(site);
          }
        }
      }
    }
    // Otherwise, check if occupied by motor head
    else if (site->motor_head_ != nullptr) {
      Kinesin *motor = site->motor_head_->motor_;
      // Motors can only exert a force if they are tethered
      if (motor->tethered_) {
        AssociatedProtein *xlink = motor->xlink_;
        // singly-bound xlinks must be on other MTs to exert a force
        if (xlink->heads_active_ == 1) {
          Tubulin *xlink_site = xlink->GetActiveHeadSite();
          if (xlink_site->mt_ != site->mt_) {
            // With 1 head active, no danger of double counting
            if (motor->heads_active_ == 1) {
              forces_summed += motor->GetTetherForce(site);
            }
            // With 2 heads active, only get force from first
            else if (site == motor->head_one_.site_) {
              forces_summed += motor->GetTetherForce(site);
            }
          }
        }
        // If xlink is doubly-bound, get reaction force
        else if (xlink->heads_active_ == 2) {
          // With 1 head active, no danger of double counting
          if (motor->heads_active_ == 1) {
            forces_summed += motor->GetTetherForce(site);
          }
          // With 2 heads active, only get force from first
          else if (site == motor->head_one_.site_) {
            forces_summed += motor->GetTetherForce(site);
          }
        }
      }
    }
  }
  return forces_summed;
}

double Microtubule::GetNetForce_Xlinks() {

  double forces_summed = 0;
  for (int i_site = 0; i_site < n_sites_; i_site++) {
    Tubulin *site = &lattice_[i_site];
    // Check if site is occupied by an xlink
    if (site->xlink_head_ != nullptr) {
      AssociatedProtein *xlink = site->xlink_head_->xlink_;
      // Xlinks can only exert forces if they're double bound
      if (xlink->heads_active_ == 2) {
        forces_summed += xlink->GetExtensionForce(site);
      }
    }
  }
  return forces_summed;
}
