#include "cylaks/protofilament.hpp"
#include "cylaks/system_namespace.hpp"
#include "cylaks/system_parameters.hpp"

void Protofilament::SetParameters() {

  using namespace Params;
  n_sites_ = Filaments::n_sites[index_];
  pos_[0] = Filaments::x_initial[index_];
  pos_[1] = Filaments::y_initial[index_];
  orientation_[0] = 1.0;
  orientation_[1] = 0.0; // Begin aligned with x-axis
  immobile_until_.resize(2);
  immobile_until_[0] = Filaments::x_immobile_until[index_] / dt; // n_steps
  immobile_until_[1] = Filaments::y_immobile_until[index_] / dt; // n_steps
  length_ = Filaments::site_size * Filaments::n_sites[index_];   // nm
  polarity_ = Filaments::polarity[index_];
  polarity_ == 0 ? dx_ = -1 : dx_ = 1;
  dt_eff_ = dt / Filaments::n_bd_per_kmc;       // s
  double ar{length_ / (2 * Filaments::radius)}; // unitless aspect ratio
  // Make sure the denominator for gamma_[2] (gamma_rot) is greater than 0.0
  if (ar <= 0.8 / 3.0 and Filaments::rotation_enabled[index_]) {
    Sys::Log("Filament #%i aspect ratio is too small for the form of gamma_rot "
             "we use. Please increase filament length.\n",
             index_);
    Sys::ErrorExit("Protofilament::SetParameters()");
  }
  double eta_adj{eta * 1e-06};                      // pN*s/nm^2
  double pi{M_PI};                                  // literally just pi
  gamma_[0] = 2 * pi * eta_adj * length_ / log(ar); // pN*s/nm
  gamma_[1] = 2 * gamma_[0];                        // pN*s/nm
  gamma_[2] = pi * eta_adj * Cube(length_) / (3 * (log(ar) - 0.8)); // pN*s*nm
  for (int i_dim{0}; i_dim < sigma_.size(); i_dim++) {
    sigma_[i_dim] = sqrt(2 * kbT * dt_eff_ / gamma_[i_dim]); // nm or rad
  }
}

void Protofilament::SetParametersNucleated(Protofilament *parent) {

  using namespace Params;
  // using namespace Filaments;
  n_sites_ = 100; // Filaments::n_sites[index_ - 1];
  double ran{SysRNG::GetRanProb()};
  pos_[0] = parent->pos_[0] + (ran - 0.5) * parent->length_;
  pos_[1] = parent->pos_[1] + 10.0;
  orientation_[0] = parent->orientation_[0]; // 1.0
  orientation_[1] = parent->orientation_[1]; // 0.0
  immobile_until_.resize(2);
  immobile_until_[0] =
      0; // Filaments::x_immobile_until[index_ - 1] / dt; // n_steps
  immobile_until_[1] = Filaments::y_immobile_until[0] / dt; // n_steps
  length_ = Filaments::site_size * n_sites_;                // nm
  polarity_ = parent->polarity_;
  polarity_ == 0 ? dx_ = -1 : dx_ = 1;
  dt_eff_ = dt / Filaments::n_bd_per_kmc; // s
  // Filaments::n_sites.push_back(n_sites_);
  // Filaments::x_initial.push_back(pos_[0]);
  // Filaments::y_initial.push_back(pos_[1]);
  // Filaments::x_immobile_until.push_back(immobile_until_[0]);
  // Filaments::y_immobile_until.push_back(immobile_until_[1]);
  // Filaments::polarity.push_back(polarity_);
  Filaments::rotation_enabled.push_back(false);
  double ar{length_ / (2 * Filaments::radius)}; // unitless aspect ratio
  // Make sure the denominator for gamma_[2] (gamma_rot) is greater than 0.0
  if (ar <= 0.8 / 3.0 and Filaments::rotation_enabled[index_ - 1]) {
    Sys::Log("Filament #%i aspect ratio is too small for the form of gamma_rot "
             "we use. Please increase filament length.\n",
             index_ - 1);
    Sys::ErrorExit("Protofilament::SetParameters()");
  }
  double eta_adj{eta * 1e-06};                      // pN*s/nm^2
  double pi{M_PI};                                  // literally just pi
  gamma_[0] = 2 * pi * eta_adj * length_ / log(ar); // pN*s/nm
  gamma_[1] = 2 * gamma_[0];                        // pN*s/nm
  gamma_[2] = pi * eta_adj * Cube(length_) / (3 * (log(ar) - 0.8)); // pN*s*nm
  for (int i_dim{0}; i_dim < sigma_.size(); i_dim++) {
    sigma_[i_dim] = sqrt(2 * kbT * dt_eff_ / gamma_[i_dim]); // nm or rad
  }
}

void Protofilament::GenerateSites() {

  sites_.resize(n_sites_);
  // Initialize sites
  for (int i_entry{0}; i_entry < n_sites_; i_entry++) {
    sites_[i_entry].Initialize(_id_site, Sys::n_objects_++, _r_site, i_entry,
                               this);
  }
  // Set site neighbors (immediately forward/behind; 2 max on a 1-D lattice)
  for (auto &&site : sites_) {
    int i_fwd{(int)site.index_ + 1};
    if (i_fwd < sites_.size()) {
      site.AddNeighbor(&sites_[i_fwd]);
    }
    int i_bck{(int)site.index_ - 1};
    if (i_bck >= 0) {
      site.AddNeighbor(&sites_[i_bck]);
    }
  }
  plus_end_ = &sites_[(n_sites_ - 1) * polarity_];
  minus_end_ = &sites_[(n_sites_ - 1) * (1.0 - polarity_)];
  // plus_end_->SetBindingAffinity(0.1);
  // for (int i_site{0}; i_site < 10; i_site++) {
  //   int index = plus_end_->index_ - (dx_ * i_site);
  //   sites_[index].SetBindingAffinity(0.1);
  //   printf("site %i binding affinity set\n", index);
  // }
  Sys::Log(2, "     plus_end = site %i\n", plus_end_->index_);
  Sys::Log(2, "     minus_end = site %i\n", minus_end_->index_);
  center_index_ = double(n_sites_ - 1) / 2;
}

void Protofilament::UpdateRodPosition() {

  double noise_par{SysRNG::GetGaussianNoise(sigma_[0])};
  double noise_perp{SysRNG::GetGaussianNoise(sigma_[1])};
  double noise_rot{SysRNG::GetGaussianNoise(sigma_[2])};

  noise_par = noise_perp = noise_rot = 0.0;

  // First row is a unit vector (in lab frame) along length of rod
  // Second row is a unit vector (in lab frame) perpendicular to length of rod
  Vec2D<double> rod_basis{GetOrthonormalBasis(orientation_)};
  /* c.f. Tao et al., J. Chem. Phys. (2005); doi.org/10.1063/1.1940031 */
  Vec2D<double> xi_inv(_n_dims_max, Vec<double>(_n_dims_max, 0.0));
  for (int i{0}; i < _n_dims_max; i++) {
    for (int j{i}; j < _n_dims_max; j++) {
      double uiuj{orientation_[i] * orientation_[j]};
      xi_inv[i][j] = xi_inv[j][i] = uiuj * (1.0 / gamma_[0] - 1.0 / gamma_[1]);
      if (i == j) {
        xi_inv[i][j] += 1.0 / gamma_[1];
      }
    }
  }
  // Apply translationl and rotational displacements
  Vec<double> torque_proj{Cross(torque_, orientation_)};
  double u_norm{0.0};
  for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
    // Only update position if protofilament isnt immobilized
    if (Sys::i_step_ > immobile_until_[i_dim]) {
      double vel{Dot(xi_inv[i_dim], force_)};
      velocity_[i_dim] = vel;
      // if (i_dim == 0 and vel != 50) {
      //   printf("v[%i] = %g\n", i_dim, vel);
      // }
      pos_[i_dim] += vel * dt_eff_;
      pos_[i_dim] += rod_basis[0][i_dim] * noise_par;
      pos_[i_dim] += rod_basis[1][i_dim] * noise_perp;
      // velocity_[i_dim] += rod_basis[0][i_dim] * noise_par / dt_eff_;
      // printf("%g\n", rod_basis[0][i_dim] * noise_par / dt_eff_);
      // velocity_[i_dim] += rod_basis[1][i_dim] * noise_perp / dt_eff_;
      // printf("%g\n", rod_basis[1][i_dim] * noise_perp / dt_eff_);
      // Check for NaN positions
      if (pos_[i_dim] != pos_[i_dim]) {
        Sys::Log("force = %g\n", force_[i_dim]);
        Sys::ErrorExit("Protofilament::UpdateRodPositions() [1]");
      }
    }
    // Only update orientation if rotation is enabled
    if (Params::Filaments::rotation_enabled[index_]) {
      orientation_[i_dim] += torque_proj[i_dim] / gamma_[2] * dt_eff_;
      orientation_[i_dim] += rod_basis[1][i_dim] * noise_rot;
      // Check for NaN orientations
      if (orientation_[i_dim] != orientation_[i_dim]) {
        Sys::Log("torque_proj = %g\n", torque_proj[i_dim]);
        Sys::ErrorExit("Protofilament::UpdateRodPositions() [2]");
      }
      u_norm += Square(orientation_[i_dim]);
    }
  }
  if (Params::Filaments::rotation_enabled[index_]) {
    // Re-normalize orientation vector
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      orientation_[i_dim] /= sqrt(u_norm);
    }
  }
}

void Protofilament::UpdateSitePositions() {

  // If proteins are disabled (i.e., just PFs present), update endpoints only
  if (Params::Motors::c_bulk == 0.0 and Params::Xlinks::c_bulk == 0.0) {
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      // Distance will be negative for first half of sites
      double p_dist{double(plus_end_->index_) - center_index_};
      p_dist *= Params::Filaments::site_size; // convert to nm
      // Orientation always points towards increasing site index
      plus_end_->pos_[i_dim] = pos_[i_dim] + p_dist * Dot(orientation_, i_dim);
      // Distance will be negative for first half of sites
      double m_dist{double(minus_end_->index_) - center_index_};
      m_dist *= Params::Filaments::site_size; // convert to nm
      // Orientation always points towards increasing site index
      minus_end_->pos_[i_dim] = pos_[i_dim] + m_dist * Dot(orientation_, i_dim);
    }
    return;
  }
  for (auto &&site : sites_) {
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      // Distance will be negative for first half of sites
      double dist{double(site.index_) - center_index_};
      dist *= Params::Filaments::site_size; // convert to nm
      // Orientation always points towards increasing site index
      site.pos_[i_dim] = pos_[i_dim] + dist * Dot(orientation_, i_dim);
    }
  }
  Sys::Log(3, "%zu & %zu\n", plus_end_->pos_.size(), minus_end_->pos_.size());
  Sys::Log(3, "plus-end: (%g, %g)\n", plus_end_->pos_[0], plus_end_->pos_[1]);
  Sys::Log(3, "minus_end: (%g, %g)\n", minus_end_->pos_[0],
           minus_end_->pos_[1]);
}

BindingSite *Protofilament::GetNeighb(BindingSite *site, int delta) {

  using namespace Params;
  Sys::Log(2, "i_site = %i, delta = %i\n", site->index_, delta);
  if (site->filament_ == this) {
    Sys::ErrorExit("Protofilament::GetNeighb()");
  }
  // First, we find which site best aligns vertically w/ given site
  double site_x{site->pos_[0]};
  // x-coords equal, so site_pos_x = (i_align - center_index) * site_size + pos
  int i_aligned{(int)std::round((site_x - pos_[0]) / Filaments::site_size +
                                center_index_)};
  // Scan relative to aligned site using given delta value
  int i_neighb{i_aligned + delta};
  Sys::Log(2, "i_neighb is %i\n", i_neighb);
  if (i_neighb < 0 or i_neighb > sites_.size() - 1) {
    return nullptr;
  }
  return &sites_[i_neighb];
}

void Protofilament::AddSite() {

  int i_plus = plus_end_->index_;
  n_sites_++;
  sites_.emplace_back();
  sites_.back().Initialize(_id_site, Sys::n_objects_++, _r_site, n_sites_,
                           this);
  if (i_plus == 0) {
    plus_end_ = &sites_[0];
    minus_end_ = &sites_.back();
    pos_[0] += -Params::Filaments::site_size / 2.0;
  } else {
    minus_end_ = &sites_[0];
    plus_end_ = &sites_.back();
    pos_[0] += Params::Filaments::site_size / 2.0;
  }
  center_index_ = double(n_sites_ - 1) / 2;
}

void Protofilament::RemoveSite() {
  int i_plus = plus_end_->index_;
  if (n_sites_ == 2) {
    return;
  }
  n_sites_--;
  sites_.pop_back();

  if (i_plus == 0) {
    plus_end_ = &sites_[0];
    minus_end_ = &sites_.back();
    pos_[0] += Params::Filaments::site_size / 2.0;
  } else {
    minus_end_ = &sites_[0];
    plus_end_ = &sites_.back();
    pos_[0] += -Params::Filaments::site_size / 2.0;
  }
  center_index_ = double(n_sites_ - 1) / 2;
}