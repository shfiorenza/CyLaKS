#include "cylaks/protofilament.hpp"

void Protofilament::SetParameters() {

  using namespace Params;
  // using namespace Filaments;
  pos_[0] = Filaments::x_initial[index_];
  pos_[1] = Filaments::y_initial[index_];
  orientation_[0] = 1.0;
  orientation_[1] = 0.0; // Begin aligned with x-axis
  length_ = Filaments::site_size * Filaments::n_sites[index_]; // nm
  polarity_ = Filaments::polarity[index_];
  polarity_ == 0 ? dx_ = -1 : dx_ = 1;
  immobile_until_ = Filaments::immobile_until[index_] / dt; // n_steps
  dt_eff_ = dt / Filaments::n_bd_per_kmc;                   // s
  double ar{length_ / (2 * Filaments::radius)}; // unitless aspect ratio
  // Make sure the denominator for gamma_[2] (gamma_rot) is greater than 0.0
  if (ar <= 0.8 / 3.0 and Filaments::rotation_enabled) {
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

void Protofilament::GenerateSites() {

  size_t n_sites{Params::Filaments::n_sites[index_]};
  sites_.resize(n_sites);
  // Initialize sites
  for (int i_entry{0}; i_entry < n_sites; i_entry++) {
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
  plus_end_ = &sites_[(n_sites - 1) * polarity_];
  minus_end_ = &sites_[(n_sites - 1) * (1.0 - polarity_)];
  Sys::Log(2, "     plus_end = site %i\n", plus_end_->index_);
  Sys::Log(2, "     minus_end = site %i\n", minus_end_->index_);
  center_index_ = double(n_sites - 1) / 2;
}

void Protofilament::UpdateRodPosition() {

  double noise_par{SysRNG::GetGaussianNoise(sigma_[0])};
  double noise_perp{SysRNG::GetGaussianNoise(sigma_[1])};
  double noise_rot{SysRNG::GetGaussianNoise(sigma_[2])};

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
    // Only update pos in dimensions with translational movement enabled
    if (Params::Filaments::translation_enabled[i_dim]) {
      pos_[i_dim] += Dot(xi_inv[i_dim], force_) * dt_eff_;
      pos_[i_dim] += rod_basis[0][i_dim] * noise_par;
      pos_[i_dim] += rod_basis[1][i_dim] * noise_perp;
      // Check for NaN positions
      if (pos_[i_dim] != pos_[i_dim]) {
        Sys::Log("force = %g\n", force_[i_dim]);
        Sys::ErrorExit("Protofilament::UpdateRodPositions() [1]");
      }
    }
    // Only update orientation if rotation is enabled
    if (Params::Filaments::rotation_enabled) {
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
  if (Params::Filaments::rotation_enabled) {
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
  // First, we find which site best aligns vertically w/ given site
  int site_x{(int)site->pos_[0]};
  // x-coords equal, so site_pos_x = (i_align - center_index) * site_size + pos
  int i_aligned{int((site_x - pos_[0]) / Filaments::site_size + center_index_)};
  // Scan relative to aligned site using given delta value
  int i_neighb{i_aligned + delta};
  Sys::Log(2, "i_neighb is %i\n", i_neighb);
  if (i_neighb < 0 or i_neighb > sites_.size() - 1) {
    return nullptr;
  }
  return &sites_[i_neighb];
}