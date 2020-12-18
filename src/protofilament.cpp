#include "protofilament.hpp"

void Protofilament::SetParameters() {

  using namespace Params;
  using namespace Filaments;
  pos_[0] = x_initial[index_];
  pos_[1] = y_initial[index_];
  orientation_[0] = 1.0;
  orientation_[1] = 0.0;                 // Begin aligned with x-axis
  length_ = site_size * n_sites[index_]; // nm
  polarity_ = polarity[index_];
  polarity_ == 0 ? dx_ = -1 : dx_ = 1;
  immobile_until_ = immobile_until[index_] / dt;
  dt_eff_ = dt / n_bd_per_kmc;
  double ar{length_ / (2 * radius)};                // unitless aspect ratio
  double eta_adj{eta * 1e-06};                      // pN*s/nm^2
  double pi{M_PI};                                  // literally just pi
  gamma_[0] = 2 * pi * eta_adj * length_ / log(ar); // pN*s/nm
  gamma_[1] = 2 * gamma_[0];                        // pN*s/nm
  gamma_[2] = pi * eta_adj * Cube(length_) / (3 * (log(ar) - 0.8)); // pN*s*nm
  Vec<Str> label{"par", "perp", "rot"};
  Vec<double> diffusion_const(3, 0.0);
  double dt_eff{dt / n_bd_per_kmc};
  Sys::Log("  For filament #%i:\n", index_);
  Sys::Log("    length = %g nm\n", length_);
  for (int i_dim{0}; i_dim < 3; i_dim++) {
    sigma_[i_dim] = sqrt(2 * kbT * dt_eff / gamma_[i_dim]); // nm or rad
    diffusion_const[i_dim] = kbT / gamma_[i_dim];           // nm^2/s or rad^2/s
    Sys::Log(1, "    gamma_%s = %g pN*s%s\n", label[i_dim].c_str(),
             gamma_[i_dim], i_dim < 2 ? "/nm" : "*nm");
    Sys::Log(2, "     sigma_%s = %g %s\n", label[i_dim].c_str(), sigma_[i_dim],
             i_dim < 2 ? "nm" : "rad");
    Sys::Log(2, "      D_%s = %g %s^2/s\n", label[i_dim].c_str(),
             diffusion_const[i_dim], i_dim < 2 ? "nm" : "rad");
  }
}

void Protofilament::GenerateSites() {

  size_t n_sites{Params::Filaments::n_sites[index_]};
  sites_.resize(n_sites);
  // Initialize sites
  for (int i_entry{0}; i_entry < n_sites; i_entry++) {
    sites_[i_entry].Initialize(_id_site, Sys::n_unique_objects_++, _r_site,
                               i_entry, this);
  }
  // Set site neighbors (immediately forward/behind; 2 max on a 1-D lattice)
  for (auto &&site : sites_) {
    int i_fwd{site.index_ + 1};
    if (i_fwd < sites_.size()) {
      site.AddNeighbor(&sites_[i_fwd]);
    }
    int i_bck{site.index_ - 1};
    if (i_bck > 0) {
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

  // Independent terms for rod trans/rotational diffusion
  double noise_par{SysRNG::GetGaussianNoise(sigma_[0])};
  double noise_perp{SysRNG::GetGaussianNoise(sigma_[1])};
  double noise_rot{SysRNG::GetGaussianNoise(sigma_[2])};
  // Construct body_frame matrix, which transforms rod body to lab frame
  // First row is a unit vector (in lab frame) along length of rod
  // Second row is a unit vector (in lab frame) perpendicular to length of rod
  // Third row (in 3-D only) is the same as 2nd, but also perpendicular to that
  double body_frame[2][_n_dims_max];
  body_frame[0][0] = orientation_[0];
  body_frame[0][1] = orientation_[1];
  body_frame[1][0] = orientation_[1];
  body_frame[1][1] = -orientation_[0];
  /* c.f. Tao et al., J. Chem. Phys. (2005); doi.org/10.1063/1.1940031 */
  // Construct xi tensor
  double xi[_n_dims_max][_n_dims_max];
  for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
    for (int j_dim{i_dim}; j_dim < _n_dims_max; j_dim++) {
      double uiuj{orientation_[i_dim] * orientation_[j_dim]};
      xi[i_dim][j_dim] = xi[j_dim][i_dim] = uiuj * (gamma_[0] - gamma_[1]);
      if (i_dim == j_dim) {
        xi[i_dim][j_dim] += gamma_[1];
      }
    }
  }
  // Calculate inverse of xi tensor, which will be used to find velocity
  double xi_inv[_n_dims_max][_n_dims_max];
  double det{xi[0][0] * xi[1][1] - xi[0][1] * xi[1][0]};
  xi_inv[0][0] = xi[1][1] / det;
  xi_inv[0][1] = -xi[1][0] / det;
  xi_inv[1][0] = -xi[0][1] / det;
  xi_inv[1][1] = xi[0][0] / det;
  // Apply translationl and rotational displacements
  Vec<double> torque_proj{Cross(torque_, orientation_)};
  double u_norm{0.0};
  for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
    for (int j_dim{0}; j_dim < _n_dims_max; j_dim++) {
      pos_[i_dim] += xi_inv[i_dim][j_dim] * force_[j_dim] * dt_eff_;
    }
    pos_[i_dim] += body_frame[0][i_dim] * noise_par;
    pos_[i_dim] += body_frame[1][i_dim] * noise_perp;
    orientation_[i_dim] += torque_proj[i_dim] * dt_eff_ / gamma_[2];
    orientation_[i_dim] += body_frame[1][i_dim] * noise_rot;
    u_norm += Square(orientation_[i_dim]);
  }
  // Re-normalize orientation vector
  for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
    orientation_[i_dim] /= sqrt(u_norm);
  }
}

void Protofilament::UpdateSitePositions() {

  for (auto &&site : sites_) {
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      // Distance will be negative for first half of sites
      double dist{double(site.index_) - center_index_};
      dist *= Params::Filaments::site_size; // convert to nm
      // Orientation always points towards increasing site index
      site.pos_[i_dim] = pos_[i_dim] + dist * Dot(orientation_, i_dim);
    }
  }
  /*
  Sys::Log("%zu & %zu\n", plus_end_->pos_.size(), minus_end_->pos_.size());
  Sys::Log("plus-end: (%g, %g)\n", plus_end_->pos_[0], plus_end_->pos_[1]);
  Sys::Log("minus_end: (%g, %g)\n", minus_end_->pos_[0], minus_end_->pos_[1]);
  */
}

BindingSite *Protofilament::GetNeighb(BindingSite *site, int delta) {

  double pos_x{site->pos_[0]};
  int i_delta{
      (int)std::round((pos_[0] - pos_x) / Params::Filaments::site_size)};
  int i_site{i_delta + (int)center_index_};
  if (i_site < 0 or i_site > sites_.size() - 1) {
    return nullptr;
  }
  return &sites_[i_site];
}