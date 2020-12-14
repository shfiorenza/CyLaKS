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
    Sys::Log(1, "     gamma_%s = %g pN*s%s\n", label[i_dim].c_str(),
             gamma_[i_dim], i_dim < 2 ? "/nm" : "*nm");
    Sys::Log(1, "      sigma_%s = %g %s\n", label[i_dim].c_str(), sigma_[i_dim],
             i_dim < 2 ? "nm" : "rad");
    Sys::Log(1, "       D_%s = %g %s^2/s\n", label[i_dim].c_str(),
             diffusion_const[i_dim], i_dim < 2 ? "nm" : "rad");
  }
}

void Protofilament::GenerateSites() {

  size_t n_sites{Params::Filaments::n_sites[index_]};
  sites_.resize(n_sites);
  for (int i_entry{0}; i_entry < n_sites; i_entry++) {
    sites_[i_entry].Initialize(_id_site, Sys::n_unique_objects_++, _r_site,
                               i_entry, this);
  }
  plus_end_ = &sites_[(n_sites - 1) * polarity_];
  minus_end_ = &sites_[(n_sites - 1) * (1.0 - polarity_)];
  Sys::Log("    plus_end = site %i\n", plus_end_->index_);
  Sys::Log("    minus_end = site %i\n", minus_end_->index_);
  center_index_ = double(n_sites - 1) / 2;
}

void Protofilament::UpdateSitePositions() {

  for (auto &&site : sites_) {
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      // Distance will be negative for first half of sites
      double dist{double(site.index_) - center_index_};
      dist *= Params::Filaments::site_size; // convert to nm
      // Orientation always points towards increasing site index
      site.pos_[i_dim] = pos_[i_dim] + dist * Dot(orientation_, i_dim);
      // Sys::Log("%g\n", Dot(orientation_, i_dim));
      // Sys::Log("%g\n", dist);
      // Sys::Log("%g\n", site.pos_[i_dim]);
    }
    // Sys::Log("site %i: (%g, %g)\n", site.index_, site.pos_[0], site.pos_[1]);
  }
  /*
  Sys::Log("%zu & %zu\n", plus_end_->pos_.size(), minus_end_->pos_.size());
  Sys::Log("plus-end: (%g, %g)\n", plus_end_->pos_[0], plus_end_->pos_[1]);
  Sys::Log("minus_end: (%g, %g)\n", minus_end_->pos_[0], minus_end_->pos_[1]);
  */
}