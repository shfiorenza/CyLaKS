#include "cylaks/filament_manager.hpp"
#include "cylaks/protein_manager.hpp"

void FilamentManager::SetParameters() {

  // Filaments are mobile so long as at least 1 dimension is enabled
  for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
    if (Params::Filaments::translation_enabled[i_dim]) {
      mobile_ = true;
    }
  }
  if (Params::Filaments::rotation_enabled) {
    mobile_ = true;
  }
  threshold_ = std::pow(2, 1.0 / 6.0) * sigma_;
  n_bd_iterations_ = Params::Filaments::n_bd_per_kmc;
  dt_eff_ = Params::dt / n_bd_iterations_;
}

void FilamentManager::GenerateFilaments() {

  protofilaments_.resize(Params::Filaments::count);
  for (int i_fil{0}; i_fil < protofilaments_.size(); i_fil++) {
    protofilaments_[i_fil].Initialize(_id_site, Sys::n_objects_++, i_fil);
  }
  for (auto &&pf : protofilaments_) {
    for (auto &&site : pf.sites_) {
      sites_.emplace_back(&site);
    }
  }
  if (protofilaments_.size() == 2) {
    protofilaments_[0].neighbor_ = &protofilaments_[1];
    protofilaments_[1].neighbor_ = &protofilaments_[0];
  }
  using namespace Sys;
  Log("  Filament variables calculated post-initialization:\n");
  for (auto const &pf : protofilaments_) {
    Log("   length[%i] = %g nm\n", pf.index_, pf.length_);
  }
  for (auto const &pf : protofilaments_) {
    Log("    D_par[%i] = %g nm^2/s\n", pf.index_, Params::kbT / pf.gamma_[0]);
  }
  for (auto const &pf : protofilaments_) {
    Log("     D_perp[%i] = %g nm^2/s\n", pf.index_, Params::kbT / pf.gamma_[1]);
  }
  for (auto const &pf : protofilaments_) {
    Log("     D_rot[%i] = %g nm^2/s\n", pf.index_, Params::kbT / pf.gamma_[2]);
  }
  for (auto const &pf : protofilaments_) {
    Log("      gamma_par[%i] = %g nm^2/s\n", pf.index_, pf.gamma_[0]);
  }
  for (auto const &pf : protofilaments_) {
    Log("       gamma_perp[%i] = %g nm^2/s\n", pf.index_, pf.gamma_[1]);
  }
  for (auto const &pf : protofilaments_) {
    Log("       gamma_rot[%i] = %g nm^2/s\n", pf.index_, pf.gamma_[2]);
  }
}

bool FilamentManager::AllFilamentsImmobile() {

  if (!mobile_) {
    return true;
  }
  for (auto const &pf : protofilaments_) {
    if (Sys::i_step_ >= pf.immobile_until_) {
      return false;
    }
  }
  return true;
}

void FilamentManager::UpdateForces() {

  for (auto &&pf : protofilaments_) {
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      pf.force_[i_dim] = Params::Filaments::f_applied[i_dim];
    }
    pf.torque_ = 0.0;
  }
  // TODO: add hard core potential
  if (Params::Filaments::wca_potential_enabled) {
    double r{protofilaments_[1].pos_[1] - protofilaments_[0].pos_[1]};
    if (r < threshold_) {
      double f_mag{
          48 * epsilon_ *
          (Pow(sigma_, 12) / Pow(r, 13) - 0.5 * Pow(sigma_, 6) / Pow(r, 7))};
      protofilaments_[1].force_[1] += f_mag;
      protofilaments_[0].force_[1] -= f_mag;
    }
  }
  proteins_->UpdateExtensions();
}

void FilamentManager::UpdateLattice() { proteins_->UpdateLatticeDeformation(); }