#include "filament_manager.hpp"
#include "protein_manager.hpp"

void FilamentManager::SetParameters() {

  using namespace Params;
  using namespace Filaments;
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
  weight_neighbs_bind_.resize(_n_neighbs_max + 1);
  weight_neighbs_unbind_.resize(_n_neighbs_max + 1);
  for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
    double dE{n_neighbs * -Xlinks::neighb_neighb_energy};
    weight_neighbs_bind_[n_neighbs] = exp(-(1.0 - _lambda_neighb) * dE);
    weight_neighbs_unbind_[n_neighbs] = exp(_lambda_neighb * dE);
  }
}

void FilamentManager::GenerateFilaments() {

  using namespace Sys;
  using namespace Params;
  proto_.resize(Filaments::count);
  for (int i_fil{0}; i_fil < proto_.size(); i_fil++) {
    proto_[i_fil].Initialize(_id_site, n_unique_objects_++, i_fil);
  }
  for (auto &&pf : proto_) {
    for (auto &&site : pf.sites_) {
      sites_.emplace_back(&site);
    }
  }
  if (proto_.size() == 2) {
    proto_[0].neighbor_ = &proto_[1];
    proto_[1].neighbor_ = &proto_[0];
  }
  // FIXME gamma_rot is invalid for MTs that are too short; need better
  // expression --- do not include for now
  int n_dims{2}; // 3};
  Vec<double> D(n_dims, 0.0);
  Log("  Filament variables calculated post-initialization:\n");
  for (auto const &pf : proto_) {
    Log("   length[%i] = %g nm\n", pf.index_, pf.length_);
  }
  for (auto const &pf : proto_) {
    Log("    D_par[%i] = %g nm^2/s\n", pf.index_, kbT / pf.gamma_[0]);
  }
  for (auto const &pf : proto_) {
    Log("     D_perp[%i] = %g nm^2/s\n", pf.index_, kbT / pf.gamma_[1]);
  }
  Log("\n");
  /*
  for (int i_fil{0}; i_fil < proto; i_fil++) {
    Sys::Log("mt #%i\n", proto_[i_fil].index_);
    Sys::Log("%zu & %zu\n", proto_[i_fil].plus_end_->pos_.size(),
             proto_[i_fil].minus_end_->pos_.size());
    Sys::Log("plus-end: (%g, %g)\n", proto_[i_fil].plus_end_->pos_[0],
             proto_[i_fil].plus_end_->pos_[1]);
    Sys::Log("minus_end: (%g, %g)\n", proto_[i_fil].minus_end_->pos_[0],
             proto_[i_fil].minus_end_->pos_[1]);
  }
  exit(1);
  */
}

bool FilamentManager::AllFilamentsImmobile() {

  if (!mobile_) {
    return true;
  }
  for (auto const &pf : proto_) {
    if (Sys::i_step_ >= pf.immobile_until_) {
      return false;
    }
  }
  return true;
}

void FilamentManager::UpdateForces() {
  for (auto &&pf : proto_) {
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      pf.force_[i_dim] = Params::Filaments::f_applied[i_dim];
    }
    pf.torque_ = 0.0;
  }
  if (proto_.size() == 2 and proteins_->xlinks_.active_) {
    // this some jank 1-D wca potential type jawn
    double r{proto_[1].pos_[1] - proto_[0].pos_[1]};
    // printf("threshold is %g\n", threshold_);
    if (r < threshold_) {
      double f_mag{
          48 * epsilon_ *
          (Pow(sigma_, 12) / Pow(r, 13) - 0.5 * Pow(sigma_, 6) / Pow(r, 7))};
      // printf("F_MAG IS %g\n", f_mag);
      proto_[1].force_[1] += f_mag;
      proto_[0].force_[1] -= f_mag;
    }
  }
  proteins_->UpdateExtensions();
}

void FilamentManager::UpdateLattice() { proteins_->UpdateLatticeDeformation(); }