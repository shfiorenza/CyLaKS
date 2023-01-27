#include "cylaks/filament_manager.hpp"
#include "cylaks/protein_manager.hpp"

void FilamentManager::SetParameters() {

  threshold_ = std::pow(2, 1.0 / 6.0) * sigma_;
  n_bd_iterations_ = Params::Filaments::n_bd_per_kmc;
  dt_eff_ = Params::dt / n_bd_iterations_;
}

void FilamentManager::GenerateFilaments() {

  // Create the appropriate number of protofilaments
  if (Params::Filaments::n_subfilaments <= 1) {
    protofilaments_.resize(Params::Filaments::count);
  } else {
    protofilaments_.resize(Params::Filaments::n_subfilaments);
    // int n_sub{Params::Filaments::n_subfilaments};
    using namespace Params::Filaments;
    n_sites = std::vector(n_subfilaments, n_sites[0]);
    polarity = std::vector(n_subfilaments, polarity[0]);
    x_initial = std::vector(n_subfilaments, x_initial[0]);
    y_initial = std::vector(n_subfilaments, y_initial[0]);
    // printf("\n\n");
    for (int i_sub{1}; i_sub < n_subfilaments; i_sub++) {
      y_initial[i_sub] =
          y_initial[0] + (i_sub * 2 * Params::Filaments::site_size);
      Sys::Log("n_sites[%i] = %zu\n", i_sub, n_sites[i_sub]);
      Sys::Log("polarity[%i] = %zu\n", i_sub, polarity[i_sub]);
    }
    // printf("\n\n");
    f_applied = std::vector(n_subfilaments, f_applied[0]);
    x_immobile_until = std::vector(n_subfilaments, x_immobile_until[0]);
    y_immobile_until = std::vector(n_subfilaments, y_immobile_until[0]);
  }
  // Initialize the protofilaments we implicitly created w/ the resize
  for (int i_fil{0}; i_fil < protofilaments_.size(); i_fil++) {
    protofilaments_[i_fil].Initialize(_id_site, Sys::n_objects_++, i_fil);
  }
  // w/ subfilaments: sets "perodic" boundaries to mimic 3-D microtubule tube
  // w/o subfilaments: sets pair of anti-parallel MTs to be eachother's neighb
  if (protofilaments_.size() > 1) {
    if (Params::Filaments::n_subfilaments > 1) {
      if (Params::Filaments::count > 1) {
        Sys::ErrorExit("Multiple MTs w/ explicit PFs not implemented yet! [1]");
      }
      // Use "top" and "bot" neighbors to designate adjacent PFs in array
      for (int i_sub{1}; i_sub < protofilaments_.size() - 1; i_sub++) {
        protofilaments_[i_sub].top_neighb_ = &protofilaments_[i_sub + 1];
        protofilaments_[i_sub].bot_neighb_ = &protofilaments_[i_sub - 1];
      }
      size_t i_end{protofilaments_.size() - 1};
      protofilaments_[0].top_neighb_ = &protofilaments_[1];
      protofilaments_[i_end].bot_neighb_ = &protofilaments_[i_end - 1];
      if (Params::Filaments::periodic_barrel) {
        protofilaments_[0].bot_neighb_ = &protofilaments_[i_end];
        protofilaments_[i_end].top_neighb_ = &protofilaments_[0];
      }
    } else if (Params::Filaments::count > 1) {
      if (Params::Filaments::n_subfilaments > 1) {
        Sys::ErrorExit("Multiple MTs w/ explicit PFs not implemented yet! [1]");
      }
      if (Params::Filaments::count > 2) {
        Sys::ErrorExit("MT bundles beyond 2 not implemented yet!");
      }
      // By default, since only 2 PFs exist, they are each other's neighbor
      protofilaments_[0].neighbor_ = &protofilaments_[1];
      protofilaments_[1].neighbor_ = &protofilaments_[0];
    } else {
      Sys::ErrorExit("how did we get here");
    }
  }
  // Add all sites across all PFs to a master site list
  for (auto &&pf : protofilaments_) {
    for (auto &&site : pf.sites_) {
      sites_.emplace_back(&site);
    }
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

  for (auto const &pf : protofilaments_) {
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      if (Sys::i_step_ > pf.immobile_until_[i_dim]) {
        return false;
      }
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
  // for (auto &&pf : protofilaments_) {
  //   printf("F = <%g, %g> for PF #%i\n", pf.force_[0], pf.force_[1],
  //   pf.index_);
  // }
}

void FilamentManager::UpdateLattice() { proteins_->UpdateLatticeDeformation(); }
