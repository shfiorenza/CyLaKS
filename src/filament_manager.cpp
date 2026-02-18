#include "cylaks/filament_manager.hpp"
#include "cylaks/protein_manager.hpp"
#include "cylaks/protofilament.hpp"
#include "cylaks/system_definitions.hpp"
#include "cylaks/system_namespace.hpp"
#include "cylaks/system_parameters.hpp"
#include "cylaks/system_rng.hpp"

void FilamentManager::SetParameters() {

  threshold_ = std::pow(2, 1.0 / 6.0) * sigma_;
  // threshold_ *= 5;
  // epsilon_ = 0; // 1e-6;
  n_bd_iterations_ = Params::Filaments::n_bd_per_kmc;
  dt_eff_ = Params::dt / n_bd_iterations_;
}

void FilamentManager::GenerateFilaments() {

  // Create the appropriate number of protofilaments
  if (Params::Filaments::axon_arrangement == true) {
    protofilaments_.resize(Params::Filaments::count);
    protofilaments_.reserve(n_pfs_max_);
    size_t block_size{8};
    size_t n_blocks{Params::Filaments::count / block_size};
    for (size_t i_block{1}; i_block < n_blocks; i_block++) {
      for (size_t i_fil{i_block * block_size};
           i_fil < (i_block + 1) * block_size; i_fil++) {
        Params::Filaments::x_initial[i_fil] = 5000 * i_block;

        Params::Filaments::y_initial[i_fil] =
            Params::Filaments::y_initial[i_fil - i_block * block_size];
        // if (i_block % 2 != 0) {
        Params::Filaments::y_initial[i_fil] +=
            (SysRNG::GetRanProb() - 0.5) * 32.0;
        // Params::Filaments::polarity[i_fil] = 0;
        // } else {
        // Params::Filaments::y_initial[i_fil] +=
        //     i_block * 12.0; //+= SysRNG::GetRanProb() * 8.0;
        // Params::Filaments::polarity[i_fil] = 1;
        // }
      }
    }
    double p_parallel{0.5};
    for (int i_fil{0}; i_fil < protofilaments_.size(); i_fil++) {
      if (SysRNG::GetRanProb() < p_parallel) {
        Params::Filaments::polarity[i_fil] = 0;
      } else {
        Params::Filaments::polarity[i_fil] = 1;
      }
    }
    for (int i_fil{0}; i_fil < protofilaments_.size(); i_fil++) {
      protofilaments_[i_fil].Initialize(_id_site, Sys::n_objects_++, i_fil);
    }
    UpdateNeighbors();
    // Use "top" and "bot" neighbors to designate adjacent PFs in axon
    // for (int i_fil{1}; i_fil < protofilaments_.size() - 1; i_fil++) {
    //   protofilaments_[i_fil].top_neighb_ = &protofilaments_[i_fil + 1];
    //   protofilaments_[i_fil].bot_neighb_ = &protofilaments_[i_fil - 1];
    // }
    // size_t i_end{protofilaments_.size() - 1};
    // protofilaments_[0].top_neighb_ = &protofilaments_[1];
    // protofilaments_[i_end].bot_neighb_ = &protofilaments_[i_end - 1];
    return;
  }
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

void FilamentManager::RunKMC() {

  double p_p2g = 0.005;
  double p_g2s = 0.0005;
  double p_s2p = 0.005;
  double v_grow = 120;  // nm/s
  double v_shrink = 60; // nm/s

  double soma_boundary = 150000; // nm

  double p_add_site = v_grow * Params::dt / Params::Filaments::site_size;
  double p_rmv_site = v_shrink * Params::dt / Params::Filaments::site_size;
  if (p_add_site >= 1.0 || p_rmv_site >= 1.0) {
    printf("filament problems. decrease timestep.\n");
    exit(1);
  }
  // Dynamic instability
  for (auto &&pf : protofilaments_) {
    if (pf.plus_end_->pos_[0] > soma_boundary ||
        pf.minus_end_->pos_[0] > soma_boundary) {
      pf.RemoveSite();
      // pf.state_ = shrink;
      // double ran2{SysRNG::GetRanProb()};
      // if (ran2 < p_add_site) {
      // }
      continue;
    }
    double ran{SysRNG::GetRanProb()};
    switch (pf.state_) {
    case pause: {
      if (ran < p_p2g) {
        pf.state_ = grow;
        break;
      }
      break;
    }
    case grow: {
      if (ran < p_g2s) {
        pf.state_ = shrink;
        break;
      }
      double ran2{SysRNG::GetRanProb()};
      if (ran2 < p_add_site) {
        pf.AddSite();
      }
      break;
    }
    case shrink: {
      if (ran < p_s2p) {
        pf.state_ = pause;
        break;
      }
      double ran2{SysRNG::GetRanProb()};
      if (ran2 < p_add_site) {
        pf.RemoveSite();
      }
      break;
    }
    }
  }
  for (auto &&pf : protofilaments_) {
    if (pf.n_sites_ == 2) {
      // printf("boop\n");
      size_t i_pf{pf.index_};
      protofilaments_.erase(protofilaments_.begin() + i_pf);
      for (int i_entry{0}; i_entry < protofilaments_.size(); i_entry++) {
        protofilaments_[i_entry].index_ = i_entry;
      }
      // protofilaments_[i_pf] = protofilaments_.back();
      // protofilaments_[i_pf].index_ = i_pf;
      // protofilaments_.pop_back();
      // protofilaments_.shrink_to_fit();
      // printf("bap\n");
      UpdateNeighbors();
      // AnnihilateProtofilament(pf);
    }
  }
  // protofilaments_.erase(
  //     std::remove_if(protofilaments_.begin(), protofilaments_.end(),
  //                    [](Protofilament pf) { return pf.n_sites_ == 2; }),
  //     protofilaments_.end());
  // Nucleation of new microtubules
  double p_nucleate = 3e-7 * Params::dt; // probability per micron
  double tot_nucleation{0.0};
  Vec<Protofilament *> targets;
  targets.reserve(protofilaments_.size());
  for (auto &&pf : protofilaments_) {
    tot_nucleation += pf.length_ * p_nucleate;
    targets.push_back(&pf);
  }
  int n_events = SysRNG::SamplePoisson(tot_nucleation);
  // printf("%i\n", n_events);
  for (int i_event{0}; i_event < n_events; i_event++) {
    double p_cum{0.0};
    double ran{SysRNG::GetRanProb()};
    for (int i_pf{0}; i_pf < targets.size(); i_pf++) {
      Protofilament *pf{targets[i_pf]};
      p_cum += pf->length_ * p_nucleate / tot_nucleation;
      if (ran < p_cum) {
        // shlould be handled by mgmt properly
        bool success{NucleateProtofilament(pf)};
        // bool success{targets[i_pf]->Nucleate()};
        if (success) {
          targets[i_pf] = targets.back();
          targets.pop_back();
          UpdateNeighbors();
        }
        break;
      }
    }
  }
}

void FilamentManager::UpdateForces() {

  for (auto &&pf : protofilaments_) {
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      // SF TODO FIX for nucleating microtubules
      // pf.force_[i_dim] = Params::Filaments::f_applied[i_dim];
      pf.force_[i_dim] = 0.0;
    }
    pf.torque_ = 0.0;
  }
  double F_factor{1.0e-8};
  double v0{67};               // nm/s
  double wall_location{-4500}; // nm
  double k_spring{2e-4};
  double r0{25.0};
  if (Params::Filaments::axon_arrangement) {
    for (auto &&pf : protofilaments_) {
      // printf("%zu neighbs\n", pf.neighbors_.size());
      for (auto &&neighb : pf.neighbors_) {
        // double dx{pf.plus_end_->pos_[0] - neighb->plus_end_->pos_[0]};
        double overlap_start{pf.sites_[0].pos_[0]};
        if (neighb->sites_[0].pos_[0] > overlap_start) {
          overlap_start = neighb->sites_[0].pos_[0];
        }
        double overlap_end{neighb->sites_.back().pos_[0]};
        if (pf.sites_.back().pos_[0] < overlap_end) {
          overlap_end = pf.sites_.back().pos_[0];
        }
        double O{overlap_end - overlap_start};
        if (O < 0.0) {
          O = 0.0;
        }
        double min_length{pf.length_ > neighb->length_ ? neighb->length_
                                                       : pf.length_};
        if (O > min_length) {
          O = min_length;
        }
        double dVel{pf.velocity_avg_[0] - neighb->velocity_avg_[0]};
        if (pf.polarity_ != neighb->polarity_) {
          // printf("%g\n", 1.0 - dVel / v0);
          // pf.force_[0] += pf.dx_ * (1.0 - dVel / v0) * O * F_factor;
          pf.force_[0] += pf.dx_ * O * F_factor; // * 1e3;

          // neighb->force_[0] -= pf.dx_ * O * F_factor;
          // printf("1\n");
        } else {
          pf.force_[0] += -dVel * O * F_factor;
          // printf("%g = %g * %g * %g\n", pf.force_[0], -dVel, O, F_factor);
          // pf.force_[0] += 1e-5;

          // pf.force_[0] += O * F_factor;
          // printf("2\n");
        }
      }
      // continue;
      // printf("+ %g <-> %g\n", pf.plus_end_->pos_[0],
      // pf.minus_end_->pos_[0]);
      if (pf.plus_end_->pos_[0] < pf.minus_end_->pos_[0]) {
        double r{pf.plus_end_->pos_[0] - wall_location};
        // printf("%g\n", r);
        if (r < r0) {
          double f_mag{-k_spring * (r - r0)};
          // double f_mag{48 * epsilon_ *
          //              (Pow(sigma_, 12) / Pow(r, 13) -
          //               0.5 * Pow(sigma_, 6) / Pow(r, 7))};
          pf.force_[0] += f_mag;
        }
      } else {
        double r{pf.minus_end_->pos_[0] - wall_location};
        if (r < r0) {
          double f_mag{-k_spring * (r - r0)};
          // double f_mag{48 * epsilon_ *
          //              (Pow(sigma_, 12) / Pow(r, 13) -
          //               0.5 * Pow(sigma_, 6) / Pow(r, 7))};
          pf.force_[0] += f_mag;
        }
      }
    }
  }
  /*
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
  */
  // for (auto &&pf : protofilaments_) {
  //   printf("F = <%g, %g> for PF #%i\n", pf.force_[0], pf.force_[1],
  //   pf.index_);
  // }
}

void FilamentManager::UpdateLattice() { proteins_->UpdateLatticeDeformation(); }

void FilamentManager::UpdateNeighbors() {

  double threshold{32};
  for (auto &&pf : protofilaments_) {
    pf.neighbors_.clear();
    for (auto &&neighb : protofilaments_) {
      if (fabs(pf.pos_[1] - neighb.pos_[1]) <= threshold) {
        pf.neighbors_.push_back(&neighb);
      }
    }
  }
  // for (auto &&pf : protofilaments_) {
  //   printf("neighb size is %zu\n", pf.neighbors_.size());
  // }
}

bool FilamentManager::NucleateProtofilament(Protofilament *parent) {

  if (protofilaments_.size() >= n_pfs_max_) {
    return false;
  }
  protofilaments_.emplace_back();
  size_t i_last{protofilaments_.size() - 1};
  bool success{protofilaments_.back().Nucleate(_id_site, Sys::n_objects_++,
                                               i_last, parent)};
  if (!success) {
    protofilaments_.pop_back();
    return false;
  }
  // protofilaments_[i_last - 1].top_neighb_ = &protofilaments_.back();
  // protofilaments_.back().bot_neighb_ = &protofilaments_[i_last - 1];
  // protofilaments_.back().top_neighb_ = nullptr;
  printf("added MT #%zu (t = %g)\n", i_last, Sys::i_step_ * Params::dt);
  return true;
}