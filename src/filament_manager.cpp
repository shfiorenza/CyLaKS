#include "cylaks/filament_manager.hpp"
#include "cylaks/protein_manager.hpp"
#include "cylaks/protofilament.hpp"
#include "cylaks/system_definitions.hpp"
#include "cylaks/system_namespace.hpp"
#include "cylaks/system_parameters.hpp"
#include "cylaks/system_rng.hpp"
#include <cmath>

void FilamentManager::SetParameters() {

  threshold_ = std::pow(2, 1.0 / 6.0) * sigma_;
  n_bd_iterations_ = Params::Filaments::n_bd_per_kmc;
  dt_eff_ = Params::dt / n_bd_iterations_;
}

void FilamentManager::GenerateFilaments() {

  // Create the appropriate number of protofilaments
  if (Params::Filaments::axon_arrangement == true) {
    protofilaments_.resize(Params::Filaments::count);
    protofilaments_.reserve(n_pfs_max_);
    size_t block_size{Params::Filaments::Neuron::block_size};
    size_t n_blocks{Params::Filaments::count / block_size};
    double x_offset{Params::Filaments::Neuron::x_offset};
    double y_offset{Params::Filaments::Neuron::y_offset};
    for (size_t i_block{0}; i_block < n_blocks; i_block++) {
      for (size_t i_fil{i_block * block_size};
           i_fil < (i_block + 1) * block_size; i_fil++) {
        if (i_block > 0) {
          Params::Filaments::x_initial[i_fil] = x_offset * i_block;
        }
        // Params::Filaments::y_initial[i_fil] =
        //     Params::Filaments::y_initial[i_fil - i_block * block_size];
        Params::Filaments::y_initial[i_fil] +=
            (i_fil - i_block * block_size) * y_offset;
        // if (i_block % 2 == 1) {
        //   Params::Filaments::y_initial[i_fil] += 0.5 * y_offset;
        // }
        Params::Filaments::y_initial[i_fil] +=
            (SysRNG::GetRanProb() - 0.5) * y_offset;
      }
    }
    size_t n_pfs_tot = protofilaments_.size();
    size_t n_plus_out =
        std::round(n_pfs_tot * Params::Filaments::Neuron::p_plus);
    size_t indices[n_pfs_tot];
    for (int i_fil{0}; i_fil < n_pfs_tot; i_fil++) {
      indices[i_fil] = i_fil;
    }
    SysRNG::Shuffle(indices, n_pfs_tot, sizeof(size_t));
    for (size_t i_fil{0}; i_fil < n_plus_out; i_fil++) {
      Params::Filaments::polarity[indices[i_fil]] = 0;
    }
    for (size_t i_fil{n_plus_out}; i_fil < n_pfs_tot; i_fil++) {
      Params::Filaments::polarity[indices[i_fil]] = 1;
    }
    /*
    for (int i_fil{0}; i_fil < protofilaments_.size(); i_fil++) {
      if (SysRNG::GetRanProb() < Params::Filaments::Neuron::p_plus) {
        Params::Filaments::polarity[i_fil] = 0;
      } else {
        Params::Filaments::polarity[i_fil] = 1;
      }
    }
    */
    for (int i_fil{0}; i_fil < protofilaments_.size(); i_fil++) {
      protofilaments_[i_fil].Initialize(_id_site, Sys::n_objects_++, i_fil);
    }
    UpdateNeighbors();
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
  // grow means unstable tip is growing
  // shrink means unstable tip is shrinking
  // pause means entire MT is stable; no labile tip
  // - grow->shrink: catastrophe
  // - shrink->grow: rescue
  // - grow->pause: stabilize
  // - pause->grow: growth
  double k_cata{Params::Filaments::Neuron::k_cata};
  double k_resc{Params::Filaments::Neuron::k_resc};
  double k_stab{Params::Filaments::Neuron::k_stab};
  double k_grow{Params::Filaments::Neuron::k_grow};
  double p_add_site = Params::Filaments::Neuron::v_grow * Params::dt /
                      Params::Filaments::site_size;
  double p_rmv_site = Params::Filaments::Neuron::v_shrink * Params::dt /
                      Params::Filaments::site_size;
  if (p_add_site >= 1.0 || p_rmv_site >= 1.0) {
    Sys::ErrorExit("filament problems. decrease timestep.\n");
  }
  // Dynamic instability
  for (auto &&pf : protofilaments_) {
    // printf("pf #%zu: %zu stable, %zu labile, %zu tot\n", pf.index_,
    //        pf.n_sites_stable_, pf.n_sites_labile_, pf.n_sites_);
    // if (pf.n_sites_stable_ > 100000) {
    //   exit(1);
    // }
    if (Params::Filaments::Neuron::soma_depoly) {
      if (pf.plus_end_->pos_[0] > Params::Filaments::Neuron::soma_pos) {
        pf.RemoveSite_PlusEnd();
        continue;
      }
      if (pf.minus_end_->pos_[0] > Params::Filaments::Neuron::soma_pos) {
        // pf.RemoveSite_PlusEnd();
        pf.RemoveSite_MinusEnd();
        continue;
      }
    }
    double ran{SysRNG::GetRanProb()};
    switch (pf.state_) {
    case pause: {
      // if (ran < Params::Filaments::Neuron::p_p2g) {
      if (ran < k_grow * Params::dt) {
        pf.state_ = grow;
        break;
      }
      break;
    }
    case grow: {
      if (ran < k_cata * Params::dt) {
        pf.state_ = shrink;
        break;
      } else if (ran < k_cata * Params::dt + k_stab * Params::dt) {
        pf.Stabilize();
        pf.state_ = pause;
        break;
      }
      double ran2{SysRNG::GetRanProb()};
      if (ran2 < p_add_site) {
        pf.AddSite_PlusEnd();
      }
      break;
    }
    case shrink: {
      if (ran < k_resc * Params::dt) {
        pf.state_ = grow;
        break;
      } else if (ran < k_resc * Params::dt + k_stab * Params::dt) {
        pf.Stabilize();
        pf.state_ = pause;
        break;
      }
      double ran2{SysRNG::GetRanProb()};
      if (ran2 < p_rmv_site) {
        pf.RemoveSite_PlusEnd();
      }
      break;
    }
    }
  }
  // SF TODO optimize with dynamic flagging
  for (auto &&pf : protofilaments_) {
    if (pf.n_sites_ == 2) {
      size_t i_pf{pf.index_};
      protofilaments_.erase(protofilaments_.begin() + i_pf);
      for (int i_entry{0}; i_entry < protofilaments_.size(); i_entry++) {
        protofilaments_[i_entry].index_ = i_entry;
      }
      UpdateNeighbors();
    }
  }
  bool MTs_added{false};
  // New MTs arriving from soma
  double p_spawn_soma{Params::Filaments::Neuron::k_spawn_soma * Params::dt};
  int n_spawn_soma = SysRNG::SamplePoisson(p_spawn_soma);
  for (int i_event{0}; i_event < n_spawn_soma; i_event++) {
    bool success{NucleateProtofilamentAtSoma()};
    if (success) {
      // printf("MT arrived from soma !\n");
      MTs_added = true;
    }
  }
  // New MTs spawning in cytoplasm
  double p_spawn_cyto{Params::Filaments::Neuron::k_spawn_cyto * Params::dt};
  p_spawn_cyto *= (Params::Filaments::Neuron::soma_pos -
                   Params::Filaments::Neuron::tip_pos);
  int n_spawn_cyto = SysRNG::SamplePoisson(p_spawn_cyto);
  for (int i_event{0}; i_event < n_spawn_cyto; i_event++) {
    bool success{NucleateProtofilamentInCyto()};
    if (success) {
      MTs_added = true;
    }
  }
  // New MTs nucleating from pre-existing MTs
  double p_nucleate{Params::Filaments::Neuron::k_nucleate *
                    Params::dt}; // per nm
  size_t n_max{500};
  p_nucleate *= (1.0 - double(protofilaments_.size()) / double(n_max));
  double tot_nucleation{0.0};
  Vec<Protofilament *> targets;
  targets.reserve(protofilaments_.size());
  for (auto &&pf : protofilaments_) {
    tot_nucleation += pf.length_ * p_nucleate;
    targets.push_back(&pf);
  }
  int n_events = SysRNG::SamplePoisson(tot_nucleation);
  for (int i_event{0}; i_event < n_events; i_event++) {
    double p_cum{0.0};
    double ran{SysRNG::GetRanProb()};
    for (int i_pf{0}; i_pf < targets.size(); i_pf++) {
      Protofilament *pf{targets[i_pf]};
      p_cum += pf->length_ * p_nucleate / tot_nucleation;
      if (ran < p_cum) {
        bool success{NucleateProtofilament(pf)};
        if (success) {
          targets[i_pf] = targets.back();
          targets.pop_back();
          MTs_added = true;
        }
        break;
      }
    }
  }
  size_t n_stable{0};
  Vec<Protofilament *> stable_pfs;
  stable_pfs.reserve(protofilaments_.size());
  for (auto &&pf : protofilaments_) {
    if (pf.n_sites_stable_ == pf.n_sites_) {
      stable_pfs.push_back(&pf);
      n_stable++;
    }
  }
  double avg_loss{Params::Filaments::Neuron::k_loss * n_stable * Params::dt};
  int n_lost = SysRNG::SamplePoisson(avg_loss);
  int indices[n_lost];
  SysRNG::SetRanIndices(indices, n_lost, n_stable);
  for (int i_event{0}; i_event < n_lost; i_event++) {
    Protofilament *pf = stable_pfs[indices[i_event]];
    int i_pf = pf->index_;
    protofilaments_.erase(protofilaments_.begin() + i_pf);
    for (int i_entry{0}; i_entry < protofilaments_.size(); i_entry++) {
      protofilaments_[i_entry].index_ = i_entry;
    }
    MTs_added = true;
  }
  if (MTs_added) {
    UpdateNeighbors();
  }
}

void FilamentManager::UpdateForces() {

  for (auto &&pf : protofilaments_) {
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      pf.force_[i_dim] = Params::Filaments::f_applied[i_dim];
    }
    pf.torque_ = 0.0;
  }
  double F_factor_slide{Params::Filaments::Neuron::F_factor_slide};
  double F_factor_para{Params::Filaments::Neuron::F_factor_para};
  double tip_pos{Params::Filaments::Neuron::tip_pos}; // nm
  double k_spring{Params::Filaments::Neuron::tip_k};
  double r0{Params::Filaments::Neuron::tip_r0};
  // double r0{pow(2.0, 1.0 / 6.0) * sigma_};
  if (Params::Filaments::axon_arrangement) {
    for (auto &&pf : protofilaments_) {
      for (auto &&neighb : pf.neighbors_) {
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
          // pf.force_[0] += pf.dx_ * (1.0 - dVel / v0) * O * F_factor;
          pf.force_[0] += pf.dx_ * O * F_factor_slide;
        } else {
          pf.force_[0] += -dVel * O * F_factor_para;
        }
      }
      if (pf.plus_end_->pos_[0] < pf.minus_end_->pos_[0]) {
        double r{pf.plus_end_->pos_[0] - tip_pos};
        if (r < r0) {
          double f_mag{-k_spring * (r - r0)};
          // double f_mag{48 * epsilon_ *
          //              (Pow(sigma_, 12) / Pow(r, 13) -
          //               0.5 * Pow(sigma_, 6) / Pow(r, 7))};
          pf.force_[0] += f_mag;
          // printf("PLUS: %g\n", f_mag);
          // pf.f_barrier_ = f_mag;
        }
      } else {
        double r{pf.minus_end_->pos_[0] - tip_pos};
        if (r < r0) {
          double f_mag{-k_spring * (r - r0)};
          // double f_mag{48 * epsilon_ *
          //              (Pow(sigma_, 12) / Pow(r, 13) -
          //               0.5 * Pow(sigma_, 6) / Pow(r, 7))};
          pf.force_[0] += f_mag;
          // printf("MINUS: %g\n", f_mag);
          // pf.f_barrier_ = f_mag;
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
  */
  proteins_->UpdateExtensions();
}

void FilamentManager::UpdateLattice() { proteins_->UpdateLatticeDeformation(); }

void FilamentManager::UpdateNeighbors() {

  double threshold{Params::Filaments::Neuron::neighb_threshold};
  for (auto &&pf : protofilaments_) {
    pf.neighbors_.clear();
    for (auto &&neighb : protofilaments_) {
      if (fabs(pf.pos_[1] - neighb.pos_[1]) <= threshold) {
        pf.neighbors_.push_back(&neighb);
      }
    }
  }
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
  Sys::Log("Added MT #%zu (t = %g)\n", i_last, Sys::i_step_ * Params::dt);
  return true;
}

bool FilamentManager::NucleateProtofilamentAtSoma() {

  if (protofilaments_.size() >= n_pfs_max_) {
    return false;
  }
  protofilaments_.emplace_back();
  size_t i_last{protofilaments_.size() - 1};
  bool success{
      protofilaments_.back().Nucleate(_id_site, Sys::n_objects_++, i_last, 1)};
  if (!success) {
    protofilaments_.pop_back();
    return false;
  }
  Sys::Log("Added MT (SOMA) #%zu (t = %g)\n", i_last,
           Sys::i_step_ * Params::dt);
  return true;
}

bool FilamentManager::NucleateProtofilamentInCyto() {

  if (protofilaments_.size() >= n_pfs_max_) {
    return false;
  }
  protofilaments_.emplace_back();
  size_t i_last{protofilaments_.size() - 1};
  bool success{
      protofilaments_.back().Nucleate(_id_site, Sys::n_objects_++, i_last, 2)};
  if (!success) {
    protofilaments_.pop_back();
    return false;
  }
  Sys::Log("Added MT (CYTO) #%zu (t = %g)\n", i_last,
           Sys::i_step_ * Params::dt);
  return true;
}