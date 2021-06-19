#include "cylaks/protein_manager.hpp"
#include "cylaks/filament_manager.hpp"

void ProteinManager::FlagFilamentsForUpdate() { filaments_->FlagForUpdate(); }

void ProteinManager::UpdateFilaments() { filaments_->UpdateUnoccupied(); }

void ProteinManager::GenerateReservoirs() {

  // To model an infinite reservoir, we generate a number of proteins equal to
  // the total number of binding sites in the simulation. This value is static.
  size_t reservoir_size{0};
  for (int i_mt{0}; i_mt < Params::Filaments::count; i_mt++) {
    reservoir_size += Params::Filaments::n_sites[i_mt];
  }
  // Convert t_active parameter to number of simulation timesteps
  size_t motor_step_active{size_t(Params::Motors::t_active / Params::dt)};
  if (Params::Motors::c_bulk == 0.0) {
    motor_step_active = std::numeric_limits<size_t>::max();
  }
  // Initialize the motor reservoir
  motors_.Initialize(_id_motor, reservoir_size, motor_step_active);
  // Convert t_active parameter to number of simulation timesteps
  size_t xlink_step_active{size_t(Params::Xlinks::t_active / Params::dt)};
  if (Params::Xlinks::c_bulk == 0.0) {
    xlink_step_active = std::numeric_limits<size_t>::max();
  }
  // Initialize the crosslinker reservoir
  xlinks_.Initialize(_id_xlink, reservoir_size, xlink_step_active);
}

void ProteinManager::InitializeWeights() {

  /*
    For events that result in a change in energy dE, we use Boltzmann factors
    to scale rates appropriately. Detailed balance is satisfied by the
    factors: exp{-(1.0 - lambda) * dE / kBT}, and exp{lambda * dE / kBT} for
    forward and reverse pathways, respectively, (e.g., binding and unbinding),
    where lambda is a constant that ranges from 0 to 1, and kBT is thermal
    energy of the system. 3 values of lambda demonstrate its function: Lambda
    = 0.0 means all energy dependences is in binding Lambda = 0.5 means energy
    dependence is equal for binding and unbinding Lambda = 1.0 means all
    energy dependence is in unbinding
   */
  // Neighbor stuff
  Str name{"neighbs"};
  Sys::weight_neighb_bind_.resize(_n_neighbs_max + 1);
  Sys::weight_neighb_unbind_.resize(_n_neighbs_max + 1);
  motors_.AddWeight(name, _n_neighbs_max + 1);
  xlinks_.AddWeight(name, _n_neighbs_max + 1);
  double lambda_neighb{1.0};
  double dE{0.0};
  for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
    dE = -1 * Params::Motors::neighb_neighb_energy * n_neighbs;
    // ! FIXME
    motors_.weights_[name].bind_[n_neighbs] = exp(-(1.0 - _lambda_neighb) * dE);
    motors_.weights_[name].unbind_[n_neighbs] = exp(_lambda_neighb * dE);
    Sys::weight_neighb_bind_[n_neighbs] = exp(-(1.0 - _lambda_neighb) * dE);
    Sys::weight_neighb_unbind_[n_neighbs] = exp(_lambda_neighb * dE);
    // ! FIXME
    dE = -1 * Params::Xlinks::neighb_neighb_energy * n_neighbs;
    xlinks_.weights_[name].bind_[n_neighbs] = exp(-(1.0 - _lambda_neighb) * dE);
    xlinks_.weights_[name].unbind_[n_neighbs] = exp(_lambda_neighb * dE);
  }
  // Array of MAXIMUM possible weights including neighbor interactions and
  // lattice deformation effects from MULTIPLE motor effects stacking
  double lattice_E_max_tot{-1 * Params::Motors::gaussian_ceiling_bulk};
  double wt_lattice_bind_max{exp(-(1.0 - _lambda_lattice) * lattice_E_max_tot)};
  double wt_lattice_unbind_max{exp(_lambda_lattice * lattice_E_max_tot)};
  Sys::weight_lattice_bind_max_.resize(_n_neighbs_max + 1);
  Sys::weight_lattice_unbind_max_.resize(_n_neighbs_max + 1);
  for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
    Sys::weight_lattice_bind_max_[n_neighbs] =
        wt_lattice_bind_max * Sys::weight_neighb_bind_[n_neighbs];
    Sys::weight_lattice_unbind_max_[n_neighbs] =
        wt_lattice_unbind_max * Sys::weight_neighb_unbind_[n_neighbs];
  }
  // Calculate lattice_alpha based on input range
  Sys::lattice_cutoff_ = Params::Motors::gaussian_range;
  double dx_cutoff{(Sys::lattice_cutoff_ + 1) * Params::Filaments::site_size};
  double lattice_E_0_solo{-1 * Params::Motors::gaussian_amp_solo};
  double lattice_alpha{-1 * lattice_E_0_solo / (dx_cutoff * dx_cutoff)};
  if (motors_.active_) {
    if (Sys::lattice_cutoff_ > 0) {
      Sys::Log("  Motor variables calculated post-initialization:\n");
      Sys::Log("   lattice_cutoff_ is %i\n", Sys::lattice_cutoff_);
      Sys::Log("   lattice_alpha_ is %g\n", lattice_alpha);
    } else {
      Sys::Log("  Lattice cooperativity is disabled.\n");
    }
  }
  // Array of binding/unbinding weights due to lattice deformation from a SINGLE
  // motor; multiplied together to get total weight for any arrangement
  Sys::weight_lattice_bind_.resize(Sys::lattice_cutoff_ + 1);
  Sys::weight_lattice_unbind_.resize(Sys::lattice_cutoff_ + 1);
  for (int delta{0}; delta <= Sys::lattice_cutoff_; delta++) {
    double dx{delta * Params::Filaments::site_size};
    double energy{lattice_alpha * dx * dx + lattice_E_0_solo}; // in kbT
    Sys::weight_lattice_bind_[delta] = exp(-(1.0 - _lambda_lattice) * energy);
    Sys::weight_lattice_unbind_[delta] = exp(_lambda_lattice * energy);
  }
}

void ProteinManager::SetParameters() {

  using namespace Params;
  // Bind I
  xlinks_.AddProb("bind_i", Xlinks::k_on * Xlinks::c_bulk * dt, "neighbs", 0);
  motors_.AddProb("bind_i", Motors::k_on * Motors::c_bulk * dt);
  // Bind_I_Teth

  // Bind_ATP -- motors only
  double p_bind_ATP{Motors::k_on_ATP * Motors::c_ATP * dt};
  double wt_ATP_i{1.0};
  if (Motors::applied_force > 0.0) {
    wt_ATP_i = exp(-Motors::applied_force * Motors::sigma_ATP / kbT);
  }
  motors_.AddProb("bind_ATP_i", p_bind_ATP * wt_ATP_i);
  double wt_ATP_ii{exp(-Motors::internal_force * Motors::sigma_ATP / kbT)};
  if (Motors::internal_force == 0.0) {
    wt_ATP_ii = 0.0;
  }
  if (Motors::applied_force > 0.0) {
    wt_ATP_ii *= exp(Motors::applied_force * Motors::sigma_ATP / kbT);
  }
  motors_.AddProb("bind_ATP_ii", p_bind_ATP * wt_ATP_ii);
  // Hydrolyze_ATP -- motors only
  motors_.AddProb("hydrolyze", Motors::k_hydrolyze * dt);
  // Bind_II
  xlinks_.AddProb("bind_ii", Xlinks::k_on * Xlinks::c_eff_bind * dt);
  motors_.AddProb("bind_ii", Motors::k_on * Motors::c_eff_bind * dt);
  // Bind_II_Teth -- xlinks only

  // Unbind_II
  xlinks_.AddProb("unbind_ii", Xlinks::k_off_ii * dt);
  double wt_unbind_ii{exp(Motors::internal_force * Motors::sigma_off_ii / kbT)};
  if (Motors::applied_force > 0.0) {
    wt_unbind_ii *= exp(-Motors::applied_force * Motors::sigma_off_ii / kbT);
  }
  motors_.AddProb("unbind_ii", Motors::k_off_ii * dt * wt_unbind_ii);
  // Unbind_II_Teth -- xlinks only

  // Unbind I
  xlinks_.AddProb("unbind_i", Xlinks::k_off_i * dt, "neighbs", 1);
  double wt_unbind_i{1.0};
  if (Motors::applied_force > 0.0) {
    wt_unbind_i = exp(Motors::applied_force * Motors::sigma_off_i / kbT);
  }
  motors_.AddProb("unbind_i", Motors::k_off_i * dt * wt_unbind_i);
  // Unbind_I_Teth

  // Diffusion
  double x_sq{Square(Filaments::site_size / 1000)}; // in um^2
  double tau_i{x_sq / (2 * Xlinks::d_i)};
  // Sys::Log("tau = %g\n", tau_i);
  double tau_ii{x_sq / (2 * Xlinks::d_ii)};
  double p_diffuse_i{dt / tau_i};
  // Sys::Log("p = %g\n", p_diffuse_i);
  double p_diffuse_ii{dt / tau_ii};
  // diff_fwd and diff_bck are two separate events, which effectively
  // doubles the probability to diffuse. Thus we divide p_diff by 2.
  p_diffuse_i /= 2.0;
  p_diffuse_ii /= 2.0;
  Vec3D<double> weight_diff{{Vec<double>(_n_neighbs_max + 1, p_diffuse_i)}};
  weight_diff[0][0][1] *= xlinks_.weights_.at("neighbs").unbind_[1];
  weight_diff[0][0][2] *= 0.0;
  xlinks_.AddProb("diffuse_i_fwd", weight_diff);
  xlinks_.AddProb("diffuse_i_bck", weight_diff);
  xlinks_.AddProb("diffuse_ii_to_rest", p_diffuse_ii);
  xlinks_.AddProb("diffuse_ii_fr_rest", p_diffuse_ii);
  // Tether_Free

  // Tether_Bound -- motors only

  // Untether_Free

  // Untether_Bound -- motors only
}

void ProteinManager::InitializeEvents() {

  //  Binomial probabilitiy distribution; sampled to predict most events
  auto binomial = [&](double p, int n) {
    if (n > 0) {
      return SysRNG::SampleBinomial(p, n);
    } else {
      return 0;
    }
  };
  // Poisson distribution; sampled to predict events w/ variable probabilities
  auto poisson = [&](double p, int n) {
    if (p > 0.0) {
      return SysRNG::SamplePoisson(p);
    } else {
      return 0;
    }
  };
  // Bind_I: Bind first head of a protein
  auto exe_bind_i = [&](auto *site, auto *pop, auto *fil) {
    // Do not bind from solution until species is active
    if (Sys::i_step_ < pop->step_active_) {
      return;
    }
    auto entry{pop->GetFreeEntry()};
    if (entry == nullptr) {
      return;
    }
    bool executed{entry->Bind(site, &entry->head_one_)};
    if (executed) {
      pop->AddToActive(entry);
      fil->FlagForUpdate();
    }
  };
  auto weight_bind_i = [](auto *site) { return site->GetWeight_Bind(); };
  auto is_unocc = [](Object *site) -> Vec<Object *> {
    if (!site->IsOccupied()) {
      return {site};
    }
    return {};
  };
  Vec<int> i_min{0, 0, 0};
  Vec<size_t> dim_size{1, 1, _n_neighbs_max + 1};
  auto get_n_neighbs = [](Object *entry) {
    Vec<int> indices_vec{entry->GetNumNeighborsOccupied()};
    return indices_vec;
  };
  if (xlinks_.active_) {
    // Add unoccupied site tracker for crosslinkers; segregated by n_neighbs
    filaments_->AddPop("xlinks", is_unocc, dim_size, i_min, get_n_neighbs);
    for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
      kmc_.events_.emplace_back(
          "bind_i", xlinks_.p_event_.at("bind_i").GetVal(n_neighbs),
          &filaments_->unoccupied_.at("xlinks").bin_size_[0][0][n_neighbs],
          &filaments_->unoccupied_.at("xlinks").bin_entries_[0][0][n_neighbs],
          binomial, [&](Object *base) {
            exe_bind_i(dynamic_cast<BindingSite *>(base), &xlinks_, filaments_);
          });
    }
  }
  if (motors_.active_) {
    // Add unoccupied site tracker for motors; no binning b/c it's Poisson-based
    filaments_->AddPop("motors", is_unocc);
    kmc_.events_.emplace_back(
        "bind_i", motors_.p_event_.at("bind_i").GetVal(),
        &filaments_->unoccupied_.at("motors").size_,
        &filaments_->unoccupied_.at("motors").entries_, poisson,
        [&](Object *base) {
          return weight_bind_i(dynamic_cast<BindingSite *>(base));
        },
        [&](Object *base) {
          exe_bind_i(dynamic_cast<BindingSite *>(base), &motors_, filaments_);
        });
  }
  // Bind_I_Teth
  // ! Need to add
  // Bind_II
  auto exe_bind_ii = [](auto *bound_head, auto *pop, auto *fil) {
    auto head{bound_head->GetOtherHead()};
    auto site{head->parent_->GetNeighbor_Bind_II()};
    auto executed{head->parent_->Bind(site, head)};
    if (executed) {
      bool still_attached{head->parent_->UpdateExtension()};
      pop->FlagForUpdate();
      fil->FlagForUpdate();
    }
  };
  auto weight_bind_ii = [](auto *head) {
    return head->parent_->GetWeight_Bind_II();
  };
  auto is_singly_bound = [](Object *protein) -> Vec<Object *> {
    if (protein->GetNumHeadsActive() == 1) {
      return {protein->GetActiveHead()};
    }
    return {};
  };
  if (xlinks_.crosslinking_active_) {
    xlinks_.AddPop("bind_ii", is_singly_bound);
    kmc_.events_.emplace_back(
        "bind_ii", xlinks_.p_event_.at("bind_ii").GetVal(),
        &xlinks_.sorted_.at("bind_ii").size_,
        &xlinks_.sorted_.at("bind_ii").entries_, poisson,
        [&](Object *base) {
          return weight_bind_ii(dynamic_cast<BindingHead *>(base));
        },
        [&](Object *base) {
          exe_bind_ii(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_);
        });
  }
  if (motors_.active_) {
    auto is_docked = [](auto *motor) -> Vec<Object *> {
      auto *docked_head{motor->GetDockedHead()};
      if (docked_head != nullptr) {
        return {docked_head->GetOtherHead()};
      }
      return {};
    };
    motors_.AddPop("bind_ii", [&](Object *base) {
      return is_docked(dynamic_cast<Motor *>(base));
    });
    kmc_.events_.emplace_back(
        "bind_ii", motors_.p_event_.at("bind_ii").GetVal(),
        &motors_.sorted_.at("bind_ii").size_,
        &motors_.sorted_.at("bind_ii").entries_, poisson,
        [&](Object *base) {
          return weight_bind_ii(dynamic_cast<CatalyticHead *>(base));
        },
        [&](Object *base) {
          exe_bind_ii(dynamic_cast<CatalyticHead *>(base), &motors_,
                      filaments_);
        });
  }
  // Unbind_II
  auto exe_unbind_ii = [](auto *head, auto *pop, auto *fil) {
    bool executed{head->Unbind()};
    if (executed) {
      pop->FlagForUpdate();
      fil->FlagForUpdate();
    }
  };
  auto weight_unbind_ii = [](auto *head) {
    return head->GetWeight_Unbind_II();
  };
  auto is_doubly_bound = [](Object *protein) -> Vec<Object *> {
    if (protein->GetNumHeadsActive() == 2) {
      return {protein->GetHeadOne(), protein->GetHeadTwo()};
    }
    return {};
  };
  if (xlinks_.crosslinking_active_) {
    xlinks_.AddPop("unbind_ii", is_doubly_bound);
    kmc_.events_.emplace_back(
        "unbind_ii", xlinks_.p_event_.at("unbind_ii").GetVal(),
        &xlinks_.sorted_.at("unbind_ii").size_,
        &xlinks_.sorted_.at("unbind_ii").entries_, poisson,
        [&](Object *base) {
          return weight_unbind_ii(dynamic_cast<BindingHead *>(base));
        },
        [&](Object *base) {
          exe_unbind_ii(dynamic_cast<BindingHead *>(base), &xlinks_,
                        filaments_);
        });
  }
  if (motors_.active_) {
    auto is_ADPP_ii_bound = [](auto *motor) -> Vec<Object *> {
      if (motor->n_heads_active_ == 2) {
        bool found_head{false};
        CatalyticHead *chosen_head{nullptr};
        if (motor->head_one_.ligand_ == CatalyticHead::Ligand::ADPP) {
          chosen_head = &motor->head_one_;
        }
        if (motor->head_two_.ligand_ == CatalyticHead::Ligand::ADPP) {
          if (chosen_head != nullptr) {
            Sys::ErrorExit("Protein_MGR::is_ADPP_ii_bound()");
          }
          chosen_head = &motor->head_two_;
        }
        if (chosen_head != nullptr) {
          return {chosen_head};
        }
      }
      return {};
    };
    motors_.AddPop("unbind_ii", [&](Object *base) {
      return is_ADPP_ii_bound(dynamic_cast<Motor *>(base));
    });
    kmc_.events_.emplace_back(
        "unbind_ii", motors_.p_event_.at("unbind_ii").GetVal(),
        &motors_.sorted_.at("unbind_ii").size_,
        &motors_.sorted_.at("unbind_ii").entries_, poisson,
        [&](Object *base) {
          return weight_unbind_ii(dynamic_cast<CatalyticHead *>(base));
        },
        [&](Object *base) {
          exe_unbind_ii(dynamic_cast<CatalyticHead *>(base), &motors_,
                        filaments_);
        });
  }
  // Unbind_I: Unbind first (singly bound) head of a protein
  auto exe_unbind_i = [](auto *head, auto *pop, auto *fil) {
    bool executed{head->Unbind()};
    if (executed) {
      head->UntetherSatellite();
      pop->RemoveFromActive(head->parent_);
      fil->FlagForUpdate();
    }
  };
  if (xlinks_.active_) {
    xlinks_.AddPop("bound_i", is_singly_bound, dim_size, i_min, get_n_neighbs);
    for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
      kmc_.events_.emplace_back(
          "unbind_i", xlinks_.p_event_.at("unbind_i").GetVal(n_neighbs),
          &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
          &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs],
          binomial, [&](Object *base) {
            exe_unbind_i(dynamic_cast<BindingHead *>(base), &xlinks_,
                         filaments_);
          });
    }
  }
  if (motors_.active_) {
    auto weight_unbind_i = [](auto *head) {
      return head->parent_->GetWeight_Unbind_I();
    };
    auto is_ADPP_i_bound = [](auto *motor) -> Vec<Object *> {
      if (motor->n_heads_active_ == 1) {
        if (motor->GetActiveHead()->ligand_ == CatalyticHead::Ligand::ADPP) {
          return {motor->GetActiveHead()};
        }
      }
      return {};
    };
    motors_.AddPop("bound_i_ADPP", [&](Object *base) {
      return is_ADPP_i_bound(dynamic_cast<Motor *>(base));
    });
    kmc_.events_.emplace_back(
        "unbind_i", motors_.p_event_.at("unbind_i").GetVal(),
        &motors_.sorted_.at("bound_i_ADPP").size_,
        &motors_.sorted_.at("bound_i_ADPP").entries_, poisson,
        [&](Object *base) {
          return weight_unbind_i(dynamic_cast<CatalyticHead *>(base));
        },
        [&](Object *base) {
          exe_unbind_i(dynamic_cast<CatalyticHead *>(base), &motors_,
                       filaments_);
        });
  }
  // Unbind_I_Teth
  // Tether_Free
  // Untether_Free

  // vv motors only vv
  // Bind_ATP
  if (motors_.active_) {
    auto exe_bind_ATP = [](auto *head, auto *pop) {
      // printf("boop\n");
      bool executed{head->parent_->Bind_ATP(head)};
      if (executed) {
        pop->FlagForUpdate();
      }
    };
    auto is_NULL_i_bound = [](auto *motor) -> Vec<Object *> {
      if (motor->n_heads_active_ == 1) {
        if (motor->GetActiveHead()->ligand_ == CatalyticHead::Ligand::NONE) {
          return {motor->GetActiveHead()};
        }
      }
      return {};
    };
    motors_.AddPop("bound_i_NULL", [&](Object *base) {
      return is_NULL_i_bound(dynamic_cast<Motor *>(base));
    });
    kmc_.events_.emplace_back(
        "bind_ATP_i", motors_.p_event_.at("bind_ATP_i").GetVal(),
        &motors_.sorted_.at("bound_i_NULL").size_,
        &motors_.sorted_.at("bound_i_NULL").entries_, binomial,
        [&](Object *base) {
          exe_bind_ATP(dynamic_cast<CatalyticHead *>(base), &motors_);
        });
    auto exe_bind_ATP_ii = [](auto *front_head, auto *pop, auto *fil) {
      auto *rear_head{front_head->GetOtherHead()};
      if (front_head->trailing_) {
        return;
      }
      bool unbound{rear_head->Unbind()};
      bool executed{front_head->parent_->Bind_ATP(front_head)};
      if (executed) {
        pop->FlagForUpdate();
        fil->FlagForUpdate();
      }
    };
    auto weight_bind_ATP_ii = [](auto *head) {
      return head->parent_->GetWeight_BindATP_II(head);
    };
    auto is_NULL_ii_bound = [](auto *motor) -> Vec<Object *> {
      if (motor->n_heads_active_ == 2) {
        bool found_head{false};
        CatalyticHead *chosen_head{nullptr};
        if (motor->head_one_.ligand_ == CatalyticHead::Ligand::NONE) {
          chosen_head = &motor->head_one_;
        }
        if (motor->head_two_.ligand_ == CatalyticHead::Ligand::NONE) {
          if (chosen_head != nullptr) {
            Sys::ErrorExit("Protein_MGR::is_NULL_ii_bound()");
          }
          chosen_head = &motor->head_two_;
        }
        if (chosen_head != nullptr) {
          return {chosen_head};
        }
      }
      return {};
    };
    motors_.AddPop("bound_ii_NULL", [&](Object *base) {
      return is_NULL_ii_bound(dynamic_cast<Motor *>(base));
    });
    kmc_.events_.emplace_back(
        "bind_ATP_ii", motors_.p_event_.at("bind_ATP_ii").GetVal(),
        &motors_.sorted_.at("bound_ii_NULL").size_,
        &motors_.sorted_.at("bound_ii_NULL").entries_, poisson,
        [&](Object *base) {
          return weight_bind_ATP_ii(dynamic_cast<CatalyticHead *>(base));
        },
        [&](Object *base) {
          exe_bind_ATP_ii(dynamic_cast<CatalyticHead *>(base), &motors_,
                          filaments_);
        });
  }
  // Hydrolyze_ATP
  if (motors_.active_) {
    auto exe_hydrolyze = [](auto *head, auto *pop) {
      bool executed{head->parent_->Hydrolyze(head)};
      if (executed) {
        pop->FlagForUpdate();
      }
    };
    auto is_ATP_i_bound = [](auto *motor) -> Vec<Object *> {
      if (motor->n_heads_active_ == 1) {
        if (motor->GetActiveHead()->ligand_ == CatalyticHead::Ligand::ATP) {
          return {motor->GetActiveHead()};
        }
      }
      return {};
    };
    motors_.AddPop("bound_i_ATP", [&](Object *base) {
      return is_ATP_i_bound(dynamic_cast<Motor *>(base));
    });
    kmc_.events_.emplace_back(
        "hydrolyze", motors_.p_event_.at("hydrolyze").GetVal(),
        &motors_.sorted_.at("bound_i_ATP").size_,
        &motors_.sorted_.at("bound_i_ATP").entries_, binomial,
        [&](Object *base) {
          exe_hydrolyze(dynamic_cast<CatalyticHead *>(base), &motors_);
        });
  }
  // Tether_Bound
  // Untether_Bound

  // vv xlinks only vv
  // Diffusion
  auto exe_diff = [](auto *head, auto *pop, auto *fil, int dir) {
    bool executed{head->Diffuse(dir)};
    if (executed) {
      bool still_attached{head->parent_->UpdateExtension()};
      if (!still_attached) {
        // printf("what\n");
      }
      pop->FlagForUpdate();
      fil->FlagForUpdate();
    }
  };
  if (xlinks_.active_) {
    for (int n_neighbs{0}; n_neighbs < _n_neighbs_max; n_neighbs++) {
      kmc_.events_.emplace_back(
          "diffuse_i_fwd",
          xlinks_.p_event_.at("diffuse_i_fwd").GetVal(n_neighbs),
          &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
          &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs],
          binomial, [&](Object *base) {
            exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_,
                     1);
          });
      kmc_.events_.emplace_back(
          "diffuse_i_bck",
          xlinks_.p_event_.at("diffuse_i_bck").GetVal(n_neighbs),
          &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
          &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs],
          binomial, [&](Object *base) {
            exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_,
                     -1);
          });
    }
  }
  if (xlinks_.crosslinking_active_) {
    auto weight_diff_ii = [](auto *head, int dir) {
      return head->GetWeight_Diffuse(dir);
    };
    xlinks_.AddPop("diffuse_ii_to_rest", is_doubly_bound);
    xlinks_.AddPop("diffuse_ii_fr_rest", is_doubly_bound);
    kmc_.events_.emplace_back(
        "diffuse_ii_to_rest",
        xlinks_.p_event_.at("diffuse_ii_to_rest").GetVal(),
        &xlinks_.sorted_.at("diffuse_ii_to_rest").size_,
        &xlinks_.sorted_.at("diffuse_ii_to_rest").entries_, poisson,
        [&](Object *base) {
          return weight_diff_ii(dynamic_cast<BindingHead *>(base), 1);
        },
        [&](Object *base) {
          exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_, 1);
        });
    kmc_.events_.emplace_back(
        "diffuse_ii_fr_rest",
        xlinks_.p_event_.at("diffuse_ii_fr_rest").GetVal(),
        &xlinks_.sorted_.at("diffuse_ii_fr_rest").size_,
        &xlinks_.sorted_.at("diffuse_ii_fr_rest").entries_, poisson,
        [&](Object *base) {
          return weight_diff_ii(dynamic_cast<BindingHead *>(base), -1);
        },
        [&](Object *base) {
          exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_, -1);
        });
  }
  // Bind_II_Teth
  // ! Need to add
  // Unbind_II_Teth
  // ! Need to add
  // ! FIXME export this to a function that is also called in test modes
  // Check that all event probabilities are valid
  for (auto const &event : kmc_.events_) {
    if (event.p_occur_ < 0.0 or event.p_occur_ > 1.0) {
      Sys::Log("Invalid probabilitiy: p_%s = %g\n", event.name_.c_str(),
               event.p_occur_);
      Sys::ErrorExit("ProteinManager::InitializeEvents()");
    } else if (event.p_occur_ > 0.1) {
      Sys::Log("WARNING: p_%s = %g. High probability may reduce accuracy\n",
               event.name_.c_str(), event.p_occur_);
    }
  }
}