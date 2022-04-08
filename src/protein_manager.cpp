#include "cylaks/protein_manager.hpp"
#include "cylaks/filament_manager.hpp"

void ProteinManager::FlagFilamentsForUpdate() { filaments_->FlagForUpdate(); }

void ProteinManager::UpdateFilaments() { filaments_->UpdateUnoccupied(); }

void ProteinManager::GenerateReservoirs() {

  // ! FIXME hella temp
  double E_max{std::log(_max_weight) * Params::kbT};
  double teth_r_min{Params::Motors::r_0 -
                    sqrt(2 * E_max / Params::Motors::k_slack)};
  double teth_r_max{Params::Motors::r_0 +
                    sqrt(2 * E_max / Params::Motors::k_spring)};
  Sys::teth_x_min_ = (int)std::floor(teth_r_min / Params::Filaments::site_size);
  Sys::teth_x_max_ = (int)std::ceil(teth_r_max / Params::Filaments::site_size);
  printf("%g < %g < %g\n", teth_r_min, Params::Motors::r_0, teth_r_max);
  printf("x_min: %i\n", Sys::teth_x_min_);
  printf("x_max: %i\n", Sys::teth_x_max_);

  // To model an infinite reservoir, we generate a number of proteins equal to
  // the total number of binding sites in the simulation. This value is static.
  size_t reservoir_size{0};
  for (int i_mt{0}; i_mt < Params::Filaments::count; i_mt++) {
    reservoir_size += Params::Filaments::n_sites[i_mt];
  }
  // Convert t_active parameter from seconds to number of simulation timesteps
  size_t motor_step_active{size_t(Params::Motors::t_active / Params::dt)};
  if (Params::Motors::c_bulk == 0.0) {
    motor_step_active = std::numeric_limits<size_t>::max();
  }
  // Initialize the motor reservoir
  motors_.Initialize(_id_motor, reservoir_size, motor_step_active);
  // Convert t_active parameter from seconds to number of simulation timesteps
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
    energy of the system. 3 values of lambda demonstrate to its function:
        Lambda = 0.0 means all energy dependences is in binding
        Lambda = 0.5 means energy dependence is equal for binding and unbinding
        Lambda = 1.0 means all energy dependence is in unbinding
   */
  // Neighbor stuff
  Str name{"neighbs"};
  Sys::weight_neighb_bind_.resize(_n_neighbs_max + 1);
  Sys::weight_neighb_unbind_.resize(_n_neighbs_max + 1);
  motors_.AddWeight("neighbs", _n_neighbs_max + 1);
  xlinks_.AddWeight("neighbs", _n_neighbs_max + 1);
  double lambda_neighb{1.0};
  double dE{0.0};
  for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
    dE = -1 * Params::Motors::neighb_neighb_energy * n_neighbs;
    // ! FIXME generalize w/o use of Sys namespace
    motors_.weights_[name].bind_[n_neighbs] = exp(-(1.0 - _lambda_neighb) * dE);
    motors_.weights_[name].unbind_[n_neighbs] = exp(_lambda_neighb * dE);
    Sys::weight_neighb_bind_[n_neighbs] = exp(-(1.0 - _lambda_neighb) * dE);
    Sys::weight_neighb_unbind_[n_neighbs] = exp(_lambda_neighb * dE);
    // ! FIXME generalize w/o use of Sys namespace
    dE = -1 * Params::Xlinks::neighb_neighb_energy * n_neighbs;
    xlinks_.weights_[name].bind_[n_neighbs] = exp(-(1.0 - _lambda_neighb) * dE);
    xlinks_.weights_[name].unbind_[n_neighbs] = exp(_lambda_neighb * dE);
  }
  // Array of MAXIMUM possible weights including neighbor interactions and
  // lattice deformation effects from MULTIPLE motor effects stacking
  double lattice_E_max_tot{-1 * Params::Motors::gaussian_ceiling_bulk};
  double wt_lattice_bind_max{exp(-(1.0 - _lambda_lattice) * lattice_E_max_tot)};
  double wt_lattice_unbind_max{exp(_lambda_lattice * lattice_E_max_tot)};
  // ! FIXME generalize w/o use of Sys namespace
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
  // ! FIXME generalize w/o use of Sys namespace
  Sys::weight_lattice_bind_.resize(Sys::lattice_cutoff_ + 1);
  Sys::weight_lattice_unbind_.resize(Sys::lattice_cutoff_ + 1);
  for (int delta{0}; delta <= Sys::lattice_cutoff_; delta++) {
    double dx{delta * Params::Filaments::site_size};
    double energy{lattice_alpha * dx * dx + lattice_E_0_solo}; // in kbT
    // printf("E = %g\n", energy);
    Sys::weight_lattice_bind_[delta] = exp(-(1.0 - _lambda_lattice) * energy);
    Sys::weight_lattice_unbind_[delta] = exp(_lambda_lattice * energy);
  }
}

void ProteinManager::SetParameters() {

  using namespace Params;
  // Bind I
  motors_.AddProb("bind_i", Motors::k_on * Motors::c_bulk * dt);
  xlinks_.AddProb("bind_i", Xlinks::k_on * Xlinks::c_bulk * dt, "neighbs", 0);
  // Bind_I_Teth
  motors_.AddProb("bind_i_teth", Motors::k_on * Motors::c_eff_tether * dt);
  xlinks_.AddProb("bind_i_teth", Xlinks::k_on * Motors::c_eff_tether * dt);
  // Bind_ATP_I -- motors only
  double p_bind_ATP{Motors::k_on_ATP * Motors::c_ATP * dt};
  double wt_ATP_i{1.0};
  if (Motors::applied_force > 0.0) {
    wt_ATP_i = exp(-Motors::applied_force * Motors::sigma_ATP / kbT);
  }
  motors_.AddProb("bind_ATP_i", p_bind_ATP * wt_ATP_i);
  motors_.AddProb("bind_ATP_i_teth", p_bind_ATP * wt_ATP_i);
  // Bind_ATP_II -- motors only
  double wt_ATP_ii{exp(-Motors::internal_force * Motors::sigma_ATP / kbT)};
  if (Motors::internal_force == 0.0 or Motors::gaussian_range == 0 or
      !Motors::gaussian_stepping_coop) {
    wt_ATP_ii = 0.0;
  }
  if (Motors::applied_force > 0.0) {
    wt_ATP_ii *= exp(Motors::applied_force * Motors::sigma_ATP / kbT);
  }
  motors_.AddProb("bind_ATP_ii", p_bind_ATP * wt_ATP_ii);
  // Hydrolyze_ATP -- motors only
  motors_.AddProb("hydrolyze", Motors::k_hydrolyze * dt);
  // Bind_II
  motors_.AddProb("bind_ii", Motors::k_on * Motors::c_eff_bind * dt);
  xlinks_.AddProb("bind_ii", Xlinks::k_on * Xlinks::c_eff_bind * dt);
  // Bind_II_Teth -- xlinks only
  xlinks_.AddProb("bind_ii_teth", Xlinks::k_on * Xlinks::c_eff_bind * dt);
  // Unbind_II
  double wt_unbind_ii{exp(Motors::internal_force * Motors::sigma_off_ii / kbT)};
  if (Motors::applied_force > 0.0) {
    wt_unbind_ii *= exp(-Motors::applied_force * Motors::sigma_off_ii / kbT);
  }
  motors_.AddProb("unbind_ii", Motors::k_off_ii * dt * wt_unbind_ii);
  xlinks_.AddProb("unbind_ii", Xlinks::k_off_ii * dt);
  // Unbind_II_Teth -- xlinks only
  xlinks_.AddProb("unbind_ii_teth", Xlinks::k_off_ii * dt);
  // Unbind I
  double wt_unbind_i{1.0};
  if (Motors::applied_force > 0.0) {
    wt_unbind_i = exp(Motors::applied_force * Motors::sigma_off_i / kbT);
  }
  motors_.AddProb("unbind_i", Motors::k_off_i * dt * wt_unbind_i);
  xlinks_.AddProb("unbind_i", Xlinks::k_off_i * dt, "neighbs", 1);
  // Unbind_I_Teth
  motors_.AddProb("unbind_i_teth", Motors::k_off_i * dt * wt_unbind_i);
  xlinks_.AddProb("unbind_i_teth", Xlinks::k_off_i * dt);
  // Diffusion -- xlinks only
  double x_sq{Square(Filaments::site_size / 1000)}; // in um^2
  double tau_i{x_sq / (2 * Xlinks::d_i)};
  double tau_ii{x_sq / (2 * Xlinks::d_ii)};
  double p_diffuse_i{dt / tau_i};
  double p_diffuse_ii{dt / tau_ii};
  // diff_fwd and diff_bck are two separate events, which effectively
  // doubles the probability to diffuse. Thus we divide p_diff by 2.
  p_diffuse_i /= 2.0;
  p_diffuse_ii /= 2.0;
  // STOKES DRAG FOR XLINKS
  // for a single MT being driven towards the right (pos dir.),
  // hydrodynamic drag would increase p_diffuse_i_bck (neg dir.),
  // and decrease p_diffuse_i_fwd (pos dir.) by the same proportion
  double mt_force{Params::Filaments::f_applied[0]};
  if (mt_force != 0.0) {
    double vel{mt_force / filaments_->protofilaments_[0].gamma_[0]};
    printf("MT VEL: %g um/s\n", vel * 0.001);
    // Fluid flow experienced by crosslinkers is in opposite direction
    // (ETA is in units of pN*s/um^2; need to convert all nm to um to be valid)
    double f_drag{-1 * 6 * M_PI * Params::eta * (_r_xlink_head * 0.001) * vel};
    // Convert force to an energy energy by assuming it does work
    // on the xlink as it "jumps" from site to site
    double dE_drag{f_drag * Params::Filaments::site_size * 0.001};
    printf("f_drag = %g\n", f_drag);
    printf("dE_drag = %g\n", dE_drag);
    // Assume lambda = 0.5 for Boltzmann factor so we don't have to worry
    // about which "jump" is the forward or reverse reaction pathway
    double weight_drag_fwd{exp(0.5 * dE_drag / Params::kbT)};
    double weight_drag_bck{exp(-0.5 * dE_drag / Params::kbT)};
    printf("weight_fwd: %g\n", weight_drag_fwd);
    printf("weight_bck: %g\n", weight_drag_bck);
    Vec3D<double> weight_fwd{{Vec<double>(_n_neighbs_max + 1, p_diffuse_i)}};
    Vec3D<double> weight_bck{{Vec<double>(_n_neighbs_max + 1, p_diffuse_i)}};
    for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
      weight_fwd[0][0][n_neighbs] *= weight_drag_fwd;
      weight_bck[0][0][n_neighbs] *= weight_drag_bck;
    }
    weight_fwd[0][0][1] *= xlinks_.weights_.at("neighbs").unbind_[1];
    weight_bck[0][0][1] *= xlinks_.weights_.at("neighbs").unbind_[1];
    // cannot diffuse with 2 neighbs
    weight_fwd[0][0][2] *= 0.0;
    weight_bck[0][0][2] *= 0.0;
    xlinks_.AddProb("diffuse_i_fwd", weight_fwd);
    xlinks_.AddProb("diffuse_i_bck", weight_bck);
  } else {
    Vec3D<double> weight_diff{{Vec<double>(_n_neighbs_max + 1, p_diffuse_i)}};
    weight_diff[0][0][1] *= xlinks_.weights_.at("neighbs").unbind_[1];
    weight_diff[0][0][2] *= 0.0;
    xlinks_.AddProb("diffuse_i_fwd", weight_diff);
    xlinks_.AddProb("diffuse_i_bck", weight_diff);
  }
  xlinks_.AddProb("diffuse_ii_to_rest", p_diffuse_ii);
  xlinks_.AddProb("diffuse_ii_fr_rest", p_diffuse_ii);
  // Tether_Free
  motors_.AddProb("tether_free", Motors::k_tether * Motors::c_bulk * dt);
  xlinks_.AddProb("tether_free", Motors::k_tether * Xlinks::c_bulk * dt);
  // Tether_Bound -- motors only
  motors_.AddProb("tether_bound", Motors::k_tether * Motors::c_eff_tether * dt);
  // Untether_Free
  motors_.AddProb("untether_satellite", Motors::k_untether * dt);
  xlinks_.AddProb("untether_satellite", Motors::k_untether * dt);
  // Untether_Bound -- motors only
  motors_.AddProb("untether_bound", Motors::k_untether * dt);
}

void ProteinManager::InitializeEvents() {

  /*
  General strategy:
    We create flexible, generalized execution functions as Lambda Expressions.
    We re-use these functions for all protein species by passing an
    in-line lambda expression that effectively binds the appropriate arguments
    to the function input.
    Since the events are constructed in-place (via emplace_back), the lambda
    functons get optimized by the compiler and we do not see the typical
    performance drawback associated with passing functions/ptrs-to-functs/etc.
    (also, the use of 'auto' facilitates this modular lambda useage)
  See documentation (TODO) for further elaboration.
  */
  /* *** Probability distributions *** */
  auto binomial = [&](double p, int n) { // For events w/ constant probs.
    if (n > 0) {
      return SysRNG::SampleBinomial(p, n);
    } else {
      return 0;
    }
  };
  auto poisson = [&](double p, int n) { // For events w/ time-variable probs.
    if (p > 0.0) {
      return SysRNG::SamplePoisson(p);
    } else {
      return 0;
    }
  };
  /* *** Sorting functions - general *** */
  auto is_unocc = [](Object *base) -> Vec<Object *> {
    BindingSite *site{dynamic_cast<BindingSite *>(base)};
    if (!site->IsOccupied()) {
      return {site};
    }
    return {};
  };
  auto is_bound_i = [](Object *base) -> Vec<Object *> {
    Protein *protein{dynamic_cast<Protein *>(base)};
    if (protein->GetNumHeadsActive() == 1) {
      if (!protein->IsTethered() or protein->HasSatellite()) {
        return {protein->GetActiveHead()};
      }
    }
    return {};
  };
  auto is_bound_ii = [](Object *base) -> Vec<Object *> {
    Protein *protein{dynamic_cast<Protein *>(base)};
    if (protein->GetNumHeadsActive() == 2) {
      if (!protein->IsTethered() or protein->HasSatellite()) {
        return {protein->GetHeadOne(), protein->GetHeadTwo()};
      }
    }
    return {};
  };
  /* *** Sorting functions - specific to ATP hydrolysis cycle in motors ** */
  auto is_bound_i_NULL = [](Object *base) -> Vec<Object *> {
    Motor *motor = dynamic_cast<Motor *>(base);
    if (motor->GetNumHeadsActive() == 1) {
      if (!motor->IsTethered() or motor->HasSatellite()) {
        if (motor->GetActiveHead()->ligand_ == Ligand::NONE) {
          return {motor->GetActiveHead()};
        }
      }
    }
    return {};
  };
  auto is_bound_i_ATP = [](Object *base) -> Vec<Object *> {
    Motor *motor = dynamic_cast<Motor *>(base);
    if (motor->GetNumHeadsActive() == 1) {
      if (motor->GetActiveHead()->ligand_ == Ligand::ATP) {
        return {motor->GetActiveHead()};
      }
    }
    return {};
  };
  auto is_bound_i_ADPP = [](Object *base) -> Vec<Object *> {
    Motor *motor = dynamic_cast<Motor *>(base);
    if (motor->GetNumHeadsActive() == 1) {
      if (!motor->IsTethered() or motor->HasSatellite()) {
        if (motor->GetActiveHead()->ligand_ == Ligand::ADPP) {
          return {motor->GetActiveHead()};
        }
      }
    }
    return {};
  };
  auto is_docked = [](Object *base) -> Vec<Object *> {
    Motor *motor = dynamic_cast<Motor *>(base);
    auto *docked_head{motor->GetDockedHead()};
    if (docked_head != nullptr) {
      return {docked_head->GetOtherHead()};
    }
    return {};
  };
  auto is_bound_ii_ADPP = [](Object *base) -> Vec<Object *> {
    Motor *motor = dynamic_cast<Motor *>(base);
    if (motor->GetNumHeadsActive() == 2) {
      Object *chosen_head{nullptr};
      if (motor->GetHeadOne()->ligand_ == Ligand::ADPP) {
        chosen_head = motor->GetHeadOne();
      }
      if (motor->GetHeadTwo()->ligand_ == Ligand::ADPP) {
        if (chosen_head != nullptr) {
          Sys::ErrorExit("Protein_MNGR::is_bound_ii_ADPP()");
        }
        chosen_head = motor->GetHeadTwo();
      }
      if (chosen_head != nullptr) {
        return {chosen_head};
      }
    }
    return {};
  };
  auto is_bound_ii_NULL = [](Object *base) -> Vec<Object *> {
    Motor *motor = dynamic_cast<Motor *>(base);
    if (motor->GetNumHeadsActive() == 2) {
      Object *chosen_head{nullptr};
      if (motor->GetHeadOne()->ligand_ == Ligand::NONE) {
        chosen_head = motor->GetHeadOne();
      }
      if (motor->GetHeadTwo()->ligand_ == Ligand::NONE) {
        if (chosen_head != nullptr) {
          Sys::ErrorExit("Protein_MGR::is_bound_ii_NULL()");
        }
        chosen_head = motor->GetHeadTwo();
      }
      if (chosen_head != nullptr) {
        return {chosen_head};
      }
    }
    return {};
  };
  /* *** Sorting functions - specific to motor-xlink tethering *** */
  // ! TODO: return protein instead of arbitrary head
  auto is_a_satellite = [](Object *base) -> Vec<Object *> {
    Protein *protein{dynamic_cast<Protein *>(base)};
    if (protein->GetNumHeadsActive() == 0 and protein->IsTethered() and
        !protein->HasSatellite()) {
      return {protein->GetHeadOne()};
    }
    return {};
  };
  auto is_bound_i_teth = [](Object *base) -> Vec<Object *> {
    Protein *protein{dynamic_cast<Protein *>(base)};
    if (protein->GetNumHeadsActive() == 1 and protein->IsTethered() and
        !protein->HasSatellite()) {
      return {protein->GetActiveHead()};
    }
    return {};
  };
  auto is_bound_unteth = [](Object *base) -> Vec<Object *> {
    Protein *protein{dynamic_cast<Protein *>(base)};
    if (protein->GetNumHeadsActive() > 0 and !protein->IsTethered()) {
      if (protein->GetNumHeadsActive() == 2) {
        return {protein->GetHeadOne()};
      } else {
        return {protein->GetActiveHead()};
      }
    }
    return {};
  };
  auto is_bound_teth = [](Object *base) -> Vec<Object *> {
    Protein *protein{dynamic_cast<Protein *>(base)};
    if (protein->GetNumHeadsActive() > 0 and protein->IsTethered() and
        !protein->HasSatellite()) {
      if (protein->GetNumHeadsActive() == 2) {
        return {protein->GetHeadOne()};
      } else {
        return {protein->GetActiveHead()};
      }
    }
    return {};
  };
  auto is_bound_i_NULL_teth = [](Object *base) -> Vec<Object *> {
    Motor *motor = dynamic_cast<Motor *>(base);
    if (motor->GetNumHeadsActive() == 1 and motor->IsTethered() and
        !motor->HasSatellite()) {
      if (motor->GetActiveHead()->ligand_ == Ligand::NONE) {
        return {motor->GetActiveHead()};
      }
    }
    return {};
  };
  auto is_bound_i_ADPP_teth = [](Object *base) -> Vec<Object *> {
    Motor *motor = dynamic_cast<Motor *>(base);
    if (motor->GetNumHeadsActive() == 1 and motor->IsTethered() and
        !motor->HasSatellite()) {
      if (motor->GetActiveHead()->ligand_ == Ligand::ADPP) {
        return {motor->GetActiveHead()};
      }
    }
    return {};
  };

  /* *** Binning functions - used to generate multi-dim. sorted lists *** */
  // Sorted array dimension {i,j,k}; 1-D for neighbs so we use k and pad i & j
  Vec<size_t> dim_size{1, 1, _n_neighbs_max + 1};
  // Starting indices {i, j, k} of array for neighb coop, only use k dimension
  Vec<int> i_min{0, 0, 0};
  auto get_n_neighbs = [](Object *entry) {
    Vec<int> indices_vec{entry->GetNumNeighborsOccupied()};
    return indices_vec;
  };
  /* *** Bind_I *** */
  auto weight_bind_i = [](auto *site) { return site->GetWeight_Bind(); };
  auto exe_bind_i = [&](auto *site, auto *pop, auto *fil) {
    if (Sys::i_step_ < pop->step_active_) {
      return false;
    }
    auto entry{pop->GetFreeEntry()};
    if (entry == nullptr) {
      return false;
    }
    bool executed{entry->Bind(site, entry->GetHeadOne())};
    if (executed) {
      pop->AddToActive(entry);
    }
    return executed;
  };
  if (motors_.active_) {
    // Add unoccupied site tracker for motors; no binning b/c it's Poisson-based
    filaments_->AddPop("unbinned", is_unocc);
    // Create a single poisson event for all binding possibilities
    kmc_.events_.emplace_back(
        "bind_i (motors)", motors_.p_event_.at("bind_i").GetVal(),
        &filaments_->unoccupied_.at("unbinned").size_,
        &filaments_->unoccupied_.at("unbinned").entries_, poisson,
        [&](Object *base) {
          return weight_bind_i(dynamic_cast<BindingSite *>(base));
        },
        [&](Object *base) {
          return exe_bind_i(dynamic_cast<BindingSite *>(base), &motors_,
                            filaments_);
        });
  }
  if (xlinks_.active_) {
    // Add unoccupied site tracker for crosslinkers; segregated by n_neighbs
    filaments_->AddPop("neighbs", is_unocc, dim_size, i_min, get_n_neighbs);
    // Create a binomial event for each n_neighb possibility
    for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
      kmc_.events_.emplace_back(
          "bind_i_" + std::to_string(n_neighbs) + " (xlinks)",
          xlinks_.p_event_.at("bind_i").GetVal(n_neighbs),
          &filaments_->unoccupied_.at("neighbs").bin_size_[0][0][n_neighbs],
          &filaments_->unoccupied_.at("neighbs").bin_entries_[0][0][n_neighbs],
          binomial, [&](Object *base) {
            return exe_bind_i(dynamic_cast<BindingSite *>(base), &xlinks_,
                              filaments_);
          });
    }
  }

  /* *** Bind_I_Teth *** */
  if (motors_.active_ and xlinks_.active_ and Params::Motors::tethers_active) {
    motors_.AddPop("satellite", is_a_satellite);
    xlinks_.AddPop("satellite", is_a_satellite);
    auto weight_bind_i_teth = [](auto *satellite_head) {
      return satellite_head->parent_->GetWeight_Bind_I_Teth();
    };
    auto exe_bind_i_teth = [](auto *satellite_head, auto *pop, auto *altpop,
                              auto *fil) {
      auto site{satellite_head->parent_->GetNeighbor_Bind_I_Teth()};
      if (site == nullptr) {
        return false;
      }
      bool executed{satellite_head->parent_->Bind(site, satellite_head)};
      return executed;
    };
    kmc_.events_.emplace_back(
        "bind_i_teth (motors)", motors_.p_event_.at("bind_i_teth").GetVal(),
        &motors_.sorted_.at("satellite").size_,
        &motors_.sorted_.at("satellite").entries_, poisson,
        [&](Object *base) {
          return weight_bind_i_teth(dynamic_cast<CatalyticHead *>(base));
        },
        [&](Object *base) {
          return exe_bind_i_teth(dynamic_cast<CatalyticHead *>(base), &motors_,
                                 &xlinks_, filaments_);
        });
    kmc_.events_.emplace_back(
        "bind_i_teth (xlinks)", xlinks_.p_event_.at("bind_i_teth").GetVal(),
        &xlinks_.sorted_.at("satellite").size_,
        &xlinks_.sorted_.at("satellite").entries_, poisson,
        [&](Object *base) {
          return weight_bind_i_teth(dynamic_cast<BindingHead *>(base));
        },
        [&](Object *base) {
          return exe_bind_i_teth(dynamic_cast<BindingHead *>(base), &xlinks_,
                                 &motors_, filaments_);
        });
  }

  /* *** Bind_II *** */
  auto weight_bind_ii = [](auto *head) {
    return head->parent_->GetWeight_Bind_II();
  };
  auto exe_bind_ii = [](auto *bound_head, auto *pop, auto *fil) {
    auto head{bound_head->GetOtherHead()};
    auto site{head->parent_->GetNeighbor_Bind_II()};
    auto executed{head->parent_->Bind(site, head)};
    if (executed) {
      // FIXME check if this check is necessary / what to do
      bool still_attached{head->parent_->UpdateExtension()};
    }
    return executed;
  };
  if (xlinks_.crosslinking_active_) {
    xlinks_.AddPop("singly_bound", is_bound_i);
    kmc_.events_.emplace_back(
        "bind_ii (xlinks)", xlinks_.p_event_.at("bind_ii").GetVal(),
        &xlinks_.sorted_.at("singly_bound").size_,
        &xlinks_.sorted_.at("singly_bound").entries_, poisson,
        [&](Object *base) {
          return weight_bind_ii(dynamic_cast<BindingHead *>(base));
        },
        [&](Object *base) {
          return exe_bind_ii(dynamic_cast<BindingHead *>(base), &xlinks_,
                             filaments_);
        });
  }
  if (motors_.active_) {
    motors_.AddPop("docked", is_docked);
    kmc_.events_.emplace_back(
        "bind_ii (motors)", motors_.p_event_.at("bind_ii").GetVal(),
        &motors_.sorted_.at("docked").size_,
        &motors_.sorted_.at("docked").entries_, poisson,
        [&](Object *base) {
          return weight_bind_ii(dynamic_cast<CatalyticHead *>(base));
        },
        [&](Object *base) {
          return exe_bind_ii(dynamic_cast<CatalyticHead *>(base), &motors_,
                             filaments_);
        });
  }
  /* *** Unbind_II *** */
  auto weight_unbind_ii = [](auto *head) {
    return head->GetWeight_Unbind_II();
  };
  auto exe_unbind_ii = [](auto *head, auto *pop, auto *fil) {
    bool executed{head->Unbind()};
    return executed;
  };
  if (xlinks_.crosslinking_active_) {
    xlinks_.AddPop("doubly_bound", is_bound_ii);
    kmc_.events_.emplace_back(
        "unbind_ii (xlinks)", xlinks_.p_event_.at("unbind_ii").GetVal(),
        &xlinks_.sorted_.at("doubly_bound").size_,
        &xlinks_.sorted_.at("doubly_bound").entries_, poisson,
        [&](Object *base) {
          return weight_unbind_ii(dynamic_cast<BindingHead *>(base));
        },
        [&](Object *base) {
          return exe_unbind_ii(dynamic_cast<BindingHead *>(base), &xlinks_,
                               filaments_);
        });
  }
  if (motors_.active_) {
    motors_.AddPop("doubly_bound_ADPP", is_bound_ii_ADPP);
    kmc_.events_.emplace_back(
        "unbind_ii (motors)", motors_.p_event_.at("unbind_ii").GetVal(),
        &motors_.sorted_.at("doubly_bound_ADPP").size_,
        &motors_.sorted_.at("doubly_bound_ADPP").entries_, poisson,
        [&](Object *base) {
          return weight_unbind_ii(dynamic_cast<CatalyticHead *>(base));
        },
        [&](Object *base) {
          return exe_unbind_ii(dynamic_cast<CatalyticHead *>(base), &motors_,
                               filaments_);
        });
  }
  /* *** Unbind_I *** */
  auto exe_unbind_i = [](auto *head, auto *pop, auto *alt_pop, auto *fil) {
    bool executed{head->Unbind()};
    if (executed) {
      if (head->parent_->HasSatellite()) {
        alt_pop->RemoveFromActive(head->parent_->teth_partner_);
        bool untethered_sat{head->UntetherSatellite()};
        if (!untethered_sat) {
          Sys::ErrorExit("exe_unbind_i");
        }
      }
      if (!head->parent_->IsTethered()) {
        pop->RemoveFromActive(head->parent_);
      }
    }
    return executed;
  };
  if (xlinks_.active_) {
    xlinks_.AddPop("bound_i", is_bound_i, dim_size, i_min, get_n_neighbs);
    for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
      kmc_.events_.emplace_back(
          "unbind_i_" + std::to_string(n_neighbs) + " (xlinks)",
          xlinks_.p_event_.at("unbind_i").GetVal(n_neighbs),
          &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
          &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs],
          binomial, [&](Object *base) {
            return exe_unbind_i(dynamic_cast<BindingHead *>(base), &xlinks_,
                                &motors_, filaments_);
          });
    }
  }
  if (motors_.active_) {
    auto weight_unbind_i = [](auto *head) {
      return head->parent_->GetWeight_Unbind_I();
    };
    motors_.AddPop("bound_i_ADPP", is_bound_i_ADPP);
    kmc_.events_.emplace_back(
        "unbind_i (motors)", motors_.p_event_.at("unbind_i").GetVal(),
        &motors_.sorted_.at("bound_i_ADPP").size_,
        &motors_.sorted_.at("bound_i_ADPP").entries_, poisson,
        [&](Object *base) {
          return weight_unbind_i(dynamic_cast<CatalyticHead *>(base));
        },
        [&](Object *base) {
          return exe_unbind_i(dynamic_cast<CatalyticHead *>(base), &motors_,
                              &xlinks_, filaments_);
        });
  }
  /* *** Unbind_I_Teth *** */
  if (motors_.active_ and xlinks_.active_ and Params::Motors::tethers_active) {
    // force-dependent unbinding for motors
    auto weight_unbind_i_f = [](Object *base) {
      CatalyticHead *head{dynamic_cast<CatalyticHead *>(base)};
      double anchor_x{head->parent_->GetAnchorCoordinate(0)};
      double site_x{head->parent_->teth_partner_->GetAnchorCoordinate(0)};
      double r_x{anchor_x - site_x};
      // printf("r_x is %g\n", r_x);
      double dr{Params::Motors::r_0 - std::fabs(r_x)};
      // printf("dr is %g\n", dr);
      double k{dr < 0.0 ? Params::Motors::k_spring : Params::Motors::k_slack};
      double f{k * dr}; // positive force means spring is compressed
      // get direction of force; first orient using coordinates
      int f_dir{r_x > 0.0 ? 1 : -1};
      // if spring is compressed, flip signs since it's pushing not pulling
      if (dr > 0.0) {
        f_dir *= -1;
      }
      int dx{head->site_->filament_->dx_};
      // if opposite direction as stepping, treat as hindering load
      // printf("f is %g\n", f);
      double weight{0.0};
      if (f_dir == -dx) {
        weight =
            exp(std::fabs(k * dr) * Params::Motors::sigma_off_i / Params::kbT);
        // printf("HINDERING\n");
      }
      // otherwise, it is an assisting load (10x unbinding rate)
      else if (f_dir == dx) {
        weight = 10 * exp(std::fabs(k * dr) *
                          (0.5 * Params::Motors::sigma_off_i) / Params::kbT);
        // printf("ASSISTING \n");
      } else {
        printf("wtf\n");
      }
      // printf("weight is %g\n", weight);
      // printf(" ----- \n\n");
      // double weight{
      //     exp(-std::fabs(k * dr) * Params::Motors::sigma_off_i /
      //     Params::kbT)};
      if (weight > _max_weight) {
        printf("huh 1\n");
        return 0.0;
      }
      return weight;
    };
    motors_.AddPop("bound_i_ADPP_teth", is_bound_i_ADPP_teth);
    kmc_.events_.emplace_back(
        "unbind_i_teth (motors)", motors_.p_event_.at("unbind_i_teth").GetVal(),
        &motors_.sorted_.at("bound_i_ADPP_teth").size_,
        &motors_.sorted_.at("bound_i_ADPP_teth").entries_, poisson,
        weight_unbind_i_f, [&](Object *base) {
          return exe_unbind_i(dynamic_cast<CatalyticHead *>(base), &motors_,
                              &xlinks_, filaments_);
        });
    // energy-dependent unbinding for crosslinkers
    auto weight_unbind_i_u = [](auto *base) {
      BindingHead *head{dynamic_cast<BindingHead *>(base)};
      double site_x{head->parent_->GetAnchorCoordinate(0)};
      double anchor_x{head->parent_->teth_partner_->GetAnchorCoordinate(0)};
      double r_x{anchor_x - site_x};
      // FIXME make sure this ABS is correct
      double dr{Params::Motors::r_0 - std::fabs(r_x)};
      double k{dr < 0.0 ? Params::Motors::k_spring : Params::Motors::k_slack};
      double weight{exp(_lambda_spring * 0.5 * k * Square(dr) / Params::kbT)};
      if (weight > _max_weight) {
        return 0.0;
        printf("huh 2\n");
      }
      return weight;
    };
    xlinks_.AddPop("bound_i_teth", is_bound_i_teth);
    kmc_.events_.emplace_back(
        "unbind_i_teth (xlinks)", xlinks_.p_event_.at("unbind_i_teth").GetVal(),
        &xlinks_.sorted_.at("bound_i_teth").size_,
        &xlinks_.sorted_.at("bound_i_teth").entries_, poisson,
        weight_unbind_i_u, [&](Object *base) {
          return exe_unbind_i(dynamic_cast<BindingHead *>(base), &xlinks_,
                              &motors_, filaments_);
        });
  }
  /* *** Tether_Free *** */
  if (motors_.active_ and xlinks_.active_ and Params::Motors::tethers_active) {
    motors_.AddPop("bound_unteth", is_bound_unteth);
    xlinks_.AddPop("bound_unteth", is_bound_unteth);
    // Here, the head belongs to the bound untethered protein (or motor)
    auto exe_tether_free = [](auto *head, auto *pop, auto *altpop, auto *fil) {
      if (head->parent_->IsTethered()) {
        printf("what the FUCL @ %zu\n", head->GetID());
      }
      auto satellite{altpop->GetFreeEntry()};
      if (satellite == nullptr) {
        return false;
      }
      bool executed{head->parent_->Tether(satellite)};
      if (executed) {
        altpop->AddToActive(satellite);
        Sys::Log(1, "Added satellite %zu to protein %zu\n", satellite->GetID(),
                 head->GetID());
      }
      return executed;
    };
    kmc_.events_.emplace_back(
        "tether_free (motors)", motors_.p_event_.at("tether_free").GetVal(),
        &motors_.sorted_.at("bound_unteth").size_,
        &motors_.sorted_.at("bound_unteth").entries_, binomial,
        [&](Object *base) {
          return exe_tether_free(dynamic_cast<CatalyticHead *>(base), &motors_,
                                 &xlinks_, filaments_);
        });
    kmc_.events_.emplace_back(
        "tether_free (xlinks)", xlinks_.p_event_.at("tether_free").GetVal(),
        &xlinks_.sorted_.at("bound_unteth").size_,
        &xlinks_.sorted_.at("bound_unteth").entries_, binomial,
        [&](Object *base) {
          return exe_tether_free(dynamic_cast<BindingHead *>(base), &xlinks_,
                                 &motors_, filaments_);
        });
  }

  /* *** Untether_Free *** */
  if (motors_.active_ and xlinks_.active_ and Params::Motors::tethers_active) {
    auto exe_untether_satellite = [](auto *head, auto *pop, auto *alt_pop) {
      if (!head->parent_->teth_partner_->HasSatellite()) {
        Sys::ErrorExit("untether_satellite");
        return false;
      }
      Sys::Log(1, "untethering satellite %zu from protein %zu\n", head->GetID(),
               head->parent_->teth_partner_->GetID());
      bool executed{head->parent_->teth_partner_->UntetherSatellite()};
      if (executed) {
        pop->RemoveFromActive(head->parent_);
      }
      return executed;
    };
    kmc_.events_.emplace_back(
        "untether_satellite (motors)",
        motors_.p_event_.at("untether_satellite").GetVal(),
        &motors_.sorted_.at("satellite").size_,
        &motors_.sorted_.at("satellite").entries_, binomial, [&](Object *base) {
          return exe_untether_satellite(dynamic_cast<CatalyticHead *>(base),
                                        &motors_, &xlinks_);
        });
    kmc_.events_.emplace_back(
        "untether_satellite (xlinks)",
        xlinks_.p_event_.at("untether_satellite").GetVal(),
        &xlinks_.sorted_.at("satellite").size_,
        &xlinks_.sorted_.at("satellite").entries_, binomial, [&](Object *base) {
          return exe_untether_satellite(dynamic_cast<BindingHead *>(base),
                                        &xlinks_, &motors_);
        });
  }
  // vv motors only vv
  /* *** Bind_ATP_I *** */
  if (motors_.active_) {
    auto exe_bind_ATP = [](auto *head, auto *pop) {
      bool executed{head->parent_->Bind_ATP(head)};
      return executed;
    };
    motors_.AddPop("bound_i_NULL", is_bound_i_NULL);
    kmc_.events_.emplace_back(
        "bind_ATP_i", motors_.p_event_.at("bind_ATP_i").GetVal(),
        &motors_.sorted_.at("bound_i_NULL").size_,
        &motors_.sorted_.at("bound_i_NULL").entries_, binomial,
        [&](Object *base) {
          return exe_bind_ATP(dynamic_cast<CatalyticHead *>(base), &motors_);
        });

    /* *** Bind_ATP_I_Teth *** */
    if (xlinks_.active_ and motors_.tethering_active_) {
      motors_.AddPop("bound_i_NULL_teth", is_bound_i_NULL_teth);
      auto weight_bind_ATP_i_teth = [](Object *base) {
        CatalyticHead *head{dynamic_cast<CatalyticHead *>(base)};
        double anchor_x{head->parent_->GetAnchorCoordinate(0)};
        double site_x{head->parent_->teth_partner_->GetAnchorCoordinate(0)};
        double r_x{anchor_x - site_x};
        // FIXME check to make sure ABS is correct
        double dr{Params::Motors::r_0 - std::fabs(r_x)};
        double k{dr < 0.0 ? Params::Motors::k_spring : Params::Motors::k_slack};
        double weight{
            exp(-std::fabs(k * dr) * Params::Motors::sigma_ATP / Params::kbT)};
        // printf("WEIGHT: %g\n", weight);
        // if (std::fabs(r_x) > 135) {
        //   printf("whut\n");
        // }
        if (weight < (1.0 / _max_weight)) {
          printf("NOPE\n");
          return 0.0;
        }
        return weight;
      };
      kmc_.events_.emplace_back(
          "bind_ATP_i_teth", motors_.p_event_.at("bind_ATP_i_teth").GetVal(),
          &motors_.sorted_.at("bound_i_NULL_teth").size_,
          &motors_.sorted_.at("bound_i_NULL_teth").entries_, poisson,
          weight_bind_ATP_i_teth, [&](Object *base) {
            return exe_bind_ATP(dynamic_cast<CatalyticHead *>(base), &motors_);
          });
    }
    /* *** Bind_ATP_II *** */
    if (Params::Motors::gaussian_stepping_coop and
        Params::Motors::gaussian_range > 0) {
      auto exe_bind_ATP_ii = [](auto *front_head, auto *pop, auto *fil) {
        auto *rear_head{front_head->GetOtherHead()};
        if (front_head->trailing_) {
          return false;
        }
        bool unbound{rear_head->Unbind()};
        bool executed{front_head->parent_->Bind_ATP(front_head)};
        return executed;
      };
      auto weight_bind_ATP_ii = [](auto *head) {
        return head->parent_->GetWeight_BindATP_II(head);
      };
      // motors_.AddPop("bound_ii_NULL", is_bound_ii_NULL);
      // kmc_.events_.emplace_back(
      //     "bind_ATP_ii", motors_.p_event_.at("bind_ATP_ii").GetVal(),
      //     &motors_.sorted_.at("bound_ii_NULL").size_,
      //     &motors_.sorted_.at("bound_ii_NULL").entries_, poisson,
      //     [&](Object *base) {
      //       return weight_bind_ATP_ii(dynamic_cast<CatalyticHead *>(base));
      //     },
      //     [&](Object *base) {
      //       return exe_bind_ATP_ii(dynamic_cast<CatalyticHead *>(base),
      //                              &motors_, filaments_);
      //     });
    }
    // *** Bind_ATP_II_Teth ***
    // ! TODO: add
  }
  /* *** Hydrolyze_ATP *** */
  if (motors_.active_) {
    auto exe_hydrolyze = [](auto *head, auto *pop) {
      bool executed{head->parent_->Hydrolyze(head)};
      return executed;
    };
    motors_.AddPop("bound_i_ATP", is_bound_i_ATP);
    kmc_.events_.emplace_back(
        "hydrolyze", motors_.p_event_.at("hydrolyze").GetVal(),
        &motors_.sorted_.at("bound_i_ATP").size_,
        &motors_.sorted_.at("bound_i_ATP").entries_, binomial,
        [&](Object *base) {
          return exe_hydrolyze(dynamic_cast<CatalyticHead *>(base), &motors_);
        });
  }
  /* *** Tether_Bound *** */
  if (motors_.active_ and xlinks_.active_ and Params::Motors::tethers_active) {
    motors_.AddPop("bound_unteth", is_bound_unteth);
    auto weight_tether = [](auto *base) {
      CatalyticHead *head{dynamic_cast<CatalyticHead *>(base)};
      return head->parent_->GetWeight_Tether();
    };
    auto exe_tether = [](auto *head, auto *pop, auto *alt_pop) {
      auto xlink{head->parent_->GetNeighbor_Tether()};
      if (xlink == nullptr) {
        return false;
      }
      bool executed{head->parent_->Tether(xlink)};
      if (executed) {
        Sys::Log(1, "Tethered motor %zu to protein %zu\n", head->GetID(),
                 xlink->GetID());
      }
      return executed;
    };
    kmc_.events_.emplace_back(
        "tether_bound (motors)", motors_.p_event_.at("tether_bound").GetVal(),
        &motors_.sorted_.at("bound_unteth").size_,
        &motors_.sorted_.at("bound_unteth").entries_, poisson, weight_tether,
        [&](Object *base) {
          return exe_tether(dynamic_cast<CatalyticHead *>(base), &motors_,
                            &xlinks_);
        });
  }
  /* *** Untether_Bound *** */
  if (motors_.active_ and xlinks_.active_ and Params::Motors::tethers_active) {
    motors_.AddPop("bound_teth", is_bound_teth);
    // energy-dependent untethering
    auto weight_untether = [](auto *base) {
      CatalyticHead *head{dynamic_cast<CatalyticHead *>(base)};
      double anchor_x{head->parent_->GetAnchorCoordinate(0)};
      double site_x{head->parent_->teth_partner_->GetAnchorCoordinate(0)};
      double r_x{anchor_x - site_x};
      // FIXME make sure this ABS is correct
      double dr{Params::Motors::r_0 - std::fabs(r_x)};
      double k{dr < 0.0 ? Params::Motors::k_spring : Params::Motors::k_slack};
      double weight{exp(_lambda_spring * 0.5 * k * Square(dr) / Params::kbT)};
      if (weight > _max_weight) {
        return 0.0;
        printf("huh 4\n");
      }
      return weight;
    };
    auto exe_untether = [](auto *head, auto *pop, auto *alt_pop) {
      Sys::Log(1, "Untethering motor %zu from protein %zu\n", head->GetID(),
               head->parent_->teth_partner_->GetID());
      bool executed{head->parent_->Untether()};
      return executed;
    };
    kmc_.events_.emplace_back("untether_bound (motors)",
                              motors_.p_event_.at("untether_bound").GetVal(),
                              &motors_.sorted_.at("bound_teth").size_,
                              &motors_.sorted_.at("bound_teth").entries_,
                              poisson, weight_untether, [&](Object *base) {
                                return exe_untether(
                                    dynamic_cast<CatalyticHead *>(base),
                                    &motors_, &xlinks_);
                              });
  }
  // vv xlinks only vv
  // Diffusion
  auto exe_diff = [](auto *head, auto *pop, auto *fil, int dir) {
    // if (head->parent_->IsTethered() and !head->parent_->HasSatellite()) {
    //   double anchor_x{head->parent_->teth_partner_->GetAnchorCoordinate(0)};
    //   // horizontal tethers for now
    //   double r_y{0.0};
    //   double r{sqrt(Square(r_y) + Square(anchor_x - head->site_->pos_[0]))};
    //   printf("r is %g\n", r);
    //   // if (r > 140) {
    //   //   exit(1);
    //   // }
    // }
    bool executed{head->Diffuse(dir)};
    // if (head->parent_->IsTethered() and !head->parent_->HasSatellite()) {
    //   double anchor_x{head->parent_->teth_partner_->GetAnchorCoordinate(0)};
    //   // horizontal tethers for now
    //   double r_y{0.0};
    //   double r{sqrt(Square(r_y) + Square(anchor_x - head->site_->pos_[0]))};
    //   printf("r is now %g\n\n", r);
    //   if (r > 140) {
    //     exit(1);
    //   }
    // }
    if (executed) {
      bool still_attached{head->parent_->UpdateExtension()};
      // TODO do I need this check??
      if (!still_attached) {
        // printf("what\n");
      }
    }
    return executed;
  };
  if (xlinks_.active_) {
    for (int n_neighbs{0}; n_neighbs < _n_neighbs_max; n_neighbs++) {
      kmc_.events_.emplace_back(
          "diffuse_i_fwd",
          xlinks_.p_event_.at("diffuse_i_fwd").GetVal(n_neighbs),
          &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
          &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs],
          binomial, [&](Object *base) {
            return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_,
                            filaments_, 1);
          });
      kmc_.events_.emplace_back(
          "diffuse_i_bck",
          xlinks_.p_event_.at("diffuse_i_bck").GetVal(n_neighbs),
          &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
          &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs],
          binomial, [&](Object *base) {
            return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_,
                            filaments_, -1);
          });
    }
    if (motors_.active_ and motors_.tethering_active_) {
      xlinks_.AddPop("bound_i_teth", is_bound_i_teth);
      auto weight_diffuse_teth = [](auto head, int dir) {
        BindingSite *old_site{head->site_};
        BindingSite *new_site{old_site->GetNeighbor(dir)};
        if (new_site == nullptr) {
          return 0.0;
        }
        double anchor_x{head->parent_->teth_partner_->GetAnchorCoordinate(0)};
        // horizontal tethers for now
        double r_y{0.0};
        double r_old{sqrt(Square(r_y) + Square(anchor_x - old_site->pos_[0]))};
        double dr_old{r_old - Params::Motors::r_0};
        double energy_old{dr_old > 0.0
                              ? 0.5 * Params::Motors::k_spring * Square(dr_old)
                              : 0.5 * Params::Motors::k_slack * Square(dr_old)};
        double r_new{sqrt(Square(r_y) + Square(anchor_x - new_site->pos_[0]))};
        double dr_new{r_new - Params::Motors::r_0};
        double energy_new{dr_new > 0.0
                              ? 0.5 * Params::Motors::k_spring * Square(dr_new)
                              : 0.5 * Params::Motors::k_slack * Square(dr_new)};
        double dE{energy_new - energy_old};
        double weight{0.0};
        if (dE < 0.0) {
          weight = exp(_lambda_spring * fabs(dE) / Params::kbT);
        } else {
          weight = exp(-(1.0 - _lambda_spring) * fabs(dE) / Params::kbT);
        }
        if (r_new > 141) {
          return 0.0;
          // printf("MOTHER @ %zu - weight: %g\n", head->GetID(), weight);
          // exit(1);
        }
        // if (weight > _max_weight) {
        //   printf("r_old: %g\n", r_old);
        //   printf("r_new: %g\n", r_new);
        //   printf("dE: %g - %g = %g\n", energy_new, energy_old, dE);
        //   printf("huh 5 - %g\n", weight);
        //   Sys::ErrorExit("U NO");
        //   // return 0.0;
        // }
        if (weight < (1.0 / _max_weight)) {
          printf("HYUK\n\n");
          // Sys::ErrorExit("U NOTT");
          return 0.0;
        }
        return weight;
      };
      kmc_.events_.emplace_back(
          "diffuse_i_fwd_teth", xlinks_.p_event_.at("diffuse_i_fwd").GetVal(0),
          &xlinks_.sorted_.at("bound_i_teth").size_,
          &xlinks_.sorted_.at("bound_i_teth").entries_, poisson,
          [&](Object *base) {
            return weight_diffuse_teth(dynamic_cast<BindingHead *>(base), 1);
          },
          [&](Object *base) {
            return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_,
                            filaments_, 1);
          });
      kmc_.events_.emplace_back(
          "diffuse_i_bck_teth", xlinks_.p_event_.at("diffuse_i_bck").GetVal(0),
          &xlinks_.sorted_.at("bound_i_teth").size_,
          &xlinks_.sorted_.at("bound_i_teth").entries_, poisson,
          [&](Object *base) {
            return weight_diffuse_teth(dynamic_cast<BindingHead *>(base), -1);
          },
          [&](Object *base) {
            return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_,
                            filaments_, -1);
          });
    }
  }
  if (xlinks_.crosslinking_active_) {
    auto weight_diff_ii = [](auto *head, int dir) {
      return head->GetWeight_Diffuse(dir);
    };
    xlinks_.AddPop("diffuse_ii_to_rest", is_bound_ii);
    xlinks_.AddPop("diffuse_ii_fr_rest", is_bound_ii);
    kmc_.events_.emplace_back(
        "diffuse_ii_to_rest",
        xlinks_.p_event_.at("diffuse_ii_to_rest").GetVal(),
        &xlinks_.sorted_.at("diffuse_ii_to_rest").size_,
        &xlinks_.sorted_.at("diffuse_ii_to_rest").entries_, poisson,
        [&](Object *base) {
          return weight_diff_ii(dynamic_cast<BindingHead *>(base), 1);
        },
        [&](Object *base) {
          return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_,
                          filaments_, 1);
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
          return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_,
                          filaments_, -1);
        });
  }
  // Diffuse_II_Teth
  // ! TODO: need to add
  // Bind_II_Teth
  // ! TODO: Need to add
  // Unbind_II_Teth
  // ! TODO: Need to add
}

void ProteinManager::RunKMC() {

  // if (Sys::i_step_ >= 1910000) {
  //   Sys::verbosity_ = 1;
  // }
  UpdateFilaments();
  motors_.PrepForKMC();
  xlinks_.PrepForKMC();
  // Sys::Log(1, "BEGIN KMC EVENTS\n");
  bool event_executed{kmc_.ExecuteEvents()};
  if (event_executed) {
    motors_.FlagForUpdate();
    xlinks_.FlagForUpdate();
    filaments_->FlagForUpdate();
  }
  // UpdateExtensions();
  // Sys::Log(1, "END KMC EVENTS\n");
}

void ProteinManager::CheckProbabilities() {
  for (auto const &event : kmc_.events_) {
    Sys::Log(1, "p_%s = %g\n", event.name_.c_str(), event.p_occur_);
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