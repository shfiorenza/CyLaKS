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

  /* ! FIXME
  switch (_n_neighbs_max) {
  case 0:
    Sys::Log("  Nearest-neighbor interactions turned off.\n");
    break;
  case 2:
    if (Params::Filaments::n_subfilaments > 1) {
      Sys::Log("  Neighbors on adjacent protofilaments are not counted.\n");
    }
    break;
  case 4:
    if (Params::Filaments::n_subfilaments == 1) {
      Sys::Log("Error. '_n_neighbs_max' must be either 0 or 2 for 1-D "
               "protofilaments.\n");
      Sys::ErrorExit("ProteinManager::InitializeWeights()");
    } else {
      Sys::Log("  Neighbors on adjacent protofilaments are counted.\n");
    }
    break;
  default:
    Sys::Log("Error; '_n_neighbs_max' must be eiither 0, 2, or 4 (if "
             "subfilaments enabled)\n");
    Sys::ErrorExit("ProteinManager::InitializeWeights()");
  }
  */
  Str name{"neighbs"};
  Sys::weight_neighb_bind_.resize(_n_neighbs_max + 1);
  Sys::weight_neighb_unbind_.resize(_n_neighbs_max + 1);
  motors_.AddWeight("neighbs", _n_neighbs_max + 1);
  xlinks_.AddWeight("neighbs", _n_neighbs_max + 1);
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
  // for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
  //   printf("bind: %g (%i neighbs)\n",
  //   xlinks_.weights_[name].bind_[n_neighbs],
  //          n_neighbs);
  //   printf("unbind: %g (%i neighbs)\n",
  //          xlinks_.weights_[name].unbind_[n_neighbs], n_neighbs);
  // }
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
  // Unbind_II
  double wt_unbind_ii{exp(Motors::internal_force * Motors::sigma_off_ii / kbT)};
  if (Motors::applied_force > 0.0) {
    wt_unbind_ii *= exp(-Motors::applied_force * Motors::sigma_off_ii / kbT);
  }
  // printf("wt_unbind_ii: %g\n", wt_unbind_ii);
  motors_.AddProb("unbind_ii", Motors::k_off_ii * dt * wt_unbind_ii);
  xlinks_.AddProb("unbind_ii", Xlinks::k_off_ii * dt);
  // Unbind I
  double wt_unbind_i{1.0};
  if (Motors::applied_force > 0.0) {
    wt_unbind_i = exp(Motors::applied_force * Motors::sigma_off_i / kbT);
  }
  motors_.AddProb("unbind_i", Motors::k_off_i * dt * wt_unbind_i);
  xlinks_.AddProb("unbind_i", Xlinks::k_off_i * dt, "neighbs", 1);
  if (Params::Filaments::n_subfilaments > 1) {
    double p_bind_i{Xlinks::k_on * Xlinks::c_bulk * dt};
    double p_unbind_i{Xlinks::k_off_i * dt};
    Vec3D<double> p_bind{{Vec<double>(2 * _n_neighbs_max + 1, p_bind_i)}};
    Vec3D<double> p_unbind{{Vec<double>(2 * _n_neighbs_max + 1, p_unbind_i)}};
    for (int n_neighbs{0}; n_neighbs <= 2 * _n_neighbs_max; n_neighbs++) {
      double dE{-1 * Params::Xlinks::neighb_neighb_energy * n_neighbs};
      p_bind[0][0][n_neighbs] *= exp(-(1.0 - _lambda_neighb) * dE);
      p_unbind[0][0][n_neighbs] *= exp(_lambda_neighb * dE);
    }
    xlinks_.AddProb("bind_i_multi", p_bind);
    xlinks_.AddProb("unbind_i_multi", p_unbind);
  }
  // Diffusion -- xlinks only
  double x_sq{Square(Filaments::site_size / 1000)}; // in um^2
  double tau_i{x_sq / (2 * Xlinks::d_i)};
  double tau_ii{x_sq / (2 * Xlinks::d_ii)};
  // ! FIXME this should use a different x_sq for side-to-side spacing ( nm ?)
  double tau_side{x_sq / (2 * Xlinks::d_side)};
  double p_diffuse_i{dt / tau_i};
  double p_diffuse_ii{dt / tau_ii};
  double p_diffuse_side{dt / tau_side};
  // diff_fwd and diff_bck are two separate events, which effectively
  // doubles the probability to diffuse. Thus we divide p_diff by 2.
  p_diffuse_i /= 2.0;
  p_diffuse_ii /= 2.0;
  p_diffuse_side /= 2.0;
  // STOKES DRAG FOR XLINKS
  // for a single MT being driven towards the right (pos dir.),
  // hydrodynamic drag would increase p_diffuse_i_bck (neg dir.),
  // and decrease p_diffuse_i_fwd (pos dir.) by the same proportion
  double mt_force{Params::Filaments::f_applied[0]};
  if (mt_force != 0.0) {
    if (Params::Filaments::n_subfilaments != 1) {
      Sys::Log("Microtubule gliding w/ multiple protofilaments is not "
               "currently implemented; exiting.\n");
      Sys::ErrorExit("ProteinManager::SetParameters()");
    }
    double vel{mt_force / filaments_->protofilaments_[0].gamma_[0]}; // nm/s
    printf("MT VEL: %g um/s\n", vel * 0.001);
    // Fluid flow experienced by crosslinkers is in opposite direction
    // (ETA is in units of pN*s/um^2; need to convert all nm to um to be valid)
    // double r_xlink_head{16.0}; // nm
    double f_drag{-1 * 6 * M_PI * Params::eta * (_r_xlink_head * 0.001) *
                  (vel * 0.001)};
    // Convert force to an energy energy by assuming it does work
    // on the xlink as it "jumps" from site to site
    double dE_drag{f_drag * Params::Filaments::site_size}; // pN*nm
    printf("f_drag = %g\n", f_drag);
    printf("dE_drag = %g\n", dE_drag);
    // Assume lambda = 0.5 for Boltzmann factor so we don't have to worry
    // about which "jump" is the forward or reverse reaction pathway
    double weight_drag_fwd{exp(0.5 * dE_drag / Params::kbT)};
    double weight_drag_bck{exp(-0.5 * dE_drag / Params::kbT)};
    printf("weight_fwd: %g\n", weight_drag_fwd);
    printf("weight_bck: %g\n", weight_drag_bck);
    Vec3D<double> p_fwd{{Vec<double>(_n_neighbs_max + 1, p_diffuse_i)}};
    Vec3D<double> p_bck{{Vec<double>(_n_neighbs_max + 1, p_diffuse_i)}};
    for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
      p_fwd[0][0][n_neighbs] *= weight_drag_fwd;
      p_bck[0][0][n_neighbs] *= weight_drag_bck;
    }
    p_fwd[0][0][1] *= xlinks_.weights_.at("neighbs").unbind_[1];
    p_bck[0][0][1] *= xlinks_.weights_.at("neighbs").unbind_[1];
    // cannot diffuse with 2 neighbs
    p_fwd[0][0][2] *= 0.0;
    p_bck[0][0][2] *= 0.0;
    xlinks_.AddProb("diffuse_i_fwd", p_fwd);
    xlinks_.AddProb("diffuse_i_bck", p_bck);

  } else {
    if (Params::Filaments::n_subfilaments == 1) {
      Vec3D<double> p_diff{{Vec<double>(_n_neighbs_max + 1, p_diffuse_i)}};
      p_diff[0][0][1] *= xlinks_.weights_.at("neighbs").unbind_[1];
      p_diff[0][0][2] *= 0.0;
      xlinks_.AddProb("diffuse_i_fwd", p_diff);
      xlinks_.AddProb("diffuse_i_bck", p_diff);
    } else {
      size_t vec_size{_n_neighbs_max + 1};
      Vec3D<double> p_diff_same{
          {Vec2D<double>(vec_size, Vec<double>(vec_size, p_diffuse_i))}};
      Vec3D<double> p_diff_side{
          {Vec2D<double>(vec_size, Vec<double>(vec_size, p_diffuse_side))}};
      // Scan over n_neighbors along same protofilament
      for (int n_same{0}; n_same <= _n_neighbs_max; n_same++) {
        // Scan over n_neighbors along adjacent protofilaments
        for (int n_side{0}; n_side <= _n_neighbs_max; n_side++) {
          double tot_weight{xlinks_.weights_.at("neighbs").unbind_[n_same] *
                            xlinks_.weights_.at("neighbs").unbind_[n_side]};
          p_diff_same[0][n_side][n_same] *= tot_weight;
          p_diff_side[0][n_side][n_same] *= tot_weight;
          // With max neighbs in respective dimensions, cannot diffuse
          if (n_same == _n_neighbs_max) {
            p_diff_same[0][n_side][n_same] *= 0.0;
          }
          if (n_side == _n_neighbs_max) {
            p_diff_side[0][n_side][n_same] *= 0.0;
          }
        }
      }
      xlinks_.AddProb("diffuse_i_fwd", p_diff_same);
      xlinks_.AddProb("diffuse_i_bck", p_diff_same);
      xlinks_.AddProb("diffuse_side_up", p_diff_side);
      xlinks_.AddProb("diffuse_side_down", p_diff_side);
    }
  }
  xlinks_.AddProb("diffuse_ii_to_rest", p_diffuse_ii);
  xlinks_.AddProb("diffuse_ii_fr_rest", p_diffuse_ii);
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
  /* *** Binning functions - used to generate multi-dim. sorted lists *** */
  // Sorted array dimension {i,j,k}; 1-D for neighbs so we use k and pad i & j
  Vec<size_t> dim_size{1, 1, _n_neighbs_max + 1};
  Vec<size_t> dim_size_tot{1, 1, 2 * _n_neighbs_max + 1};
  // Starting indices {i, j, k} of array for neighb coop, only use k dimension
  Vec<int> i_min{0, 0, 0};
  auto get_n_neighbs = [](Object *entry) {
    Vec<int> indices_vec{entry->GetNumNeighborsOccupied()};
    return indices_vec;
  };
  auto get_n_neighbs_tot = [](Object *entry) {
    Vec<int> indices_vec{entry->GetNumNeighborsOccupied() +
                         entry->GetNumNeighborsOccupied_Side()};
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
    if (Params::Filaments::n_subfilaments == 1) {
      // Add unoccupied site tracker for crosslinkers; segregated by n_neighbs
      filaments_->AddPop("neighbs", is_unocc, dim_size, i_min, get_n_neighbs);
      // Create a binomial event for each n_neighb possibility
      for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
        kmc_.events_.emplace_back(
            "bind_i_" + std::to_string(n_neighbs) + " (xlinks)",
            xlinks_.p_event_.at("bind_i").GetVal(n_neighbs),
            &filaments_->unoccupied_.at("neighbs").bin_size_[0][0][n_neighbs],
            &filaments_->unoccupied_.at("neighbs")
                 .bin_entries_[0][0][n_neighbs],
            binomial, [&](Object *base) {
              return exe_bind_i(dynamic_cast<BindingSite *>(base), &xlinks_,
                                filaments_);
            });
      }
    } else {
      filaments_->AddPop("neighbs_tot", is_unocc, dim_size_tot, i_min,
                         get_n_neighbs_tot);
      for (int n_neighbs{0}; n_neighbs <= 2 * _n_neighbs_max; n_neighbs++) {
        kmc_.events_.emplace_back(
            "bind_i_" + std::to_string(n_neighbs) + " (xlinks)",
            xlinks_.p_event_.at("bind_i_multi").GetVal(n_neighbs),
            &filaments_->unoccupied_.at("neighbs_tot")
                 .bin_size_[0][0][n_neighbs],
            &filaments_->unoccupied_.at("neighbs_tot")
                 .bin_entries_[0][0][n_neighbs],
            binomial, [&](Object *base) {
              return exe_bind_i(dynamic_cast<BindingSite *>(base), &xlinks_,
                                filaments_);
            });
      }
    }
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
    if (Params::Filaments::n_subfilaments == 1) {
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
    } else {
      xlinks_.AddPop("bound_i", is_bound_i, dim_size_tot, i_min,
                     get_n_neighbs_tot);
      for (int n_neighbs{0}; n_neighbs <= 2 * _n_neighbs_max; n_neighbs++) {
        kmc_.events_.emplace_back(
            "unbind_i_" + std::to_string(n_neighbs) + " (xlinks)",
            xlinks_.p_event_.at("unbind_i_multi").GetVal(n_neighbs),
            &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
            &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs],
            binomial, [&](Object *base) {
              return exe_unbind_i(dynamic_cast<BindingHead *>(base), &xlinks_,
                                  &motors_, filaments_);
            });
      }
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
      motors_.AddPop("bound_ii_NULL", is_bound_ii_NULL);
      kmc_.events_.emplace_back(
          "bind_ATP_ii", motors_.p_event_.at("bind_ATP_ii").GetVal(),
          &motors_.sorted_.at("bound_ii_NULL").size_,
          &motors_.sorted_.at("bound_ii_NULL").entries_, poisson,
          [&](Object *base) {
            return weight_bind_ATP_ii(dynamic_cast<CatalyticHead *>(base));
          },
          [&](Object *base) {
            return exe_bind_ATP_ii(dynamic_cast<CatalyticHead *>(base),
                                   &motors_, filaments_);
          });
    }
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

  // vv xlinks only vv
  if (xlinks_.active_) {
    // Diffusion
    auto exe_diff = [](auto *head, auto *pop, auto *fil, int dir) {
      bool executed{head->Diffuse(dir)};
      if (executed) {
        pop->FlagForUpdate();
        fil->FlagForUpdate();
        bool still_attached{head->parent_->UpdateExtension()};
        // TODO do I need this check??
        if (!still_attached) {
          // printf("what\n");
        }
      }
      return executed;
    };
    // Simple 1-D binomial for single-PF diffusion
    if (Params::Filaments::n_subfilaments == 1) {
      for (int n_neighbs{0}; n_neighbs < _n_neighbs_max; n_neighbs++) {
        kmc_.events_.emplace_back(
            "diffuse_i_" + std::to_string(n_neighbs) + "_fwd",
            xlinks_.p_event_.at("diffuse_i_fwd").GetVal(n_neighbs),
            &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
            &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs],
            binomial, [&](Object *base) {
              return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_,
                              filaments_, 1);
            });
        kmc_.events_.emplace_back(
            "diffuse_i_" + std::to_string(n_neighbs) + "_bck",
            xlinks_.p_event_.at("diffuse_i_bck").GetVal(n_neighbs),
            &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
            &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs],
            binomial, [&](Object *base) {
              return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_,
                              filaments_, -1);
            });
        // SIDE STEP DIFFUSION
      }
    }
    // If we have multiple PFs, need to account for neighbs in different dims.
    else {
      // For multi-PF diffusion, need to segregate neighbors on same PF versus
      // adjacent PFs, so we just pad i
      Vec<size_t> dim_size_multi{1, _n_neighbs_max + 1, _n_neighbs_max + 1};
      auto get_n_neighbs_multi = [](Object *entry) {
        Vec<int> indices_vec{entry->GetNumNeighborsOccupied(),
                             entry->GetNumNeighborsOccupied_Side()};
        return indices_vec;
      };
      auto exe_diff_side = [](auto *head, auto *pop, auto *fil, int dir) {
        bool executed{head->Diffuse_Side(dir)};
        if (executed) {
          pop->FlagForUpdate();
          fil->FlagForUpdate();
        }
        return executed;
      };
      xlinks_.AddPop("bound_i_multi", is_bound_i, dim_size_multi, i_min,
                     get_n_neighbs_multi);
      // FIXME make 2-d loop
      for (int n_same{0}; n_same <= _n_neighbs_max; n_same++) {
        for (int n_side{0}; n_side <= _n_neighbs_max; n_side++) {
          Str neighb_label{std::to_string(n_side) + "_" +
                           std::to_string(n_same)};
          if (n_same != _n_neighbs_max) {
            kmc_.events_.emplace_back(
                "diffuse_i_fwd_" + neighb_label,
                xlinks_.p_event_.at("diffuse_i_fwd").GetVal(n_side, n_same),
                &xlinks_.sorted_.at("bound_i_multi")
                     .bin_size_[0][n_side][n_same],
                &xlinks_.sorted_.at("bound_i_multi")
                     .bin_entries_[0][n_side][n_same],
                binomial, [&](Object *base) {
                  return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_,
                                  filaments_, 1);
                });
            kmc_.events_.emplace_back(
                "diffuse_i_bck_" + neighb_label,
                xlinks_.p_event_.at("diffuse_i_bck").GetVal(n_side, n_same),
                &xlinks_.sorted_.at("bound_i_multi")
                     .bin_size_[0][n_side][n_same],
                &xlinks_.sorted_.at("bound_i_multi")
                     .bin_entries_[0][n_side][n_same],
                binomial, [&](Object *base) {
                  return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_,
                                  filaments_, -1);
                });
          }
          // SIDE STEP DIFFUSION
          if (n_side != _n_neighbs_max) {
            kmc_.events_.emplace_back(
                "diffuse_i_up_" + neighb_label,
                xlinks_.p_event_.at("diffuse_side_up").GetVal(n_side, n_same),
                &xlinks_.sorted_.at("bound_i_multi")
                     .bin_size_[0][n_side][n_same],
                &xlinks_.sorted_.at("bound_i_multi")
                     .bin_entries_[0][n_side][n_same],
                binomial, [&](Object *base) {
                  return exe_diff_side(dynamic_cast<BindingHead *>(base),
                                       &xlinks_, filaments_, 1);
                });
            kmc_.events_.emplace_back(
                "diffuse_i_down_" + neighb_label,
                xlinks_.p_event_.at("diffuse_side_down").GetVal(n_side, n_same),
                &xlinks_.sorted_.at("bound_i_multi")
                     .bin_size_[0][n_side][n_same],
                &xlinks_.sorted_.at("bound_i_multi")
                     .bin_entries_[0][n_side][n_same],
                binomial, [&](Object *base) {
                  return exe_diff_side(dynamic_cast<BindingHead *>(base),
                                       &xlinks_, filaments_, -1);
                });
          }
        }
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
  }
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