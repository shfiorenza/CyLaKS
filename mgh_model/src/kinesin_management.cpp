#include "kinesin_management.h"
#include "master_header.h"

KinesinManagement::KinesinManagement() {}

void KinesinManagement::Initialize(system_parameters *parameters,
                                   system_properties *properties) {

  // verbosity_ = 1;
  wally_ = &properties->wallace;
  gsl_ = &properties->gsl;
  parameters_ = parameters;
  properties_ = properties;
  SetParameters();
  GenerateMotors();
  InitializeLists();
  if (wally_->test_mode_ != nullptr) {
    InitializeTestEnvironment();
  }
  InitializeEvents();
}

void KinesinManagement::SetParameters() {

  double kbT{parameters_->kbT};
  double delta_t{parameters_->delta_t};
  step_active_ = (size_t)std::ceil(parameters_->motors.t_active / delta_t);
  /*
    For events that result in a change in energy, we use Boltzmann factors to
    scale rates appropriately. Detailed balance is satisfied with the factors:
                  exp{-(1 - lambda)*(delta_E)/(kB*T)}, and
                  exp{-(lambda)*(delta_E)/(kB*T)}
    for an event and its complement (e.g., binding and unbinding), where lambda
    is a constant that ranges from 0 to 1, delta_E is the change in energy that
    results from the event, kB is Boltzmann's constant, and T is the temperature
 */
  // Lambda = 0.5 means energy dependence is equal for binding and unbinding
  double lambda_lattice{0.5}; // Lambda for lattice interactions
  // Lambda = 1.0 means all energy dependence is in unbinding
  double lambda_neighb{1.0}; // Lambda for neighbor cooperativity

  // Array of Boltzmann weights due to neighbor-neighbor interactions
  double energy_per_neighb{-1 * parameters_->motors.interaction_energy}; // kbT
  weight_neighbs_bind_.resize(max_neighbs_ + 1);
  weight_neighbs_unbind_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    double energy{n_neighbs * energy_per_neighb}; // ! in kbT
    weight_neighbs_bind_[n_neighbs] = exp(-(1.0 - lambda_neighb) * energy);
    weight_neighbs_unbind_[n_neighbs] = exp(lambda_neighb * energy);
  }
  // Array of maximum possible weights including neighbor interactions and
  // lattice deformation effects from MULTIPLE motors stacking
  double lattice_E_max_tot{-1 * parameters_->motors.lattice_coop_Emax_bulk};
  double wt_lattice_bind_max{exp(-(1.0 - lambda_lattice) * lattice_E_max_tot)};
  double wt_lattice_unbind_max{exp(lambda_lattice * lattice_E_max_tot)};
  weight_lattice_bind_max_.resize(max_neighbs_ + 1);
  weight_lattice_unbind_max_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    weight_lattice_bind_max_[n_neighbs] =
        wt_lattice_bind_max * weight_neighbs_bind_[n_neighbs];
    weight_lattice_unbind_max_[n_neighbs] =
        wt_lattice_unbind_max * weight_neighbs_unbind_[n_neighbs];
  }
  // Calculate lattice_alpha based on input range
  lattice_cutoff_ = parameters_->motors.lattice_coop_range;
  // Set lattice_alpha_ so that E = 0 at site immediately after cutoff
  double site_size{parameters_->microtubules.site_size};
  double dx_cutoff{(lattice_cutoff_ + 1) * site_size};
  double lattice_E_0_solo{-1 * parameters_->motors.lattice_coop_Emax_solo};
  double lattice_alpha{-1 * lattice_E_0_solo / (dx_cutoff * dx_cutoff)};
  wally_->Log("\nFor motors:\n");
  if (lattice_cutoff_ > 0) {
    lattice_coop_active_ = true;
    wally_->Log("  lattice_cutoff_ is %i\n", lattice_cutoff_);
    wally_->Log("  lattice_alpha_ is %g\n", lattice_alpha);
  } else {
    wally_->Log("  Lattice cooperativity is disabled.\n");
  }
  // Array of binding/unbinding weights due to lattice deformation from a SINGLE
  // motor; multiplied together to get total weight for any arrangement
  weight_lattice_bind_.resize(lattice_cutoff_ + 1);
  weight_lattice_unbind_.resize(lattice_cutoff_ + 1);
  for (int delta{0}; delta <= lattice_cutoff_; delta++) {
    double dx{delta * site_size};
    double energy{lattice_alpha * dx * dx + lattice_E_0_solo}; // in kbT
    weight_lattice_bind_[delta] = exp(-(1.0 - lambda_lattice) * energy);
    weight_lattice_unbind_[delta] = exp(lambda_lattice * energy);
  }
  // Calculate probabilities for each possible KMC event
  double k_on{parameters_->motors.k_on};
  double c_bulk{parameters_->motors.c_bulk};
  p_avg_bind_i_ = k_on * c_bulk * delta_t;
  double c_eff_bind{parameters_->motors.c_eff_bind};
  p_avg_bind_ii_ = k_on * c_eff_bind * delta_t;
  double k_on_ATP{parameters_->motors.k_on_ATP};
  double c_ATP{parameters_->motors.c_ATP};
  p_bind_ATP_i_ = k_on_ATP * c_ATP * delta_t;
  double internal_force{parameters_->motors.internal_force};
  double sigma_ATP{parameters_->motors.sigma_ATP};
  double weight_internal_ATP{exp(-internal_force * sigma_ATP / kbT)};
  p_avg_bind_ATP_ii_ = p_bind_ATP_i_ * weight_internal_ATP;
  double k_hydrolyze{parameters_->motors.k_hydrolyze};
  p_hydrolyze_ = k_hydrolyze * delta_t;
  double sigma_off_ii{parameters_->motors.sigma_off_ii};
  double weight_internal_off_ii{exp(internal_force * sigma_off_ii / kbT)};
  double k_off_ii{parameters_->motors.k_off_ii};
  p_avg_unbind_ii_ = k_off_ii * delta_t * weight_internal_off_ii;
  double k_off_i{parameters_->motors.k_off_i};
  p_avg_unbind_i_ = k_off_i * delta_t;
  // Force perpetually applied to motors, e.g., by optical trapping
  if (parameters_->motors.applied_force > 0.0) {
    double kbT{parameters_->kbT};
    double force{parameters_->motors.applied_force};
    wally_->Log("\n Motor applied force is %g\n", force);
    double p_unbind_i_old{p_avg_unbind_i_};
    double p_unbind_ii_old{p_avg_unbind_ii_};
    double p_bind_ATP_i_old{p_bind_ATP_i_};
    double p_bind_ATP_ii_old{p_avg_bind_ATP_ii_};
    double sigma_off_i{parameters_->motors.sigma_off_i};
    double sigma_off_ii{parameters_->motors.sigma_off_ii};
    double sigma_ATP{parameters_->motors.sigma_ATP};
    // Motor singly-bound unbinding rate increases with increasing force
    p_avg_unbind_i_ *= exp(force * sigma_off_i / kbT);
    wally_->Log("p_avg_unbind_i scaled from %g to %g \n", p_unbind_i_old,
                p_avg_unbind_i_);
    // Motor doubly-bound unbinding rate decreases with increasing force
    // because it partially cancels out the internal necklinker forces
    p_avg_unbind_ii_ *= exp(-1 * force * sigma_off_ii / kbT);
    wally_->Log("p_avg_unbind_ii scaled from %g to %g \n", p_unbind_ii_old,
                p_avg_unbind_ii_);
    // ATP binding rate decreases with increasing force
    p_bind_ATP_i_ *= exp(-1 * force * sigma_ATP / kbT);
    wally_->Log("p_bind_ATP_i scaled from %g to %g \n", p_bind_ATP_i_old,
                p_bind_ATP_i_);
    p_avg_bind_ATP_ii_ *= exp(-1 * force * sigma_ATP / kbT);
    wally_->Log("p_avg_bind_ATP_ii scaled from %g to %g \n", p_bind_ATP_ii_old,
                p_avg_bind_ATP_ii_);
  }
  p_theory_["bind_i"] = p_avg_bind_i_;
  p_theory_["bind_ATP_i"] = p_bind_ATP_i_;
  p_theory_["bind_ATP_ii"] = p_avg_bind_ATP_ii_;
  p_theory_["hydrolyze"] = p_hydrolyze_;
  p_theory_["bind_ii"] = p_avg_bind_ii_;
  p_theory_["unbind_ii"] = p_avg_unbind_ii_;
  p_theory_["unbind_i"] = p_avg_unbind_i_;
  // READ EM OUT
  for (const auto &entry : p_theory_) {
    auto label = entry.first;
    std::string name{"p_" + label};
    auto value = entry.second;
    if (value > 1.0) {
      wally_->Log("Error! %s = %g\n", name.c_str(), value);
      wally_->ErrorExit("Kin_MGMT:SetParameters()");
    }
    wally_->Log("%s = %g\n", name.c_str(), value);
  }
}

void KinesinManagement::GenerateMotors() {

  int n_mts{parameters_->microtubules.count};
  // Calculate total number of motors to have in reservoir
  for (int i_mt = 0; i_mt < n_mts; i_mt++) {
    n_motors_ += parameters_->microtubules.length[i_mt];
  }
  // Since only one head has to be bound, the most that will ever
  // be needed (all single-bound) is the total number of sites
  motors_.resize(n_motors_);
  for (int id = 0; id < n_motors_; id++)
    motors_[id].Initialize(parameters_, properties_, id);
}

void KinesinManagement::InitializeLists() {

  // 1-D stuff -- indexed by i_motor
  active_.resize(n_motors_);
  bound_NULL_i_.resize(n_motors_);
  bound_NULL_ii_.resize(n_motors_);
  bound_ATP_.resize(n_motors_);
  docked_.resize(n_motors_);
  bound_ADPP_ii_.resize(n_motors_);
  bound_ADPP_i_.resize(n_motors_);
}

void KinesinManagement::InitializeTestEnvironment() {

  if (strcmp(wally_->test_mode_, "motor_lattice_bind") == 0) {
    wally_->Log("Initializing test %s\n", wally_->test_mode_);
    MicrotubuleManagement *mts = &properties_->microtubules;
    mts->InitializeTestEnvironment();
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      weight_neighbs_bind_[n_neighbs] = 1.0;
      weight_neighbs_unbind_[n_neighbs] = 1.0;
    }
    // Set p_avg_bind_i_ to a large value for better statistics
    double c_bulk{10.0};
    double k_on{1.0};
    p_avg_bind_i_ = k_on * c_bulk * parameters_->delta_t;
    // Choose middle MT site and bind motor to it
    int i_site{lattice_cutoff_};
    Tubulin *site{&properties_->microtubules.mt_list_[0].lattice_[i_site]};
    Bind_I(site);
    properties_->microtubules.FlagForUpdate();
    properties_->microtubules.UpdateUnoccupied();
    // Initialize stat tracker
    Vec<std::pair<size_t, size_t>> zeros(lattice_cutoff_ + 1, {0, 0});
    test_stats_.emplace("lattice_bind", zeros);
    Vec<double> p_theory(lattice_cutoff_ + 1, p_avg_bind_i_);
    for (int delta{0}; delta <= lattice_cutoff_; delta++) {
      p_theory[delta] *= weight_lattice_bind_[delta];
    }
    test_ref_.emplace("lattice_bind", p_theory);
  }
  if (strcmp(wally_->test_mode_, "motor_lattice_step") == 0) {
    wally_->Log("Initializing test %s\n", wally_->test_mode_);
    properties_->microtubules.InitializeTestEnvironment();
    // If test_delta is left uninitialized, test checks against self-coop
    double test_weight_bind{1.0};
    double test_weight_unbind{1.0};
    int test_delta{5};
    if (test_delta > 0) {
      test_weight_bind = weight_lattice_bind_[test_delta];
      test_weight_unbind = weight_lattice_unbind_[test_delta];
      for (int delta{0}; delta <= lattice_cutoff_; delta++) {
        weight_lattice_bind_[delta] = test_weight_bind;
        weight_lattice_unbind_[delta] = test_weight_unbind;
      }
    }
    // Initialize stat trackers
    Vec<std::pair<size_t, size_t>> zeros(1, {0, 0});
    Vec<double> p_theory_bind_ii(1, p_avg_bind_ii_ * test_weight_bind);
    test_stats_.emplace("lattice_step_bind_ii", zeros);
    test_ref_.emplace("lattice_step_bind_ii", p_theory_bind_ii);
    Vec<double> p_theory_unbind_ii(1, p_avg_unbind_ii_ * test_weight_unbind);
    test_stats_.emplace("lattice_step_unbind_ii", zeros);
    test_ref_.emplace("lattice_step_unbind_ii", p_theory_unbind_ii);
    Vec<double> p_theory_unbind_i(1, p_avg_unbind_i_ * test_weight_unbind);
    test_stats_.emplace("lattice_step_unbind_i", zeros);
    test_ref_.emplace("lattice_step_unbind_i", p_theory_unbind_i);
    // Place motor on minus end of microtubule
    Microtubule *mt = &properties_->microtubules.mt_list_[0];
    Tubulin *site{&mt->lattice_[mt->minus_end_]};
    Bind_I(site);
    properties_->microtubules.FlagForUpdate();
    properties_->microtubules.UpdateUnoccupied();
    lattice_coop_active_ = false;
  }
}

void KinesinManagement::InitializeEvents() {

  //  Binomial probabilitiy distribution; sampled to predict most events
  auto binomial = [&](double p, int n) {
    if (n > 0) {
      return properties_->gsl.SampleBinomialDist(p, n);
    } else {
      return 0;
    }
  };
  // Poisson distribution; sampled to predict events w/ variable probabilities
  auto poisson = [&](double p, int n) {
    if (p > 0.0) {
      return gsl_->SamplePoissonDist(p);
    } else {
      return 0;
    }
  };
  // Function that returns a random probability in the range [0.0, 1.0)
  auto get_ran_prob = [&](void) { return gsl_->GetRanProb(); };
  //  Function that sets n random indices from the range [0, m)
  auto set_ran_indices = [&](int *indices, int n, int m) {
    if (m > 1) {
      properties_->gsl.SetRanIndices(indices, n, m);
    } else {
      indices[0] = 0;
    }
  };
  //  Event entries
  std::string event_name; // scratch space to construct each event name
  // If we are testing a mechanism, initialize test events only
  if (wally_->test_mode_ != nullptr) {
    if (strcmp(wally_->test_mode_, "motor_lattice_bind") == 0) {
      // Bind_I: bind first motor head to MT and release ADP
      event_name = "bind_i";
      auto exe_bind_i = [&](ENTRY_T target) {
        Tubulin *site{std::get<SITE_T *>(target)};
        int i_site{site->index_};
        // 'main' kinesin motor will always be at index_ = lattice_cutoff_
        int delta{abs(i_site - lattice_cutoff_)};
        // Just count the event; do not actually bind a motor to this site
        test_stats_["lattice_bind"][delta].first++;
        properties_->microtubules.FlagForUpdate();
      };
      auto weight_bind_i = [&](ENTRY_T target) {
        return std::get<SITE_T *>(target)->weight_bind_;
      };
      auto poisson_bind_i = [&](double p, int n) {
        // Each delta distance has 2 sites available to it each timestep
        // Do not count delta = 0, where 'main' motor is permanently bound to
        for (int delta{1}; delta <= lattice_cutoff_; delta++) {
          test_stats_["lattice_bind"][delta].second += 2;
        }
        if (p > 0.0) {
          return gsl_->SamplePoissonDist(p);
        }
        return 0;
      };
      events_.emplace_back(
          event_name, p_avg_bind_i_, &properties_->microtubules.n_unocc_motor_,
          &properties_->microtubules.unocc_motor_, exe_bind_i, weight_bind_i,
          poisson_bind_i, get_ran_prob, set_ran_indices);
    }
    if (strcmp(wally_->test_mode_, "motor_lattice_step") == 0) {
      // Bind_ATP: bind ATP to motor heads that have released their ADP
      auto exe_bind_ATP = [&](ENTRY_T target) {
        Bind_ATP(std::get<POP_T *>(target));
      };
      event_name = "bind_ATP_i";
      events_.emplace_back(event_name, p_bind_ATP_i_, &n_bound_NULL_i_,
                           &bound_NULL_i_, exe_bind_ATP, binomial,
                           set_ran_indices);
      // Bind_ATP_II: bind ATP to doubly-bound motor heads; rate modified by
      // lattice
      event_name = "bind_ATP_ii";
      auto exe_bind_ATP_ii = [&](ENTRY_T target) {
        try {
          POP_T *front_head{std::get<POP_T *>(target)};
          POP_T *rear_head{front_head->GetOtherHead()};
          wally_->Log(2, "starting unbind_ii\n");
          Unbind_II(rear_head);
          wally_->Log(2, "starting bind_ATP\n");
          Bind_ATP(front_head);
        } catch (...) {
          wally_->ErrorExit("K_MGMT::exe_bind_ATP_II");
        }
      };
      auto weight_bind_ATP_ii = [&](ENTRY_T target) {
        try {
          return std::get<POP_T *>(target)->motor_->GetWeight_BindATP_II();
        } catch (...) {
          wally_->ErrorExit("K_MGMT::weight_bind_ATP_ii");
        }
      };
      events_.emplace_back(event_name, p_avg_bind_ATP_ii_, &n_bound_NULL_ii_,
                           &bound_NULL_ii_, exe_bind_ATP_ii, weight_bind_ATP_ii,
                           poisson, get_ran_prob, set_ran_indices);
      // Hydrolyze: convert ATP to ADPP
      event_name = "hydrolyze";
      auto exe_hydrolyze = [&](ENTRY_T target) {
        Hydrolyze(std::get<POP_T *>(target));
      };
      events_.emplace_back(event_name, p_hydrolyze_, &n_bound_ATP_, &bound_ATP_,
                           exe_hydrolyze, binomial, set_ran_indices);
      // Bind_II: binds docked motor head to MT and releases its ADP
      event_name = "bind_ii";
      auto exe_bind_ii = [&](ENTRY_T target) {
        test_stats_["lattice_step_bind_ii"][0].first++;
        // If dock site is plus end, unbind motor and place it on minus end
        POP_T *docked_head{std::get<POP_T *>(target)};
        SITE_T *dock_site{docked_head->motor_->GetDockSite()};
        if (dock_site->index_ == dock_site->mt_->plus_end_) {
          Unbind_I(docked_head->motor_->GetActiveHead());
          SITE_T *new_site{
              &dock_site->mt_->lattice_[dock_site->mt_->minus_end_]};
          Bind_I(new_site);
          Bind_ATP(new_site->motor_head_);
          Hydrolyze(new_site->motor_head_);
          docked_head = new_site->motor_head_->motor_->GetDockedHead();
        }
        Bind_II(docked_head);
      };
      auto weight_bind_ii = [&](ENTRY_T target) {
        return std::get<POP_T *>(target)->motor_->GetWeight_Bind_II();
      };
      auto poisson_bind_ii = [&](double p, int n) {
        test_stats_["lattice_step_bind_ii"][0].second += n_docked_;
        if (p > 0.0) {
          return gsl_->SamplePoissonDist(p);
        }
        return 0;
      };
      events_.emplace_back(event_name, p_avg_bind_ii_, &n_docked_, &docked_,
                           exe_bind_ii, weight_bind_ii, poisson_bind_ii,
                           get_ran_prob, set_ran_indices);
      // Unbind_II: Converts ADPP to ADP and unbinds a doubly-bound head
      event_name = "unbind_ii";
      auto exe_unbind_ii = [&](ENTRY_T target) {
        test_stats_["lattice_step_unbind_ii"][0].first++;
        Unbind_II(std::get<POP_T *>(target));
      };
      auto weight_unbind_ii = [&](ENTRY_T target) {
        double weight{std::get<POP_T *>(target)->motor_->GetWeight_Unbind_II()};
        return weight;
      };
      auto poisson_unbind_ii = [&](double p, int n) {
        test_stats_["lattice_step_unbind_ii"][0].second += n_bound_ADPP_ii_;
        if (p > 0.0) {
          return gsl_->SamplePoissonDist(p);
        }
        return 0;
      };
      events_.emplace_back(event_name, p_avg_unbind_ii_, &n_bound_ADPP_ii_,
                           &bound_ADPP_ii_, exe_unbind_ii, weight_unbind_ii,
                           poisson_unbind_ii, get_ran_prob, set_ran_indices);
      // Unbind_I: Converts ADPP to ADP and unbinds a singly-bound head
      event_name = "unbind_i";
      auto exe_unbind_i = [&](ENTRY_T target) {
        // Count stats for unbind_i but do not actually execute it
        test_stats_["lattice_step_unbind_i"][0].first++;
      };
      auto weight_unbind_i = [&](ENTRY_T target) {
        return std::get<POP_T *>(target)->site_->weight_unbind_;
      };
      auto poisson_unbind_i = [&](double p, int n) {
        test_stats_["lattice_step_unbind_i"][0].second += n_bound_ADPP_i_;
        if (p > 0.0) {
          return gsl_->SamplePoissonDist(p);
        }
        return 0;
      };
      events_.emplace_back(event_name, p_avg_unbind_i_, &n_bound_ADPP_i_,
                           &bound_ADPP_i_, exe_unbind_i, weight_unbind_i,
                           poisson_unbind_i, get_ran_prob, set_ran_indices);
    }
    return;
  }
  // Bind_I: bind first motor head to MT and release ADP
  event_name = "bind_i";
  auto exe_bind_i = [&](ENTRY_T target) {
    try {
      Bind_I(std::get<SITE_T *>(target));
    } catch (...) {
      wally_->ErrorExit("K_MGMT::exe_bind_i");
    }
  };
  auto weight_bind_i = [&](ENTRY_T target) {
    try {
      return std::get<SITE_T *>(target)->weight_bind_;
    } catch (...) {
      wally_->ErrorExit("K_MGMT::weight_bind_i");
    }
  };
  events_.emplace_back(event_name, p_avg_bind_i_,
                       &properties_->microtubules.n_unocc_motor_,
                       &properties_->microtubules.unocc_motor_, exe_bind_i,
                       weight_bind_i, poisson, get_ran_prob, set_ran_indices);
  // Bind_ATP: bind ATP to motor heads that have released their ADP
  event_name = "bind_ATP_i";
  auto exe_bind_ATP = [&](ENTRY_T target) {
    try {
      Bind_ATP(std::get<POP_T *>(target));
    } catch (...) {
      wally_->ErrorExit("K_MGMT::exe_bind_ATP");
    }
  };
  events_.emplace_back(event_name, p_bind_ATP_i_, &n_bound_NULL_i_,
                       &bound_NULL_i_, exe_bind_ATP, binomial, set_ran_indices);
  // Bind_ATP_II: bind ATP to doubly-bound motor heads; rate modified by lattice
  event_name = "bind_ATP_ii";
  auto exe_bind_ATP_ii = [&](ENTRY_T target) {
    try {
      POP_T *front_head{std::get<POP_T *>(target)};
      POP_T *rear_head{front_head->GetOtherHead()};
      // If front head is trailing, this motor is end-paused -- skip
      if (front_head->trailing_) {
        return;
      }
      Unbind_II(rear_head);
      Bind_ATP(front_head);
    } catch (...) {
      wally_->ErrorExit("K_MGMT::exe_bind_ATP_II");
    }
  };
  auto weight_bind_ATP_ii = [&](ENTRY_T target) {
    try {
      return std::get<POP_T *>(target)->motor_->GetWeight_BindATP_II();
    } catch (...) {
      wally_->ErrorExit("K_MGMT::weight_bind_ATP_ii");
    }
  };
  events_.emplace_back(event_name, p_avg_bind_ATP_ii_, &n_bound_NULL_ii_,
                       &bound_NULL_ii_, exe_bind_ATP_ii, weight_bind_ATP_ii,
                       poisson, get_ran_prob, set_ran_indices);
  // Hydrolyze: convert ATP to ADPP
  event_name = "hydrolyze";
  auto exe_hydrolyze = [&](ENTRY_T target) {
    try {
      Hydrolyze(std::get<POP_T *>(target));
    } catch (...) {
      wally_->ErrorExit("K_MGMT::exe_hydrolyze");
    }
  };
  events_.emplace_back(event_name, p_hydrolyze_, &n_bound_ATP_, &bound_ATP_,
                       exe_hydrolyze, binomial, set_ran_indices);
  // Bind_II: binds docked motor head to MT and releases its ADP
  event_name = "bind_ii";
  auto exe_bind_ii = [&](ENTRY_T target) {
    try {
      Bind_II(std::get<POP_T *>(target));
    } catch (...) {
      wally_->ErrorExit("K_MGMT::exe_bind_ii");
    }
  };
  auto weight_bind_ii = [&](ENTRY_T target) {
    try {
      return std::get<POP_T *>(target)->motor_->GetWeight_Bind_II();
    } catch (...) {
      wally_->ErrorExit("K_MGMT::weight_bind_ii");
    }
  };
  events_.emplace_back(event_name, p_avg_bind_ii_, &n_docked_, &docked_,
                       exe_bind_ii, weight_bind_ii, poisson, get_ran_prob,
                       set_ran_indices);

  // Unbind_II: Converts ADPP to ADP and unbinds a doubly-bound head
  event_name = "unbind_ii";
  auto exe_unbind_ii = [&](ENTRY_T target) {
    try {
      Unbind_II(std::get<POP_T *>(target));
    } catch (...) {
      wally_->ErrorExit("K_MGMT::exe_unbind_ii");
    }
  };
  auto weight_unbind_ii = [&](ENTRY_T target) {
    try {
      return std::get<POP_T *>(target)->motor_->GetWeight_Unbind_II();
    } catch (...) {
      wally_->ErrorExit("K_MGMT::weight_unbind_ii");
    }
  };
  events_.emplace_back(event_name, p_avg_unbind_ii_, &n_bound_ADPP_ii_,
                       &bound_ADPP_ii_, exe_unbind_ii, weight_unbind_ii,
                       poisson, get_ran_prob, set_ran_indices);
  // Unbind_I: Converts ADPP to ADP and unbinds a singly-bound head
  event_name = "unbind_i";
  auto exe_unbind_i = [&](ENTRY_T target) {
    try {
      Unbind_I(std::get<POP_T *>(target));
    } catch (...) {
      wally_->ErrorExit("K_MGMT::exe_unbind_i");
    }
  };
  auto weight_unbind_i = [&](ENTRY_T target) {
    try {
      return std::get<POP_T *>(target)->site_->weight_unbind_;
    } catch (...) {
      wally_->ErrorExit("K_MGMT::weight_unbind_i");
    }
  };
  events_.emplace_back(event_name, p_avg_unbind_i_, &n_bound_ADPP_i_,
                       &bound_ADPP_i_, exe_unbind_i, weight_unbind_i, poisson,
                       get_ran_prob, set_ran_indices);
  if (wally_->verbosity_ >= 3) {
    wally_->Log("\nMotor events: \n");
    for (const auto &event : events_) {
      wally_->Log("   %s\n", event.name_.c_str());
    }
  }
}

void KinesinManagement::ReportProbabilities() {

  // Baseline function -- readout of all event probabilities
  for (const auto &entry : p_theory_) {
    auto name = entry.first;
    auto value = entry.second;
    auto event = std::find_if(
        events_.begin(), events_.end(),
        [&name](EVENT_T const &event) { return event.name_ == name; });
    if (event->name_ != name) {
      wally_->Log("Couldn't find %s\n", name.c_str());
      continue;
    }
    if (event->n_opportunities_tot_ == 0) {
      wally_->Log("No statistics for %s\n", name.c_str());
      continue;
    }
    size_t n_exe{event->n_executed_tot_};
    size_t n_opp{event->n_opportunities_tot_};
    wally_->Log("Statistics for kinesin event '%s':\n", name.c_str());
    wally_->Log("   p_theory = %g\n", value);
    wally_->Log("   p_actual = %g", double(n_exe) / n_opp);
    wally_->Log(" (n_exe = %lu | n_opp = %lu)\n", n_exe, n_opp);
  }
  // Test-specific event probabilities
  if (wally_->test_mode_ != nullptr) {
    wally_->Log("Statistics for test %s:\n", wally_->test_mode_);
  }
  for (const auto &entry : test_stats_) {
    auto name = entry.first;
    auto p_theory =
        std::find_if(test_ref_.begin(), test_ref_.end(),
                     [&name](std::pair<std::string, Vec<double>> parameter) {
                       return parameter.first == name;
                     });
    if (p_theory->first != name) {
      wally_->Log("ERROR! '%s' not found\n", name.c_str());
      continue;
    }
    wally_->Log(" Event '%s':\n", name.c_str());
    for (int index{0}; index < entry.second.size(); index++) {
      wally_->Log("  p_theory_[%i] = %g\n", index, p_theory->second[index]);
      size_t n_exe{entry.second[index].first};
      size_t n_opp{entry.second[index].second};
      wally_->Log("  p_actual_[%i] = %g ", index, double(n_exe) / n_opp);
      wally_->Log("(n_exe = %i | n_opp = %i)\n", n_exe, n_opp);
    }
  }
}

Kinesin *KinesinManagement::GetFreeMotor() {

  // Randomly pick a motor from the reservoir
  int i_motor = properties_->gsl.GetRanInt(n_motors_);
  Kinesin *motor = &motors_[i_motor];
  int attempts = 0;
  while (motor->heads_active_ > 0) {
    i_motor++;
    if (i_motor == n_motors_) {
      i_motor = 0;
    }
    motor = &motors_[i_motor];
    attempts++;
    if (attempts > n_motors_) {
      wally_->ErrorExit("Kin_MGMT::GetFreeMotor()");
    }
  }
  return motor;
}

void KinesinManagement::UpdateLatticeWeights() {

  if (!lattice_coop_active_) {
    return;
  }
  // Add up all lattice deformation weights on top of baseline neighb weights
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    Kinesin *motor{active_[i_entry]};
    if (motor->heads_active_ == 0) {
      continue;
    }
    Tubulin *epicenter{nullptr};
    if (motor->heads_active_ == 1) {
      epicenter = motor->GetActiveHead()->site_;
    } else if (motor->head_one_.trailing_) {
      epicenter = motor->head_one_.site_;
    } else {
      epicenter = motor->head_two_.site_;
    }
    int i_epicenter{epicenter->index_};
    for (int delta{1}; delta <= lattice_cutoff_; delta++) {
      for (int dir{-1}; dir <= 1; dir += 2) {
        int i_scan{i_epicenter + dir * delta};
        // Only access sites that exist
        if (i_scan < 0 or i_scan >= epicenter->mt_->n_sites_) {
          continue;
        }
        Tubulin *site{&epicenter->mt_->lattice_[i_scan]};
        int n_neighbs{site->GetKif4ANeighborCount()};
        if (site->weight_bind_ > weight_lattice_bind_max_[n_neighbs]) {
          site->weight_bind_ = weight_lattice_bind_max_[n_neighbs];
          site->weight_unbind_ = weight_lattice_unbind_max_[n_neighbs];
          continue;
        }
        site->weight_bind_ *= weight_lattice_bind_[delta];
        site->weight_unbind_ *= weight_lattice_unbind_[delta];
      }
    }
  }
}

void KinesinManagement::AddToActive(Kinesin *motor) {

  active_[n_active_] = motor;
  motor->active_index_ = n_active_;
  n_active_++;
}

void KinesinManagement::RemoveFromActive(Kinesin *motor) {

  Kinesin *last_entry{active_[n_active_ - 1]};
  int this_index{motor->active_index_};
  active_[this_index] = last_entry;
  last_entry->active_index_ = this_index;
  n_active_--;
}

void KinesinManagement::RunKMC() {

  if (properties_->current_step_ < step_active_) {
    return;
  }
  // if (properties_->current_step_ >= 16910100) {
  //   properties_->wallace.verbosity_ = 3;
  // }
  UpdateLists();
  SampleEventStatistics();
  GenerateExecutionSequence();
  ExecuteEvents();
}

void KinesinManagement::UpdateLists() {

  wally_->Log(1, "Starting K_MGMT::UpdateLists()\n");
  if (lists_up_to_date_) {
    wally_->Log(1, " -- already up-to-date --\n");
    return;
  }
  lists_up_to_date_ = true;
  // Update MT's unoccupied list -- also calls UpdateLatticeWeights();
  properties_->microtubules.UpdateUnoccupied();
  // Reset all population sizes to zero before rebuilding lists
  n_bound_NULL_i_ = 0;
  n_bound_NULL_ii_ = 0;
  n_bound_ATP_ = 0;
  n_docked_ = 0;
  n_bound_ADPP_ii_ = 0;
  n_bound_ADPP_i_ = 0;
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    Kinesin *motor = active_[i_entry];
    if (motor->heads_active_ == 1) {
      POP_T *head{motor->GetActiveHead()};
      if (head->ligand_ == "NULL") {
        bound_NULL_i_[n_bound_NULL_i_++] = head;
        wally_->Log(2, "motor #%i added to bound_NULL_i_[%i]\n",
                    head->motor_->id_, n_bound_NULL_i_ - 1);
      } else if (head->ligand_ == "ATP") {
        bound_ATP_[n_bound_ATP_++] = head;
        wally_->Log(2, "motor #%i added to bound_ATP_[%i]\n", head->motor_->id_,
                    n_bound_ATP_ - 1);
      } else if (head->ligand_ == "ADPP") {
        bound_ADPP_i_[n_bound_ADPP_i_++] = head;
        wally_->Log(2, "motor #%i added to bound_ADPP_i_[%i]\n",
                    head->motor_->id_, n_bound_ADPP_i_ - 1);
        Tubulin *dock_site{motor->GetDockSite()};
        if (!dock_site->occupied_) {
          docked_[n_docked_++] = motor->GetDockedHead();
          wally_->Log(2, "motor #%i added to docked_[%i]\n", head->motor_->id_,
                      n_docked_ - 1);
        }
      } else {
        wally_->ErrorExit("K_MGMT::UpdateLists() [singly-bound]");
      }
    } else if (motor->heads_active_ == 2) {
      POP_T *head_one{&motor->head_one_};
      POP_T *head_two{&motor->head_two_};
      bool found_NULL{false};
      if (head_one->ligand_ == "NULL") {
        bound_NULL_ii_[n_bound_NULL_ii_++] = head_one;
        found_NULL = true;
        wally_->Log(2, "motor #%i (head_one) added to bound_NULL_ii_[%i]\n",
                    head_one->motor_->id_, n_bound_NULL_ii_ - 1);
      } else if (head_two->ligand_ == "NULL") {
        if (found_NULL) {
          wally_->ErrorExit("[doubly-bound -- NULL]");
        }
        bound_NULL_ii_[n_bound_NULL_ii_++] = head_two;
        wally_->Log(2, "motor #%i (head_two) added to bound_NULL_ii_[%i]\n",
                    head_two->motor_->id_, n_bound_NULL_ii_ - 1);
      }
      bool found_ADPP{false};
      if (head_one->ligand_ == "ADPP") {
        bound_ADPP_ii_[n_bound_ADPP_ii_++] = head_one;
        wally_->Log(2, "motor #%i (head one) added to bound_ADPP_ii_[%i]\n",
                    head_one->motor_->id_, n_bound_ADPP_ii_ - 1);
        found_ADPP = true;
      } else if (head_two->ligand_ == "ADPP") {
        if (found_ADPP) {
          wally_->ErrorExit("[doubly-bound -- ADPP]");
        }
        bound_ADPP_ii_[n_bound_ADPP_ii_++] = head_two;
        wally_->Log(2, "motor #%i (head_two) added to bound_ADPP_ii_[%i]\n",
                    head_two->motor_->id_, n_bound_ADPP_ii_ - 1);
      }
    } else {
      wally_->ErrorExit("K_MGMT::UpdateLists()");
    }
  }
  wally_->Log(1, " -- finished --\n");
}

void KinesinManagement::SampleEventStatistics() {

  wally_->Log(1, "Starting K_MGMT::SampleEventStatistics()\n");
  // Scan through all events & get expected number of occurrences
  n_events_to_exe_ = 0;
  for (auto &&event : events_) {
    n_events_to_exe_ += event.SampleStatistics();
    // Ensure no funny business occurs with statistics
    if (event.n_expected_ > *event.n_avail_) {
      wally_->Log("Error; %i events expected but only %i available for %s\n",
                  event.n_expected_, *event.n_avail_, event.name_.c_str());
      wally_->ErrorExit("Kin_MGMT::SampleEventStatistics()");
    }
    if (wally_->verbosity_ >= 2 and event.n_expected_ > 0) {
      wally_->Log(" %i events expected for %s (%i avail)\n", event.n_expected_,
                  event.name_.c_str(), *event.n_avail_);
      for (int i_entry{0}; i_entry < event.n_expected_; i_entry++) {
        POP_T *head{nullptr};
        try {
          head = std::get<POP_T *>(event.targets_[i_entry]);
        } catch (...) {
          unsigned long index{event.targets_[i_entry].index()};
          wally_->Log("    -target %i: INDEX %lu ?!\n", i_entry, index);
          continue;
        }
        wally_->Log("    -target %i: motor %i\n", i_entry, head->motor_->id_);
      }
    }
  }
  if (n_events_to_exe_ <= 1) {
    return;
  }
  wally_->Log(2, " %i events expected pre-correction\n", n_events_to_exe_);
  // Put all events with >1 target into an active_ array
  std::pair<EVENT_T *, ENTRY_T> active_events[n_events_to_exe_];
  int i_active{0};
  for (auto &&event : events_) {
    // Add a ptr to the event for each target it has
    for (int i_tar{0}; i_tar < event.n_expected_; i_tar++) {
      active_events[i_active++] = std::make_pair(&event, event.targets_[i_tar]);
      wally_->Log(2, " Added %s to active_events\n", event.name_.c_str());
    }
  }
  // Scan through all active events to ensure that no two target the same
  // motor
  for (int i_entry{0}; i_entry < n_events_to_exe_; i_entry++) {
    EVENT_T *event_i{active_events[i_entry].first};
    wally_->Log(2, "event_i = %s\n", event_i->name_.c_str());
    ENTRY_T tar_i{active_events[i_entry].second};
    for (int j_entry{i_entry + 1}; j_entry < n_events_to_exe_; j_entry++) {
      EVENT_T *event_j{active_events[j_entry].first};
      wally_->Log(2, "   event_j = %s\n", event_j->name_.c_str());
      ENTRY_T tar_j{active_events[j_entry].second};
      bool motors_are_identical{false};
      if (tar_i.index() == 2 and tar_j.index() == 2) {
        if (std::get<2>(tar_i)->motor_ == std::get<2>(tar_j)->motor_) {
          motors_are_identical = true;
        }
      }
      // If event_i and event_j target different motors, continue
      if (tar_i == tar_j or motors_are_identical) {
        double p_one{event_i->p_occur_};
        double p_two{event_j->p_occur_};
        double ran{gsl_->GetRanProb()};
        if (ran < p_one / (p_one + p_two)) {
          wally_->Log(2, " Removed %s from active_events\n",
                      event_i->name_.c_str());
          event_i->RemoveTarget(active_events[i_entry].second);
          active_events[i_entry] = active_events[n_events_to_exe_ - 1];
          i_entry--;
          n_events_to_exe_--;
          break;
        } else {
          wally_->Log(2, " Removed %s from active_events\n",
                      event_j->name_.c_str());
          event_j->RemoveTarget(active_events[j_entry].second);
          active_events[j_entry] = active_events[n_events_to_exe_ - 1];
          j_entry--;
          n_events_to_exe_--;
        }
      }
    }
  }
  wally_->Log(1, " -- finished --\n");
}

void KinesinManagement::GenerateExecutionSequence() {

  if (n_events_to_exe_ == 0) {
    return;
  }
  wally_->Log(1, "Starting K_MGMT::GenerateExecutionSequence()\n");
  int i_array{0};
  EVENT_T *pre_array[n_events_to_exe_];
  for (auto &&event : events_) {
    for (int i_entry{0}; i_entry < event.n_expected_; i_entry++) {
      pre_array[i_array++] = &event;
    }
  }
  if (i_array != n_events_to_exe_) {
    wally_->ErrorExit("Kin_MGMT::GenerateExecutionSequence()");
  }
  if (n_events_to_exe_ > 1) {
    gsl_->Shuffle(pre_array, n_events_to_exe_, sizeof(EVENT_T *));
  }
  if (n_events_to_exe_ > events_to_exe_.size()) {
    events_to_exe_.resize(n_events_to_exe_);
  }
  for (int i_event{0}; i_event < n_events_to_exe_; i_event++) {
    events_to_exe_[i_event] = pre_array[i_event];
  }
  wally_->Log(1, " -- finished --\n");
}

void KinesinManagement::ExecuteEvents() {

  for (int i_event{0}; i_event < n_events_to_exe_; i_event++) {
    wally_->Log(1, "Executing %s\n", events_to_exe_[i_event]->name_.c_str());
    events_to_exe_[i_event]->Execute();
    wally_->Log(1, " -- finished -- \n");
    lists_up_to_date_ = false;
  }
}

void KinesinManagement::Bind_I(SITE_T *site) {

  if (site->occupied_) {
    wally_->Log("??? in K_MGMT::Bind_I()");
    return;
  }
  // Get random free motor
  Kinesin *motor{GetFreeMotor()};
  // Update site details
  site->motor_head_ = &motor->head_one_;
  site->occupied_ = true;
  // Update motor details
  motor->mt_ = site->mt_;
  motor->head_one_.site_ = site;
  motor->head_one_.ligand_ = "NULL";
  motor->head_one_.trailing_ = false;
  motor->head_two_.trailing_ = true;
  motor->heads_active_++;
  // Update active_ list
  AddToActive(motor);
  properties_->microtubules.FlagForUpdate();
}

void KinesinManagement::Bind_ATP(POP_T *head) {

  if (head->ligand_ != "NULL") {
    wally_->ErrorExit("K_MGMT::Bind_ATP()");
  }
  // Update motor head
  head->ligand_ = "ATP";
  // Do not change conformation of trailing heads (for end-pausing)
  if (!head->trailing_) {
    head->motor_->ChangeConformation();
  }
}

void KinesinManagement::Hydrolyze(POP_T *head) {

  if (head->ligand_ != "ATP") {
    wally_->ErrorExit("K_MGMT::Hydrolyze");
  }
  head->ligand_ = "ADPP";
}

void KinesinManagement::Bind_II(POP_T *head) {

  // Verify that proposed site is unoccupied
  Kinesin *motor{head->motor_};
  Tubulin *dock_site{motor->GetDockSite()};
  if (dock_site->occupied_) {
    return;
  }
  // Update site
  dock_site->motor_head_ = head;
  dock_site->occupied_ = true;
  // Update motor
  head->site_ = dock_site;
  head->ligand_ = "NULL";
  motor->heads_active_++;
  // Flag microtubules for update
  properties_->microtubules.FlagForUpdate();
}

void KinesinManagement::Unbind_II(POP_T *head) {

  Kinesin *motor{head->motor_};
  if (motor->heads_active_ != 2) {
    wally_->ErrorExit("Kin_MGMT::Unbind_II()");
  }
  // Update site
  head->site_->occupied_ = false;
  head->site_->motor_head_ = nullptr;
  // Update motor
  head->site_ = nullptr;
  head->ligand_ = "ADP";
  motor->heads_active_--;
  // Flag microtubules for update
  properties_->microtubules.FlagForUpdate();
}

void KinesinManagement::Unbind_I(POP_T *head) {

  Kinesin *motor{head->motor_};
  if (motor->heads_active_ != 1) {
    wally_->ErrorExit("Kin_MGMT::Unbind_I()");
  }
  // Update site
  head->site_->occupied_ = false;
  head->site_->motor_head_ = nullptr;
  // Update motor
  head->site_ = nullptr;
  head->ligand_ = "ADP";
  motor->heads_active_--;
  motor->mt_ = nullptr;
  RemoveFromActive(motor);
  // Flag microtubules for update
  properties_->microtubules.FlagForUpdate();
  // Update number of kinesin runs
  n_runs_executed_++;
  if (n_runs_executed_ >= n_runs_desired_) {
    n_events_to_exe_ = 0;
    properties_->sim_running_ = false;
    wally_->Log("Sim '%s' terminated after %i successful Kinesin runs\n",
                wally_->sim_name_, n_runs_executed_);
    size_t steps_recorded{properties_->current_step_ - wally_->data_threshold_};
    size_t n_datapoints{steps_recorded / wally_->n_steps_per_output_};
    wally_->Log("N_STEPS = %zu\n", properties_->current_step_);
    wally_->Log("N_DATAPOINTS = %zu\n", n_datapoints);
  }
}
