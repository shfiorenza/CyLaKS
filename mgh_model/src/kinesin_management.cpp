#include "kinesin_management.h"
#include "master_header.h"

KinesinManagement::KinesinManagement() {}

void KinesinManagement::Initialize(system_parameters *parameters,
                                   system_properties *properties) {
  // verbose_ = true;
  parameters_ = parameters;
  properties_ = properties;
  wally_ = &properties_->wallace;
  GenerateMotors();
  SetParameters();
  InitializeLists();
  InitializeEvents();
}

void KinesinManagement::GenerateMotors() {

  int n_mts = parameters_->microtubules.count;
  n_motors_ = 0;
  for (int i_mt = 0; i_mt < n_mts; i_mt++) {
    n_motors_ += parameters_->microtubules.length[i_mt];
  }
  // Since only one head has to be bound, the most that will ever
  // be needed (all single-bound) is the total number of sites
  motors_.resize(n_motors_);
  for (int id = 0; id < n_motors_; id++)
    motors_[id].Initialize(parameters_, properties_, id);
}

void KinesinManagement::SetParameters() {

  // Simulation constants
  double delta_t = parameters_->delta_t;
  double site_size = parameters_->microtubules.site_size;
  double kbT = parameters_->kbT;
  double r_0 = motors_[0].r_0_;
  double k_spring = motors_[0].k_spring_;
  double k_slack = motors_[0].k_slack_;
  double r_y = parameters_->microtubules.y_dist / 2;

  // lattice cooperativity jawn - temp
  n_affinities_ = 1;
  std::vector<int> affinity_weights(n_affinities_, 0);
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    double weight = exp((double)i_aff / 5);
    affinity_weights[i_aff] = weight;
  }

  // Non-tethered statistics
  double k_on = parameters_->motors.k_on;
  double c_motor = parameters_->motors.c_bulk;
  double p_bind_i_base = k_on * c_motor * delta_t;
  double k_on_ATP = parameters_->motors.k_on_ATP;
  double c_ATP = parameters_->motors.c_ATP;
  p_bind_ATP_ = k_on_ATP * c_ATP * delta_t;
  double k_hydrolyze = parameters_->motors.k_hydrolyze;
  p_hydrolyze_ = k_hydrolyze * delta_t;
  double k_hydrolyze_stalled = parameters_->motors.k_hydrolyze_stalled;
  p_hydrolyze_stalled_ = k_hydrolyze_stalled * delta_t;
  double c_eff = parameters_->motors.c_eff_bind;
  double p_bind_ii_base = k_on * c_eff * delta_t;
  double k_off_ii = parameters_->motors.k_off_ii;
  double p_unbind_ii_base = k_off_ii * delta_t;
  double k_off_i = parameters_->motors.k_off_i;
  double p_unbind_i_base = k_off_i * delta_t;
  double k_off_i_st = parameters_->motors.k_off_i_stalled;
  double p_unbind_i_stalled_base = k_off_i_st * delta_t;
  // Sound the alarm if our timestep is too large
  if (p_bind_ATP_ > 1)
    printf("WARNING: p_bind_ATP=%g for motors\n", p_bind_ATP_);
  if (p_hydrolyze_ > 1)
    printf("WARNING: p_phos=%g for motors\n", p_hydrolyze_);
  /* artifical force perpetually applied to motors */
  double app_force = parameters_->motors.applied_force;
  if (app_force > 0) {
    double sigma_off = 1.2;
    double weight = exp(app_force * sigma_off / kbT);
    p_unbind_i_base = k_off_i * weight * delta_t;
    printf("p_unbind_i_ scaled from %g to %g \n", k_off_i * delta_t,
           k_off_i * weight * delta_t);
  }
  p_bind_i_.resize(n_affinities_);
  p_bind_ii_.resize(n_affinities_);
  p_unbind_ii_.resize(n_affinities_);
  p_unbind_i_.resize(n_affinities_);
  p_unbind_i_stalled_.resize(n_affinities_);
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    p_bind_i_[i_aff] = (p_bind_i_base * affinity_weights[i_aff]);
    p_bind_ii_[i_aff] = (p_bind_ii_base * affinity_weights[i_aff]);
    p_unbind_ii_[i_aff] = (p_unbind_ii_base / affinity_weights[i_aff]);
    p_unbind_i_[i_aff] = (p_unbind_i_base / affinity_weights[i_aff]);
    p_unbind_i_stalled_[i_aff] =
        (p_unbind_i_stalled_base / affinity_weights[i_aff]);
    if (p_bind_i_[i_aff] > 1)
      printf("WARNING: p_bind_i=%g for motors\n", p_bind_i_[i_aff]);
    if (p_bind_ii_[i_aff] > 1)
      printf("WARNING: p_bind_ii=%g for motors\n", p_bind_ii_[i_aff]);
    if (p_unbind_ii_[i_aff] > 1)
      printf("WARNING: p_unbind_ii=%g for motors\n", p_unbind_ii_[i_aff]);
    if (p_unbind_i_[i_aff] > 1)
      printf("WARNING: p_unbind_i=%g for motors\n", p_unbind_i_[i_aff]);
    if (p_unbind_i_stalled_[i_aff] > 1)
      printf("WARNING: p_unbind_i_stalled=%g for motors\n",
             p_unbind_i_stalled_[i_aff]);
  }

  // Get tether information from motors
  rest_dist_ = motors_[0].rest_dist_;
  comp_cutoff_ = motors_[0].comp_cutoff_;
  dist_cutoff_ = motors_[0].dist_cutoff_;
  if (parameters_->motors.tethers_active) {
    properties_->wallace.Log("\nFor motors:\n");
    properties_->wallace.Log("  rest_dist is %g\n", rest_dist_);
    properties_->wallace.Log("  comp_cutoff is %i\n", comp_cutoff_);
    properties_->wallace.Log("  dist_cutoff is %i\n", dist_cutoff_);
  }
  // Tethered statistics
  double k_tether = parameters_->motors.k_tether;
  double k_untether = parameters_->motors.k_untether;
  double c_eff_teth = parameters_->motors.c_eff_tether;
  if (!parameters_->motors.tethers_active) {
    k_tether = 0;
    k_untether = 0;
    c_eff_teth = 0;
  }
  p_bind_i_tethered_ = k_on * c_eff_teth * delta_t;
  p_tether_free_ = k_tether * c_motor * delta_t;
  p_tether_bound_ = k_tether * c_eff_teth * delta_t;
  p_untether_free_ = k_untether * delta_t;
  // Sound the alarm if our timestep is too large
  if (p_bind_i_tethered_ > 1)
    printf("WARNING: p_bind_i_teth=%g for mots\n", p_bind_i_tethered_);
  if (p_tether_free_ > 1)
    printf("WARNING: p_teth_free=%g for mots\n", p_tether_free_);
  if (p_tether_bound_ > 1)
    printf("WARNING: p_teth_bound=%g for mots\n", p_tether_bound_);
  if (p_untether_free_ > 1)
    printf("WARNING: p_unteth_free=%g for mots\n", p_untether_free_);
  // For events that depend on tether stretch, each different extension
  // has its own rate; "base" refers to when the tether is unstretched
  double p_bind_base = k_on * c_eff_teth * delta_t;
  double p_unbind_base = p_unbind_i_base;
  double p_unbind_st_base = p_unbind_i_stalled_base;
  double p_teth_base = k_tether * c_eff_teth * delta_t;
  double p_unteth_base = p_untether_free_;
  p_bind_ATP_tethered_.resize(2 * dist_cutoff_ + 1);
  p_bind_ii_tethered_.resize(2 * dist_cutoff_ + 1);
  p_unbind_i_tethered_.resize(2 * dist_cutoff_ + 1);
  p_unbind_i_teth_st_.resize(2 * dist_cutoff_ + 1);
  p_untether_bound_.resize(2 * dist_cutoff_ + 1);
  for (int x_dub = 0; x_dub <= 2 * dist_cutoff_; x_dub++) {
    // Calculate tether length for this x_dist
    double r_x = (double)x_dub * site_size / 2;
    double r = sqrt(r_y * r_y + r_x * r_x);
    double cosine = r_x / r;
    // Calculate extension of tether for given x_dub
    double dr = r - r_0;
    // Get appropriate spring constant
    double k = 0;
    if (dr < 0)
      k = k_slack;
    else
      k = k_spring;
    // Calculate spring force and potential energy of this extension
    double f_x = fabs(k * dr * cosine);
    double U_teth = k * dr * dr / 2;
    // Calculate weight
    double weight = exp(U_teth / (2 * kbT));
    double sigma_off = 1.5;
    double weight_alt = exp(f_x * sigma_off / kbT);
    // If tethering is disabled, all weights are automatically zero
    if (!parameters_->motors.tethers_active) {
      weight = 0;
      weight_alt = 0;
    }
    // Otherwise, only weights below comp_cutoff_ are zero
    else if (x_dub < 2 * comp_cutoff_) {
      weight = 0;
      weight_alt = 0;
    }
    // Calculate appropriately-weighted probabilities
    p_bind_ATP_tethered_[x_dub] = p_bind_ATP_;
    p_bind_ii_tethered_[x_dub] = p_bind_ii_base;
    p_unbind_i_tethered_[x_dub] = weight_alt * p_unbind_base;
    p_unbind_i_teth_st_[x_dub] = weight_alt * p_unbind_st_base;
    p_untether_bound_[x_dub] = weight * p_unteth_base;
    // Sound the alarm if our timestep is too large
    if (p_bind_ATP_tethered_[x_dub] > 1)
      printf("WaRNING: p_bind_ATP_teth=%g for 2x=%i\n",
             p_bind_ATP_tethered_[x_dub], x_dub);
    if (p_bind_ii_tethered_[x_dub] > 1)
      printf("WARNING: p_bind_ii_teth=%g for 2x=%i\n",
             p_bind_ii_tethered_[x_dub], x_dub);
    if (p_unbind_i_tethered_[x_dub] > 1)
      printf("WARNING: p_unbind_i_teth=%g for 2x=%i\n",
             p_unbind_i_tethered_[x_dub], x_dub);
    if (p_unbind_i_teth_st_[x_dub] > 1)
      printf("WARNING: p_unbind_i_teth_st=%g for 2x=%i\n",
             p_unbind_i_teth_st_[x_dub], x_dub);
    if (p_untether_bound_[x_dub] > 1)
      printf("WARNING: p_unteth_bound=%g for 2x=%i\n", p_untether_bound_[x_dub],
             x_dub);
  }
}

void KinesinManagement::InitializeLists() {

  scratch_.resize(n_motors_);
  // One dimensional stuff
  active_.resize(n_motors_);
  free_tethered_.resize(n_motors_);
  bound_untethered_.resize(n_motors_);
  bound_NULL_.resize(n_motors_);
  bound_ATP_.resize(n_motors_);
  bound_ATP_stalled_.resize(n_motors_);
  // Two dimensional stuff - coop
  n_docked_.resize(n_affinities_);
  n_bound_ADPP_i_.resize(n_affinities_);
  n_bound_ADPP_i_stalled_.resize(n_affinities_);
  n_bound_ADPP_ii_.resize(n_affinities_);
  docked_.resize(n_affinities_);
  bound_ADPP_i_.resize(n_affinities_);
  bound_ADPP_i_stalled_.resize(n_affinities_);
  bound_ADPP_ii_.resize(n_affinities_);
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    n_docked_[i_aff] = 0;
    n_bound_ADPP_i_[i_aff] = 0;
    n_bound_ADPP_i_stalled_[i_aff] = 0;
    n_bound_ADPP_ii_[i_aff] = 0;
    docked_[i_aff].resize(n_motors_);
    bound_ADPP_i_[i_aff].resize(n_motors_);
    bound_ADPP_i_stalled_[i_aff].resize(n_motors_);
    bound_ADPP_ii_[i_aff].resize(n_motors_);
  }
  // Two dimensional stuff - tethers
  n_docked_tethered_.resize(2 * dist_cutoff_ + 1);
  n_bound_NULL_tethered_.resize(2 * dist_cutoff_ + 1);
  n_bound_ADPP_i_tethered_.resize(2 * dist_cutoff_ + 1);
  n_bound_ADPP_i_teth_st_.resize(2 * dist_cutoff_ + 1);
  n_bound_tethered_.resize(2 * dist_cutoff_ + 1);
  docked_tethered_.resize(2 * dist_cutoff_ + 1);
  bound_NULL_tethered_.resize(2 * dist_cutoff_ + 1);
  bound_ADPP_i_tethered_.resize(2 * dist_cutoff_ + 1);
  bound_ADPP_i_teth_st_.resize(2 * dist_cutoff_ + 1);
  bound_tethered_.resize(2 * dist_cutoff_ + 1);
  for (int x_dub = 0; x_dub <= 2 * dist_cutoff_; x_dub++) {
    n_docked_tethered_[x_dub] = 0;
    n_bound_NULL_tethered_[x_dub] = 0;
    n_bound_ADPP_i_tethered_[x_dub] = 0;
    n_bound_ADPP_i_teth_st_[x_dub] = 0;
    n_bound_tethered_[x_dub] = 0;
    docked_tethered_[x_dub].resize(n_motors_);
    bound_NULL_tethered_[x_dub].resize(n_motors_);
    bound_ADPP_i_tethered_[x_dub].resize(n_motors_);
    bound_ADPP_i_teth_st_[x_dub].resize(n_motors_);
    bound_tethered_[x_dub].resize(n_motors_);
  }
}

void KinesinManagement::InitializeEvents() {

  /* * Serialized & unique index of each KMC event * */
  int id{0};
  /* * Basic random integer generator; passed to all events * */
  auto ran_int = [&](int n) {
    if (n > 0) {
      return properties_->gsl.GetRanInt(n);
    } else {
      return 0;
    }
  };
  /* * Binomial probabilitiy distribution; sampled to predict most events * */
  auto binomial = [&](double p, int n) {
    if (n > 0) {
      return properties_->gsl.SampleBinomialDist(p, n);
    } else {
      return 0;
    }
  };
  /* *** Event entries *** */
  // Bind_I: bind first motor head to MT and release ADP
  auto update_unnoc = [&]() { properties_->microtubules.UpdateUnoccupied(); };
  auto exe_bind_i = [&](ENTRY_T target) {
    SITE_T *site = std::get<SITE_T *>(target);
    if (site != nullptr) {
      ReportExecutionOf("bind_i");
      KMC_Bind_I(site);
    } else {
      ReportFailureOf("bind_i");
    }
  };
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    std::string affinity = std::to_string(i_aff);
    std::string name = "bind_i_" + affinity;
    std::string tar = "unnoc_" + affinity;
    events_.emplace_back(this, id++, name, Vec<Str>{tar}, p_bind_i_[i_aff],
                         &properties_->microtubules.n_unoccupied_mot_[i_aff],
                         &properties_->microtubules.unoccupied_list_mot_[i_aff],
                         update_unnoc, exe_bind_i, binomial, ran_int);
  }
  // Bind_ATP: bind ATP to motor heads that have released their ADP
  auto update_NULL = [&]() { Update_Bound_NULL(); };
  auto exe_bind_ATP = [&](ENTRY_T target) {
    POP_T *head = std::get<POP_T *>(target);
    if (head != nullptr) {
      ReportExecutionOf("bind_ATP");
      KMC_Bind_ATP(head);
    } else {
      ReportFailureOf("bind_ATP");
    }
  };
  events_.emplace_back(this, id++, "bind_ATP", Vec<Str>{"bound_NULL"},
                       p_bind_ATP_, &n_bound_NULL_, &bound_NULL_, update_NULL,
                       exe_bind_ATP, binomial, ran_int);
  // Hydrolyze: convert ATP to ADPP
  auto update_ATP = [&]() { Update_Bound_ATP(); };
  auto exe_hydrolyze = [&](ENTRY_T target) {
    POP_T *head = std::get<POP_T *>(target);
    if (head != nullptr) {
      ReportExecutionOf("hydrolyze");
      KMC_Hydrolyze(head);
    } else {
      ReportFailureOf("hydrolyze");
    }
  };
  events_.emplace_back(this, id++, "hydrolyze", Vec<Str>{"bound_ATP"},
                       p_hydrolyze_, &n_bound_ATP_, &bound_ATP_, update_ATP,
                       exe_hydrolyze, binomial, ran_int);
  // Hydrolyze_stalled: same as Hydrolyze but only targets stalled motors
  auto update_ATP_st = [&]() { Update_Bound_ATP_Stalled(); };
  events_.emplace_back(this, id++, "hydrolyze_st", Vec<Str>{"bound_ATP_st"},
                       p_hydrolyze_stalled_, &n_bound_ATP_stalled_,
                       &bound_ATP_stalled_, update_ATP_st, exe_hydrolyze,
                       binomial, ran_int);
  // Bind_II: binds docked motor head to MT and releases its ADP
  auto update_docked = [&]() { Update_Docked(); };
  auto exe_bind_ii = [&](ENTRY_T target) {
    POP_T *head = std::get<POP_T *>(target);
    if (head != nullptr) {
      ReportExecutionOf("bind_ii");
      KMC_Bind_II(head);
    } else {
      ReportFailureOf("bind_ii");
    }
  };
  /*
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    std::string affinity = std::to_string(i_aff);
    std::string name = "bind_ii_" + affinity;
    std::string tar1 = "docked_" + affinity;
    std::string tar2 = "bound_ADPP_i_" + affinity;
    std::string tar3 = "bound_ADPP_i_st_" + affinity;
    events_.emplace_back(this, id++, name, Vec<Str>{tar1, tar2, tar3}, false,
                         p_bind_ii_[i_aff], &n_docked_[i_aff], &docked_[i_aff],
                         update_docked, exe_bind_ii, binomial, ran_int);
  }
  */
  // Unbind_II: Converts ADPP to ADP and unbinds a doubly-bound head
  auto update_ADPP_ii = [&]() { Update_Bound_ADPP_II(); };
  auto exe_unbind_ii = [&](ENTRY_T target) {
    POP_T *head = std::get<POP_T *>(target);
    if (head != nullptr) {
      ReportExecutionOf("unbind_ii");
      KMC_Unbind_II(head);
    } else {
      ReportFailureOf("unbind_ii");
    }
  };
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    std::string affinity = std::to_string(i_aff);
    std::string name = "unbind_ii_" + affinity;
    std::string tar = "bound_ADPP_ii_" + affinity;
    events_.emplace_back(this, id++, name, Vec<Str>{tar}, p_unbind_ii_[i_aff],
                         &n_bound_ADPP_ii_[i_aff], &bound_ADPP_ii_[i_aff],
                         update_ADPP_ii, exe_unbind_ii, binomial, ran_int);
  }
  // Unbind_I: Converts ADPP to ADP and unbinds a singly-bound head
  auto update_ADPP_i = [&]() { Update_Bound_ADPP_I(); };
  auto exe_unbind_i = [&](ENTRY_T target) {
    POP_T *head = std::get<POP_T *>(target);
    if (head != nullptr) {
      ReportExecutionOf("unbind_i");
      KMC_Unbind_I(head);
    } else {
      ReportFailureOf("unbind_i");
    }
  };
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    std::string affinity = std::to_string(i_aff);
    std::string name = "unbind_i_" + affinity;
    std::string tar = "bound_ADPP_i_" + affinity;
    events_.emplace_back(this, id++, name, Vec<Str>{tar}, p_unbind_i_[i_aff],
                         &n_bound_ADPP_i_[i_aff], &bound_ADPP_i_[i_aff],
                         update_ADPP_i, exe_unbind_i, binomial, ran_int);
  }
  // Unbind_I_Stalled: Same as unbind_I but only targets stalled motors
  auto update_ADPP_i_st = [&]() { Update_Bound_ADPP_I_Stalled(); };
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    std::string affinity = std::to_string(i_aff);
    std::string name = "unbind_i_st_" + affinity;
    std::string tar = "bound_ADPP_i_st_" + affinity;
    events_.emplace_back(
        this, id++, name, Vec<Str>{tar}, p_unbind_i_stalled_[i_aff],
        &n_bound_ADPP_i_stalled_[i_aff], &bound_ADPP_i_stalled_[i_aff],
        update_ADPP_i_st, exe_unbind_i, binomial, ran_int);
  }
  /*
  // If tethers ARE active, serialize all extension-based events
  if (parameters_->motors.tethers_active) {
  // Poisson dist. is used w/ partition function for E-dependent binding
  auto poisson_bind = [&](double p, int n) {
    if (n > 0) {
      double n_avg = p * GetWeight_BindTethered();
      return properties_->gsl.SamplePoissonDist(n_avg);
    } else
      return 0;
  };
  // Poisson dist. is used w/ partition function for E-dependent binding
  auto poisson_teth = [&](double p, int n) {
    if (n > 0) {
      double n_avg = p * GetWeight_TetherBound();
      return properties_->gsl.SamplePoissonDist(n_avg);
    } else
      return 0;
  };
  events_.emplace_back(event(index++, 11, "bind_i_teth", "free_teth",
                             poisson_bind, &n_free_tethered_,
                             p_bind_i_tethered_));
  for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * dist_cutoff_; x_dub++) {
    events_.emplace_back(event(
        index++, 200 + x_dub, "bind_ATP_teth", "bound_NULL_teth", binomial,
        &n_bound_NULL_tethered_[x_dub], p_bind_ATP_tethered_[x_dub]));
    }
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * dist_cutoff_; x_dub++) {
      events_.emplace_back(event(
          index++, 400 + x_dub, "bind_ii_teth", "bound_ADPP_i_teth",
  binomial, &n_docked_tethered_[x_dub], p_bind_ii_tethered_[x_dub]));
    }
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * dist_cutoff_; x_dub++) {
      events_.emplace_back(event(
          index++, 600 + x_dub, "unbind_i_teth", "bound_ADPP_i_teth",
  binomial, &n_bound_ADPP_i_tethered_[x_dub], p_unbind_i_tethered_[x_dub]));
    }
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * dist_cutoff_; x_dub++) {
      events_.emplace_back(event(index++, 700 + x_dub, "unbind_i_teth_st",
                                 "bound_ADPP_i_teth_st", binomial,
                                 &n_bound_ADPP_i_teth_st_[x_dub],
                                 p_unbind_i_teth_st_[x_dub]));
    }
    events_.emplace_back(event(index++, 70, "tether_free",
  "untethered_xlinks", binomial, &properties_->prc1.n_bound_unteth_,
                               p_tether_free_));
    events_.emplace_back(event(index++, 71, "tether_bound",
  "untethered_xlinks", poisson_teth, &n_bound_untethered_,
                               p_tether_bound_));
    events_.emplace_back(event(index++, 80, "untether_free",
  "free_tethered", binomial, &n_free_tethered_, p_untether_free_)); for (int
  x_dub(2 * comp_cutoff_); x_dub <= 2 * dist_cutoff_; x_dub++) {
      events_.emplace_back(
          event(index++, 800 + x_dub, "untether_bound", "bound_tethered",
                binomial, &n_bound_tethered_[x_dub],
  p_untether_bound_[x_dub]));
    }
  }
  */
  /* ** Segregate events_ into events_by_pop_ based on target pop. ** */
  // Scan over all events in simulation; segregate based on MAIN target first
  for (int i_event{0}; i_event < events_.size(); i_event++) {
    std::string main_target = events_[i_event].targets_[0];
    // Assume target population has not been seen before
    bool new_pop{true};
    // If assumption is wrong, store index of recorded pop
    int pop_index{0};
    // Scan over all previously-seen populations to check our assumption
    for (int i_pop{0}; i_pop < events_by_pop_.size(); i_pop++) {
      std::string recorded_target = events_by_pop_[i_pop][0]->targets_[0];
      // If target pop matches record, this population has been seen before
      if (main_target == recorded_target) {
        new_pop = false;
        pop_index = i_pop;
        break;
      }
    }
    // If this is indeed a new pop, start a new row in entries array
    if (new_pop) {
      std::vector<EVENT_T *> new_row = {&events_[i_event]};
      events_by_pop_.push_back(new_row);
    }
    // Otherwise, add this event to row of previously-recorded population
    else {
      events_by_pop_[pop_index].push_back(&events_[i_event]);
    }
  }
  int n_primary_pops = events_by_pop_.size();
  // Scan over all events again; segregate based on SECONDARY targets
  for (int i_event{0}; i_event < events_.size(); i_event++) {
    if (events_[i_event].targets_.size() > 1) {
      std::vector<EVENT_T *> new_row = {&events_[i_event]};
      for (int i_tar{1}; i_tar < events_[i_event].targets_.size(); i_tar++) {
        // Get name of secondary target
        std::string secondary_target = events_[i_event].targets_[i_tar];
        // Create new row for this population; will be a secondary
        // stat-correction
        // Scan over events_by_pop_, add all entries whose target matches this
        for (int i_pop{0}; i_pop < n_primary_pops; i_pop++) {
          std::string this_target = events_by_pop_[i_pop][0]->targets_[0];
          if (this_target == secondary_target) {
            for (int i_entry{0}; i_entry < events_by_pop_[i_pop].size();
                 i_entry++) {
              new_row.push_back(events_by_pop_[i_pop][i_entry]);
            }
          }
        }
      }
      events_by_pop_.push_back(new_row);
    }
  }
  /*
  // Transfer info to events_by_pop_ structure
  events_by_pop_.resize(n_distinct_pops);
  for (int i_pop{0}; i_pop < n_distinct_pops; i_pop++) {
    std::string main_target = entries[i_pop][0]->targets_[0];
    // 2nd index corresponds to entry of events that target this pop.
    events_by_pop_[i_pop].resize(n_entries[i_pop]);
    for (int i_entry{0}; i_entry < n_entries[i_pop]; i_entry++) {
      events_by_pop_[i_pop][i_entry] = entries[i_pop][i_entry];
    }
  }
  // Total number of unique general target populations in simulation
  int n_gen_pops = 0;
  // Number of events that target each distinct population
  int n_gen_entries[events_.size()];
  // id of each event for each distinct population:
  int gen_entry_IDs[events_.size()]
                   [2 * (2 * dist_cutoff_ + 1) * (dist_cutoff_ + 1)];
  // Scan over all events in simulation
  for (int i_entry(0); i_entry < events_.size(); i_entry++) {
    std::string gen_target = events_[i_entry].general_target_;
    if (gen_target != "NULL") {
      // Assume target population has not been seen before
      bool new_pop{true};
      // If assumption is wrong, store index of recorded pop
      int pop_index{0};
      // Scan over all previously-seen populations to check our assumption
      for (int i_pop(0); i_pop < n_gen_pops; i_pop++) {
        // Get id of first-stored event that targets this population
        int recorded_ID = gen_entry_IDs[i_pop][0];
        // Use this id to get name of recorded population
        std::string record = events_[gen_entry_IDs[i_pop][0]].general_target_;
        // If target pop matches record, this population has been seen before
        if (gen_target == record) {
          new_pop = false;
          pop_index = i_pop;
          break;
        }
      }
      // If this is indeed a new pop, start a new row in entry_IDs array
      if (new_pop) {
        gen_entry_IDs[n_gen_pops][0] = events_[i_entry].id_;
        n_gen_entries[n_gen_pops] = 1;
        n_gen_pops++;
      }
      // Otherwise, add this event to row of previously-recorded population
      else {
        int entry_no = n_gen_entries[pop_index];
        gen_entry_IDs[pop_index][entry_no] = events_[i_entry].id_;
        n_gen_entries[pop_index]++;
      }
    }
  }
  // Transfer information to IDs_by_pop_ structure
  IDs_by_pop_.resize(n_distinct_pops + n_gen_pops);
  for (int i_gen(0); i_gen < n_gen_pops; i_gen++) {
    std::string target_pop = events_[gen_entry_IDs[i_gen][0]].general_target_;
    printf("KMC GENERAL Pop #%i - ", i_gen);
    std::cout << target_pop;
    printf(":");
    int i_pop = i_gen + n_distinct_pops;
    // 2nd index corresponds to entry of events that target this pop.
    IDs_by_pop_[i_pop].resize(n_gen_entries[i_gen]);
    for (int i_entry(0); i_entry < n_gen_entries[i_gen]; i_entry++) {
      IDs_by_pop_[i_pop][i_entry] = gen_entry_IDs[i_gen][i_entry];
      printf(" %i,", gen_entry_IDs[i_gen][i_entry]);
    }
    printf("\n");
    // Ensure that the event with general_tar == specific_tar is first
    for (int i_entry(0); i_entry < n_gen_entries[i_gen]; i_entry++) {
      std::string spec_tar =
          events_[gen_entry_IDs[i_gen][i_entry]].specific_target_;
      if (spec_tar == target_pop) {
        int first_id = IDs_by_pop_[i_pop][0];
        int this_id = IDs_by_pop_[i_pop][i_entry];
        IDs_by_pop_[i_pop][0] = this_id;
        IDs_by_pop_[i_pop][i_entry] = first_id;
        printf("Swapped position of ID %i and %i in IDs_by_pop\n", first_id,
               this_id);
      }
    }
  }
  */
  for (int i_pop(0); i_pop < events_by_pop_.size(); i_pop++) {
    printf("KMC event partition #%i (size %lu): ", i_pop,
           events_by_pop_[i_pop].size());
    for (int i_entry{0}; i_entry < events_by_pop_[i_pop].size(); i_entry++) {
      std::cout << events_by_pop_[i_pop][i_entry]->name_ << " (targets "
                << events_by_pop_[i_pop][i_entry]->targets_[0] << "), ";
    }
    printf("\n");
  }
}

Kinesin *KinesinManagement::GetFreeMotor() {

  // Randomly pick a motor from the reservoir
  int i_motor = properties_->gsl.GetRanInt(n_motors_);
  Kinesin *motor = &motors_[i_motor];
  int attempts = 0;
  while (motor->heads_active_ > 0 || motor->tethered_ == true) {
    i_motor++;
    if (i_motor == n_motors_)
      i_motor = 0;
    motor = &motors_[i_motor];
    attempts++;
    if (attempts > n_motors_) {
      printf("error in get_free_motor\n");
      exit(1);
    }
  }
  return motor;
}

/*
Kinesin *KinesinManagement::GetBoundUntetheredMotor() {

  Update_Bound_Unteth();
  int i_motor = properties_->gsl.GetRanInt(n_bound_untethered_);
  return std::get<POP_T *>(bound_untethered_[i_motor])->motor_;
}

int KinesinManagement::GetNumBoundUntethered() {

  Update_Bound_Unteth();
  return n_bound_untethered_;
}
*/

/*
void KinesinManagement::Update_Relay(std::string event, std::string taret_pop) {

  if (event == std::string{"bind_i"}) {
    properties_->microtubules.UpdateUnoccupied();
  } else if (event == std::string{"bind_ATP"}) {
    Update_Bound_NULL();
  } else if (event == std::string{"hydrolyze"}) {
    Update_Bound_ATP();
  } else if (event == std::string{"hydrolyze_st"}) {
    Update_Bound_ATP_Stalled();
  } else if (event == std::string{"bind_ii"}) {
    properties_->microtubules.UpdateUnoccupied();
    Update_Docked();
  } else if (event == std::string{"unbind_ii"}) {
    Update_Bound_ADPP_II();
  } else if (event == std::string{"unbind_i"}) {
    Update_Bound_ADPP_I();
  } else if (event == std::string{"unbind_i_st"}) {
    Update_Bound_ADPP_I_Stalled();
  }
}
*/

void KinesinManagement::Update_All_Lists() {

  if (verbose_) {
    printf("starting UPDATE ALL\n");
  }
  properties_->microtubules.UpdateUnoccupied();
  Update_Docked();
  Update_Bound_NULL();
  Update_Bound_ATP();
  Update_Bound_ATP_Stalled();
  Update_Bound_ADPP_I();
  Update_Bound_ADPP_I_Stalled();
  Update_Bound_ADPP_II();
  if (parameters_->motors.tethers_active) {
    properties_->prc1.Update_Bound_Unteth();
    Update_Free_Teth();
    Update_Docked_Teth();
    Update_Bound_NULL_Teth();
    Update_Bound_ADPP_I_Tethered();
    Update_Bound_ADPP_I_Tethered_Stalled();
    Update_Bound_Teth();
    Update_Bound_Unteth();
  }
}

void KinesinManagement::Update_Free_Teth() {

  /*
    n_free_tethered_ = 0;
    for (int i_motor = 0; i_motor < n_active_; i_motor++) {
      Kinesin *motor = active_[i_motor];
      if (motor->heads_active_ == 0 && motor->tethered_ == true) {
        free_tethered_[n_free_tethered_] = &motor->head_one_;
        n_free_tethered_++;
      }
    }
    */
}

void KinesinManagement::Update_Docked() {

  if (verbose_) {
    printf("starting UPDATE DOCKED\n");
  }
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    n_docked_[i_aff] = 0;
  }
  for (int i_motor = 0; i_motor < n_active_; i_motor++) {
    Kinesin *motor = active_[i_motor];
    bool eligible{false};
    if (motor->heads_active_ == 1) { // and !motor->is_outdated_) {
      if (!motor->tethered_) {
        eligible = true;
      }
      // Motors with satellite xlinks behave as if untethered
      else if (motor->xlink_->heads_active_ == 0) {
        eligible = true;
      }
    }
    if (eligible) {
      if (motor->GetActiveHead()->ligand_ == "ADPP") {
        double site_coord = motor->GetDockedCoordinate();
        int i_site = site_coord - motor->mt_->coord_;
        if (i_site >= 0 and i_site <= motor->mt_->n_sites_ - 1) {
          // Ensure site isn't occupied
          if (!motor->mt_->lattice_[i_site].occupied_) {
            int aff = motor->mt_->lattice_[i_site].affinity_;
            POP_T *docked_head = motor->GetDockedHead();
            docked_[aff][n_docked_[aff]] = docked_head;
            n_docked_[aff]++;
          }
        }
      }
    }
  }
}

void KinesinManagement::Update_Docked_Teth() {

  /*
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * dist_cutoff_; x_dub++) {
      n_docked_tethered_[x_dub] = 0;
    }
    for (int i_motor = 0; i_motor < n_active_; i_motor++) {
      Kinesin *motor = active_[i_motor];
      if (motor->heads_active_ == 1 && motor->tethered_) {
        // Don't count motors w/ just satellite xlinks
        if (motor->xlink_->heads_active_ > 0) {
          motor->UpdateExtension();
          // Make sure we don't force an untether event
          if (motor->tethered_) {
            if (motor->GetActiveHead()->ligand_ == "ADPP") {
              double site_coord = motor->GetDockedCoordinate();
              int i_site = site_coord - motor->mt_->coord_;
              if (i_site >= 0 && i_site <= motor->mt_->n_sites_ - 1) {
                if (!motor->mt_->lattice_[i_site].occupied_) {
                  int x_dub = motor->x_dist_doubled_;
                  int i = n_docked_tethered_[x_dub];
                  docked_tethered_[x_dub][i] = motor->GetDockedHead();
                  n_docked_tethered_[x_dub]++;
                }
              }
            }
          }
        }
      }
    }
    */
}

void KinesinManagement::Update_Bound_NULL() {

  n_bound_NULL_ = 0;
  for (int i_motor = 0; i_motor < n_active_; i_motor++) {
    Kinesin *motor = active_[i_motor];
    bool eligible{false};
    //    if (!motor->is_outdated_) {
    if (!motor->tethered_) {
      eligible = true;
    } else if (motor->xlink_->heads_active_ == 0) {
      eligible = true;
    }
    //    }
    if (eligible) {
      if (motor->head_one_.site_ != nullptr and
          motor->head_one_.ligand_ == "NULL") {
        bound_NULL_[n_bound_NULL_] = &motor->head_one_;
        n_bound_NULL_++;
      }
      if (motor->head_two_.site_ != nullptr and
          motor->head_two_.ligand_ == "NULL") {
        bound_NULL_[n_bound_NULL_] = &motor->head_two_;
        n_bound_NULL_++;
      }
    }
  }
}

void KinesinManagement::Update_Bound_NULL_Teth() {

  /*
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * dist_cutoff_; x_dub++) {
      n_bound_NULL_tethered_[x_dub] = 0;
    }
    for (int i_motor = 0; i_motor < n_active_; i_motor++) {
      Kinesin *motor = active_[i_motor];
      if (motor->tethered_) {
        if (motor->xlink_->heads_active_ > 0) {
          motor->UpdateExtension();
          // Make sure we don't force an untether event
          if (motor->tethered_) {
            int x_dub = motor->x_dist_doubled_;
            int i = n_bound_NULL_tethered_[x_dub];
            bool counted = false;
            if (motor->head_one_.site_ != nullptr &&
                motor->head_one_.ligand_ == "NULL") {
              bound_NULL_tethered_[x_dub][i] = &motor->head_one_;
              n_bound_NULL_tethered_[x_dub]++;
              counted = true;
            }
            if (motor->head_two_.site_ != nullptr &&
                motor->head_two_.ligand_ == "NULL") {
              bound_NULL_tethered_[x_dub][i] = &motor->head_two_;
              n_bound_NULL_tethered_[x_dub]++;
              if (counted) {
                printf("why - Update_Bound_NULL\n");
                exit(1);
              }
            }
          }
        }
      }
    }
    */
}

void KinesinManagement::Update_Bound_ATP() {

  n_bound_ATP_ = 0;
  for (int i_motor = 0; i_motor < n_active_; i_motor++) {
    Kinesin *motor = active_[i_motor];
    if (!motor->IsStalled()) { // and !motor->is_outdated_) {
      if (motor->head_one_.site_ != nullptr and
          motor->head_one_.ligand_ == "ATP") {
        bound_ATP_[n_bound_ATP_] = &motor->head_one_;
        n_bound_ATP_++;
      }
      if (motor->head_two_.site_ != nullptr and
          motor->head_two_.ligand_ == "ATP") {
        bound_ATP_[n_bound_ATP_] = &motor->head_two_;
        n_bound_ATP_++;
      }
    }
  }
}

void KinesinManagement::Update_Bound_ATP_Stalled() {

  n_bound_ATP_stalled_ = 0;
  for (int i_motor = 0; i_motor < n_active_; i_motor++) {
    Kinesin *motor = active_[i_motor];
    if (motor->IsStalled()) { // and !motor->is_outdated_) {
      if (motor->head_one_.site_ != nullptr and
          motor->head_one_.ligand_ == "ATP") {
        bound_ATP_stalled_[n_bound_ATP_stalled_] = &motor->head_one_;
        n_bound_ATP_stalled_++;
      }
      if (motor->head_two_.site_ != nullptr and
          motor->head_two_.ligand_ == "ATP") {
        bound_ATP_stalled_[n_bound_ATP_stalled_] = &motor->head_two_;
        n_bound_ATP_stalled_++;
      }
    }
  }
}

void KinesinManagement::Update_Bound_ADPP_I() {

  if (verbose_) {
    printf("starting UPDATE ADPP_I\n");
  }
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    n_bound_ADPP_i_[i_aff] = 0;
  }
  for (int i_motor = 0; i_motor < n_active_; i_motor++) {
    Kinesin *motor = active_[i_motor];
    bool eligible{false};
    if (motor->heads_active_ == 1 and !motor->IsStalled()) {
      //    and !motor->is_outdated_) {
      // Only count untethered motors
      if (!motor->tethered_) {
        eligible = true;
      }
      // Motors with satellite xlinks behave as if untethered
      else if (motor->xlink_->heads_active_ == 0) {
        eligible = true;
      }
    }
    if (eligible) {
      if (motor->head_one_.site_ != nullptr and
          motor->head_one_.ligand_ == "ADPP") {
        int aff = motor->head_one_.site_->affinity_;
        bound_ADPP_i_[aff][n_bound_ADPP_i_[aff]] = &motor->head_one_;
        n_bound_ADPP_i_[aff]++;
      }
      if (motor->head_two_.site_ != nullptr and
          motor->head_two_.ligand_ == "ADPP") {
        int aff = motor->head_two_.site_->affinity_;
        bound_ADPP_i_[aff][n_bound_ADPP_i_[aff]] = &motor->head_two_;
        n_bound_ADPP_i_[aff]++;
      }
    }
  }
}

void KinesinManagement::Update_Bound_ADPP_I_Stalled() {

  if (verbose_) {
    printf("starting UPDATE ADPP_I_ST\n");
  }
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    n_bound_ADPP_i_stalled_[i_aff] = 0;
  }
  for (int i_motor = 0; i_motor < n_active_; i_motor++) {
    Kinesin *motor = active_[i_motor];
    bool eligible{false};
    if (motor->heads_active_ == 1 and motor->IsStalled()) {
      //    and !motor->is_outdated_) {
      // Only count untethered motors
      if (!motor->tethered_) {
        eligible = true;
      }
      // Motors with satellite xlinks behave as if untethered
      else if (motor->xlink_->heads_active_ == 0) {
        eligible = true;
      }
    }
    if (eligible) {
      if (motor->head_one_.site_ != nullptr and
          motor->head_one_.ligand_ == "ADPP") {
        int aff = motor->head_one_.site_->affinity_;
        bound_ADPP_i_stalled_[aff][n_bound_ADPP_i_stalled_[aff]] =
            &motor->head_one_;
        n_bound_ADPP_i_stalled_[aff]++;
      }
      if (motor->head_two_.site_ != nullptr and
          motor->head_two_.ligand_ == "ADPP") {
        int aff = motor->head_two_.site_->affinity_;
        bound_ADPP_i_stalled_[aff][n_bound_ADPP_i_stalled_[aff]] =
            &motor->head_two_;
        n_bound_ADPP_i_stalled_[aff]++;
      }
    }
  }
}

void KinesinManagement::Update_Bound_ADPP_I_Tethered() {

  /*
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * dist_cutoff_; x_dub++) {
      n_bound_ADPP_i_tethered_[x_dub] = 0;
    }
    for (int i_motor = 0; i_motor < n_active_; i_motor++) {
      Kinesin *motor = active_[i_motor];
      if (motor->heads_active_ == 1 && motor->tethered_ && !motor->IsStalled())
    {
        // Don't count motors with just satellite xlinks
        if (motor->xlink_->heads_active_ > 0) {
          motor->UpdateExtension();
          // Make sure we didn't force an untether event
          if (motor->tethered_) {
            int x_dub = motor->x_dist_doubled_;
            int index = n_bound_ADPP_i_tethered_[x_dub];
            if (motor->head_one_.site_ != nullptr &&
                motor->head_one_.ligand_ == "ADPP") {
              bound_ADPP_i_tethered_[x_dub][index] = &motor->head_one_;
              n_bound_ADPP_i_tethered_[x_dub]++;
              index++;
            }
            if (motor->head_two_.site_ != nullptr &&
                motor->head_two_.ligand_ == "ADPP") {
              bound_ADPP_i_tethered_[x_dub][index] = &motor->head_two_;
              n_bound_ADPP_i_tethered_[x_dub]++;
              index++;
            }
          }
        }
      }
    }
    */
}

void KinesinManagement::Update_Bound_ADPP_I_Tethered_Stalled() {

  /*
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * dist_cutoff_; x_dub++) {
      n_bound_ADPP_i_teth_st_[x_dub] = 0;
    }
    for (int i_motor = 0; i_motor < n_active_; i_motor++) {
      Kinesin *motor = active_[i_motor];
      if (motor->heads_active_ == 1 && motor->tethered_ && motor->IsStalled()) {
        // Don't count motors with just satellite xlinks
        if (motor->xlink_->heads_active_ > 0) {
          motor->UpdateExtension();
          // Make sure we didn't force an untether event
          if (motor->tethered_) {
            int x_dub = motor->x_dist_doubled_;
            int index = n_bound_ADPP_i_teth_st_[x_dub];
            if (motor->head_one_.site_ != nullptr &&
                motor->head_one_.ligand_ == "ADPP") {
              bound_ADPP_i_teth_st_[x_dub][index] = &motor->head_one_;
              n_bound_ADPP_i_teth_st_[x_dub]++;
              index++;
            }
            if (motor->head_two_.site_ != nullptr &&
                motor->head_two_.ligand_ == "ADPP") {
              bound_ADPP_i_teth_st_[x_dub][index] = &motor->head_two_;
              n_bound_ADPP_i_teth_st_[x_dub]++;
              index++;
            }
          }
        }
      }
    }
    */
}

void KinesinManagement::Update_Bound_ADPP_II() {

  if (verbose_) {
    printf("starting UPDATE ADPP_II\n");
  }
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    n_bound_ADPP_ii_[i_aff] = 0;
  }
  for (int i_motor = 0; i_motor < n_active_; i_motor++) {
    Kinesin *motor = active_[i_motor];
    if (motor->heads_active_ == 2) { // && !motor->is_outdated_) {
      // If both heads are ADPP-bound, pick either front or rear
      // based on weight corresponding to unbinding rates
      // (rear heads k_off_ii; front heads k_off_i)
      if (motor->head_one_.ligand_ == "ADPP" and
          motor->head_two_.ligand_ == "ADPP") {
        double ran = properties_->gsl.GetRanProb();
        double p_tot = p_unbind_ii_[0] + p_unbind_i_[0];
        // Pick rear head if roll in unbind_ii range
        if (ran < p_unbind_ii_[0] / p_tot) {
          if (motor->head_one_.trailing_) {
            int aff = motor->head_one_.site_->affinity_;
            bound_ADPP_ii_[aff][n_bound_ADPP_ii_[aff]] = &motor->head_one_;
            n_bound_ADPP_ii_[aff]++;
          } else {
            int aff = motor->head_two_.site_->affinity_;
            bound_ADPP_ii_[aff][n_bound_ADPP_ii_[aff]] = &motor->head_two_;
            n_bound_ADPP_ii_[aff]++;
          }
        }
        // Otherwise, pick front head
        else {
          if (!motor->head_one_.trailing_) {
            int aff = motor->head_one_.site_->affinity_;
            bound_ADPP_ii_[aff][n_bound_ADPP_ii_[aff]] = &motor->head_one_;
            n_bound_ADPP_ii_[aff]++;
          } else {
            int aff = motor->head_two_.site_->affinity_;
            bound_ADPP_ii_[aff][n_bound_ADPP_ii_[aff]] = &motor->head_two_;
            n_bound_ADPP_ii_[aff]++;
          }
        }
      } else {
        if (motor->head_one_.site_ != nullptr and
            motor->head_one_.ligand_ == "ADPP") {
          int aff = motor->head_one_.site_->affinity_;
          bound_ADPP_ii_[aff][n_bound_ADPP_ii_[aff]] = &motor->head_one_;
          n_bound_ADPP_ii_[aff]++;
        }
        if (motor->head_two_.site_ != nullptr and
            motor->head_two_.ligand_ == "ADPP") {
          int aff = motor->head_two_.site_->affinity_;
          bound_ADPP_ii_[aff][n_bound_ADPP_ii_[aff]] = &motor->head_two_;
          n_bound_ADPP_ii_[aff]++;
        }
      }
    }
  }
}

void KinesinManagement::Update_Bound_Teth() {

  /*
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * dist_cutoff_; x_dub++) {
      n_bound_tethered_[x_dub] = 0;
    }
    for (int i_motor = 0; i_motor < n_active_; i_motor++) {
      Kinesin *motor = active_[i_motor];
      if (motor->heads_active_ > 0 && motor->tethered_ == true) {
        if (motor->xlink_->heads_active_ > 0) {
          motor->UpdateExtension();
          if (motor->tethered_ == true) {
            int x_dub = motor->x_dist_doubled_;
            int index = n_bound_tethered_[x_dub];
            if (motor->heads_active_ == 2) {
              bound_tethered_[x_dub][index] = &motor->head_one_;
            } else {
              bound_tethered_[x_dub][index] = motor->GetActiveHead();
            }
            n_bound_tethered_[x_dub]++;
          }
        }
      }
    }
    */
}

void KinesinManagement::Update_Bound_Unteth() {

  /*
    n_bound_untethered_ = 0;
    for (int i_motor = 0; i_motor < n_active_; i_motor++) {
      Kinesin *motor = active_[i_motor];
      if (!motor->tethered_) {
        if (motor->heads_active_ == 1) {
          bound_untethered_[n_bound_untethered_] = motor->GetActiveHead();
          n_bound_untethered_++;
        } else if (motor->heads_active_ == 2) {
          bound_untethered_[n_bound_untethered_] = &motor->head_one_;
          n_bound_untethered_++;
        }
      }
    }
    */
}

void KinesinManagement::Run_KMC() {

  if (parameters_->motors.c_bulk == 0) {
    return;
  }
  sys_time start = sys_clock::now();
  Update_All_Lists();
  sys_time finish_list = sys_clock::now();
  properties_->wallace.t_motors_[1] += (finish_list - start).count();
  Refresh_Populations();
  sys_time finish_pops = sys_clock::now();
  properties_->wallace.t_motors_[2] += (finish_pops - finish_list).count();
  Generate_Execution_Sequence();
  sys_time finish_seq = sys_clock::now();
  properties_->wallace.t_motors_[3] += (finish_seq - finish_pops).count();
  if (verbose_) {
    printf("Start of motor KMC cycle: %lu events to execute\n",
           events_to_exe_.size());
  }
  for (int i_event = 0; i_event < events_to_exe_.size(); i_event++) {
    events_to_exe_[i_event]->Execute();
  }
  sys_time finish_all = sys_clock::now();
  properties_->wallace.t_motors_[4] += (finish_all - finish_seq).count();
  properties_->wallace.t_motors_[0] += (finish_all - start).count();
}

void KinesinManagement::Refresh_Populations() {

  // If tethers are active, update all tether extensions
  if (parameters_->motors.tethers_active) {
    for (int i_entry{0}; i_entry < n_active_; i_entry++) {
      active_[i_entry]->UpdateExtension();
    }
  }
  // Define our 'refresh' function to be called on any outdated motors
  auto refresh = [&](POP_T *head) {
    Kinesin *entry = head->motor_;
    entry->head_one_.in_scratch_ = false;
    entry->head_two_.in_scratch_ = false;
    std::string state;
    std::string root, ligand;
    if (entry->heads_active_ == 0) {
      // If tethered, automatically classified as "free_teth"
      if (entry->tethered_)
        state = std::string("free_teth");
      // Otherwise, homeboy was iced and needs to be removed
      else
        state = std::string("unbound");
      entry->head_one_.state_ = state;
      entry->head_two_.state_ = state;
    } else if (entry->heads_active_ == 1) {
      std::string inactive_state;
      POP_T *active_head = entry->GetActiveHead();
      if (active_head->ligand_ != std::string{"ADPP"}) {
        state = std::string{"bound_"} + active_head->ligand_;
        inactive_state = std::string{"unbound"};
      } else {
        state = std::string{"bound_ADPP_i"};
        int i_dock = entry->GetDockedCoordinate() - entry->mt_->coord_;
        if (!entry->mt_->lattice_[i_dock].occupied_) {
          inactive_state = std::string{"docked"};
        } else {
          inactive_state = std::string{"unbound"};
        }
      }
      if (entry->IsStalled()) {
        state += std::string{"_st"};
      }
      if (active_head == &entry->head_one_) {
        entry->head_one_.state_ = state;
        entry->head_two_.state_ = inactive_state;
      } else {
        entry->head_one_.state_ = inactive_state;
        entry->head_two_.state_ = state;
      }
    }
    // Doubly-bound motors
    else {
      // Do head one first
      if (entry->head_one_.ligand_ != std::string{"ADPP"}) {
        state = std::string{"bound_"} + entry->head_one_.ligand_;
        if (entry->IsStalled()) {
          state += std::string{"_st"};
        }
      } else {
        state = std::string{"bound_ADPP_ii"};
      }
      entry->head_one_.state_ = state;
      // Head two next
      if (entry->head_two_.ligand_ != std::string{"ADPP"}) {
        state = std::string{"bound_"} + entry->head_two_.ligand_;
        if (entry->IsStalled()) {
          state += std::string{"_st"};
        }
      } else {
        state = std::string{"bound_ADPP_ii"};
      }
      entry->head_two_.state_ = state;
    }
  };
  // Update all current entries in scratch_
  for (int i_entry{0}; i_entry < n_scratched_; i_entry++) {
    if (scratch_[i_entry]->in_scratch_) {
      refresh(scratch_[i_entry]);
    }
  }
  // Clear scratch_ list
  n_scratched_ = 0;
}

void KinesinManagement::Generate_Execution_Sequence() {

  int n_events = Sample_Event_Statistics();
  if (n_events > 0) {
    EVENT_T *pre_array[n_events];
    int i_array{0};
    for (int i_event{0}; i_event < events_.size(); i_event++) {
      int n_expected = events_[i_event].n_expected_;
      for (int i_entry(0); i_entry < n_expected; i_entry++) {
        // Make sure we don't bamboozle ourselves
        if (i_array >= n_events) {
          printf("WOAH BUDDY. check MOT genKMC\n");
          exit(1);
        }
        // Store & pass ptr to this event
        pre_array[i_array] = &events_[i_event];
        i_array++;
      }
    }
    if (i_array != n_events) {
      printf("NOT SURE in MOTUR GEN KMC \n");
      exit(1);
    }
    // Shuffle for some guuuud random exe order
    if (n_events > 1) {
      gsl_ran_shuffle(properties_->gsl.rng_, pre_array, n_events,
                      sizeof(EVENT_T *));
    }
    // Transfer info to permanent vector owned by MGMT
    events_to_exe_.resize(n_events);
    for (int i_entry(0); i_entry < n_events; i_entry++) {
      events_to_exe_[i_entry] = pre_array[i_entry];
    }
  } else {
    events_to_exe_.clear();
  }
  /*
  sys_time start = sys_clock::now();
  // Track the time it takes to update lists and sample statistics
  sys_time finish = sys_clock::now();
  auto elapsed = std::chrono::duration_cast<t_unit>(finish - start);
  properties_->wallace.t_motors_[1] += elapsed.count();
  start = sys_clock::now();
  // Stat correcton; ensure there aren't more events than population size
  // Bind_ii and unbind_i compete for n_bound_ADPP_i_
  while (events_[4].n_expected_ + events_[6].n_expected_ > n_docked_ &&
         n_docked_ > 0) {
    //	> n_bound_ADPP_i_){
    double ran = properties_->gsl.GetRanProb();
    double p_tot = events_[4].p_occur_ + events_[6].p_occur_;
    if (ran < events_[4].p_occur_ / p_tot && events_[4].n_expected_ > 0)
      events_[4].n_expected_--;
    else if (events_[6].n_expected_ > 0)
      events_[6].n_expected_--;
    else
      events_[4].n_expected_--;
  }
  // Stat correction for tether-based events
  if (parameters_->motors.tethers_active) {
    int list_size = 2 * (dist_cutoff_ - comp_cutoff_) + 1;
    int offset = 2 * comp_cutoff_;
    // Scan over different possible extensions first
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * dist_cutoff_; x_dub++) {
      // Bind_ATP_teth and untether_bound compete
      int i_ATP_teth = 9 + x_dub - offset;
      int i_unteth = 12 + 4 * list_size + x_dub - offset;
      while (events_[i_ATP_teth].n_expected_ + events_[i_unteth].n_expected_ >
                 n_bound_NULL_tethered_[x_dub] &&
             n_bound_NULL_tethered_[x_dub] > 0) {
        //	> n_bound_tethered_[x_dub]){
        double ran = properties_->gsl.GetRanProb();
        double p_tot =
            events_[i_ATP_teth].p_occur_ + events_[i_unteth].p_occur_;
        if (ran < events_[i_ATP_teth].p_occur_ / p_tot &&
            events_[i_ATP_teth].n_expected_ > 0)
          events_[i_ATP_teth].n_expected_--;
        else if (events_[i_unteth].n_expected_ > 0)
          events_[i_unteth].n_expected_--;
      }
      // Bind_ii_teth, unbind_i_teth(_st), & unteth_bound compete
      int i_bind_ii = 9 + list_size + x_dub - offset;
      int i_unbind_i_teth = 9 + 2 * list_size + x_dub - offset;
      int i_unbind_i_teth_st = 9 + 3 * list_size + x_dub - offset;
      // Randomize order that stall vs mobile stats are corrected
      int i_0 = properties_->gsl.GetRanInt(2);
      for (int i_step = 0; i_step < 2; i_step++) {
        int index = (i_0 + i_step) % 2;
        if (index == 0) {
          // Non-stalled stats
          while (events_[i_bind_ii].n_expected_ +
                         events_[i_unbind_i_teth].n_expected_ +
                         events_[i_unteth].n_expected_ >
                     n_bound_ADPP_i_tethered_[x_dub] &&
                 n_bound_ADPP_i_tethered_[x_dub] > 0) {
            double ran = properties_->gsl.GetRanProb();
            double p_tot = events_[i_bind_ii].p_occur_ +
                           events_[i_unbind_i_teth].p_occur_ +
                           events_[i_unteth].p_occur_;
            double p_mid =
                events_[i_bind_ii].p_occur_ +
  events_[i_unbind_i_teth].p_occur_; if (ran < events_[i_bind_ii].p_occur_ /
  p_tot && events_[i_bind_ii].n_expected_ > 0)
              events_[i_bind_ii].n_expected_--;
            else if (ran < p_mid / p_tot &&
                     events_[i_unbind_i_teth].n_expected_ > 0)
              events_[i_unbind_i_teth].n_expected_--;
            else if (events_[i_unteth].n_expected_ > 0)
              events_[i_unteth].n_expected_--;
          }
        } else {
          // Stalled stats
          while (events_[i_bind_ii].n_expected_ +
                         events_[i_unbind_i_teth_st].n_expected_ +
                         events_[i_unteth].n_expected_ >
                     n_bound_ADPP_i_teth_st_[x_dub] &&
                 n_bound_ADPP_i_teth_st_[x_dub] > 0) {
            double ran = properties_->gsl.GetRanProb();
            double p_tot = events_[i_bind_ii].p_occur_ +
                           events_[i_unbind_i_teth_st].p_occur_ +
                           events_[i_unteth].p_occur_;
            double p_mid = events_[i_bind_ii].p_occur_ +
                           events_[i_unbind_i_teth_st].p_occur_;
            if (ran < events_[i_bind_ii].p_occur_ / p_tot &&
                events_[i_bind_ii].n_expected_ > 0)
              events_[i_bind_ii].n_expected_--;
            else if (ran < p_mid / p_tot &&
                     events_[i_unbind_i_teth_st].n_expected_ > 0)
              events_[i_unbind_i_teth_st].n_expected_--;
            else if (events_[i_unteth].n_expected_ > 0)
              events_[i_unteth].n_expected_--;
          }
        }
      }
    }
    // Tether_free & tether_bound compete for untethered xlinks
    int i_teth_free = 9 + 4 * list_size;
    int i_teth_bound = 10 + 4 * list_size;
    while (events_[i_teth_free].n_expected_ +
               events_[i_teth_bound].n_expected_ >
           properties_->prc1.n_bound_unteth_) {
      double ran = properties_->gsl.GetRanProb();
      double p_tot =
          events_[i_teth_free].p_occur_ + events_[i_teth_bound].p_occur_;
      if (ran < events_[i_teth_free].p_occur_ / p_tot &&
          events_[i_teth_free].n_expected_ > 0)
        events_[i_teth_free].n_expected_--;
      else if (events_[i_teth_bound].n_expected_ > 0)
        events_[i_teth_bound].n_expected_--;
    }
    // Bind_i_teth & untether_free compete for free_tethered motors
    int i_bind_i_teth = 8;
    int i_untether = 11 + 4 * list_size;
    while (events_[i_bind_i_teth].n_expected_ +
               events_[i_untether].n_expected_ >
           n_free_tethered_) {
      double ran = properties_->gsl.GetRanProb();
      double p_tot =
          events_[i_bind_i_teth].p_occur_ + events_[i_untether].p_occur_;
      if (ran < events_[i_bind_i_teth].p_occur_ / p_tot &&
          events_[i_bind_i_teth].n_expected_ > 0)
        events_[i_bind_i_teth].n_expected_--;
      else if (events_[i_untether].n_expected_ > 0)
        events_[i_untether].n_expected_--;
    }
  }
  int n_events = 0;
  int pre_array[100 * events_.size()];
  // Scan over all events; record those with >0 expected in this timestep
  for (int i_entry = 0; i_entry < events_.size(); i_entry++) {
    for (int i = 0; i < events_[i_entry].n_expected_; i++) {
      // Make sure we don't bamboozle ourselves here
      if (n_events > 100 * events_.size()) {
        printf("Error in GenerateKMCList for motors!!\n");
        exit(1);
      }
      pre_array[n_events] = events_[i_entry].kmc_code_;
      n_events++;
    }
  }
  // If total expected events is greater than 0, construct IDs_to_exe_
  if (n_events > 0) {
    // Trim array to appropriate size (only non-zero data)
    int reduced_array[n_events];
    for (int i_entry = 0; i_entry < n_events; i_entry++) {
      reduced_array[i_entry] = pre_array[i_entry];
    }
    // If there are more than 1 expected events, shuffle their order
    if (n_events > 1)
      gsl_ran_shuffle(properties_->gsl.rng_, reduced_array, n_events,
                      sizeof(int));
    // Transfer shuffled array into IDs_to_exe_ structure
    IDs_to_exe_.resize(n_events);
    for (int i_entry = 0; i_entry < n_events; i_entry++) {
      IDs_to_exe_[i_entry] = reduced_array[i_entry];
    }
  }
  // Otherwise, simply clear IDs_to_exe_
  else {
    IDs_to_exe_.clear();
  }
  // Track the time it takes to construct KMC list
  finish = sys_clock::now();
  elapsed = std::chrono::duration_cast<t_unit>(finish - start);
  properties_->wallace.t_motors_[2] += elapsed.count();
  */
}
int KinesinManagement::Sample_Event_Statistics() {

  // Scan through all events & get expected number of occurrences
  int n_events_tot{0};
  for (int i_event{0}; i_event < events_.size(); i_event++) {
    events_[i_event].SampleStatistics();
    n_events_tot += events_[i_event].n_expected_;
    if (verbose_) {
      printf("expect %i for event ", events_[i_event].n_expected_);
      std::cout << events_[i_event].name_;
      printf(" (%i available)\n", *events_[i_event].n_avail_);
    }
  }
  /* Scan through all target pops. & ensure none will become negative */
  for (int i_pop{0}; i_pop < events_by_pop_.size(); i_pop++) {
    int n_competitors = events_by_pop_[i_pop].size();
    if (verbose_) {
      printf("INSPECTIN' POP #%i\n", i_pop);
    }
    if (n_competitors > 1) {
      int n_active_competitors{0};
      int n_events_loc{0};
      int n_avail_loc{0};
      double p_tot{0.0};
      if (verbose_) {
        printf("checking that ");
        std::cout << events_by_pop_[i_pop][1]->targets_[0];
        printf(" doesn't go negative\n");
      }
      for (int i_entry{0}; i_entry < n_competitors; i_entry++) {
        if (verbose_) {
          std::cout << events_by_pop_[i_pop][i_entry]->name_;
          printf(" is a competitor\n");
        }
        n_events_loc += events_by_pop_[i_pop][i_entry]->n_expected_;
        p_tot += (events_by_pop_[i_pop][i_entry]->n_expected_ *
                  events_by_pop_[i_pop][i_entry]->p_occur_);
        if (events_by_pop_[i_pop][i_entry]->n_expected_ > 0 and
            events_by_pop_[i_pop][i_entry]->is_not_redundant_) {
          n_avail_loc += *events_by_pop_[i_pop][i_entry]->n_avail_;
          n_active_competitors++;
        }
      }
      if (verbose_) {
        printf("%i events loc\n", n_events_loc);
        printf("%i tot available\n", n_avail_loc);
      }
      if (n_active_competitors > 0) {
        if (n_avail_loc == 0 and n_events_loc > 0) {
          printf("uhhhh ?? - %i\n", n_events_loc);
          std::cout << events_by_pop_[i_pop][0]->name_ << std::endl;
          std::cout << events_by_pop_[i_pop][0]->targets_[0] << std::endl;
          exit(1);
        }
        while (n_events_loc > n_avail_loc) {
          double p_cum{0};
          double ran{properties_->gsl.GetRanProb()};
          for (int i_entry(0); i_entry < n_competitors; i_entry++) {
            p_cum += (events_by_pop_[i_pop][i_entry]->n_expected_ *
                      events_by_pop_[i_pop][i_entry]->p_occur_ / p_tot);
            if (p_cum >= ran and
                events_by_pop_[i_pop][i_entry]->n_expected_ > 0) {
              events_by_pop_[i_pop][i_entry]->n_expected_--;
              n_events_tot--;
              n_events_loc--;
              p_tot -= events_by_pop_[i_pop][i_entry]->p_occur_;
              break;
            }
          }
        }
      }
    }
  }
  return n_events_tot;
}

double KinesinManagement::GetWeight_BindTethered() {

  double weight = 0;
  for (int i_entry = 0; i_entry < n_free_tethered_; i_entry++) {
    POP_T *head = std::get<POP_T *>(free_tethered_[i_entry]);
    weight += head->motor_->GetTotalBindingWeight();
  }
  return weight;
}

double KinesinManagement::GetWeight_TetherBound() {

  double weight = 0;
  for (int i_entry = 0; i_entry < n_bound_untethered_; i_entry++) {
    POP_T *head = std::get<POP_T *>(bound_untethered_[i_entry]);
    weight += head->motor_->GetTotalTetheringWeight();
  }
  return weight;
}

Kinesin::head *KinesinManagement::CheckScratchFor(std::string pop) {

  if (verbose_) {
    printf("Scratch was checked for ");
    std::cout << pop << std::endl;
  }
  POP_T *head(nullptr);
  for (int i_entry(0); i_entry < n_scratched_; i_entry++) {
    if (verbose_) {
      printf("   entry #%i: ", i_entry);
      std::cout << scratch_[i_entry]->state_ << std::endl;
    }
    if (scratch_[i_entry]->state_ == pop) {
      head = scratch_[i_entry];
      if (verbose_) {
        printf("successfully found!\n");
      }
      break;
    }
  }
  return head;
}

void KinesinManagement::SaveToScratch(POP_T *head) {

  if (head->in_scratch_) {
    if (verbose_) {
      printf("no sir for ");
      std::cout << head->state_ << std::endl;
    }
    return;
  }
  head->in_scratch_ = true;
  scratch_[n_scratched_] = head;
  n_scratched_++;
  if (verbose_) {
    std::cout << head->state_;
    printf(" was saved to scratch.\n");
  }
}

void KinesinManagement::ReportExecutionOf(std::string event_name) {

  if (not verbose_) {
    return;
  }
  printf("Executing ");
  std::cout << event_name << std::endl;
}

void KinesinManagement::ReportFailureOf(std::string event_name) {

  printf("yo we failed to execute ");
  std::cout << event_name << std::endl;
  //  properties_->wallace.PauseSim(5);
}

/*
void KinesinManagement::Execution_Relay(ENTRY_T target, int code) {

  bool verbose{false};
  POP_T *head(nullptr);
  SITE_T *site(nullptr);
  ALT_T *xlink_head(nullptr);
  try {
    head = std::get<POP_T *>(target);
  } catch (...) {
    try {
      site = std::get<SITE_T *>(target);
    } catch (...) {
      xlink_head = std::get<ALT_T *>(target);
    }
  }
  if (head == nullptr && site == nullptr && xlink_head == nullptr) {
    printf("What the FUCK!! RUN TO THE CENTER!!! (... motor kmc_relay)\n");
    printf("code: %i\n", code);
    return;
    // exit(1);
  }
  switch (code) {
  case 10:
    if (verbose) {
      printf("bind_i on MT %i, site %i\n", site->index_, site->mt_->index_);
    }
    KMC_Bind_I(site);
    break;
  case 11:
      if (verbose) {
        printf("bind_i_teth on ");
        std::cout << head->state_ << std::endl;
      }
      Bind_I_Teth(head);
    break;
  case 20:
    if (verbose) {
      printf("bind_ATP on ");
      std::cout << head->state_ << std::endl;
    }
    KMC_Bind_ATP(head);
    break;
  case 21:
      if (verbose) {
        printf("bind_ii_teth on ");
        std::cout << head->state_ << std::endl;
      }
      Bind_II(head);
    break;
  case 30:
    if (verbose) {
      printf("hydrolysis on ");
      std::cout << head->state_ << std::endl;
    }
    KMC_Hydrolyze(head);
    break;
  case 31:
    if (verbose) {
      printf("hydrolysis_st on ");
      std::cout << head->state_ << std::endl;
    }
    KMC_Hydrolyze(head);
    break;
  case 40:
    if (verbose) {
      printf("bind_ii on ");
      std::cout << head->state_ << std::endl;
    }
    KMC_Bind_II(head);
    break;
  case 41:
      if (verbose) {
        printf("unbind_ii_to_teth on ");
        std::cout << head->state_ << std::endl;
      }
      Unbind_II(head);
    break;
  case 42:
      if (verbose) {
        printf("unbind_ii_fr_teth on ");
        std::cout << head->state_ << std::endl;
      }
      Unbind_II(head);
    break;
  case 50:
    if (verbose) {
      printf("unbind_ii on ");
      std::cout << head->state_ << std::endl;
    }
    KMC_Unbind_II(head);
    break;
  case 60:
    if (verbose) {
      printf("unbind_i on ");
      std::cout << head->state_ << std::endl;
    }
    KMC_Unbind_I(head);
    break;
  case 61:
    if (verbose) {
      printf("unbind_i_st on ");
      std::cout << head->state_ << std::endl;
    }
    KMC_Unbind_I(head);
  }
}
*/

void KinesinManagement::KMC_Bind_I(SITE_T *site) {

  // Get random free motor
  Kinesin *motor = GetFreeMotor();
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
  active_[n_active_] = motor;
  motor->active_index_ = n_active_;
  n_active_++;
  // For binding (and only binding), save to scratch AFTER event
  SaveToScratch(&motor->head_one_);
}
/*
void KinesinManagement::KMC_Bind_I_Tethered() {

  Update_Free_Teth();
  properties_->microtubules.UpdateUnoccupied();
  if (n_free_tethered_ > 0 && properties_->microtubules.n_unoccupied_ > 0) {
    // Pick a random tethered_free motor to bind
    int i_motor = properties_->gsl.GetRanInt(n_free_tethered_);
    Kinesin *motor = free_tethered_[i_motor];
    // Make sure we get a motor that has valid neighbors
    motor->UpdateNeighborSites();
    int attempts(0);
    bool neighbors_exist(true);
    while (motor->n_neighbor_sites_ == 0) {
      if (attempts > 10 * n_free_tethered_) {
        neighbors_exist = false;
        break;
      }
      i_motor = properties_->gsl.GetRanInt(n_free_tethered_);
      motor = free_tethered_[i_motor];
      motor->UpdateNeighborSites();
      attempts++;
    }
    if (neighbors_exist) {
      Tubulin *site = motor->GetWeightedNeighborSite();
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
      motor->UpdateExtension();
    } else {
      printf("Failed to bind_i_teth\n");
      //			exit(1);
    }
  } else if (properties_->microtubules.n_unoccupied_ == 0) {
    printf("Error in Bind_I_Tethered: no unoccupied sites\n");
    //      exit(1);
  } else {
    printf("Error in Bind_I_Tethered: no tethered free motors\n");
    //		exit(1);
  }
}
*/
void KinesinManagement::KMC_Bind_ATP(POP_T *head) {

  // Save this entry to scratch before modifying it
  SaveToScratch(head);
  // Update motor head
  head->ligand_ = "ATP";
  // If head is leading, change conformation
  if (!head->trailing_) {
    int i_plus_end = head->site_->mt_->plus_end_;
    int i_site = head->site_->index_;
    // If endpausing is active, don't step off MT boundary sites
    if (parameters_->motors.endpausing_active && i_site != i_plus_end) {
      head->motor_->ChangeConformation();
    } else {
      head->motor_->ChangeConformation();
    }
  }
}

/* FIXME delete this jawn once u verify it isn't needed
void KinesinManagement::KMC_Bind_ATP_Tethered(int x_dub) {

  Update_Bound_NULL_Teth();
  //	printf("%i bound_NULL\n", n_bound_NULL_);
  int n_bound = n_bound_NULL_tethered_[x_dub];
  if (n_bound > 0) {
    // Get a random bound_NULL motor
    int i_entry = properties_->gsl.GetRanInt(n_bound);
    Kinesin::head *head = bound_NULL_tethered_[x_dub][i_entry];
    // Verify correctness
    if (head->ligand_ != "NULL" || head->site_ == nullptr) {
      printf("Error in KMC_Bind_ATP_Teth()\n");
      exit(1);
    }
    // Update motor head
    head->ligand_ = "ATP";
    // If head is leading, change conformation
    if (!head->trailing_) {
      // If endpausing is active, don't step off MT boundary sites
      if (parameters_->motors.endpausing_active) {
        Microtubule *mt = head->site_->mt_;
        int i_site = head->site_->index_;
        int dx = mt->delta_x_;
        int mt_length = mt->n_sites_ - 1;
        if (!(i_site == 0 && dx == -1) && !(i_site == mt_length && dx == 1)) {
          head->motor_->ChangeConformation();
        }
      } else
        head->motor_->ChangeConformation();
    }
    head->motor_->UpdateExtension();
  } else {
    printf("Failed to Bind_ATP_Teth: no bound_NULL motors.\n");
    //		exit(1);
  }
}
*/

void KinesinManagement::KMC_Hydrolyze(POP_T *head) {

  // Save this entry to scratch before modifying it
  SaveToScratch(head);
  head->ligand_ = "ADPP";
}

void KinesinManagement::KMC_Bind_II(POP_T *head) {

  if (head->motor_->heads_active_ == 1) {
    // Save active head to scratch before binding
    SaveToScratch(head->motor_->GetActiveHead());
    // Verify that proposed site is unoccupied
    int i_dock =
        head->motor_->GetDockedCoordinate() - head->motor_->mt_->coord_;
    Tubulin *dock_site = &head->motor_->mt_->lattice_[i_dock];
    if (dock_site->occupied_) {
      printf("Error in KMC_Bind_II(): Dock is occupied!!\n");
      exit(1);
    }
    // Update site
    dock_site->motor_head_ = head;
    dock_site->occupied_ = true;
    // Update motor
    head->site_ = dock_site;
    head->state_ = "old_docked_boi";
    head->ligand_ = "NULL";
    head->motor_->heads_active_++;
    // Save newly-bound head to scratch
    SaveToScratch(head);
  } else if (head->stored_dock_site_ != nullptr) {
    printf("noooo\n");
    Tubulin *dock_site = head->stored_dock_site_;
    if (dock_site->occupied_) {
      printf("Error in Bind_II: option 2\n");
      exit(1);
    }
    // Update site
    dock_site->motor_head_ = head;
    dock_site->occupied_ = true;
    // Update motor
    head->site_ = dock_site;
    head->ligand_ = "NULL";
    head->motor_->heads_active_++;
    // Save newly-bound head to scratch
    SaveToScratch(head);
    head->stored_dock_site_ = nullptr;
  } else {
    printf("failed to Bind_II\n");
    exit(1);
  }
}
/*
void KinesinManagement::KMC_Bind_II_Tethered(int x_dub) {

  Update_Docked_Teth();
  if (n_docked_tethered_[x_dub] > 0) {
    // Get a random docked motor
    int i_entry = properties_->gsl.GetRanInt(n_docked_tethered_[x_dub]);
    Kinesin::head *docked_head = docked_tethered_[x_dub][i_entry];
    Kinesin *motor = docked_head->motor_;
    // Verify correctness
    if (docked_head->ligand_ != "ADP" ||
        motor->GetActiveHead()->ligand_ != "ADPP" || motor->frustrated_) {
      printf("Error in KMC_Bind_II()\n");
      exit(1);
    }
    // Verify that proposed site is unoccupied
    int i_dock = motor->GetDockedCoordinate() - motor->mt_->coord_;
    Tubulin *dock_site = &motor->mt_->lattice_[i_dock];
    if (dock_site->occupied_) {
      printf("Error in KMC_Bind_II(): Dock is occupied!!\n");
      exit(1);
    }
    // Update site
    dock_site->motor_head_ = docked_head;
    dock_site->occupied_ = true;
    // Update motor
    docked_head->site_ = dock_site;
    docked_head->ligand_ = "NULL";
    motor->heads_active_++;
  } else {
    printf("Failed to Bind_II_Tethered(%i): no docked motors.\n", x_dub);
    //		exit(1);
  }
}
*/
void KinesinManagement::KMC_Unbind_II(POP_T *head) {

  if (head->motor_->heads_active_ == 2) {
    // Save entry to scratch before we modify it
    SaveToScratch(head);
    SaveToScratch(head->GetOtherHead());
    // Update site
    head->site_->occupied_ = false;
    head->site_->motor_head_ = nullptr;
    // Update motor
    head->site_ = nullptr;
    head->ligand_ = "ADP";
    head->motor_->heads_active_--;
    if (head->motor_->frustrated_) {
      head->motor_->frustrated_ = false;
      // Only step if rear head was just unbound
      if (head->trailing_) {
        // If endpausing is active, don't step off boundary sites
        if (parameters_->motors.endpausing_active) {
          Tubulin *site = head->motor_->GetActiveHead()->site_;
          int i_site = site->index_;
          int i_plus = site->mt_->plus_end_;
          if (i_site != i_plus) {
            head->motor_->ChangeConformation();
          } else {
            head->motor_->frustrated_ = false;
          }
          // If no endpausing, step regardless
        } else {
          head->motor_->ChangeConformation();
        }
        // If front head was unbound, no stepping & hydrolysis was futile
      } else {
        head->motor_->frustrated_ = false;
        // printf("FUTILE HYDROLYSIS!!!\n");
      }
    }
    if (head->motor_->frustrated_) {
      printf("Error THREE in motor KMC_Unbind(). - THIS\n");
      printf("heads active: %i\n", head->motor_->heads_active_);
      exit(1);
    }
  } else {
    printf("WACK in unbind_ii\n");
    exit(1);
  }
}

void KinesinManagement::KMC_Unbind_I(POP_T *head) {

  if (head->motor_->heads_active_ == 1) {
    SaveToScratch(head);
    POP_T *docked_head = head->motor_->StoreDockSite();
    if (docked_head != nullptr) {
      SaveToScratch(docked_head);
    }
    // Update site
    head->site_->occupied_ = false;
    head->site_->motor_head_ = nullptr;
    // Update motor
    head->site_ = nullptr;
    head->ligand_ = "ADP";
    head->motor_->heads_active_--;
    head->motor_->mt_ = nullptr;
    // If this motor has a satellite xlink, untether it
    if (head->motor_->xlink_ != nullptr) {
      if (head->motor_->xlink_->heads_active_ == 0) {
        head->motor_->UntetherSatellite();
      } else {
        printf("Error in KMC_Unbind_I: tethed to bound xlink?\n");
        exit(1);
      }
    }
    if (head->motor_->heads_active_ == 0) {
      // Remove this motor from active_, replace with last entry
      int this_index = head->motor_->active_index_;
      if (this_index != n_active_ - 1) {
        Kinesin *last_entry = active_[n_active_ - 1];
        active_[this_index] = last_entry;
        last_entry->active_index_ = this_index;
      }
      n_active_--;
    }
  }
  /*
  else {
    printf("WHUT\n");
    exit(1);
  }
  */
}

/*
void KinesinManagement::KMC_Unbind_I_Stalled() {

  Update_Bound_ADPP_I_Stalled();
  if (n_bound_ADPP_i_stalled_ > 0) {
    // Get a random bound_ADPP motor head
    int i_entry = properties_->gsl.GetRanInt(n_bound_ADPP_i_stalled_);
    Kinesin::head *head = bound_ADPP_i_stalled_[i_entry];
    if (head->site_ == nullptr || head->ligand_ != "ADPP") {
      printf("Error in KMC_Unbind_I_Stalled (motors): ");
      std::cout << head->ligand_;
      printf(" bound to head\n");
      exit(1);
    }
    if (head->motor_->heads_active_ == 2) {
      printf("Error TWO in KMC_Unbind_I_Stalled (motors).\n");
      exit(1);
    }
    // Update site
    head->site_->occupied_ = false;
    head->site_->motor_head_ = nullptr;
    // Update motor
    head->site_ = nullptr;
    head->ligand_ = "ADP";
    head->motor_->heads_active_--;
    head->motor_->mt_ = nullptr;
    // If this motor has a satellite xlink, untether it
    if (head->motor_->tethered_) {
      if (head->motor_->xlink_->heads_active_ == 0) {
        head->motor_->UntetherSatellite();
      }
    }
    if (!head->motor_->tethered_) {
      // Remove this motor from active_, replace with last entry
      int this_index = head->motor_->active_index_;
      if (this_index != n_active_ - 1) {
        Kinesin *last_entry = active_[n_active_ - 1];
        active_[this_index] = last_entry;
        last_entry->active_index_ = this_index;
      }
      n_active_--;
    }
  } else {
    printf("Failed to KMC_Unbind_I: no bound_ADPP motors.\n");
    //		exit(1);
  }
}
void KinesinManagement::KMC_Unbind_I_Tethered(int x_dub) {

  Update_Bound_ADPP_I_Tethered();
  int n_bound = n_bound_ADPP_i_tethered_[x_dub];
  if (n_bound > 0) {
    // Randomly pick a motor
    int i_entry = properties_->gsl.GetRanInt(n_bound);
    Kinesin::head *head = bound_ADPP_i_tethered_[x_dub][i_entry];
    if (head->motor_->x_dist_doubled_ != x_dub || head->motor_->IsStalled()) {
      printf("\"not a fan\" -- unbind_i_teth (motors)\n");
      exit(1);
    }
    // Update site
    head->site_->motor_head_ = nullptr;
    head->site_->occupied_ = false;
    // Update motor details
    head->site_ = nullptr;
    head->ligand_ = "ADP";
    head->motor_->heads_active_--;
    head->motor_->mt_ = nullptr;
  } else {
    printf("Error in Unbind_I_Tethered: no pseudo bound motors!\n");
    //		exit(1);
  }
}

void KinesinManagement::KMC_Unbind_I_Tethered_Stalled(int x_dub) {

  Update_Bound_ADPP_I_Tethered_Stalled();
  int n_bound = n_bound_ADPP_i_teth_st_[x_dub];
  if (n_bound > 0) {
    // Randomly pick a motor
    int i_entry = properties_->gsl.GetRanInt(n_bound);
    Kinesin::head *head = bound_ADPP_i_teth_st_[x_dub][i_entry];
    if (head->motor_->x_dist_doubled_ != x_dub || !head->motor_->IsStalled())
{ printf("\"not a fan\" -- unbind_i_teth_STALL (motors)\n"); exit(1);
    }
    // Update site
    head->site_->motor_head_ = nullptr;
    head->site_->occupied_ = false;
    // Update motor details
    head->site_ = nullptr;
    head->ligand_ = "ADP";
    head->motor_->heads_active_--;
    head->motor_->mt_ = nullptr;
  } else {
    printf("Error in Unbind_I_Teth_ST: no stage-1 STALLED mots\n");
  }
}

void KinesinManagement::KMC_Tether_Free() {

  properties_->prc1.Update_Bound_Unteth();
  if (properties_->prc1.n_bound_unteth_ > 0) {
    Kinesin *motor = GetFreeMotor();
    // Randomly pick an xlink
    AssociatedProtein *xlink = properties_->prc1.GetBoundUntetheredXlink();
    // Update motor and xlink details
    motor->xlink_ = xlink;
    motor->tethered_ = true;
    xlink->tethered_ = true;
    xlink->motor_ = motor;
    // Update active_ list
    active_[n_active_] = motor;
    motor->active_index_ = n_active_;
    n_active_++;
  } else {
    printf("Error in Tether_Free: no bound untethered xlinks\n");
    //		exit(1);
  }
}

void KinesinManagement::KMC_Tether_Bound() {

  Update_Bound();
  if (n_bound_untethered_ > 0) {
    // Randomly pick a motor
    int i_motor = properties_->gsl.GetRanInt(n_bound_untethered_);
    Kinesin *motor = std::get<POP_T *>(bound_untethered_[i_motor])->motor_;
    // Update motor's neighbor xlinks (i.e., those it can teth to)
    motor->UpdateNeighborXlinks();
    int attempts = 0;
    bool neighbors_exist = true;
    // Ensure we get a motor that has eligible neighbors
    while (motor->n_neighbor_xlinks_ == 0) {
      if (attempts > 10 * n_bound_untethered_) {
        neighbors_exist = false;
        break;
      }
      i_motor = properties_->gsl.GetRanInt(n_bound_untethered_);
      motor = std::get<POP_T *>(bound_untethered_[i_motor])->motor_;
      motor->UpdateNeighborXlinks();
      attempts++;
    }
    if (neighbors_exist) {
      AssociatedProtein *xlink = motor->GetWeightedNeighborXlink();
      // Update motor and xlink details
      motor->xlink_ = xlink;
      motor->tethered_ = true;
      xlink->tethered_ = true;
      xlink->motor_ = motor;
      motor->UpdateExtension();
      if (motor->tethered_ == false) {
        printf("what in tether_boundnation\n");
        exit(1);
      }
    } else {
      printf("Failed to tether bound motor\n");
    }
  } else {
    printf("Error in Tether_Bound: no bound untethered motors!\n");
    //		exit(1);
  }
}

void KinesinManagement::KMC_Untether_Free() {

  Update_Free_Teth();
  if (n_free_tethered_ > 0) {
    // Randomly pick a motor
    int i_entry = properties_->gsl.GetRanInt(n_free_tethered_);
    Kinesin *motor = free_tethered_[i_entry];
    AssociatedProtein *xlink = motor->xlink_;
    // Update motor and xlink detail
    xlink->motor_ = nullptr;
    xlink->tethered_ = false;
    motor->xlink_ = nullptr;
    motor->tethered_ = false;
    // Remove this motor from active_, replace with last entry
    // (only if there are more than 1 active motor)
    if (n_active_ > 1) {
      int this_index = motor->active_index_;
      Kinesin *last_entry = active_[n_active_ - 1];
      last_entry->active_index_ = this_index;
      active_[this_index] = last_entry;
    }
    n_active_--;
  } else {
    printf("Error in Untether_Free: no free tethered motors!\n");
  }
}

void KinesinManagement::KMC_Untether_Bound(int x_dub) {

  Update_Bound_Teth();
  int n_bound_tethered = n_bound_tethered_[x_dub];
  if (n_bound_tethered > 0) {
    // Randomly pick a motor
    int i_entry = properties_->gsl.GetRanInt(n_bound_tethered);
    Kinesin *motor = bound_tethered_[x_dub][i_entry];
    if (motor->x_dist_doubled_ != x_dub) {
      printf("error in Untether_Bound (motor) x_dub: ");
      printf("%i received; %i in motor\n", x_dub, motor->x_dist_doubled_);
      exit(1);
    }
    AssociatedProtein *xlink = motor->xlink_;
    // Update motor and xlink details
    xlink->motor_ = nullptr;
    xlink->tethered_ = false;
    motor->xlink_ = nullptr;
    motor->tethered_ = false;
    // Update statistics
    motor->UpdateExtension();
  } else {
    printf("Error in Untether_Bound: no bound tethered motors!\n");
    //		exit(1);
  }
}
*/