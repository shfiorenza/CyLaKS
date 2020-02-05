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
  // Calculate total number of motors to have in reservoir
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

  // Lambda = 0.0 means all energy dependence is in binding
  // Lambda = 0.5 means energy dependence is equal for binding and unbinding
  // Lambda = 1.0 means all energy dependence is in unbinding
  double lambda_teth{0.5};    // For Boltzmann factor of tethering mechanisms
  double lambda_neighb{1.0};  // For Boltzmann factor of neighbor cooperativity
  double lambda_lattice{0.5}; // For Boltzmann factor of lattice interactions

  // Get general sim params
  double delta_t = parameters_->delta_t;
  double kbT{parameters_->kbT};
  double r_y{parameters_->microtubules.y_dist / 2};
  double site_size{parameters_->microtubules.site_size};
  // Get base tether parameters
  rest_dist_ = motors_[0].rest_dist_;
  comp_cutoff_ = motors_[0].comp_cutoff_;
  teth_cutoff_ = motors_[0].teth_cutoff_;
  double r_0{motors_[0].r_0_};
  double k_teth{motors_[0].k_spring_};
  double k_slack{motors_[0].k_slack_};
  // Report tether parameters if active
  if (parameters_->motors.tethers_active and parameters_->xlinks.c_bulk > 0.0) {
    wally_->Log("\nFor motors:\n");
    wally_->Log("  rest_dist is %g\n", rest_dist_);
    wally_->Log("  comp_cutoff is %i\n", comp_cutoff_);
    wally_->Log("  dist_cutoff is %i\n", teth_cutoff_);
  } else {
    wally_->Log("Tethering is disabled for motors.\n");
  }
  // Array of tether energies (in kbT) for any given extension
  std::vector<double> tether_energy(2 * teth_cutoff_ + 1, 0.0);
  for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
    double r_x{(double)x_dub * site_size / 2};
    double r{sqrt(r_x * r_x + r_y * r_y)};
    double dr{r - r_0};
    double k_spring{0.0};
    if (dr < 0) {
      k_spring = k_slack;
    } else {
      k_spring = k_teth;
    }
    tether_energy[x_dub] = (0.5 * k_spring * dr * dr) / kbT;
  }
  // Array of interaction energies (in kbT) due to neighbor cooperativity
  double int_energy{-1 * parameters_->motors.interaction_energy}; // kbT
  std::vector<double> neighb_energy(max_neighbs_ + 1, 0.0);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    neighb_energy[n_neighbs] = n_neighbs * int_energy;
  }

  // Get lattice cooperativity jawn - preliminary
  int range{(int)parameters_->motors.lattice_coop_range}; // in n_sites
  double amp{parameters_->motors.lattice_coop_amp}; // amplitude of Gaussian
  // Set the range of our gaussian to correspond to +/- 3 sigma
  double sigma{double(range) / 3};
  double sites_per_bin{0};
  if (n_affinities_ > 1) {
    sites_per_bin = (double)range / (n_affinities_ - 1);
  }
  // Check the no. of sites per bin is an integer
  if (sites_per_bin != (int)sites_per_bin) {
    printf("fix your stupid lattice coop ya dingus\n");
    exit(1);
  }
  // Array of affinity weights (unitless) due to lattice deformations
  std::vector<double> wt_lattice_bind(n_affinities_, 0.0);
  std::vector<double> wt_lattice_unbind(n_affinities_, 0.0);
  for (int i_aff{0}; i_aff < n_affinities_ - 1; i_aff++) {
    double dist{sites_per_bin * i_aff};
    printf("dist is %g\n", dist);
    double gauss{amp * exp(-1 * dist * dist / (2 * sigma * sigma))};
    wt_lattice_bind[i_aff] = 1.0 + gauss;
    wt_lattice_unbind[i_aff] = 1.0 / (1.0 + gauss);
    printf("wt_lattice_bind[%i] is %g\n", i_aff, wt_lattice_bind[i_aff]);
    printf("wt_lattice_unbind[%i] is %g\n", i_aff, wt_lattice_unbind[i_aff]);
  }
  wt_lattice_bind[n_affinities_ - 1] = 1.0;
  wt_lattice_unbind[n_affinities_ - 1] = 1.0;
  printf("wt_lattice_bind[%i] is %g\n", n_affinities_ - 1,
         wt_lattice_bind[n_affinities_ - 1]);
  printf("wt_lattice_unbind[%i] is %g\n", n_affinities_ - 1,
         wt_lattice_unbind[n_affinities_ - 1]);

  // Base motor-MT binding rate for first head without any other effects
  double k_on{parameters_->motors.k_on};
  double c_motor{parameters_->motors.c_bulk};
  double p_bind_i{k_on * c_motor * delta_t};
  // Motor-ATP binding rate for NULL-bound motor heads
  double k_on_ATP{parameters_->motors.k_on_ATP};
  double c_ATP{parameters_->motors.c_ATP};
  p_bind_ATP_ = k_on_ATP * c_ATP * delta_t;
  // ATP hydrolysis rate for ATP-bound motor heads
  double k_hydrolyze{parameters_->motors.k_hydrolyze};
  p_hydrolyze_ = k_hydrolyze * delta_t;
  // Same as above but for *stalled* motors
  double k_hydrolyze_stalled{parameters_->motors.k_hydrolyze_stalled};
  p_hydrolyze_stalled_ = k_hydrolyze_stalled * delta_t;
  // Base motor-MT binding rate for second head without any other effects
  double c_eff{parameters_->motors.c_eff_bind};
  double p_bind_ii{k_on * c_eff * delta_t};
  // Base MT unbinding rate for doubly-bound motor heads w/o any other effects
  double k_off_ii{parameters_->motors.k_off_ii};
  double p_unbind_ii{k_off_ii * delta_t};
  // Base MT unbinding rate for singly-bound motor heads w/o any other effects
  double k_off_i{parameters_->motors.k_off_i};
  double p_unbind_i{k_off_i * delta_t};
  // Same as above but for *stalled* singly-bound motors
  double k_off_i_st{parameters_->motors.k_off_i_stalled};
  double p_unbind_i_st{k_off_i_st * delta_t};
  // Sound the alarm if our timestep is too large
  if (p_bind_ATP_ > 1)
    printf("WARNING: p_bind_ATP=%g for motors\n", p_bind_ATP_);
  if (p_hydrolyze_ > 1)
    printf("WARNING: p_hydro=%g for motors\n", p_hydrolyze_);
  // Force perpetually applied to motors, e.g., by optical trapping
  if (parameters_->motors.applied_force > 0.0) {
    double old_prob{p_unbind_i};
    double sigma_off = 1.2;
    double app_force = parameters_->motors.applied_force;
    double weight = exp(app_force * sigma_off / kbT);
    p_unbind_i = k_off_i * weight * delta_t;
    printf("p_unbind_i scaled from %g to %g \n", old_prob, p_unbind_i);
  }
  p_bind_i_.resize(n_affinities_);
  p_bind_ii_.resize(n_affinities_);
  p_unbind_ii_.resize(n_affinities_);
  p_unbind_i_.resize(n_affinities_);
  p_unbind_i_stalled_.resize(n_affinities_);
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    p_bind_i_[i_aff].resize(max_neighbs_ + 1);
    p_bind_ii_[i_aff].resize(max_neighbs_ + 1);
    p_unbind_ii_[i_aff].resize(max_neighbs_ + 1);
    p_unbind_i_[i_aff].resize(max_neighbs_ + 1);
    p_unbind_i_stalled_[i_aff].resize(max_neighbs_ + 1);
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      // dU = U_f - U_i
      double dU_bind{neighb_energy[n_neighbs] - 0.0};
      double dU_unbind{0.0 - neighb_energy[n_neighbs]};
      // Weights for neighbor interaction energy via Boltzmann factor
      double wt_neighb_bind{exp(-(1.0 - lambda_neighb) * dU_bind)};
      double wt_neighb_unbind{exp(-lambda_neighb * dU_unbind)};
      // Total weights that include lattice deformation
      double wt_tot_bind{wt_neighb_bind * wt_lattice_bind[i_aff]};
      double wt_tot_unbind{wt_neighb_unbind * wt_lattice_unbind[i_aff]};
      p_bind_i_[i_aff][n_neighbs] = p_bind_i * wt_tot_bind;
      if (p_bind_i_[i_aff][n_neighbs] > 1) {
        printf("WARNING: p_bind_i_%i_%i = %g\n", i_aff, n_neighbs,
               p_bind_i_[i_aff][n_neighbs]);
      }
      p_bind_ii_[i_aff][n_neighbs] = p_bind_ii * wt_tot_bind;
      if (p_bind_ii_[i_aff][n_neighbs] > 1) {
        printf("WARNING: p_bind_ii_%i_%i = %g\n", i_aff, n_neighbs,
               p_bind_ii_[i_aff][n_neighbs]);
      }
      p_unbind_ii_[i_aff][n_neighbs] = p_unbind_ii * wt_tot_unbind;
      if (p_unbind_ii_[i_aff][n_neighbs] > 1) {
        printf("WARNING: p_unbind_ii_%i_%i = %g\n", i_aff, n_neighbs,
               p_unbind_ii_[i_aff][n_neighbs]);
      }
      p_unbind_i_[i_aff][n_neighbs] = p_unbind_i * wt_tot_unbind;
      if (p_unbind_i_[i_aff][n_neighbs] > 1) {
        printf("WARNING: p_unbind_i_%i_%i = %g\n", i_aff, n_neighbs,
               p_unbind_i_[i_aff][n_neighbs]);
      }
      p_unbind_i_stalled_[i_aff][n_neighbs] = p_unbind_i_st * wt_tot_unbind;
      if (p_unbind_i_stalled_[i_aff][n_neighbs] > 1) {
        printf("WARNING: p_unbind_ii_st_%i_%i = %g\n", i_aff, n_neighbs,
               p_unbind_i_stalled_[i_aff][n_neighbs]);
      }
    }
  }

  double k_tether{parameters_->motors.k_tether};
  double k_untether{parameters_->motors.k_untether};
  double c_eff_teth{parameters_->motors.c_eff_tether};
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
  double p_teth_base = k_tether * c_eff_teth * delta_t;
  double p_unteth_base = p_untether_free_;
  p_bind_ATP_tethered_.resize(2 * teth_cutoff_ + 1);
  p_bind_ii_tethered_.resize(2 * teth_cutoff_ + 1);
  p_unbind_i_tethered_.resize(2 * teth_cutoff_ + 1);
  p_unbind_i_teth_st_.resize(2 * teth_cutoff_ + 1);
  p_untether_bound_.resize(2 * teth_cutoff_ + 1);
  for (int x_dub = 0; x_dub <= 2 * teth_cutoff_; x_dub++) {
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
      k = k_teth;
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
    p_bind_ii_tethered_[x_dub] = p_bind_ii;
    p_unbind_i_tethered_[x_dub] = weight_alt * p_unbind_i;
    p_unbind_i_teth_st_[x_dub] = weight_alt * p_unbind_i_st;
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
    n_docked_[i_aff].resize(max_neighbs_ + 1);
    n_bound_ADPP_i_[i_aff].resize(max_neighbs_ + 1);
    n_bound_ADPP_i_stalled_[i_aff].resize(max_neighbs_ + 1);
    n_bound_ADPP_ii_[i_aff].resize(max_neighbs_ + 1);
    docked_[i_aff].resize(max_neighbs_ + 1);
    bound_ADPP_i_[i_aff].resize(max_neighbs_ + 1);
    bound_ADPP_i_stalled_[i_aff].resize(max_neighbs_ + 1);
    bound_ADPP_ii_[i_aff].resize(max_neighbs_ + 1);
    for (int i_neighb{0}; i_neighb <= max_neighbs_; i_neighb++) {
      n_docked_[i_aff][i_neighb] = 0;
      n_bound_ADPP_i_[i_aff][i_neighb] = 0;
      n_bound_ADPP_i_stalled_[i_aff][i_neighb] = 0;
      n_bound_ADPP_ii_[i_aff][i_neighb] = 0;
      docked_[i_aff][i_neighb].resize(n_motors_);
      bound_ADPP_i_[i_aff][i_neighb].resize(n_motors_);
      bound_ADPP_i_stalled_[i_aff][i_neighb].resize(n_motors_);
      bound_ADPP_ii_[i_aff][i_neighb].resize(n_motors_);
    }
  }
  // Two dimensional stuff - tethers
  n_docked_tethered_.resize(2 * teth_cutoff_ + 1);
  n_bound_NULL_tethered_.resize(2 * teth_cutoff_ + 1);
  n_bound_ADPP_i_tethered_.resize(2 * teth_cutoff_ + 1);
  n_bound_ADPP_i_teth_st_.resize(2 * teth_cutoff_ + 1);
  n_bound_tethered_.resize(2 * teth_cutoff_ + 1);
  docked_tethered_.resize(2 * teth_cutoff_ + 1);
  bound_NULL_tethered_.resize(2 * teth_cutoff_ + 1);
  bound_ADPP_i_tethered_.resize(2 * teth_cutoff_ + 1);
  bound_ADPP_i_teth_st_.resize(2 * teth_cutoff_ + 1);
  bound_tethered_.resize(2 * teth_cutoff_ + 1);
  for (int x_dub = 0; x_dub <= 2 * teth_cutoff_; x_dub++) {
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
    for (int i_neighb{0}; i_neighb <= max_neighbs_; i_neighb++) {
      std::string affinity = std::to_string(i_aff);
      std::string n_neighbs = std::to_string(i_neighb);
      std::string name = "bind_i_" + affinity + "_" + n_neighbs;
      std::string tar = "unnoc_" + affinity + "_" + n_neighbs;
      events_.emplace_back(
          this, id++, name, Vec<Str>{tar}, p_bind_i_[i_aff][i_neighb],
          &properties_->microtubules.n_unoccupied_mot_[i_aff][i_neighb],
          &properties_->microtubules.unoccupied_list_mot_[i_aff][i_neighb],
          update_unnoc, exe_bind_i, binomial, ran_int);
    }
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
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    for (int i_neighb{0}; i_neighb <= max_neighbs_; i_neighb++) {
      std::string affinity = std::to_string(i_aff);
      std::string n_neighbs = std::to_string(i_neighb);
      std::string name = "bind_ii_" + affinity + "_" + n_neighbs;
      std::string tar1 = "docked_" + affinity + "_" + n_neighbs;
      std::string tar2 = "bound_ADPP_i_" + affinity + "_" + n_neighbs;
      events_.emplace_back(this, id++, name, Vec<Str>{tar1, tar2}, false,
                           p_bind_ii_[i_aff][i_neighb],
                           &n_docked_[i_aff][i_neighb],
                           &docked_[i_aff][i_neighb], update_docked,
                           exe_bind_ii, binomial, ran_int);
    }
  }
  // Unbind_II: Converts ADPP to ADP and unbinds a doubly-bound head
  auto update_ADPP_ii = [&]() { Update_Bound_ADPP_II(); };
  auto exe_unbind_ii = [&](ENTRY_T target) {
    POP_T *head = std::get<POP_T *>(target);
    if (head != nullptr) {
      ReportExecutionOf("unbind_ii");
      KMC_Unbind_II(head);
    } else {
      ReportFailureOf("unbind_ii");
      printf("step #%i!!\n", properties_->current_step_);
      exit(1);
    }
  };
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    for (int i_neighb{0}; i_neighb <= max_neighbs_; i_neighb++) {
      std::string affinity = std::to_string(i_aff);
      std::string n_neighbs = std::to_string(i_neighb);
      std::string name = "unbind_ii_" + affinity + "_" + n_neighbs;
      std::string tar = "bound_ADPP_ii_" + affinity + "_" + n_neighbs;
      events_.emplace_back(
          this, id++, name, Vec<Str>{tar}, p_unbind_ii_[i_aff][i_neighb],
          &n_bound_ADPP_ii_[i_aff][i_neighb], &bound_ADPP_ii_[i_aff][i_neighb],
          update_ADPP_ii, exe_unbind_ii, binomial, ran_int);
    }
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
    for (int i_neighb{0}; i_neighb <= max_neighbs_; i_neighb++) {
      std::string affinity = std::to_string(i_aff);
      std::string n_neighbs = std::to_string(i_neighb);
      std::string name = "unbind_i_" + affinity + "_" + n_neighbs;
      std::string tar = "bound_ADPP_i_" + affinity + "_" + n_neighbs;
      events_.emplace_back(
          this, id++, name, Vec<Str>{tar}, p_unbind_i_[i_aff][i_neighb],
          &n_bound_ADPP_i_[i_aff][i_neighb], &bound_ADPP_i_[i_aff][i_neighb],
          update_ADPP_i, exe_unbind_i, binomial, ran_int);
    }
  }
  // Unbind_I_Stalled: Same as unbind_I but only targets stalled motors
  auto update_ADPP_i_st = [&]() { Update_Bound_ADPP_I_Stalled(); };
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    for (int i_neighb{0}; i_neighb <= max_neighbs_; i_neighb++) {
      std::string affinity = std::to_string(i_aff);
      std::string n_neighbs = std::to_string(i_neighb);
      std::string name = "unbind_i_st_" + affinity + "_" + n_neighbs;
      std::string tar = "bound_ADPP_i_st_" + affinity + "_" + n_neighbs;
      events_.emplace_back(this, id++, name, Vec<Str>{tar},
                           p_unbind_i_stalled_[i_aff][i_neighb],
                           &n_bound_ADPP_i_stalled_[i_aff][i_neighb],
                           &bound_ADPP_i_stalled_[i_aff][i_neighb],
                           update_ADPP_i_st, exe_unbind_i, binomial, ran_int);
    }
  }
  /*
  if tethers ARE active, serialize all extension-based events
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
  for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
    events_.emplace_back(event(
        index++, 200 + x_dub, "bind_ATP_teth", "bound_NULL_teth", binomial,
        &n_bound_NULL_tethered_[x_dub], p_bind_ATP_tethered_[x_dub]));
    }
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
      events_.emplace_back(event(
          index++, 400 + x_dub, "bind_ii_teth", "bound_ADPP_i_teth",
  binomial, &n_docked_tethered_[x_dub], p_bind_ii_tethered_[x_dub]));
    }
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
      events_.emplace_back(event(
          index++, 600 + x_dub, "unbind_i_teth", "bound_ADPP_i_teth",
  binomial, &n_bound_ADPP_i_tethered_[x_dub], p_unbind_i_tethered_[x_dub]));
    }
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
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
  x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
      events_.emplace_back(
          event(index++, 800 + x_dub, "untether_bound", "bound_tethered",
                binomial, &n_bound_tethered_[x_dub],
  p_untether_bound_[x_dub]));
    }
  }
  */
  // Scan through all events and segregate them based on target population
  for (int i_event{0}; i_event < events_.size(); i_event++) {
    std::vector<EVENT_T *> competitors = {&events_[i_event]};
    // For each event, compare all target populations to those of other events
    for (const auto &tar_pop_main : events_[i_event].targets_) {
      // Only compare to events with higher index to avoid double-counting
      for (int j_event{i_event + 1}; j_event < events_.size(); j_event++) {
        // Look at each target population of this alterate event
        for (const auto &tar_pop_alt : events_[j_event].targets_) {
          // If target populations are equal, add this event to partition
          if (tar_pop_main == tar_pop_alt) {
            competitors.push_back(&events_[j_event]);
          }
        }
      }
    }
    // Only add populations that have more than one competitor
    if (competitors.size() > 1) {
      events_by_pop_.push_back(competitors);
    }
  }

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

  if (lists_up_to_date_) {
    return;
  }
  lists_up_to_date_ = true;
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
    for (int i_neighb{0}; i_neighb <= max_neighbs_; i_neighb++) {
      n_docked_[i_aff][i_neighb] = 0;
    }
  }
  for (int i_motor = 0; i_motor < n_active_; i_motor++) {
    Kinesin *motor = active_[i_motor];
    bool eligible{false};
    if (motor->heads_active_ == 1) {
      if (motor->GetActiveHead()->ligand_ == "ADPP") {
        if (!motor->tethered_) {
          eligible = true;
        }
        // Motors with satellite xlinks behave as if untethered
        else if (motor->xlink_->heads_active_ == 0) {
          eligible = true;
        }
      }
    }
    if (eligible) {
      Tubulin *dock_site = motor->GetDockSite();
      if (dock_site != nullptr) {
        int aff{motor->GetActiveHead()->GetAffinity()};
        // Discount active head as a site neighbor to prevent self-coop
        int neighbs_eff{dock_site->GetKIF4ANeighborCount() - 1};
        int i_entry{n_docked_[aff][neighbs_eff]++};
        docked_[aff][neighbs_eff][i_entry] = motor->GetDockedHead();
      }
    }
  }
}

void KinesinManagement::Update_Docked_Teth() {

  /*
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
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
    if (!motor->tethered_) {
      eligible = true;
    } else if (motor->xlink_->heads_active_ == 0) {
      eligible = true;
    }
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
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
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
    if (!motor->IsStalled()) {
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
    if (motor->IsStalled()) {
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
    for (int i_neighb{0}; i_neighb <= max_neighbs_; i_neighb++) {
      n_bound_ADPP_i_[i_aff][i_neighb] = 0;
    }
  }
  for (int i_motor = 0; i_motor < n_active_; i_motor++) {
    Kinesin *motor = active_[i_motor];
    bool eligible{false};
    if (motor->heads_active_ == 1 and !motor->IsStalled()) {
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
      POP_T *active_head = motor->GetActiveHead();
      if (active_head->ligand_ == "ADPP") {
        int aff = active_head->GetAffinity();
        int neighbs = active_head->GetKIF4ANeighbCount();
        int index = n_bound_ADPP_i_[aff][neighbs]++;
        bound_ADPP_i_[aff][neighbs][index] = active_head;
      }
    }
  }
}

void KinesinManagement::Update_Bound_ADPP_I_Stalled() {

  if (verbose_) {
    printf("starting UPDATE ADPP_I_ST\n");
  }
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    for (int i_neighb{0}; i_neighb <= max_neighbs_; i_neighb++) {
      n_bound_ADPP_i_stalled_[i_aff][i_neighb] = 0;
    }
  }
  for (int i_motor = 0; i_motor < n_active_; i_motor++) {
    Kinesin *motor = active_[i_motor];
    bool eligible{false};
    if (motor->heads_active_ == 1 and motor->IsStalled()) {
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
      POP_T *active_head = motor->GetActiveHead();
      if (active_head->ligand_ == "ADPP") {
        int aff = active_head->GetAffinity();
        int neighbs = active_head->GetKIF4ANeighbCount();
        int index = n_bound_ADPP_i_stalled_[aff][neighbs]++;
        bound_ADPP_i_stalled_[aff][neighbs][index] = active_head;
      }
    }
  }
}

void KinesinManagement::Update_Bound_ADPP_I_Tethered() {

  /*
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
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
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
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
    for (int i_neighb{0}; i_neighb <= max_neighbs_; i_neighb++) {
      n_bound_ADPP_ii_[i_aff][i_neighb] = 0;
    }
  }
  for (int i_motor = 0; i_motor < n_active_; i_motor++) {
    Kinesin *motor = active_[i_motor];
    if (motor->heads_active_ == 2) {
      Kinesin::head *chosen_head{nullptr};
      // If both heads are ADPP-bound, pick rear head to unbind
      // FIXME when attempting to incorporate back-stepping cycle
      if (motor->head_one_.ligand_ == "ADPP" and
          motor->head_two_.ligand_ == "ADPP") {
        if (motor->head_one_.trailing_) {
          chosen_head = &motor->head_one_;
        } else {
          chosen_head = &motor->head_two_;
        }
      } else if (motor->head_one_.ligand_ == "ADPP") {
        chosen_head = &motor->head_one_;
      } else if (motor->head_two_.ligand_ == "ADPP") {
        chosen_head = &motor->head_two_;
      }
      if (chosen_head != nullptr) {
        int aff{chosen_head->GetAffinity()};
        int neighbs{chosen_head->GetKIF4ANeighbCount()};
        int i_entry{n_bound_ADPP_ii_[aff][neighbs]++};
        bound_ADPP_ii_[aff][neighbs][i_entry] = chosen_head;
      }
    }
  }
}

void KinesinManagement::Update_Bound_Teth() {

  /*
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
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

  if (parameters_->motors.c_bulk == 0.0) {
    return;
  }
  if (properties_->current_step_ >= 1417783) {
    // if (properties_->current_step_ >= 35809720) {
    // verbose_ = true;
  }

  // sys_time start = sys_clock::now();
  if (verbose_) {
    printf("Starting Update_All_Lists()\n");
  }
  Update_All_Lists();
  if (verbose_) {
    printf("Finished Update_All_Lists()\n");
  }
  // sys_time finish_list = sys_clock::now();
  // properties_->wallace.t_motors_[1] += (finish_list - start).count();
  if (verbose_) {
    printf("Starting Refresh_Populations()\n");
  }
  Refresh_Populations();
  if (verbose_) {
    printf("Finished Refresh_Populations()\n");
  }
  // sys_time finish_pops = sys_clock::now();
  // properties_->wallace.t_motors_[2] += (finish_pops - finish_list).count();
  if (verbose_) {
    printf("Starting Generate_Execution_Sequence()\n");
  }
  Generate_Execution_Sequence();
  if (verbose_) {
    printf("Finished Generate_Execution_Sequence()\n");
  }
  // sys_time finish_seq = sys_clock::now();
  // properties_->wallace.t_motors_[3] += (finish_seq - finish_pops).count();
  if (verbose_) {
    printf("Starting motor KMC cycle: %lu events to execute\n",
           events_to_exe_.size());
  }
  for (int i_event = 0; i_event < events_to_exe_.size(); i_event++) {
    events_to_exe_[i_event]->Execute();
    lists_up_to_date_ = false;
  }
  if (verbose_) {
    printf("Finished motor KMC cycle\n\n");
  }
  // sys_time finish_all = sys_clock::now();
  // properties_->wallace.t_motors_[4] += (finish_all - finish_seq).count();
  // properties_->wallace.t_motors_[0] += (finish_all - start).count();
}

void KinesinManagement::Refresh_Populations() {

  // If tethers are active, update all tether extensions
  if (parameters_->motors.tethers_active and parameters_->xlinks.c_bulk > 0.0) {
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
    if (entry->heads_active_ == 0) {
      // If tethered, automatically classified as "free_teth"
      if (entry->tethered_) {
        state = std::string("free_teth");
      }
      // Otherwise, homeboy was iced and needs to be removed
      else {
        state = std::string("unbound");
      }
      entry->head_one_.state_ = state;
      entry->head_two_.state_ = state;
    }
    // Singly-bound motors
    else if (entry->heads_active_ == 1) {
      std::string inactive_state;
      POP_T *active_head = entry->GetActiveHead();
      // If active head is not bound to ADPP, inactive head is simply "unbound"
      if (active_head->ligand_ != std::string{"ADPP"}) {
        state = std::string{"bound_"} + active_head->ligand_;
        if (active_head->ligand_ == "ATP" and entry->IsStalled()) {
          state += std::string{"_st"};
        }
        inactive_state = std::string{"unbound"};
      }
      // Otherwise if active head is bound to ADPP, inactive head is docked
      else {
        state = std::string{"bound_ADPP_i"};
        if (entry->IsStalled()) {
          state += std::string{"_st"};
        }
        state += "_";
        state += std::to_string(active_head->GetAffinity());
        state += "_";
        state += std::to_string(active_head->GetKIF4ANeighbCount());
        Tubulin *dock_site = entry->GetDockSite();
        if (dock_site != nullptr) {
          inactive_state = std::string{"docked_"};
          // Use active head's affinity to prevent self-cooperativity
          inactive_state += std::to_string(active_head->GetAffinity());
          inactive_state += "_";
          // Discount active head as a site neighbor to prevent self-coop
          int neighbs_eff = dock_site->GetKIF4ANeighborCount() - 1;
          inactive_state += std::to_string(neighbs_eff);
        } else {
          inactive_state = std::string{"unbound"};
        }
      }
      active_head->state_ = state;
      active_head->GetOtherHead()->state_ = inactive_state;
    }
    // Doubly-bound motors
    else {
      // Do head one first
      if (entry->head_one_.ligand_ != std::string{"ADPP"}) {
        state = std::string{"bound_"} + entry->head_one_.ligand_;
        if (entry->head_one_.ligand_ == "ATP" and entry->IsStalled()) {
          state += std::string{"_st"};
        }
      } else {
        state = std::string{"bound_ADPP_ii"};
        state += "_";
        state += std::to_string(entry->head_one_.GetAffinity());
        state += "_";
        state += std::to_string(entry->head_one_.GetKIF4ANeighbCount());
      }
      entry->head_one_.state_ = state;
      // Head two next
      if (entry->head_two_.ligand_ != std::string{"ADPP"}) {
        state = std::string{"bound_"} + entry->head_two_.ligand_;
        if (entry->head_two_.ligand_ == "ATP" and entry->IsStalled()) {
          state += std::string{"_st"};
        }
      } else {
        state = std::string{"bound_ADPP_ii"};
        state = std::string{"bound_ADPP_ii"};
        state += "_";
        state += std::to_string(entry->head_two_.GetAffinity());
        state += "_";
        state += std::to_string(entry->head_two_.GetKIF4ANeighbCount());
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
      // if (properties_->current_step_ == 249056) {
      //   printf("%i expected for %s\n", n_expected,
      //          events_[i_event].name_.c_str());
      // }
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
      printf("i_array: %i | n_events: %i\n", i_array, n_events);
      printf("step: %i\n", properties_->current_step_);
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
}

int KinesinManagement::Sample_Event_Statistics() {

  // Scan through all events & get expected number of occurrences
  int n_events_tot{0};
  for (auto &&event : events_) {
    n_events_tot += event.SampleStatistics();
    // Ensure no funny business occurs with statistics
    if (event.n_expected_ > *event.n_avail_) {
      printf("Error; %i events expected but only %i available for %s\n",
             event.n_expected_, *event.n_avail_, event.name_.c_str());
    }
    // If verbose flag is active, report statistics
    if (verbose_) {
      printf("Expect %i events for %s (%i available)\n", event.n_expected_,
             event.name_.c_str(), *event.n_avail_);
    }
  }
  /*
  // Ensure that competing events will not turn any populations negative
  for (auto &&competitors : events_by_pop_) {
    if (verbose_) {
      printf("Checking the following competitors:\n");
    }
    int n_active_competitors{0};
    int n_events_loc{0};
    int n_avail_loc{0};
    double p_tot{0.0};
    for (auto &&event_ptr : competitors) {
      if (verbose_) {
        printf("  %s: %i expected\n", event_ptr->name_.c_str(),
               event_ptr->n_expected_);
      }
      n_events_loc += event_ptr->n_expected_;
      p_tot += event_ptr->n_expected_ * event_ptr->p_occur_;
      if (event_ptr->n_expected_ > 0) {
        n_active_competitors++;
      }
      if (event_ptr->is_not_redundant_) {
        n_avail_loc += *event_ptr->n_avail_;
      }
    }
    if (verbose_) {
      printf("%i total events; %i targets available\n", n_events_loc,
             n_avail_loc);
    }
    if (n_avail_loc == 0 and n_events_loc > 0) {
      printf("Error; %i total events, but only %i targets available\n",
             n_events_loc, n_avail_loc);
      for (const auto &event_ptr : competitors) {
        printf("  %s: %i expected (%i avail)\n", event_ptr->name_.c_str(),
               event_ptr->n_expected_, *event_ptr->n_avail_);
      }
      printf("(@ step #%i)\n", properties_->current_step_);
      exit(1);
    }
    while (n_events_loc > n_avail_loc) {
      double p_cum{0.0};
      double ran{properties_->gsl.GetRanProb()};
      for (auto &&event_ptr : competitors) {
        p_cum += event_ptr->n_expected_ * event_ptr->p_occur_ / p_tot;
        if (ran < p_cum) {
          event_ptr->n_expected_--;
          n_events_tot--;
          n_events_loc--;
          p_tot -= event_ptr->p_occur_;
          break;
        }
      }
    }
  }
  */
  // Scan through all target pops. & ensure none will become negative
  for (int i_pop{0}; i_pop < events_by_pop_.size(); i_pop++) {
    int n_competitors = events_by_pop_[i_pop].size();
    // Skip entries that only have 1 competitor; they shouldn't go negative
    if (n_competitors == 1) {
      continue;
    }
    int n_active_competitors{0};
    int n_events_loc{0};
    int n_avail_loc{0};
    double p_tot{0.0};
    if (verbose_) {
      printf("Checking partition #%i: \n", i_pop);
    }
    for (int i_entry{0}; i_entry < n_competitors; i_entry++) {
      if (verbose_) {
        printf("    ");
        std::cout << events_by_pop_[i_pop][i_entry]->name_;
        printf(" is a competitor");
      }
      n_events_loc += events_by_pop_[i_pop][i_entry]->n_expected_;
      if (verbose_) {
        printf(" (%i expected)", events_by_pop_[i_pop][i_entry]->n_expected_);
      }
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

  if (!head->in_scratch_) {
    head->in_scratch_ = true;
    scratch_[n_scratched_++] = head;
    if (verbose_) {
      std::cout << head->state_;
      printf(" was saved to scratch.\n");
    }
  }
  int n_neighbs{0};
  if (head->site_ != nullptr) {
    n_neighbs = head->site_->GetKIF4ANeighborCount();
  }
  if (n_neighbs == 1) {
    POP_T *neighb;
    int i_site = head->site_->index_;
    int mt_length = head->site_->mt_->n_sites_;
    if (i_site == 0) {
      neighb = head->site_->mt_->lattice_[i_site + 1].motor_head_;
    } else if (i_site == mt_length - 1) {
      neighb = head->site_->mt_->lattice_[i_site - 1].motor_head_;
    } else if (head->site_->mt_->lattice_[i_site + 1].motor_head_ != nullptr) {
      neighb = head->site_->mt_->lattice_[i_site + 1].motor_head_;
    } else if (head->site_->mt_->lattice_[i_site - 1].motor_head_ != nullptr) {
      neighb = head->site_->mt_->lattice_[i_site - 1].motor_head_;
    }
    if (!neighb->in_scratch_) {
      neighb->in_scratch_ = true;
      scratch_[n_scratched_++] = neighb;
    }
  } else if (n_neighbs == 2) {
    POP_T *neighb_one, *neighb_two;
    int i_site = head->site_->index_;
    neighb_one = head->site_->mt_->lattice_[i_site + 1].motor_head_;
    if (!neighb_one->in_scratch_) {
      neighb_one->in_scratch_ = true;
      scratch_[n_scratched_++] = neighb_one;
    }
    neighb_two = head->site_->mt_->lattice_[i_site - 1].motor_head_;
    if (!neighb_two->in_scratch_) {
      neighb_two->in_scratch_ = true;
      scratch_[n_scratched_++] = neighb_two;
    }
  }
}

void KinesinManagement::SaveNeighbsToScratch(SITE_T *site) {

  int n_neighbs = site->GetKIF4ANeighborCount();
  if (n_neighbs == 1) {
    POP_T *neighb;
    int i_site = site->index_;
    int mt_length = site->mt_->n_sites_;
    if (i_site == 0) {
      neighb = site->mt_->lattice_[i_site + 1].motor_head_;
    } else if (i_site == mt_length - 1) {
      neighb = site->mt_->lattice_[i_site - 1].motor_head_;
    } else if (site->mt_->lattice_[i_site + 1].motor_head_ != nullptr) {
      neighb = site->mt_->lattice_[i_site + 1].motor_head_;
    } else if (site->mt_->lattice_[i_site - 1].motor_head_ != nullptr) {
      neighb = site->mt_->lattice_[i_site - 1].motor_head_;
    }
    if (!neighb->in_scratch_) {
      neighb->in_scratch_ = true;
      scratch_[n_scratched_++] = neighb;
    }
  } else if (n_neighbs == 2) {
    POP_T *neighb_one, *neighb_two;
    int i_site = site->index_;
    neighb_one = site->mt_->lattice_[i_site + 1].motor_head_;
    if (!neighb_one->in_scratch_) {
      neighb_one->in_scratch_ = true;
      scratch_[n_scratched_++] = neighb_one;
    }
    neighb_two = site->mt_->lattice_[i_site - 1].motor_head_;
    if (!neighb_two->in_scratch_) {
      neighb_two->in_scratch_ = true;
      scratch_[n_scratched_++] = neighb_two;
    }
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
}

void KinesinManagement::KMC_Bind_I(SITE_T *site) {

  // Save neighbors of unbound site to scratch_ before modifying them
  SaveNeighbsToScratch(site);
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
  properties_->microtubules.FlagForUpdate();
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
    head->motor_->ChangeConformation();
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
    SaveToScratch(head->GetOtherHead());
    // Verify that proposed site is unoccupied
    Tubulin *dock_site = head->motor_->GetDockSite();
    // Save neighbors of dock_site to scratch_ before modifying them
    SaveNeighbsToScratch(dock_site);
    if (dock_site == nullptr) {
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
    properties_->microtubules.FlagForUpdate();
  } else if (head->stored_dock_site_ != nullptr) {
    printf("YO YO YO THE STORED SITE WORKS\n");
    // Save active head to scratch before binding
    SaveToScratch(head->GetOtherHead());
    Tubulin *dock_site = head->stored_dock_site_;
    if (dock_site->occupied_) {
      printf("Error in Bind_II: option 2\n");
      exit(1);
    }
    SaveNeighbsToScratch(dock_site);
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
    properties_->microtubules.FlagForUpdate();
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
    SaveToScratch(head->GetOtherHead());
    SaveNeighbsToScratch(head->site_);
    // Update site
    head->site_->occupied_ = false;
    head->site_->motor_head_ = nullptr;
    // Update motor
    head->site_ = nullptr;
    head->ligand_ = "ADP";
    head->motor_->heads_active_--;
    properties_->microtubules.FlagForUpdate();
    if (head->motor_->frustrated_) {
      if (head->trailing_) {
        head->motor_->ChangeConformation();
      } else {
        head->motor_->frustrated_ = false;
      }
    }
    if (head->motor_->frustrated_) {
      printf("somethin wicked in motor unbind_ii\n");
      exit(1);
    }
  } else {
    printf("WACK in unbind_ii\n");
    // exit(1);
  }
}

void KinesinManagement::KMC_Unbind_I(POP_T *head) {

  if (head->motor_->heads_active_ == 1) {
    POP_T *docked_head = head->motor_->StoreDockSite();
    if (docked_head != nullptr) {
      SaveToScratch(docked_head);
    }
    SaveNeighbsToScratch(head->site_);
    // Update site
    head->site_->occupied_ = false;
    head->site_->motor_head_ = nullptr;
    // Update motor
    head->site_ = nullptr;
    head->ligand_ = "ADP";
    head->motor_->heads_active_--;
    head->motor_->mt_ = nullptr;
    properties_->microtubules.FlagForUpdate();
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
