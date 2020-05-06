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
  CalculateCutoffs();
  SetParameters();
  GenerateMotors();
  InitializeLists();
  if (wally_->test_mode_ == nullptr) {
    InitializeEvents();
  } else {
    InitializeTestEnvironment();
    InitializeTestEvents();
  }
}

void KinesinManagement::CalculateCutoffs() {

  double kbT{parameters_->kbT};
  double r_0{parameters_->motors.r_0};
  double r_y{parameters_->microtubules.y_dist / 2};
  double k_spring{parameters_->motors.k_spring};
  double k_slack{parameters_->motors.k_slack};
  double site_size{parameters_->microtubules.site_size};
  lattice_E_0_solo_ = -1 * parameters_->motors.lattice_coop_Emax_solo;
  lattice_cutoff_ = parameters_->motors.lattice_coop_range;
  // Set lattice_alpha_ so that E = 0 at site immediately after cutoff
  double dx_post_cutoff{(lattice_cutoff_ + 1) * site_size};
  lattice_alpha_ = -1 * lattice_E_0_solo_ / (dx_post_cutoff * dx_post_cutoff);
  // Calculate rest_dist_ in number of sites
  int approx_rest{int(sqrt(r_0 * r_0 - r_y * r_y) / site_size)};
  double rest_scan[3];
  double scan_force[3];
  for (int i_scan = -1; i_scan <= 1; i_scan++) {
    rest_scan[i_scan + 1] = approx_rest + ((double)i_scan * 0.5);
    double rest_scan_length = rest_scan[i_scan + 1] * site_size;
    double r_scan = sqrt(r_y * r_y + rest_scan_length * rest_scan_length);
    if (r_scan >= r_0)
      scan_force[i_scan + 1] = (r_scan - r_0) * k_spring;
    else
      scan_force[i_scan + 1] = (r_scan - r_0) * k_slack;
  }
  double min_force{1000};
  for (int i_scan = -1; i_scan <= 1; i_scan++) {
    double force = fabs(scan_force[i_scan + 1]);
    if (force < min_force) {
      min_force = force;
      rest_dist_ = rest_scan[i_scan + 1];
    }
  }
  // Next, calculate compression distance cutoff
  comp_cutoff_ = 0;
  for (int x_dub{int(2 * rest_dist_)}; x_dub >= 0; x_dub--) {
    // increment by 0.5x
    double r_x{x_dub * site_size / 2};
    double r{sqrt(r_y * r_y + r_x * r_x)};
    double dr{r - r_0};
    double spring_energy{0.0};
    if (dr < 0) {
      spring_energy = 0.5 * k_slack * dr * dr;
    } else {
      spring_energy = 0.5 * k_spring * dr * dr;
    }
    double boltzmann_weight = exp(0.5 * spring_energy / kbT);
    if (boltzmann_weight > 100) {
      comp_cutoff_ = x_dub / 2;
      break;
    }
  }
  // Finally, calculate extension distance cutoff
  for (int x_dub{int(2 * rest_dist_)}; x_dub < 1000; x_dub++) {
    // increment by 0.5x
    double r_x{x_dub * site_size / 2};
    double r{sqrt(r_y * r_y + r_x * r_x)};
    double dr{r - r_0};
    double spring_energy{0.0};
    if (dr < 0) {
      spring_energy = 0.5 * k_slack * dr * dr;
    } else {
      spring_energy = 0.5 * k_spring * dr * dr;
    }
    double boltzmann_weight = exp(0.5 * spring_energy / kbT);
    if (boltzmann_weight > 100) {
      teth_cutoff_ = x_dub / 2;
      break;
    }
  }
  wally_->Log("\nFor motors:\n");
  if (parameters_->motors.tethers_active and parameters_->xlinks.c_bulk > 0.0) {
    tethering_active_ = true;
    wally_->Log("  rest_dist is %g\n", rest_dist_);
    wally_->Log("  comp_cutoff is %i\n", comp_cutoff_);
    wally_->Log("  dist_cutoff is %i\n", teth_cutoff_);
  } else {
    wally_->Log("  Tethering is disabled.\n");
  }
  if (lattice_cutoff_ > 0) {
    lattice_coop_active_ = true;
    wally_->Log("  lattice_cutoff_ is %i\n", lattice_cutoff_);
    wally_->Log("  lattice_alpha_ is %g\n", lattice_alpha_);
  } else {
    wally_->Log("  Lattice cooperativity is disabled.\n");
  }
}

void KinesinManagement::SetParameters() {

  t_active_ = parameters_->motors.t_active;
  /*
    For events that result in a change in energy, we use Boltzmann factors to
    scale rates appropriately. Detailed balance is satisfied with the factors:
                  exp{-(1 - lambda)*(delta_E)/(kB*T)}, and
                  exp{-(lambda)*(delta_E)/(kB*T)}
    for an event and its complement (e.g., binding and unbinding), where lambda
    is a constant that ranges from 0 to 1, delta_E is the change in energy that
    results from the event, kB is Boltzmann's constant, and T is the temperature
 */
  /*
    For tethered motor binding & unbinding, we use force-dependence instead:
                exp{-(1 - lambda)*(f_teth * delta_off)/(kB*T)}, and
                exp{lambda*(f_teth * delta_off)/(kB*T)}
    for binding and unbinding, respectively.
  */
  /*
     For stepping rates (ATP binding), we use a linear scaling based on the
     motor's so-called "stall force," the force at which it's velocity is 0:
                p_step -> p_step * (1 - f_x / f_stall)
  */
  // Lambda = 0.5 means energy dependence is equal for binding and unbinding
  double lambda_lattice{0.5}; // Lambda for lattice interactions
  // Lambda = 1.0 means all energy dependence is in unbinding
  double lambda_neighb{1.0}; // Lambda for neighbor cooperativity
  // Lambda = 0.5 means energy dependence is equal for binding and unbinding
  double lambda_teth{0.5}; // Lambda for tethering mechanisms

  // Array Boltzmann weights due to neighbor-neighbor interactions
  double energy_per_neighb{-1 * parameters_->motors.interaction_energy}; // kbT
  weight_neighbs_bind_.resize(max_neighbs_ + 1);
  weight_neighbs_unbind_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    double energy{n_neighbs * energy_per_neighb}; // ! in kbT
    weight_neighbs_bind_[n_neighbs] = exp(-(1.0 - lambda_neighb) * energy);
    weight_neighbs_unbind_[n_neighbs] = exp(lambda_neighb * energy);
  }
  // Array of binding/unbinding weights due to lattice deformation from a SINGLE
  // motor; multiplied together to get total weight for any arrangement
  double site_size{parameters_->microtubules.site_size};
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
  weight_lattice_bind_.resize(lattice_cutoff_ + 1);
  weight_lattice_unbind_.resize(lattice_cutoff_ + 1);
  for (int delta{0}; delta <= lattice_cutoff_; delta++) {
    double dx{delta * site_size};
    double energy{lattice_alpha_ * dx * dx + lattice_E_0_solo_}; // in kbT
    weight_lattice_bind_[delta] = exp(-(1 - lambda_lattice) * energy);
    weight_lattice_unbind_[delta] = exp(lambda_lattice * energy);
  }
  /*
  // Array of tether energies (in kbT) for any given extension
  double k_slack{parameters_->motors.k_slack};
  double k_teth{parameters_->motors.k_spring};
  double r_y{parameters_->microtubules.y_dist / 2};
  double r_0{parameters_->motors.r_0};
  double kbT{parameters_->kbT};
  std::vector<double> tether_energy(2 * teth_cutoff_ + 1, 0.0);
  std::vector<double> tether_force(2 * teth_cutoff_ + 1, 0.0);
  std::vector<double> tether_cosine(2 * teth_cutoff_ + 1, 0.0);
  for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
    double r_x{(double)x_dub * site_size / 2};
    double r{sqrt(r_x * r_x + r_y * r_y)};
    double dr{r - r_0};
    double k_spring{0.0};
    if (dr < 0.0) {
      k_spring = k_slack;
    } else {
      k_spring = k_teth;
    }
    tether_energy[x_dub] = 0.5 * k_spring * dr * dr / kbT; // ! in kbT
    tether_force[x_dub] = fabs(k_spring * dr);
    tether_cosine[x_dub] = r_x / r;
  }
  */
  double delta_t{parameters_->delta_t};
  double k_on{parameters_->motors.k_on};
  double c_bulk{parameters_->motors.c_bulk};
  p_avg_bind_i_ = k_on * c_bulk * delta_t;
  double c_eff_bind{parameters_->motors.c_eff_bind};
  p_avg_bind_ii_ = k_on * c_eff_bind * delta_t;
  double k_on_ATP{parameters_->motors.k_on_ATP};
  double c_ATP{parameters_->motors.c_ATP};
  p_bind_ATP_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    double weight{weight_neighbs_unbind_[n_neighbs]};
    if (n_neighbs == 2) {
      weight = weight_neighbs_unbind_[1];
    }
    p_bind_ATP_[n_neighbs] = k_on_ATP * c_ATP * delta_t * weight;
  }
  double k_hydrolyze{parameters_->motors.k_hydrolyze};
  p_hydrolyze_ = k_hydrolyze * delta_t;
  double k_off_ii{parameters_->motors.k_off_ii};
  p_avg_unbind_ii_ = k_off_ii * delta_t;
  double k_off_i{parameters_->motors.k_off_i};
  p_avg_unbind_i_ = k_off_i * delta_t;
  /*
  double c_eff_teth{parameters_->motors.c_eff_tether};
  double k_tether{parameters_->motors.k_tether};
  double k_untether{parameters_->motors.k_untether};
  p_tether_free_ = k_tether * c_bulk * delta_t;
  p_untether_satellite_ = k_untether * delta_t;
  p_avg_bind_i_teth_ = k_on * c_eff_teth * delta_t;
  p_avg_tether_bound_ = k_tether * c_eff_teth * delta_t;
  p_untether_bound_.resize(2 * teth_cutoff_ + 1);
  weight_tether_bound_.resize(2 * teth_cutoff_ + 1);
  for (int x_dub{0}; x_dub <= 2 * teth_cutoff_; x_dub++) {
    double dU_teth{tether_energy[x_dub] - 0.0};
    double dU_unteth{0.0 - tether_energy[x_dub]};
    double wt_tether{exp(-(1.0 - lambda_teth) * dU_teth)};
    double wt_untether{exp(-lambda_teth * dU_unteth)};
    if (x_dub < 2 * comp_cutoff_) {
      wt_tether = 0.0;
      wt_untether = 0.0;
    }
    p_untether_bound_[x_dub] = p_untether_satellite_ * wt_untether;
    weight_tether_bound_[x_dub] = p_avg_tether_bound_ * wt_tether;
  }
  */
  /*
  double f_stall{parameters_->motors.stall_force};
  double p_bind_ATP{k_on_ATP * c_ATP * delta_t};
  p_bind_ATP_.resize(max_neighbs_ + 1);
  p_bind_ATP_to_teth_.resize(max_neighbs_ + 1);
  p_bind_ATP_fr_teth_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    double dU_unbind{0.0 - neighb_energy[n_neighbs]};
    double wt_neighb_unbind{exp(-lambda_neighb * dU_unbind)};
    p_bind_ATP_[n_neighbs] = p_bind_ATP * wt_neighb_unbind;
    p_bind_ATP_to_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    p_bind_ATP_fr_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    for (int x_dub{0}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      double f_x{tether_force[x_dub] * tether_cosine[x_dub]};
      double wt_fr{sqrt(1 - f_x / f_stall)};
      double wt_to{1.0};
      // Set weights to zero for invalid x_dub ranges
      if (x_dub < 2 * comp_cutoff_ or x_dub > 2 * teth_cutoff_) {
        wt_to = 0.0;
        wt_fr = 0.0;
      }
      if (x_dub <= 2 * comp_cutoff_ + 1 or x_dub >= 2 * teth_cutoff_ - 1) {
        wt_fr = 0.0;
      }
      p_bind_ATP_to_teth_[n_neighbs][x_dub] = p_bind_ATP_[n_neighbs] * wt_to;
      p_bind_ATP_fr_teth_[n_neighbs][x_dub] = p_bind_ATP_[n_neighbs] * wt_fr;
    }
  }
  */
  // Force perpetually applied to motors, e.g., by optical trapping
  if (parameters_->motors.applied_force > 0.0) {
    wally_->ErrorExit("APPLIED FORCE COMMENTED OUT IN K_MGMT!!");
    /*
    double force{parameters_->motors.applied_force};
    wally_->Log("\n Motor applied force is %g\n", force);
    double p_unbind_i_old{p_avg_unbind_i_};
    double p_unbind_ii_old{p_avg_unbind_ii_};
    double p_bind_ATP_old{p_bind_ATP_};
    double sigma_off_i{2.0};
    double sigma_off_ii{0.35};
    double sigma_ATP{4.6};
    double kbT{parameters_->kbT};
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
    p_bind_ATP_ *= exp(-1 * force * sigma_ATP / kbT);
    wally_->Log("p_bind_ATP scaled from %g to %g \n", p_bind_ATP_old,
                p_bind_ATP_);
    */
  }
  // Update_Weights();
  p_theory_["hydrolyze"] = {{{p_hydrolyze_}}};
  p_theory_["hydrolyze_st"] = {{{p_hydrolyze_st_}}};
  if (tethering_active_) {
    p_theory_["tether_free"] = {{{p_tether_free_}}};
    p_theory_["untether_satellite"] = {{{p_untether_satellite_}}};
    p_theory_["bind_i_teth"] = {{{p_avg_bind_i_teth_}}};
    p_theory_["tether_bound"] = {{{p_avg_tether_bound_}}};
    p_theory_["untether_bound"] = {{p_untether_bound_}};
  }
  p_theory_["bind_ATP"] = {{p_bind_ATP_}};
  if (tethering_active_) {
    p_theory_["bind_ATP_to_teth"] = {p_bind_ATP_to_teth_};
    p_theory_["bind_ATP_fr_teth"] = {p_bind_ATP_fr_teth_};
  }
  p_theory_["bind_i"] = {{{p_avg_bind_i_}}};
  p_theory_["bind_ii"] = {{{p_avg_bind_ii_}}};
  p_theory_["unbind_ii"] = {{{p_avg_unbind_ii_}}};
  p_theory_["unbind_i"] = {{{p_avg_unbind_i_}}};
  // READ EM OUT
  for (const auto &entry : p_theory_) {
    auto label = entry.first;
    auto value = entry.second;
    for (int i{0}; i < value.size(); i++) {
      for (int j{0}; j < value[i].size(); j++) {
        for (int k{0}; k < value[i][j].size(); k++) {
          std::string name{"p_" + label + "_"};
          if (value.size() > 1) {
            name += "[" + std::to_string(i) + "]";
          }
          if (value[i].size() > 1) {
            name += "[" + std::to_string(j) + "]";
          }
          if (value[i][j].size() > 1) {
            name += "[" + std::to_string(k) + "]";
          }
          // if (verbosity_ >= 3) {
          wally_->Log("%s = %g\n", name.c_str(), value[i][j][k]);
          // }
          if (value[i][j][k] > 1.0) {
            wally_->Log("Error! %s = %g\n", name.c_str(), value[i][j][k]);
            wally_->ErrorExit("Kin_MGMT:SetParameters()");
          }
        }
      }
    }
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
  // bound_NULL_.resize(n_motors_);
  bound_ATP_.resize(n_motors_);
  bound_ATP_st_.resize(n_motors_);
  bound_unteth_.resize(n_motors_);
  satellites_.resize(n_motors_);
  bind_i_teth_candidates_.resize(n_motors_);
  bind_ii_candidates_.resize(n_motors_);
  unbind_ii_candidates_.resize(n_motors_);
  unbind_i_candidates_.resize(n_motors_);
  tether_bound_candidates_.resize(n_motors_);
  // 2-D stuff -- indexed by x_dub & i_motor
  n_bound_teth_.resize(2 * teth_cutoff_ + 1);
  bound_teth_.resize(2 * teth_cutoff_ + 1);
  for (int x_dub{0}; x_dub <= 2 * teth_cutoff_; x_dub++) {
    n_bound_teth_[x_dub] = 0;
    bound_teth_[x_dub].resize(n_motors_);
  }
  // 2-D stuff -- indexed by n_neighbs & i_motor
  n_bound_NULL_.resize(max_neighbs_ + 1);
  bound_NULL_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    n_bound_NULL_[n_neighbs] = 0;
    bound_NULL_[n_neighbs].resize(n_motors_);
  }
  // 3-D stuff -- indexed by n_neighbs, x_dub, & i_motor
  n_bound_NULL_to_teth_.resize(max_neighbs_ + 1);
  n_bound_NULL_fr_teth_.resize(max_neighbs_ + 1);
  bound_NULL_to_teth_.resize(max_neighbs_ + 1);
  bound_NULL_fr_teth_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    n_bound_NULL_to_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    n_bound_NULL_fr_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    bound_NULL_to_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    bound_NULL_fr_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    for (int x_dub{0}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      n_bound_NULL_to_teth_[n_neighbs][x_dub] = 0;
      n_bound_NULL_fr_teth_[n_neighbs][x_dub] = 0;
      bound_NULL_to_teth_[n_neighbs][x_dub].resize(n_motors_);
      bound_NULL_fr_teth_[n_neighbs][x_dub].resize(n_motors_);
    }
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
  auto exe_bind_ATP = [&](ENTRY_T target) {
    try {
      Bind_ATP(std::get<POP_T *>(target));
    } catch (...) {
      wally_->ErrorExit("K_MGMT::exe_bind_ATP");
    }
  };
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    event_name = "bind_ATP_" + std::to_string(n_neighbs);
    events_.emplace_back(event_name, p_bind_ATP_[n_neighbs],
                         &n_bound_NULL_[n_neighbs], &bound_NULL_[n_neighbs],
                         exe_bind_ATP, binomial, set_ran_indices);
  }
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
  events_.emplace_back(event_name, p_avg_bind_ii_, &n_bind_ii_candidates_,
                       &bind_ii_candidates_, exe_bind_ii, weight_bind_ii,
                       poisson, get_ran_prob, set_ran_indices);

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
  events_.emplace_back(event_name, p_avg_unbind_ii_, &n_unbind_ii_candidates_,
                       &unbind_ii_candidates_, exe_unbind_ii, weight_unbind_ii,
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
  events_.emplace_back(event_name, p_avg_unbind_i_, &n_unbind_i_candidates_,
                       &unbind_i_candidates_, exe_unbind_i, weight_unbind_i,
                       poisson, get_ran_prob, set_ran_indices);
  if (verbosity_ >= 3) {
    wally_->Log("\nMotor events: \n");
    for (const auto &event : events_) {
      wally_->Log("   %s\n", event.name_.c_str());
    }
  }
}

void KinesinManagement::InitializeTestEnvironment() {

  properties_->microtubules.InitializeTestEnvironment();
  if (strcmp(wally_->test_mode_, "motor_lattice_bind") == 0) {
    wally_->Log("Initializing test %s\n", wally_->test_mode_);
    MicrotubuleManagement *mts = &properties_->microtubules;
    for (int i_mt{0}; i_mt < mts->mt_list_.size(); i_mt++) {
      for (int i_site{0}; i_site < mts->mt_list_[i_mt].n_sites_; i_site++) {
        Tubulin *site{&mts->mt_list_[i_mt].lattice_[i_site]};
        int n_neighbs{site->GetKif4ANeighborCount()};
        site->weight_bind_ = 1.0;
        site->weight_unbind_ = 1.0;
      }
    }
    // Set p_avg_bind_i_ to a large value for better statistics
    double c_bulk{10.0};
    double k_on{1.0};
    p_avg_bind_i_ = k_on * c_bulk * parameters_->delta_t;
    // Initialize stat tracker
    lattice_bind_stats_.resize(lattice_cutoff_ + 1);
    for (int delta{0}; delta <= lattice_cutoff_; delta++) {
      lattice_bind_stats_[delta].first = 0;
      lattice_bind_stats_[delta].second = 0;
    }
    // Choose middle MT site and bind motor to it
    int i_site{lattice_cutoff_};
    Tubulin *site{&properties_->microtubules.mt_list_[0].lattice_[i_site]};
    Bind_I(site);
    properties_->microtubules.FlagForUpdate();
    properties_->microtubules.UpdateUnoccupied();
    lattice_coop_active_ = false;
  }
  if (strcmp(wally_->test_mode_, "motor_lattice_step") == 0) {
    wally_->Log("Initializing test %s\n", wally_->test_mode_);
    Microtubule *mt = &properties_->microtubules.mt_list_[0];
    // If test_delta_ is left uninitialized, test checks against self-coop
    test_delta_ = 5;
    if (test_delta_ > 0) {
      // Disable automatic weight updates
      lattice_coop_active_ = false;
      // Disable neighbor interactions
      parameters_->motors.interaction_energy = 0.0;
      for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
        weight_neighbs_bind_[n_neighbs] = 1.0;
        weight_neighbs_unbind_[n_neighbs] = 1.0;
      }
      // Set all weights to delta value we wish to test
      for (int i_site{0}; i_site < mt->n_sites_; i_site++) {
        mt->lattice_[i_site].weight_bind_ = weight_lattice_bind_[test_delta_];
        mt->lattice_[i_site].weight_unbind_ =
            weight_lattice_unbind_[test_delta_];
      }
    }
    // Initialize stat trackers
    lattice_step_bind_ii_stats_ = std::make_pair(0, 0);
    lattice_step_unbind_i_stats_ = std::make_pair(0, 0);
    lattice_step_unbind_ii_stats_ = std::make_pair(0, 0);
    // Place motor on minus end of microtubule
    int i_site{mt->minus_end_};
    Tubulin *site{&mt->lattice_[i_site]};
    Bind_I(site);
  }
}

void KinesinManagement::InitializeTestEvents() {

  std::string event_name;
  // Binomial probabilitiy distribution; sampled to predict most events
  auto binomial = [&](double p, int n) {
    if (n == 0) {
      return 0;
    }
    if (p == 1.0) {
      return n;
    }
    return properties_->gsl.SampleBinomialDist(p, n);
  };
  // Function that returns a random probability in the range [0.0, 1.0)
  auto get_ran_prob = [&](void) { return gsl_->GetRanProb(); };
  // Function that sets n random indices from the range [0, m)
  auto set_ran_indices = [&](int *indices, int n, int m) {
    if (m > 1) {
      properties_->gsl.SetRanIndices(indices, n, m);
    } else {
      indices[0] = 0;
    }
  };
  if (strcmp(wally_->test_mode_, "motor_lattice_bind") == 0) {
    // Bind_I: bind first motor head to MT and release ADP
    event_name = "bind_i";
    auto exe_bind_i = [&](ENTRY_T target) {
      Tubulin *site{std::get<SITE_T *>(target)};
      int i_site{site->index_};
      // printf("i_site: %i\n", i_site);
      // 'main' kinesin motor will always be at index_ = lattice_cutoff_
      int delta{abs(i_site - lattice_cutoff_)};
      // Just count the event; do not actually bind a motor to this site
      // printf("delta of %i += 1\n", delta);
      lattice_bind_stats_[delta].first++;
      properties_->microtubules.FlagForUpdate();
    };
    auto weight_bind_i = [&](ENTRY_T target) {
      return std::get<SITE_T *>(target)->weight_bind_;
    };
    auto poisson_bind_i = [&](double p, int n) {
      // Each delta distance has 2 sites available to it each timestep
      // Do not count delta = 0, where 'main' motor is permanently bound to
      for (int delta{1}; delta <= lattice_cutoff_; delta++) {
        lattice_bind_stats_[delta].second += 2;
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
    wally_->ErrorExit("FIX TEST\n");
    // Bind_ATP: bind ATP to motor heads that have released their ADP
    event_name = "bind_ATP";
    auto exe_bind_ATP = [&](ENTRY_T target) {
      Bind_ATP(std::get<POP_T *>(target));
    };
    // events_.emplace_back(event_name, p_bind_ATP_, &n_bound_NULL_,
    // &bound_NULL_,
    //                      exe_bind_ATP, binomial, set_ran_indices);
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
      lattice_step_bind_ii_stats_.first++;
      // If dock site is plus end, unbind motor and place it on minus end
      POP_T *docked_head{std::get<POP_T *>(target)};
      SITE_T *dock_site{docked_head->motor_->GetDockSite()};
      if (dock_site->index_ == dock_site->mt_->plus_end_) {
        Unbind_I(docked_head->motor_->GetActiveHead());
        SITE_T *new_site{&dock_site->mt_->lattice_[dock_site->mt_->minus_end_]};
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
      if (n_bind_ii_candidates_ > 0) {
        lattice_step_bind_ii_stats_.second++;
      }
      if (p > 0.0) {
        return gsl_->SamplePoissonDist(p);
      }
      return 0;
    };
    events_.emplace_back(event_name, p_avg_bind_ii_, &n_bind_ii_candidates_,
                         &bind_ii_candidates_, exe_bind_ii, weight_bind_ii,
                         poisson_bind_ii, get_ran_prob, set_ran_indices);
    // Unbind_II: Converts ADPP to ADP and unbinds a doubly-bound head
    event_name = "unbind_ii";
    auto exe_unbind_ii = [&](ENTRY_T target) {
      lattice_step_unbind_ii_stats_.first++;
      Unbind_II(std::get<POP_T *>(target));
    };
    auto weight_unbind_ii = [&](ENTRY_T target) {
      double weight{std::get<POP_T *>(target)->motor_->GetWeight_Unbind_II()};
      // printf("wt_unbind_ii is %g\n", weight);
      return weight;
    };
    auto poisson_unbind_ii = [&](double p, int n) {
      if (n_unbind_ii_candidates_ > 0) {
        lattice_step_unbind_ii_stats_.second++;
      }
      if (p > 0.0) {
        return gsl_->SamplePoissonDist(p);
      }
      return 0;
    };
    events_.emplace_back(event_name, p_avg_unbind_ii_, &n_unbind_ii_candidates_,
                         &unbind_ii_candidates_, exe_unbind_ii,
                         weight_unbind_ii, poisson_unbind_ii, get_ran_prob,
                         set_ran_indices);
    // Unbind_I: Converts ADPP to ADP and unbinds a singly-bound head
    event_name = "unbind_i";
    auto exe_unbind_i = [&](ENTRY_T target) {
      // Count stats for unbind_i but do not actually execute it
      lattice_step_unbind_i_stats_.first++;
    };
    auto weight_unbind_i = [&](ENTRY_T target) {
      return std::get<POP_T *>(target)->site_->weight_unbind_;
    };
    auto poisson_unbind_i = [&](double p, int n) {
      if (n_unbind_i_candidates_ > 0) {
        lattice_step_unbind_i_stats_.second++;
      }
      if (p > 0.0) {
        return gsl_->SamplePoissonDist(p);
      }
      return 0;
    };
    events_.emplace_back(event_name, p_avg_unbind_i_, &n_unbind_i_candidates_,
                         &unbind_i_candidates_, exe_unbind_i, weight_unbind_i,
                         poisson_unbind_i, get_ran_prob, set_ran_indices);
  }
}

void KinesinManagement::ReportProbabilities() {

  for (const auto &entry : p_theory_) {
    auto label = entry.first;
    auto value = entry.second;
    for (int i{0}; i < value.size(); i++) {
      for (int j{0}; j < value[i].size(); j++) {
        for (int k{0}; k < value[i][j].size(); k++) {
          std::string name{label};
          if (value.size() > 1) {
            name += "_" + std::to_string(i);
          }
          if (value[i].size() > 1) {
            name += "_" + std::to_string(j);
          }
          if (value[i][j].size() > 1) {
            name += "_" + std::to_string(k);
          }
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
          double n_exe_tot{(double)event->n_executed_tot_};
          unsigned long n_opp_tot{event->n_opportunities_tot_};
          if (n_opp_tot < 0.0) {
            wally_->Log("WHAT?? n_opp = %g\n", n_opp_tot);
            exit(1);
          }
          wally_->Log("For Kin event %s:\n", event->name_.c_str());
          wally_->Log("   p_theory = %g\n", value[i][j][k]);
          wally_->Log("   p_actual = %g", n_exe_tot / n_opp_tot);
          wally_->Log(" (n_exe = %lu)", event->n_executed_tot_);
          wally_->Log(" (n_opp = %lu)\n", event->n_opportunities_tot_);
        }
      }
    }
  }
  if (wally_->test_mode_ != nullptr) {
    if (strcmp(wally_->test_mode_, "motor_lattice_bind") == 0) {
      printf("For motor_lattice_bind:\n");
      for (int delta{0}; delta <= lattice_cutoff_; delta++) {
        int n_exe{lattice_bind_stats_[delta].first};
        int n_opp{lattice_bind_stats_[delta].second};
        double p{double(n_exe) / n_opp};
        printf("p_theory_[%i] = %g\n", delta,
               p_avg_bind_i_ * weight_lattice_bind_[delta]);
        printf("p_actual_[%i] = %g (n_exe = %i)\n", delta, p,
               lattice_bind_stats_[delta].first);
      }
    }
    if (strcmp(wally_->test_mode_, "motor_lattice_step") == 0) {
      wally_->Log("For motor_lattice_step (delta = %i):\n", test_delta_);
      // Bind_ii
      double p_bind_ii{p_avg_bind_ii_};
      if (test_delta_ >= 0) {
        p_bind_ii *= weight_lattice_bind_[test_delta_];
      }
      wally_->Log("p_theory_bind_ii = %g\n", p_bind_ii);
      wally_->Log("p_actual_bind_ii = %g (n_exe = %i)\n",
                  double(lattice_step_bind_ii_stats_.first) /
                      lattice_step_bind_ii_stats_.second,
                  lattice_step_bind_ii_stats_.first);
      // Unbind_ii
      double p_unbind_ii{p_avg_unbind_ii_};
      if (test_delta_ >= 0) {
        p_unbind_ii *= weight_lattice_unbind_[test_delta_];
      }
      wally_->Log("p_theory_unbind_ii = %g\n", p_unbind_ii);
      wally_->Log("p_actual_unbind_ii = %g (n_exe = %i)\n",
                  double(lattice_step_unbind_ii_stats_.first) /
                      lattice_step_unbind_ii_stats_.second,
                  lattice_step_unbind_ii_stats_.first);
      // Unbind_i
      double p_unbind_i{p_avg_unbind_i_};
      if (test_delta_ >= 0) {
        p_unbind_i *= weight_lattice_unbind_[test_delta_];
      }
      wally_->Log("p_theory_unbind_i = %g\n", p_unbind_i);
      wally_->Log("p_actual_unbind_i = %g (n_exe = %i)\n",
                  double(lattice_step_unbind_i_stats_.first) /
                      lattice_step_unbind_i_stats_.second,
                  lattice_step_unbind_i_stats_.first);
    }
  }
}

Kinesin *KinesinManagement::GetFreeMotor() {

  // Randomly pick a motor from the reservoir
  int i_motor = properties_->gsl.GetRanInt(n_motors_);
  Kinesin *motor = &motors_[i_motor];
  int attempts = 0;
  while (motor->heads_active_ > 0 || motor->tethered_ == true) {
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

void KinesinManagement::ReportExecutionOf(std::string event_name) {

  if (verbosity_ < 1) {
    return;
  }
  wally_->Log("Executing ");
  std::cout << event_name << std::endl;
}

void KinesinManagement::ReportFailureOf(std::string event_name) {

  wally_->Log("yo we failed to execute ");
  std::cout << event_name << std::endl;
}

void KinesinManagement::Update_Weights() {

  /*
  if (!dynamic_weights_) {
    return;
  }
  MicrotubuleManagement *mts = &properties_->microtubules;
  // 'Reset' all tubulin weights by setting them to neighbor weights only
  for (int i_mt{0}; i_mt < mts->mt_list_.size(); i_mt++) {
    for (int i_site{0}; i_site < mts->mt_list_[i_mt].n_sites_; i_site++) {
      Tubulin *site{&mts->mt_list_[i_mt].lattice_[i_site]};
      int n_neighbs{site->GetKif4ANeighborCount()};
      site->weight_bind_ = weight_neighbs_bind_[n_neighbs];
      site->weight_unbind_ = weight_neighbs_unbind_[n_neighbs];
    }
  }
  */
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

void KinesinManagement::Update_Extensions() {

  if (!tethering_active_) {
    return;
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    active_[i_entry]->UpdateExtension();
  }
}

void KinesinManagement::Update_Docked() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting K_MGMT::Update_Docked()");
  }
  n_bind_ii_candidates_ = 0;
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    Kinesin *motor = active_[i_entry];
    if (motor->heads_active_ != 1) {
      continue;
    }
    POP_T *head{motor->GetActiveHead()};
    if (head->ligand_ == "ADPP") {
      Tubulin *dock_site{motor->GetDockSite()};
      if (!dock_site->occupied_) {
        bind_ii_candidates_[n_bind_ii_candidates_++] = motor->GetDockedHead();
      }
    }
  }
  if (verbosity_ >= 1) {
    wally_->Log(" n_bind_ii_candidates = %i\n", n_bind_ii_candidates_);
    for (int i_entry{0}; i_entry < n_bind_ii_candidates_; i_entry++) {
      POP_T *head{std::get<POP_T *>(bind_ii_candidates_[i_entry])};
      wally_->Log("   - motor %i ", head->motor_->id_);
      wally_->Log("(weight = %g)\n", head->motor_->GetWeight_Bind_II());
    }
    wally_->Log(" -- finished\n");
  }
}

void KinesinManagement::Update_Bound_NULL() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting K_MGMT::Update_Bound_NULL()");
  }
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    n_bound_NULL_[n_neighbs] = 0;
  }
  // n_bound_NULL_ = 0;
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    Kinesin *motor = active_[i_entry];
    if (motor->heads_active_ != 1) {
      continue;
    }
    if (motor->tethered_ and !motor->HasSatellite()) {
      continue;
    }
    POP_T *head{motor->GetActiveHead()};
    if (head->ligand_ == "NULL") {
      int n_neighbs{head->GetKif4ANeighborCount()};
      bound_NULL_[n_neighbs][n_bound_NULL_[n_neighbs]++] = head;
    }
    /*
    if (motor->head_one_.ligand_ == "NULL") {
      int n_neighbs{motor->head_one_.GetKif4ANeighborCount_Step()};
      int index{n_bound_NULL_[n_neighbs]++};
      bound_NULL_[n_neighbs][index] = &motor->head_one_;
    }
    if (motor->head_two_.ligand_ == "NULL") {
      int n_neighbs{motor->head_two_.GetKif4ANeighborCount_Step()};
      int index{n_bound_NULL_[n_neighbs]++};
      bound_NULL_[n_neighbs][index] = &motor->head_two_;
    }
    */
  }
  if (verbosity_ > 1) {
    wally_->Log(" -- finished\n");
  }
}

void KinesinManagement::Update_Bound_NULL_Teth() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting K_MGMT::Update_Bound_NULL_Teth()");
  }
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      n_bound_NULL_to_teth_[n_neighbs][x_dub] = 0;
      n_bound_NULL_fr_teth_[n_neighbs][x_dub] = 0;
    }
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    Kinesin *motor = active_[i_entry];
    if (motor->tethered_ and !motor->HasSatellite()) {
      int x_dub{motor->x_dist_doubled_};
      if (motor->head_one_.ligand_ == "NULL") {
        int n_neighbs{motor->head_one_.GetKif4ANeighborCount_Step()};
        if (motor->GetDirectionTowardRest() == motor->mt_->delta_x_) {
          int index{n_bound_NULL_to_teth_[n_neighbs][x_dub]++};
          bound_NULL_to_teth_[n_neighbs][x_dub][index] = &motor->head_one_;
        } else {
          int index{n_bound_NULL_fr_teth_[n_neighbs][x_dub]++};
          bound_NULL_fr_teth_[n_neighbs][x_dub][index] = &motor->head_one_;
        }
      }
      if (motor->head_two_.ligand_ == "NULL") {
        int n_neighbs{motor->head_two_.GetKif4ANeighborCount_Step()};
        if (motor->GetDirectionTowardRest() == motor->mt_->delta_x_) {
          int index{n_bound_NULL_to_teth_[n_neighbs][x_dub]++};
          bound_NULL_to_teth_[n_neighbs][x_dub][index] = &motor->head_two_;
        } else {
          int index{n_bound_NULL_fr_teth_[n_neighbs][x_dub]++};
          bound_NULL_fr_teth_[n_neighbs][x_dub][index] = &motor->head_two_;
        }
      }
    }
  }
  if (verbosity_ > 1) {
    wally_->Log(" -- finished\n");
  }
}

void KinesinManagement::Update_Bound_ATP() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting K_MGMT::Update_Bound_ATP() ");
  }
  n_bound_ATP_ = 0;
  for (int i_motor = 0; i_motor < n_active_; i_motor++) {
    Kinesin *motor = active_[i_motor];
    // if (!motor->IsStalled()) {
    if (motor->head_one_.ligand_ == "ATP") {
      bound_ATP_[n_bound_ATP_++] = &motor->head_one_;
    }
    if (motor->head_two_.ligand_ == "ATP") {
      bound_ATP_[n_bound_ATP_++] = &motor->head_two_;
    }
    // }
  }
  if (verbosity_ > 1) {
    wally_->Log("-- finished\n");
  }
}

void KinesinManagement::Update_Bound_ATP_Stalled() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting K_MGMT::Update_Bound_ATP_Stalled() ");
  }
  n_bound_ATP_st_ = 0;
  for (int i_motor = 0; i_motor < n_active_; i_motor++) {
    Kinesin *motor = active_[i_motor];
    if (motor->IsStalled()) {
      if (motor->head_one_.ligand_ == "ATP") {
        bound_ATP_st_[n_bound_ATP_st_++] = &motor->head_one_;
      }
      if (motor->head_two_.ligand_ == "ATP") {
        bound_ATP_st_[n_bound_ATP_st_++] = &motor->head_two_;
      }
    }
  }
  if (verbosity_ > 1) {
    wally_->Log("-- finished\n");
  }
}

void KinesinManagement::Update_Bound_ADPP_I() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting K_MGMT::Update_Bound_ADPP_I()");
  }
  /*
  for (int i_aff{0}; i_aff < n_affinities_tot_; i_aff++) {
    for (int i_neighb{0}; i_neighb <= max_neighbs_; i_neighb++) {
      n_bound_ADPP_i_[i_aff][i_neighb] = 0;
    }
  }
  */
  n_unbind_i_candidates_ = 0;
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    Kinesin *motor = active_[i_entry];
    if (motor->heads_active_ != 1) { // or motor->IsStalled()) {
      continue;
    }
    // if (!motor->tethered_ or motor->HasSatellite()) {
    POP_T *head = motor->GetActiveHead();
    if (head->ligand_ == "ADPP") {
      unbind_i_candidates_[n_unbind_i_candidates_++] = head;
      /*
      int aff = head->GetAffinity();
      int neighbs = head->GetKif4ANeighborCount();
      int index = n_bound_ADPP_i_[aff][neighbs]++;
      bound_ADPP_i_[aff][neighbs][index] = head;
      */
    }
    // }
  }
  if (verbosity_ > 1) {
    wally_->Log(" -- finished\n");
  }
}

/*
void KinesinManagement::Update_Bound_ADPP_I_Stalled() {

  if (verbosity_ >= 1) {
    wally_->Log("starting UPDATE ADPP_I_ST\n");
  }
  for (int i_aff{0}; i_aff < n_affinities_tot_; i_aff++) {
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      n_bound_ADPP_i_st_[i_aff][n_neighbs] = 0;
    }
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    Kinesin *motor = active_[i_entry];
    if (motor->heads_active_ != 1 or !motor->IsStalled()) {
      continue;
    }
    if (!motor->tethered_ or motor->HasSatellite()) {
      POP_T *active_head = motor->GetActiveHead();
      if (active_head->ligand_ == "ADPP") {
        int aff = active_head->GetAffinity();
        int neighbs = active_head->GetKif4ANeighborCount();
        int index = n_bound_ADPP_i_st_[aff][neighbs]++;
        bound_ADPP_i_st_[aff][neighbs][index] = active_head;
      }
    }
  }
}
*/
/*
void KinesinManagement::Update_Bound_ADPP_I_Teth() {

  for (int i_aff{0}; i_aff < n_affinities_tot_; i_aff++) {
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
        n_bound_ADPP_i_teth_[i_aff][n_neighbs][x_dub] = 0;
      }
    }
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    Kinesin *motor = active_[i_entry];
    if (motor->heads_active_ != 1 or motor->IsStalled()) {
      continue;
    }
    if (motor->tethered_ and !motor->HasSatellite()) {
      POP_T *active_head{motor->GetActiveHead()};
      if (active_head->ligand_ == "ADPP") {
        int aff{active_head->GetAffinity()};
        int n_neighbs{active_head->GetKif4ANeighborCount()};
        int x_dub{motor->x_dist_doubled_};
        int index{n_bound_ADPP_i_teth_[aff][n_neighbs][x_dub]++};
        bound_ADPP_i_teth_[aff][n_neighbs][x_dub][index] = active_head;
      }
    }
  }
}
*/
/*
void KinesinManagement::Update_Bound_ADPP_I_Teth_Stalled() {

  for (int i_aff{0}; i_aff < n_affinities_tot_; i_aff++) {
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
        n_bound_ADPP_i_teth_st_[i_aff][n_neighbs][x_dub] = 0;
      }
    }
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    Kinesin *motor = active_[i_entry];
    if (motor->heads_active_ != 1 or !motor->IsStalled()) {
      continue;
    }
    if (motor->tethered_ and !motor->HasSatellite()) {
      POP_T *active_head{motor->GetActiveHead()};
      if (active_head->ligand_ == "ADPP") {
        int aff{active_head->GetAffinity()};
        int n_neighbs{active_head->GetKif4ANeighborCount()};
        int x_dub{motor->x_dist_doubled_};
        int index{n_bound_ADPP_i_teth_st_[aff][n_neighbs][x_dub]++};
        bound_ADPP_i_teth_st_[aff][n_neighbs][x_dub][index] = active_head;
      }
    }
  }
}
*/
void KinesinManagement::Update_Bound_ADPP_II() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting K_MGMT::Update_Bound_ADPP_II()");
  }
  n_unbind_ii_candidates_ = 0;
  for (int i_motor = 0; i_motor < n_active_; i_motor++) {
    Kinesin *motor = active_[i_motor];
    if (motor->heads_active_ != 2) {
      continue;
    }
    POP_T *chosen_head{nullptr};
    // If both heads are ADPP-bound, pick rear head to unbind
    // FIXME when attempting to incorporate back-stepping cycle
    if (motor->head_one_.ligand_ == "ADPP" and
        motor->head_two_.ligand_ == "ADPP") {
      wally_->ErrorExit("K_MGMT::Update_Bound_ADPP_II()");
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
      unbind_ii_candidates_[n_unbind_ii_candidates_++] = chosen_head;
      /*
      int aff{chosen_head->GetAffinity()};
      int neighbs{chosen_head->GetKif4ANeighborCount()};
      int i_entry{n_bound_ADPP_ii_[aff][neighbs]++};
      bound_ADPP_ii_[aff][neighbs][i_entry] = chosen_head;
      */
    } else {
      wally_->ErrorExit("KinMGMT::Update_Bound_ADPP_II");
    }
  }
  // if (n_unbind_ii_candidates_ > 0) {
  //   printf("%i unbind_ii candidates\n", n_unbind_ii_candidates_);
  // }
  if (verbosity_ > 1) {
    wally_->Log(" -- finished\n");
  }
}

void KinesinManagement::Update_Bound_Unteth() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting K_MGMT::Update_Bound_Unteth()");
  }
  n_bound_unteth_ = 0;
  n_tether_bound_candidates_ = 0;
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    Kinesin *motor = active_[i_entry];
    if (motor->tethered_) {
      continue;
    }
    if (motor->heads_active_ == 1) {
      bound_unteth_[n_bound_unteth_++] = motor->GetActiveHead();
      tether_bound_candidates_[n_tether_bound_candidates_++] =
          motor->GetActiveHead();
    } else {
      bound_unteth_[n_bound_unteth_++] = &motor->head_one_;
      tether_bound_candidates_[n_tether_bound_candidates_++] =
          &motor->head_one_;
    }
  }
  if (verbosity_ > 1) {
    wally_->Log(" -- finished\n");
  }
}

void KinesinManagement::Update_Bound_Teth() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting K_MGMT::Update_Bound_Teth()");
  }
  for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
    n_bound_teth_[x_dub] = 0;
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    Kinesin *motor = active_[i_entry];
    if (motor->heads_active_ == 0) {
      continue;
    }
    if (motor->tethered_ and !motor->HasSatellite()) {
      int x_dub{motor->x_dist_doubled_};
      if (motor->heads_active_ == 2) {
        bound_teth_[x_dub][n_bound_teth_[x_dub]++] = &motor->head_one_;
      } else {
        bound_teth_[x_dub][n_bound_teth_[x_dub]++] = motor->GetActiveHead();
      }
    }
  }
  if (verbosity_ > 1) {
    wally_->Log(" -- finished\n");
  }
}

void KinesinManagement::Update_Free_Teth() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting K_MGMT::Update_Bound_Unteth()");
  }
  n_satellites_ = 0;
  n_bind_i_teth_candidates_ = 0;
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    Kinesin *motor = active_[i_entry];
    if (motor->heads_active_ == 0 and motor->tethered_) {
      satellites_[n_satellites_++] = &motor->head_one_;
      bind_i_teth_candidates_[n_bind_i_teth_candidates_++] = &motor->head_one_;
    }
  }
  if (verbosity_ > 1) {
    wally_->Log(" -- finished\n");
  }
}

void KinesinManagement::RunKMC() {

  // if (properties_->current_step_ == 1384055) {
  //   verbosity_ = 3;
  // }
  if (properties_->current_step_ * parameters_->delta_t >= t_active_) {
    population_active_ = true;
  }
  if (!population_active_) {
    return;
  }
  UpdateLists();
  SampleEventStatistics();
  GenerateExecutionSequence();
  ExecuteEvents();
}

void KinesinManagement::UpdateLists() {

  if (lists_up_to_date_) {
    return;
  }
  properties_->microtubules.UpdateUnoccupied();
  Update_Docked();
  Update_Bound_NULL();
  Update_Bound_ATP();
  Update_Bound_ADPP_I();
  Update_Bound_ADPP_II();
  if (parameters_->motors.tethers_active) {
    Update_Extensions();
    Update_Free_Teth();
    Update_Bound_NULL_Teth();
    // Update_Bound_ADPP_I_Teth();
    // Update_Bound_ADPP_I_Teth_Stalled();
    Update_Bound_Teth();
    Update_Bound_Unteth();
    properties_->prc1.Update_Bound_Unteth();
  }
}

void KinesinManagement::SampleEventStatistics() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting K_MGMT::SampleEventStatistics()\n");
  }
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
    if (verbosity_ >= 2 and event.n_expected_ > 0) {
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
  if (verbosity_ >= 2) {
    wally_->Log(" %i events expected pre-correction\n", n_events_to_exe_);
  }
  if (n_events_to_exe_ <= 1) {
    return;
  }
  // Put all events with >1 target into an active_ array
  std::pair<EVENT_T *, ENTRY_T> active_events[n_events_to_exe_];
  int i_active{0};
  for (auto &&event : events_) {
    // Add a ptr to the event for each target it has
    for (int i_tar{0}; i_tar < event.n_expected_; i_tar++) {
      active_events[i_active++] = std::make_pair(&event, event.targets_[i_tar]);
      if (verbosity_ >= 2) {
        wally_->Log(" Added %s to active_events\n", event.name_.c_str());
      }
    }
  }
  // Scan through all active events to ensure that no two target the same
  // motor
  for (int i_entry{0}; i_entry < n_events_to_exe_; i_entry++) {
    EVENT_T *event_i = active_events[i_entry].first;
    Kinesin *motor_i{nullptr};
    try {
      motor_i = std::get<POP_T *>(active_events[i_entry].second)->motor_;
    } catch (...) {
      continue;
    }
    for (int j_entry{i_entry + 1}; j_entry < n_events_to_exe_; j_entry++) {
      EVENT_T *event_j = active_events[j_entry].first;
      Kinesin *motor_j{nullptr};
      try {
        motor_j = std::get<POP_T *>(active_events[j_entry].second)->motor_;
      } catch (...) {
        continue;
      }
      // If event_i and event_j target different motors, continue
      if (motor_i != motor_j) {
        continue;
      }
      double p_one{event_i->p_occur_};
      double p_two{event_j->p_occur_};
      double ran{gsl_->GetRanProb()};
      if (ran < p_one / (p_one + p_two)) {
        if (verbosity_ >= 2) {
          wally_->Log("Removed xlink #%i from %s's targets\n", motor_i->id_,
                      event_i->name_.c_str());
        }
        event_i->RemoveTarget(active_events[i_entry].second);
        active_events[i_entry] = active_events[n_events_to_exe_ - 1];
        i_entry--;
        n_events_to_exe_--;
        break;
      } else {
        if (verbosity_ >= 2) {
          wally_->Log("Removed xlink #%i from %s's targets\n", motor_j->id_,
                      event_j->name_.c_str());
        }
        event_j->RemoveTarget(active_events[j_entry].second);
        active_events[j_entry] = active_events[n_events_to_exe_ - 1];
        j_entry--;
        n_events_to_exe_--;
      }
    }
  }
}

void KinesinManagement::GenerateExecutionSequence() {

  if (n_events_to_exe_ == 0) {
    return;
  }
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
}

void KinesinManagement::ExecuteEvents() {

  for (int i_event{0}; i_event < n_events_to_exe_; i_event++) {
    if (verbosity_ > 1) {
      wally_->Log("Executing %s\n", events_to_exe_[i_event]->name_.c_str());
    }
    events_to_exe_[i_event]->Execute();
    if (events_to_exe_[i_event]->name_ != "bind_ATP_0") {
      lists_up_to_date_ = false;
    }
    if (verbosity_ > 1) {
      wally_->Log(" -- Execution successful -- \n");
    }
  }
}

void KinesinManagement::Bind_I(SITE_T *site) {

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

void KinesinManagement::Bind_I_Teth(POP_T *satellite_head) {

  Kinesin *motor{satellite_head->motor_};
  Tubulin *site{motor->GetWeightedSite_Bind_I_Teth()};
  if (site == nullptr) {
    wally_->Log("failed to Bind_I_Free_Teth (MOTOR)\n");
    return;
  }
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
  properties_->microtubules.FlagForUpdate();
}

void KinesinManagement::Bind_ATP(POP_T *head) {

  // Update motor head
  head->ligand_ = "ATP";
  // If head is leading, change conformation
  if (!head->trailing_) {
    head->motor_->ChangeConformation();
  } else {
    // Not an error; perfectly valid @ plus end due to end pausing
    // wally_->ErrorExit("K_MGMT::Bind_ATP()");
  }
}

void KinesinManagement::Hydrolyze(POP_T *head) { head->ligand_ = "ADPP"; }

void KinesinManagement::Bind_II(POP_T *head) {

  // Verify that proposed site is unoccupied
  Kinesin *motor{head->motor_};
  Tubulin *dock_site{motor->GetDockSite()};
  if (dock_site->occupied_) {
    printf("woah.\n");
    return;
  }
  // Update site
  dock_site->motor_head_ = head;
  dock_site->occupied_ = true;
  // Update motor
  head->site_ = dock_site;
  head->ligand_ = "NULL";
  motor->heads_active_++;
  properties_->microtubules.FlagForUpdate();
}

void KinesinManagement::Unbind_II(POP_T *head) {

  if (head->motor_->heads_active_ != 2) {
    wally_->ErrorExit("Kin_MGMT::Unbind_II()");
  }
  // Update site
  head->site_->occupied_ = false;
  head->site_->motor_head_ = nullptr;
  if (verbosity_ >= 2) {
    printf("Unbound site #%i\n", head->site_->index_);
  }
  // Update motor
  head->site_ = nullptr;
  head->ligand_ = "ADP";
  head->motor_->heads_active_--;
  head->RelieveFrustration();
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
  motor->UntetherSatellite();
  if (!motor->tethered_) {
    RemoveFromActive(motor);
  }
  properties_->microtubules.FlagForUpdate();
}

void KinesinManagement::Tether_Free(ALT_T *untethered_head) {

  Kinesin *motor{GetFreeMotor()};
  AssociatedProtein *xlink{untethered_head->xlink_};
  motor->xlink_ = xlink;
  motor->tethered_ = true;
  xlink->tethered_ = true;
  xlink->motor_ = motor;
  AddToActive(motor);
  properties_->prc1.FlagForUpdate();
}

void KinesinManagement::Tether_Bound(POP_T *bound_head) {

  Kinesin *motor{bound_head->motor_};
  AssociatedProtein *xlink{motor->GetWeightedXlink_Tether_Bound()};
  // Update motor and xlink details
  motor->xlink_ = xlink;
  motor->tethered_ = true;
  xlink->tethered_ = true;
  xlink->motor_ = motor;
  if (!motor->tethered_) {
    wally_->ErrorExit("Kin_MGMT::Tether_Bound()");
  }
  properties_->prc1.FlagForUpdate();
}

void KinesinManagement::Untether(POP_T *head) {

  Kinesin *motor{head->motor_};
  AssociatedProtein *xlink{motor->xlink_};
  // Update motor and xlink details
  xlink->motor_ = nullptr;
  xlink->tethered_ = false;
  motor->xlink_ = nullptr;
  motor->tethered_ = false;
  if (motor->heads_active_ == 0) {
    RemoveFromActive(motor);
  }
  properties_->prc1.FlagForUpdate();
}
