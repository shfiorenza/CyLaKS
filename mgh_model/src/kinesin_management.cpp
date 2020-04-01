#include "kinesin_management.h"
#include "master_header.h"

KinesinManagement::KinesinManagement() {}

void KinesinManagement::Initialize(system_parameters *parameters,
                                   system_properties *properties) {

  // verbosity_ = 1;
  wally_ = &properties->wallace;
  parameters_ = parameters;
  properties_ = properties;
  CalculateCutoffs();
  SetParameters();
  GenerateMotors();
  InitializeLists();
  if (wally_->test_mode_ == nullptr) {
    InitializeEvents();
  } else {
    t_active_ = std::numeric_limits<double>::max();
  }
}

void KinesinManagement::CalculateCutoffs() {

  double kbT{parameters_->kbT};
  double r_0{parameters_->motors.r_0};
  double r_y{parameters_->microtubules.y_dist / 2};
  double k_spring{parameters_->motors.k_spring};
  double k_slack{parameters_->motors.k_slack};
  double site_size{parameters_->microtubules.site_size};
  // First, calculate rest_dist_ in number of sites
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
  if (parameters_->xlinks.c_bulk > 0.0 and parameters_->motors.c_bulk > 0.0 and
      parameters_->motors.tethers_active) {
    tethering_active_ = true;
    wally_->Log("\nFor motors:\n");
    wally_->Log("  rest_dist is %g\n", rest_dist_);
    wally_->Log("  comp_cutoff is %i\n", comp_cutoff_);
    wally_->Log("  dist_cutoff is %i\n", teth_cutoff_);
  } else {
    wally_->Log("\nTethering is disabled for motors.\n");
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
    double gauss{amp * exp(-1 * dist * dist / (2 * sigma * sigma))};
    wt_lattice_bind[i_aff] = 1.0 + gauss;
    wt_lattice_unbind[i_aff] = 1.0 / (1.0 + gauss);
  }
  wt_lattice_bind[n_affinities_ - 1] = 1.0;
  wt_lattice_unbind[n_affinities_ - 1] = 1.0;
  // Array of interaction energies (in kbT) due to neighbor cooperativity
  double int_energy{-1 * parameters_->motors.interaction_energy}; // kbT
  std::vector<double> neighb_energy(max_neighbs_ + 1, 0.0);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    neighb_energy[n_neighbs] = n_neighbs * int_energy; // ! in kbT
  }
  // Array of tether energies (in kbT) for any given extension
  double site_size{parameters_->microtubules.site_size};
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
    if (dr < 0) {
      k_spring = k_slack;
    } else {
      k_spring = k_teth;
    }
    tether_energy[x_dub] = 0.5 * k_spring * dr * dr / kbT; // ! in kbT
    tether_force[x_dub] = fabs(k_spring * dr);
    tether_cosine[x_dub] = r_x / r;
  }
  double delta_t{parameters_->delta_t};
  double k_on{parameters_->motors.k_on};
  double c_motor{parameters_->motors.c_bulk};
  double k_hydrolyze{parameters_->motors.k_hydrolyze};
  double k_hydrolyze_stalled{parameters_->motors.k_hydrolyze_stalled};
  double k_tether{parameters_->motors.k_tether};
  double k_untether{parameters_->motors.k_untether};
  double c_eff_teth{parameters_->motors.c_eff_tether};
  p_hydrolyze_ = k_hydrolyze * delta_t;
  p_hydrolyze_st_ = k_hydrolyze_stalled * delta_t;
  p_tether_free_ = k_tether * c_motor * delta_t;
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
  double f_stall{parameters_->motors.stall_force};
  double k_on_ATP{parameters_->motors.k_on_ATP};
  double c_ATP{parameters_->motors.c_ATP};
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
  double c_eff{parameters_->motors.c_eff_bind};
  double k_off_ii{parameters_->motors.k_off_ii};
  double k_off_i{parameters_->motors.k_off_i};
  double k_off_i_st{parameters_->motors.k_off_i_stalled};
  double p_bind_i{k_on * c_motor * delta_t};
  double p_bind_ii{k_on * c_eff * delta_t};
  double p_unbind_ii{k_off_ii * delta_t};
  double p_unbind_i{k_off_i * delta_t};
  double p_unbind_i_st{k_off_i_st * delta_t};
  double sigma_off{1.5};
  // Force perpetually applied to motors, e.g., by optical trapping
  if (parameters_->motors.applied_force > 0.0) {
    printf("FIX APPLIED MOT FORCE\n");
    exit(1);
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
  p_unbind_i_st_.resize(n_affinities_);
  p_unbind_i_teth_.resize(n_affinities_);
  p_unbind_i_teth_st_.resize(n_affinities_);
  weight_bind_i_teth_.resize(n_affinities_);
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    p_bind_i_[i_aff].resize(max_neighbs_ + 1);
    p_bind_ii_[i_aff].resize(max_neighbs_ + 1);
    p_unbind_ii_[i_aff].resize(max_neighbs_ + 1);
    p_unbind_i_[i_aff].resize(max_neighbs_ + 1);
    p_unbind_i_st_[i_aff].resize(max_neighbs_ + 1);
    p_unbind_i_teth_[i_aff].resize(max_neighbs_ + 1);
    p_unbind_i_teth_st_[i_aff].resize(max_neighbs_ + 1);
    weight_bind_i_teth_[i_aff].resize(max_neighbs_ + 1);
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
      p_bind_ii_[i_aff][n_neighbs] = p_bind_ii * wt_tot_bind;
      p_unbind_ii_[i_aff][n_neighbs] = p_unbind_ii * wt_tot_unbind;
      p_unbind_i_[i_aff][n_neighbs] = p_unbind_i * wt_tot_unbind;
      p_unbind_i_st_[i_aff][n_neighbs] = p_unbind_i_st * wt_tot_unbind;
      p_unbind_i_teth_[i_aff][n_neighbs].resize(2 * teth_cutoff_ + 1);
      p_unbind_i_teth_st_[i_aff][n_neighbs].resize(2 * teth_cutoff_ + 1);
      weight_bind_i_teth_[i_aff][n_neighbs].resize(2 * teth_cutoff_ + 1);
      for (int x_dub{0}; x_dub <= 2 * teth_cutoff_; x_dub++) {
        double f_dependence{tether_force[x_dub] * sigma_off / kbT};
        double wt_teth{exp(-(1.0 - lambda_teth) * f_dependence)};
        double wt_unteth{exp(lambda_teth * f_dependence)};
        if (x_dub < 2 * comp_cutoff_) {
          wt_teth = 0.0;
          wt_unteth = 0.0;
        }
        weight_bind_i_teth_[i_aff][n_neighbs][x_dub] = wt_tot_bind * wt_teth;
        p_unbind_i_teth_[i_aff][n_neighbs][x_dub] =
            p_unbind_i_[i_aff][n_neighbs] * wt_unteth;
        p_unbind_i_teth_st_[i_aff][n_neighbs][x_dub] =
            p_unbind_i_st_[i_aff][n_neighbs] * wt_unteth;
      }
    }
  }
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
  p_theory_["bind_i"] = {p_bind_i_};
  p_theory_["bind_ii"] = {p_bind_ii_};
  p_theory_["unbind_ii"] = {p_unbind_ii_};
  p_theory_["unbind_i"] = {p_unbind_i_};
  p_theory_["unbind_i_st"] = {p_unbind_i_st_};
  if (tethering_active_) {
    p_theory_["unbind_i_teth"] = p_unbind_i_teth_;
    p_theory_["unbind_i_teth_st"] = p_unbind_i_teth_st_;
  }
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
          if (verbosity_ >= 3) {
            wally_->Log("%s = %g\n", name.c_str(), value[i][j][k]);
          }
          if (value[i][j][k] > 0.5) {
            wally_->Log("Error! %s = %g\n", name.c_str(), value[i][j][k]);
            // wally_->ErrorExit("Kin_MGMT:SetParameters()");
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
  bound_ATP_.resize(n_motors_);
  bound_ATP_st_.resize(n_motors_);
  bound_unteth_.resize(n_motors_);
  satellites_.resize(n_motors_);
  bind_i_teth_candidates_.resize(n_motors_);
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
  // 3-D stuff -- indexed by tubulin_affinity, n_neighbs, & i_motor
  n_docked_.resize(n_affinities_);
  n_bound_ADPP_ii_.resize(n_affinities_);
  n_bound_ADPP_i_.resize(n_affinities_);
  n_bound_ADPP_i_st_.resize(n_affinities_);
  docked_.resize(n_affinities_);
  bound_ADPP_ii_.resize(n_affinities_);
  bound_ADPP_i_.resize(n_affinities_);
  bound_ADPP_i_st_.resize(n_affinities_);
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    n_docked_[i_aff].resize(max_neighbs_ + 1);
    n_bound_ADPP_ii_[i_aff].resize(max_neighbs_ + 1);
    n_bound_ADPP_i_[i_aff].resize(max_neighbs_ + 1);
    n_bound_ADPP_i_st_[i_aff].resize(max_neighbs_ + 1);
    docked_[i_aff].resize(max_neighbs_ + 1);
    bound_ADPP_ii_[i_aff].resize(max_neighbs_ + 1);
    bound_ADPP_i_[i_aff].resize(max_neighbs_ + 1);
    bound_ADPP_i_st_[i_aff].resize(max_neighbs_ + 1);
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      n_docked_[i_aff][n_neighbs] = 0;
      n_bound_ADPP_ii_[i_aff][n_neighbs] = 0;
      n_bound_ADPP_i_[i_aff][n_neighbs] = 0;
      n_bound_ADPP_i_st_[i_aff][n_neighbs] = 0;
      docked_[i_aff][n_neighbs].resize(n_motors_);
      bound_ADPP_ii_[i_aff][n_neighbs].resize(n_motors_);
      bound_ADPP_i_[i_aff][n_neighbs].resize(n_motors_);
      bound_ADPP_i_st_[i_aff][n_neighbs].resize(n_motors_);
    }
  }
  // 4-D stuff -- indexed by tubulin_affinity, n_neighbs, x_dub, & i_motor
  n_bound_ADPP_i_teth_.resize(n_affinities_);
  n_bound_ADPP_i_teth_st_.resize(n_affinities_);
  bound_ADPP_i_teth_.resize(n_affinities_);
  bound_ADPP_i_teth_st_.resize(n_affinities_);
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    n_bound_ADPP_i_teth_[i_aff].resize(max_neighbs_ + 1);
    n_bound_ADPP_i_teth_st_[i_aff].resize(max_neighbs_ + 1);
    bound_ADPP_i_teth_[i_aff].resize(max_neighbs_ + 1);
    bound_ADPP_i_teth_st_[i_aff].resize(max_neighbs_ + 1);
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      n_bound_ADPP_i_teth_[i_aff][n_neighbs].resize(2 * teth_cutoff_ + 1);
      n_bound_ADPP_i_teth_st_[i_aff][n_neighbs].resize(2 * teth_cutoff_ + 1);
      bound_ADPP_i_teth_[i_aff][n_neighbs].resize(2 * teth_cutoff_ + 1);
      bound_ADPP_i_teth_st_[i_aff][n_neighbs].resize(2 * teth_cutoff_ + 1);
      for (int x_dub{0}; x_dub <= 2 * teth_cutoff_; x_dub++) {
        n_bound_ADPP_i_teth_[i_aff][n_neighbs][x_dub] = 0;
        n_bound_ADPP_i_teth_st_[i_aff][n_neighbs][x_dub] = 0;
        bound_ADPP_i_teth_[i_aff][n_neighbs][x_dub].resize(n_motors_);
        bound_ADPP_i_teth_st_[i_aff][n_neighbs][x_dub].resize(n_motors_);
      }
    }
  }
}

void KinesinManagement::InitializeEvents() {

  /* * Function that sets n random indices from the range [0, m) * */
  auto set_ran_indices = [&](int *indices, int n, int m) {
    if (m > 1) {
      properties_->gsl.SetRanIndices(indices, n, m);
    } else {
      indices[0] = 0;
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
  /* * Event entries * */
  std::string event_name; // scratch space to construct each event name
  // Bind_I: bind first motor head to MT and release ADP
  auto exe_bind_i = [&](ENTRY_T target) {
    SITE_T *site = std::get<SITE_T *>(target);
    if (site != nullptr) {
      ReportExecutionOf("bind_i");
      Bind_I(site);
    } else {
      ReportFailureOf("bind_i");
    }
  };
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    std::string AFFINITY{"_" + std::to_string(i_aff)};
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
      event_name = "bind_i";
      if (n_affinities_ > 1) {
        event_name += AFFINITY;
      }
      if (max_neighbs_ > 0) {
        event_name += N_NEIGHBS;
      }
      events_.emplace_back(event_name, p_bind_i_[i_aff][n_neighbs],
                           &properties_->microtubules.n_unocc_motor_[n_neighbs],
                           &properties_->microtubules.unocc_motor_[n_neighbs],
                           exe_bind_i, binomial, set_ran_indices);
    }
  }
  // Bind_I_Teth: same as bind_i but for tethered populations
  if (tethering_active_) {
    event_name = "bind_i_teth";
    auto exe_bind_i_teth = [&](ENTRY_T target) {
      Bind_I_Teth(std::get<POP_T *>(target));
    };
    auto poisson_i_teth = [&](double p, int n) {
      double weight{GetWeight_Bind_I_Teth()};
      if (weight == 0.0) {
        return 0;
      }
      int n_expected{properties_->gsl.SamplePoissonDist(p * weight)};
      if (n_expected > 0) {
        if (n_expected > n) {
          n_expected = n;
        }
        Set_Bind_I_Teth_Candidates(n_expected);
      }
      return n_expected;
    };
  }
  // Bind_ATP: bind ATP to motor heads that have released their ADP
  auto exe_bind_ATP = [&](ENTRY_T target) {
    POP_T *head = std::get<POP_T *>(target);
    if (head != nullptr) {
      ReportExecutionOf("bind_ATP");
      Bind_ATP(head);
    } else {
      ReportFailureOf("bind_ATP");
    }
  };
  int max_neighbs_step{max_neighbs_ - 1};
  if (max_neighbs_ == 0) {
    max_neighbs_step = 0;
  }
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_step; n_neighbs++) {
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    event_name = "bind_ATP";
    if (max_neighbs_ > 0) {
      event_name += N_NEIGHBS;
    }
    events_.emplace_back(event_name, p_bind_ATP_[n_neighbs],
                         &n_bound_NULL_[n_neighbs], &bound_NULL_[n_neighbs],
                         exe_bind_ATP, binomial, set_ran_indices);
    if (!tethering_active_) {
      continue;
    }
    // Bind_ATP_Teth
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      std::string X_DUB{"_" + std::to_string(x_dub)};
      event_name = "bind_ATP_to_teth" + N_NEIGHBS + X_DUB;
      events_.emplace_back(event_name, p_bind_ATP_to_teth_[n_neighbs][x_dub],
                           &n_bound_NULL_to_teth_[n_neighbs][x_dub],
                           &bound_NULL_to_teth_[n_neighbs][x_dub], exe_bind_ATP,
                           binomial, set_ran_indices);
      event_name = "bind_ATP_fr_teth" + N_NEIGHBS + X_DUB;
      events_.emplace_back(event_name, p_bind_ATP_fr_teth_[n_neighbs][x_dub],
                           &n_bound_NULL_fr_teth_[n_neighbs][x_dub],
                           &bound_NULL_fr_teth_[n_neighbs][x_dub], exe_bind_ATP,
                           binomial, set_ran_indices);
    }
  }
  // Hydrolyze: convert ATP to ADPP
  event_name = "hydrolyze";
  auto exe_hydrolyze = [&](ENTRY_T target) {
    POP_T *head = std::get<POP_T *>(target);
    if (head != nullptr) {
      ReportExecutionOf("hydrolyze");
      Hydrolyze(head);
    } else {
      ReportFailureOf("hydrolyze");
    }
  };
  events_.emplace_back(event_name, p_hydrolyze_, &n_bound_ATP_, &bound_ATP_,
                       exe_hydrolyze, binomial, set_ran_indices);
  // Hydrolyze_stalled: same as Hydrolyze but only targets stalled motors
  event_name = "hydrolyze_st";
  events_.emplace_back(event_name, p_hydrolyze_st_, &n_bound_ATP_st_,
                       &bound_ATP_st_, exe_hydrolyze, binomial,
                       set_ran_indices);
  // Bind_II: binds docked motor head to MT and releases its ADP
  auto exe_bind_ii = [&](ENTRY_T target) {
    POP_T *head = std::get<POP_T *>(target);
    if (head != nullptr) {
      ReportExecutionOf("bind_ii");
      Bind_II(head);
    } else {
      ReportFailureOf("bind_ii");
    }
  };
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    std::string AFFINITY{"_" + std::to_string(i_aff)};
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
      event_name = "bind_ii";
      if (n_affinities_ > 1) {
        event_name += AFFINITY;
      }
      if (max_neighbs_ > 0) {
        event_name += N_NEIGHBS;
      }
      events_.emplace_back(event_name, p_bind_ii_[i_aff][n_neighbs],
                           &n_docked_[i_aff][n_neighbs],
                           &docked_[i_aff][n_neighbs], exe_bind_ii, binomial,
                           set_ran_indices);
    }
  }
  // Unbind_II: Converts ADPP to ADP and unbinds a doubly-bound head
  auto exe_unbind_ii = [&](ENTRY_T target) {
    POP_T *head = std::get<POP_T *>(target);
    if (head != nullptr) {
      ReportExecutionOf("unbind_ii");
      Unbind_II(head);
    } else {
      ReportFailureOf("unbind_ii");
      printf("step #%lu!!\n", properties_->current_step_);
      exit(1);
    }
  };
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    std::string AFFINITY{"_" + std::to_string(i_aff)};
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
      event_name = "unbind_ii";
      if (n_affinities_ > 1) {
        event_name += AFFINITY;
      }
      if (max_neighbs_ > 0) {
        event_name += N_NEIGHBS;
      }
      events_.emplace_back(event_name, p_unbind_ii_[i_aff][n_neighbs],
                           &n_bound_ADPP_ii_[i_aff][n_neighbs],
                           &bound_ADPP_ii_[i_aff][n_neighbs], exe_unbind_ii,
                           binomial, set_ran_indices);
    }
  }
  // Unbind_I: Converts ADPP to ADP and unbinds a singly-bound head
  auto exe_unbind_i = [&](ENTRY_T target) {
    POP_T *head = std::get<POP_T *>(target);
    if (head != nullptr) {
      ReportExecutionOf("unbind_i");
      Unbind_I(head);
    } else {
      ReportFailureOf("unbind_i");
    }
  };
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    std::string AFFINITY{"_" + std::to_string(i_aff)};
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
      event_name = "unbind_i";
      if (n_affinities_ > 1) {
        event_name += AFFINITY;
      }
      if (max_neighbs_ > 0) {
        event_name += N_NEIGHBS;
      }
      events_.emplace_back(event_name, p_unbind_i_[i_aff][n_neighbs],
                           &n_bound_ADPP_i_[i_aff][n_neighbs],
                           &bound_ADPP_i_[i_aff][n_neighbs], exe_unbind_i,
                           binomial, set_ran_indices);
      // Unbind_I_Stalled: Same as unbind_I but only targets stalled motors
      event_name = "unbind_i_st";
      if (n_affinities_ > 1) {
        event_name += AFFINITY;
      }
      if (max_neighbs_ > 0) {
        event_name += N_NEIGHBS;
      }
      events_.emplace_back(event_name, p_unbind_i_st_[i_aff][n_neighbs],
                           &n_bound_ADPP_i_st_[i_aff][n_neighbs],
                           &bound_ADPP_i_st_[i_aff][n_neighbs], exe_unbind_i,
                           binomial, set_ran_indices);
      if (!tethering_active_) {
        continue;
      }
      for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
        std::string X_DUB{"_" + std::to_string(x_dub)};
        // Unbind_I_Tethered
        event_name = "unbind_i_teth" + AFFINITY + N_NEIGHBS + X_DUB;
        events_.emplace_back(event_name,
                             p_unbind_i_teth_[i_aff][n_neighbs][x_dub],
                             &n_bound_ADPP_i_teth_[i_aff][n_neighbs][x_dub],
                             &bound_ADPP_i_teth_[i_aff][n_neighbs][x_dub],
                             exe_unbind_i, binomial, set_ran_indices);
        // Unbind_I_Tethered_Stalled
        event_name = "unbind_i_teth_st" + AFFINITY + N_NEIGHBS + X_DUB;
        events_.emplace_back(event_name,
                             p_unbind_i_teth_st_[i_aff][n_neighbs][x_dub],
                             &n_bound_ADPP_i_teth_st_[i_aff][n_neighbs][x_dub],
                             &bound_ADPP_i_teth_st_[i_aff][n_neighbs][x_dub],
                             exe_unbind_i, binomial, set_ran_indices);
      }
    }
  }
  if (tethering_active_) {
    // Tether_bound
    auto exe_tether_bound = [&](ENTRY_T target) {
      Tether_Bound(std::get<POP_T *>(target));
    };
    auto poisson_tether_bound = [&](int n, double p) {
      double weight{GetWeight_Tether_Bound()};
      if (weight == 0.0) {
        return 0;
      }
      int n_expected{properties_->gsl.SamplePoissonDist(p * weight)};
      if (n_expected > 0) {
        if (n_expected > n) {
          n_expected = n;
        }
        Set_Tether_Bound_Candidates(n_expected);
      }
      return n_expected;
    };
    event_name = "tether_bound";
    events_.emplace_back(event_name, p_avg_tether_bound_,
                         &n_tether_bound_candidates_, &tether_bound_candidates_,
                         exe_tether_bound, poisson_tether_bound,
                         set_ran_indices);
    // Tether_free
    auto exe_tether_free = [&](ENTRY_T target) {
      Tether_Free(std::get<ALT_T *>(target));
    };
    event_name = "tether_free";
    events_.emplace_back(event_name, p_tether_free_,
                         &properties_->prc1.n_bound_unteth_,
                         &properties_->prc1.bound_unteth_, exe_tether_free,
                         binomial, set_ran_indices);
    // Untether_bound
    auto exe_untether = [&](ENTRY_T target) {
      Untether(std::get<POP_T *>(target));
    };
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      std::string X_DUB{"_" + std::to_string(x_dub)};
      event_name = "untether_bound" + X_DUB;
      events_.emplace_back(event_name, p_untether_bound_[x_dub],
                           &n_bound_teth_[x_dub], &bound_teth_[x_dub],
                           exe_untether, binomial, set_ran_indices);
    }
    // Untether_satellite
    event_name = "untether_satellite";
    events_.emplace_back(event_name, p_untether_satellite_, &n_satellites_,
                         &satellites_, exe_untether, binomial, set_ran_indices);
  }
  if (verbosity_ >= 3) {
    wally_->Log("\nMotor events: \n");
    for (const auto &event : events_) {
      wally_->Log("   %s\n", event.name_.c_str());
    }
  }
}

void KinesinManagement::ReportProbabilities() {

  if (!population_active_) {
    return;
  }
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
            printf("Couldn't find %s\n", name.c_str());
            continue;
          }
          if (event->n_opportunities_tot_ == 0) {
            printf("No statistics for %s\n", name.c_str());
            continue;
          }
          double n_exe_tot{(double)event->n_executed_tot_};
          double n_opp_tot{(double)event->n_opportunities_tot_};
          if (n_opp_tot < 0.0) {
            printf("WHAT?? n_opp = %g\n", n_opp_tot);
            exit(1);
          }
          wally_->Log("For Kin event %s:\n", event->name_.c_str());
          wally_->Log("   p_theory = %g\n", value[i][j][k]);
          wally_->Log("   p_actual = %g", n_exe_tot / n_opp_tot);
          wally_->Log(" (n_exe = %i)\n", event->n_executed_tot_);
        }
      }
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
  printf("Executing ");
  std::cout << event_name << std::endl;
}

void KinesinManagement::ReportFailureOf(std::string event_name) {

  printf("yo we failed to execute ");
  std::cout << event_name << std::endl;
}

void KinesinManagement::FlagForUpdate() { lists_up_to_date_ = false; }

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
    printf("starting UPDATE DOCKED\n");
  }
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      n_docked_[i_aff][n_neighbs] = 0;
    }
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    Kinesin *motor = active_[i_entry];
    if (motor->heads_active_ != 1) {
      continue;
    }
    if (motor->GetActiveHead()->ligand_ == "ADPP") {
      Tubulin *dock_site{motor->GetDockSite()};
      if (!dock_site->occupied_) {
        int aff{0};
        int n_neighbs_eff{0};
        if (max_neighbs_ > 0) {
          n_neighbs_eff = dock_site->GetKif4ANeighborCount() - 1;
        }
        int index{n_docked_[aff][n_neighbs_eff]++};
        docked_[aff][n_neighbs_eff][index] = motor->GetDockedHead();
      }
    }
  }
}

void KinesinManagement::Update_Bound_NULL() {

  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    n_bound_NULL_[n_neighbs] = 0;
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    Kinesin *motor = active_[i_entry];
    if (!motor->tethered_ or motor->HasSatellite()) {
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
    }
  }
}

void KinesinManagement::Update_Bound_NULL_Teth() {

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
}

void KinesinManagement::Update_Bound_ATP() {

  n_bound_ATP_ = 0;
  for (int i_motor = 0; i_motor < n_active_; i_motor++) {
    Kinesin *motor = active_[i_motor];
    if (!motor->IsStalled()) {
      if (motor->head_one_.ligand_ == "ATP") {
        bound_ATP_[n_bound_ATP_++] = &motor->head_one_;
      }
      if (motor->head_two_.ligand_ == "ATP") {
        bound_ATP_[n_bound_ATP_++] = &motor->head_two_;
      }
    }
  }
}

void KinesinManagement::Update_Bound_ATP_Stalled() {

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
}

void KinesinManagement::Update_Bound_ADPP_I() {

  if (verbosity_ >= 1) {
    printf("starting UPDATE ADPP_I\n");
  }
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    for (int i_neighb{0}; i_neighb <= max_neighbs_; i_neighb++) {
      n_bound_ADPP_i_[i_aff][i_neighb] = 0;
    }
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    Kinesin *motor = active_[i_entry];
    if (motor->heads_active_ != 1 or motor->IsStalled()) {
      continue;
    }
    if (!motor->tethered_ or motor->HasSatellite()) {
      POP_T *active_head = motor->GetActiveHead();
      if (active_head->ligand_ == "ADPP") {
        int aff = active_head->GetAffinity();
        int neighbs = active_head->GetKif4ANeighborCount();
        int index = n_bound_ADPP_i_[aff][neighbs]++;
        bound_ADPP_i_[aff][neighbs][index] = active_head;
      }
    }
  }
}

void KinesinManagement::Update_Bound_ADPP_I_Stalled() {

  if (verbosity_ >= 1) {
    printf("starting UPDATE ADPP_I_ST\n");
  }
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
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

void KinesinManagement::Update_Bound_ADPP_I_Teth() {

  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
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

void KinesinManagement::Update_Bound_ADPP_I_Teth_Stalled() {

  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
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

void KinesinManagement::Update_Bound_ADPP_II() {

  if (verbosity_ >= 1) {
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
      Kinesin::Monomer *chosen_head{nullptr};
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
        int neighbs{chosen_head->GetKif4ANeighborCount()};
        int i_entry{n_bound_ADPP_ii_[aff][neighbs]++};
        bound_ADPP_ii_[aff][neighbs][i_entry] = chosen_head;
      }
    }
  }
}

void KinesinManagement::Update_Bound_Unteth() {

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
}

void KinesinManagement::Update_Bound_Teth() {

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
}

void KinesinManagement::Update_Free_Teth() {

  n_satellites_ = 0;
  n_bind_i_teth_candidates_ = 0;
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    Kinesin *motor = active_[i_entry];
    if (motor->heads_active_ == 0 and motor->tethered_) {
      satellites_[n_satellites_++] = &motor->head_one_;
      bind_i_teth_candidates_[n_bind_i_teth_candidates_++] = &motor->head_one_;
    }
  }
}

double KinesinManagement::GetWeight_Bind_I_Teth() {
  double weight = 0;
  for (int i_entry = 0; i_entry < n_satellites_; i_entry++) {
    POP_T *head = std::get<POP_T *>(satellites_[i_entry]);
    weight += head->motor_->GetTotalWeight_Bind_I_Teth();
  }
  return weight;
}

double KinesinManagement::GetWeight_Tether_Bound() {
  double weight = 0;
  for (int i_entry = 0; i_entry < n_bound_unteth_; i_entry++) {
    POP_T *head = std::get<POP_T *>(bound_unteth_[i_entry]);
    weight += head->motor_->GetTotalWeight_Tether_Bound();
  }
  return weight;
}

void KinesinManagement::Set_Bind_I_Teth_Candidates(int n_to_set) {
  double weight[n_bind_i_teth_candidates_];
  double weight_total{0.0};
  for (int i_entry{0}; i_entry < n_bind_i_teth_candidates_; i_entry++) {
    POP_T *head = std::get<POP_T *>(bind_i_teth_candidates_[i_entry]);
    weight[i_entry] = head->motor_->GetTotalWeight_Bind_I_Teth();
    weight_total += weight[i_entry];
  }
  if (weight_total == 0.0) {
    wally_->ErrorExit("Kin_MGMT::Set_Bind_I_Teth_Candidates()");
  }
  ENTRY_T selected_candidates[n_to_set];
  for (int i_set{0}; i_set < n_to_set; i_set++) {
    double p_cum{0.0};
    double ran{properties_->gsl.GetRanProb()};
    for (int i_entry{0}; i_entry < n_bind_i_teth_candidates_; i_entry++) {
      p_cum += (weight[i_entry] / weight_total);
      if (ran < p_cum) {
        selected_candidates[i_set] = bind_i_teth_candidates_[i_entry];
        weight_total -= weight[i_entry];
        weight[i_entry] = weight[n_bind_i_teth_candidates_ - 1];
        n_bind_i_teth_candidates_--;
        break;
      }
    }
  }
  n_bind_i_teth_candidates_ = n_to_set;
  for (int i_entry{0}; i_entry < n_bind_i_teth_candidates_; i_entry++) {
    bind_i_teth_candidates_[i_entry] = selected_candidates[i_entry];
  }
}

void KinesinManagement::Set_Tether_Bound_Candidates(int n_to_set) {

  double weight[n_tether_bound_candidates_];
  double weight_total{0.0};
  for (int i_entry{0}; i_entry < n_tether_bound_candidates_; i_entry++) {
    POP_T *head = std::get<POP_T *>(tether_bound_candidates_[i_entry]);
    weight[i_entry] = head->motor_->GetTotalWeight_Tether_Bound();
    weight_total += weight[i_entry];
  }
  if (weight_total == 0.0) {
    wally_->ErrorExit("Kin_MGMT::Set_Tether_Bound_Candidates()");
  }
  ENTRY_T selected_candidates[n_to_set];
  for (int i_set{0}; i_set < n_to_set; i_set++) {
    double p_cum{0.0};
    double ran{properties_->gsl.GetRanProb()};
    for (int i_entry{0}; i_entry < n_tether_bound_candidates_; i_entry++) {
      p_cum += (weight[i_entry] / weight_total);
      if (ran < p_cum) {
        selected_candidates[i_set] = tether_bound_candidates_[i_entry];
        weight_total -= weight[i_entry];
        weight[i_entry] = weight[n_tether_bound_candidates_ - 1];
        n_tether_bound_candidates_--;
        break;
      }
    }
  }
  n_tether_bound_candidates_ = n_to_set;
  for (int i_entry{0}; i_entry < n_tether_bound_candidates_; i_entry++) {
    tether_bound_candidates_[i_entry] = selected_candidates[i_entry];
  }
}

void KinesinManagement::RunKMC() {

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

  /*
  if (lists_up_to_date_) {
    return;
  }
  lists_up_to_date_ = true;
  */

  Update_Extensions();
  properties_->microtubules.UpdateUnoccupied();
  Update_Docked();
  Update_Bound_NULL();
  Update_Bound_ATP();
  Update_Bound_ATP_Stalled();
  Update_Bound_ADPP_I();
  Update_Bound_ADPP_I_Stalled();
  Update_Bound_ADPP_II();
  if (parameters_->motors.tethers_active) {
    Update_Free_Teth();
    Update_Bound_NULL_Teth();
    Update_Bound_ADPP_I_Teth();
    Update_Bound_ADPP_I_Teth_Stalled();
    Update_Bound_Teth();
    Update_Bound_Unteth();
    properties_->prc1.Update_Bound_Unteth();
  }
}

void KinesinManagement::SampleEventStatistics() {

  // Scan through all events & get expected number of occurrences
  n_events_to_exe_ = 0;
  for (auto &&event : events_) {
    n_events_to_exe_ += event.SampleStatistics();
    // Ensure no funny business occurs with statistics
    if (event.n_expected_ > *event.n_avail_) {
      printf("Error; %i events expected but only %i available for %s\n",
             event.n_expected_, *event.n_avail_, event.name_.c_str());
      wally_->ErrorExit("Kin_MGMT::SampleEventStatistics()");
    }
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
  // Scan through all active events to ensure that no two target the same motor
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
      double ran{properties_->gsl.GetRanProb()};
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
    properties_->gsl.Shuffle(pre_array, n_events_to_exe_, sizeof(EVENT_T *));
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
    events_to_exe_[i_event]->Execute();
    FlagForUpdate();
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
    printf("failed to Bind_I_Free_Teth (MOTOR)\n");
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
  }
}

void KinesinManagement::Hydrolyze(POP_T *head) { head->ligand_ = "ADPP"; }

void KinesinManagement::Bind_II(POP_T *head) {

  // Verify that proposed site is unoccupied
  Kinesin *motor{head->motor_};
  Tubulin *dock_site{motor->GetDockSite()};
  if (dock_site == nullptr) {
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