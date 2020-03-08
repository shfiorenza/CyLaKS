#include "associated_protein_management.h"
#include "master_header.h"

AssociatedProteinManagement::AssociatedProteinManagement() {}

void AssociatedProteinManagement::Initialize(system_parameters *parameters,
                                             system_properties *properties) {

  wally_ = &properties->wallace;
  parameters_ = parameters;
  properties_ = properties;
  CalculateCutoffs();
  SetParameters();
  GenerateXLinks();
  InitializeLists();
  if (wally_->test_mode_ == nullptr) {
    InitializeEvents();
  } else {
    InitializeTestEvents();
    InitializeTestEnvironment();
  }
}

void AssociatedProteinManagement::CalculateCutoffs() {

  double kbT{parameters_->kbT};
  double r_0{parameters_->xlinks.r_0};
  double r_y{parameters_->microtubules.y_dist};
  double k_spring{parameters_->xlinks.k_spring};
  double site_size{parameters_->microtubules.site_size};
  /* First, calculate rest_dist_ in number of sites */
  // For now, it should always come out to zero
  int approx_rest{(int)(sqrt(r_0 * r_0 - r_y * r_y) / site_size)};
  if (approx_rest == 0) {
    rest_dist_ = approx_rest;
  } else {
    wally_->ErrorExit("AP_MGMT::CalculateCutoffs()");
  }
  double force_cutoff{20}; // in pN
  double f_spring{0.0};
  /* next, calculate extension distance cutoff */
  for (int x_dist{rest_dist_}; x_dist < 1000; x_dist++) {
    int r_x = x_dist * site_size;
    double r = sqrt(r_y * r_y + r_x * r_x);
    double dr = r - r_0;
    f_spring = fabs(dr * k_spring);
    /*
    if (f_spring >= force_cutoff) {
      dist_cutoff_ = x_dist;
      break;
    }
    */
    double U = (k_spring / 2) * dr * dr;
    double boltzmann_weight = exp(U / (2 * kbT));
    if (boltzmann_weight > 1e2) {
      dist_cutoff_ = x_dist;
      break;
    }
  }
  /*
  if (dist_cutoff_ > 9) {
    printf("As of now, x_dist must be less than 10");
    printf(" (currenly %i with k_spring=%g)\n", dist_cutoff_, k_spring);
    wally_->ErrorExit("AP_MGMT::CalculateCutoffs() [2]");
  }
  */
  if (parameters_->microtubules.count > 1) {
    crosslinking_active_ = true;
    wally_->Log("\nFor crosslinkers:\n");
    wally_->Log("  rest_dist is %i\n", rest_dist_);
    wally_->Log("  dist_cutoff is %i (f = %g pN)\n\n", dist_cutoff_, f_spring);
  } else {
    wally_->Log("\nCrosslinker double-binding is inactive.\n\n");
  }
}

void AssociatedProteinManagement::SetParameters() {

  if (parameters_->xlinks.c_bulk > 0.0) {
    population_active_ = true;
  }
  if (parameters_->motors.tethers_active and parameters_->motors.c_bulk > 0.0) {
    tethering_active_ = true;
  }
  /*
    For events that result in a change in energy, we use Boltzmann factors to
    scale rates appropriately. Detailed balance is satisfied with the factors:
                  exp{-(1 - lambda)*(delta_E)/(kB*T)}, and
                  exp{(lambda)*(delta_E)/(kB*T)}
    for an event and its complement (e.g., binding and unbinding), where lambda
    is a constant that ranges from 0 to 1, delta_E is the change in energy that
    results from the event, kB is Boltzmann's constant, and T is the
    temperature
  */
  // Lambda = 1.0 means all of the weight goes into unbinding
  double lambda_neighb{1.0}; // Lambda for neighbor interaction energies
  // Lambda = 0.0 means all of the weight goes into binding
  double lambda_spring{0.5}; // Lambda for spring energies
  // Lambda = 0.5 means the weight is equally split between binding & unbinding
  double lambda_teth{0.5}; // lambda for tether spring energies
  // Calculate neighbor interaction energies
  interaction_energy_ = -1 * parameters_->xlinks.interaction_energy;
  std::vector<double> neighb_energy(max_neighbs_ + 1, 0.0);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    // Neighbor energies are negative since its an attractive potential
    neighb_energy[n_neighbs] = n_neighbs * interaction_energy_; //! in kBT
  }
  if (lambda_neighb != 1.0) {
    printf("Lambda != 1.0 for neighb interactions not implemented yet!\n");
    wally_->ErrorExit("AssociatedProteinManagement::SetParameters()\n");
  }
  // Calculate spring extension energies
  double site_size{parameters_->microtubules.site_size};
  double k_spring{parameters_->xlinks.k_spring};
  double r_y{parameters_->microtubules.y_dist};
  double r_0{parameters_->xlinks.r_0};
  std::vector<double> spring_energy(dist_cutoff_ + 1, 0.0);
  for (int x{0}; x <= dist_cutoff_; x++) {
    double r_x{x * site_size};
    double r{sqrt(r_x * r_x + r_y * r_y)};
    double dr{r - r_0};
    // Spring energies are positive by definition
    spring_energy[x] = 0.5 * k_spring * dr * dr; //! in pN*nm
    printf("spring_energy[%i] = %g\n", x, spring_energy[x]);
  }
  // Calculate tether extension energies
  double k_teth{parameters_->motors.k_spring};
  double k_slack{parameters_->motors.k_slack};
  double r_0_teth{parameters_->motors.r_0};
  double r_y_teth{parameters_->microtubules.y_dist / 2};
  double rest_dist_teth{properties_->kinesin4.rest_dist_};
  teth_cutoff_ = properties_->kinesin4.teth_cutoff_;
  comp_cutoff_ = properties_->kinesin4.comp_cutoff_;
  std::vector<double> teth_energy(2 * teth_cutoff_ + 1, 0.0);
  for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
    double r_x_teth{((double)x_dub / 2) * site_size};
    double r_teth{sqrt(r_x_teth * r_x_teth + r_y_teth * r_y_teth)};
    double dr_teth{r_teth - r_0_teth};
    // Tether spring energies are positive by definition
    if (dr_teth > 0) {
      // If tether is extended past rest dist, use k_teth as spring constant
      teth_energy[x_dub] = 0.5 * k_teth * dr_teth * dr_teth; // in pN*nm
    } else {
      // Otherwise, if tether is compressed past rest dist, use k_slack
      teth_energy[x_dub] = 0.5 * k_slack * dr_teth * dr_teth; // in pN*nm
    }
  }
  // [DIFFUSION STATISTICS FOR CROSSLINKER W/O TETH BELOW] //
  double kbT{parameters_->kbT};
  double delta_t{parameters_->delta_t};
  double x_squared{(site_size / 1000) * (site_size / 1000)}; //! in um^2
  // Characteristic timescale for singly-bound xlink to diffuse 1 site
  double tau_i{x_squared / (2 * parameters_->xlinks.diffu_coeff_i)};
  double p_diffu_i{parameters_->delta_t / tau_i};
  // Characteristic timescale for doubly-bound xlink to diffuse 1 site
  double tau_ii{x_squared / (2 * parameters_->xlinks.diffu_coeff_ii)};
  double p_diffu_ii{parameters_->delta_t / tau_ii};
  // We consider diff_fwd and diff_bck as two separate events, which
  // effectively doubles the probability to diffuse. To counteract this, divide
  // p_diff by 2.
  p_diffu_i /= 2;
  p_diffu_ii /= 2;
  p_diffuse_i_fwd_.resize(max_neighbs_ + 1);
  p_diffuse_i_bck_.resize(max_neighbs_ + 1);
  p_diffuse_ii_to_rest_.resize(max_neighbs_ + 1);
  p_diffuse_ii_fr_rest_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    // For neighb. interactions, we only consider the energy penalty for
    // diffusing away from neighbors (same as unbinding) since lambda = 1.0
    //     dE = E_f - E_i
    double dE = 0.0 - neighb_energy[n_neighbs];
    // dE has units of kbT, so dividing by its numerical value isn't necessary
    double weight_neighb{exp(-lambda_neighb * dE)};
    // With two neighbors, the binding head is jammed and cannot diffuse
    if (n_neighbs == 2) {
      weight_neighb = 0.0;
    }
    p_diffuse_i_fwd_[n_neighbs] = weight_neighb * p_diffu_i;
    p_diffuse_i_bck_[n_neighbs] = weight_neighb * p_diffu_i;
    p_diffuse_ii_to_rest_[n_neighbs].resize(dist_cutoff_ + 1);
    p_diffuse_ii_fr_rest_[n_neighbs].resize(dist_cutoff_ + 1);
    for (int x{0}; x <= dist_cutoff_; x++) {
      // x_dist if head were to diffuse towards (to) spring rest
      int x_to{0};
      if (x > 0) {
        x_to = x - 1;
      }
      // x_dist if head were to diffuse away from (fr) spring rest
      int x_fr{0};
      if (x < dist_cutoff_) {
        x_fr = x + 1;
      }
      // dU (U_f - U_i) if head diffuses towards (to) spring rest
      double dU_to{spring_energy[x_to] - spring_energy[x]};
      // dU (U_f - U_i) if head diffuses away from (fr) spring rest
      double dU_fr{spring_energy[x_fr] - spring_energy[x]};
      // Diffusing towards rest is considered an unbinding-type event in
      // regards to Boltzmann factors, since both events let the spring relax
      double weight_to{exp(-lambda_spring * dU_to / kbT)};
      // Diffusing away from rest is considered a binding-type event in regards
      // to Boltzmann factors, since both events stretch the spring out
      double weight_fr{exp(-(1.0 - lambda_spring) * dU_fr / kbT)};
      if (x == 0) {
        weight_to = 0.0;
        weight_fr *= 2;
      } else if (x == dist_cutoff_) {
        weight_fr = 0.0;
      }
      p_diffuse_ii_to_rest_[n_neighbs][x] =
          weight_neighb * weight_to * p_diffu_ii;
      p_diffuse_ii_fr_rest_[n_neighbs][x] =
          weight_neighb * weight_fr * p_diffu_ii;
    }
  }
  // [DIFFUSION STATISTICS INVOLVING TETHER BELOW] //
  p_diffuse_i_to_teth_rest_.resize(max_neighbs_ + 1);
  p_diffuse_i_fr_teth_rest_.resize(max_neighbs_ + 1);
  p_diffuse_ii_to_both_.resize(max_neighbs_ + 1);
  p_diffuse_ii_fr_both_.resize(max_neighbs_ + 1);
  p_diffuse_ii_to_self_fr_teth_.resize(max_neighbs_ + 1);
  p_diffuse_ii_fr_self_to_teth_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    p_diffuse_i_to_teth_rest_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    p_diffuse_i_fr_teth_rest_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    p_diffuse_ii_to_both_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    p_diffuse_ii_to_self_fr_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    p_diffuse_ii_fr_self_to_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    p_diffuse_ii_fr_both_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      // Change in x_dub that brings teth towards rest; -1 if extended
      // (default)
      int dx_rest{-1};
      // However, if tether is compressed, increased x_dub is towards rest
      if (x_dub < 2 * rest_dist_teth) {
        dx_rest = 1;
      }
      // dU (U_f - U_i) for tether if head steps towards (to) tether rest
      double dU_to{teth_energy[x_dub + dx_rest] - teth_energy[x_dub]};
      // dU (U_f - U_i) for tether if head steps away from (fr) tether rest
      double dU_fr{teth_energy[x_dub - dx_rest] - teth_energy[x_dub]};
      // Diffusing towards rest is considered an unbinding-type event in
      // regards to Boltzmann factors, since both events let the spring relax
      double weight_to_teth{exp(-lambda_teth * dU_to / kbT)};
      // Diffusing away from rest is considered a binding-type event in regards
      // to Boltzmann factors, since both events stretch the spring out
      double weight_fr_teth{exp(-(1.0 - lambda_teth) * dU_fr / kbT)};
      if (x_dub == 2 * comp_cutoff_ or x_dub == 2 * teth_cutoff_) {
        weight_fr_teth = 0.0;
      }
      p_diffuse_i_to_teth_rest_[n_neighbs][x_dub] =
          p_diffuse_i_fwd_[n_neighbs] * weight_to_teth;
      // For singly-bound xlinks, diffusing one site changes x_dub by 2, so
      // keep probability 0.0 for 2*cutoff +/- 1
      if (x_dub > 2 * comp_cutoff_ + 1 and x_dub < 2 * teth_cutoff_ - 1) {
        p_diffuse_i_fr_teth_rest_[n_neighbs][x_dub] =
            p_diffuse_i_bck_[n_neighbs] * weight_fr_teth;
      }
      p_diffuse_ii_to_both_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      p_diffuse_ii_fr_both_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      p_diffuse_ii_to_self_fr_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      p_diffuse_ii_fr_self_to_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      for (int x{0}; x <= dist_cutoff_; x++) {
        // Diffuse towards both own and tether rest
        p_diffuse_ii_to_both_[n_neighbs][x_dub][x] =
            p_diffuse_ii_to_rest_[n_neighbs][x] * weight_to_teth;
        // Diffuse away from both own and tether rest
        p_diffuse_ii_fr_both_[n_neighbs][x_dub][x] =
            p_diffuse_ii_fr_rest_[n_neighbs][x] * weight_fr_teth;
        // Diffuse towards own rest, but away from tether rest
        p_diffuse_ii_to_self_fr_teth_[n_neighbs][x_dub][x] =
            p_diffuse_ii_to_rest_[n_neighbs][x] * weight_fr_teth;
        // Diffuse away from own rest, but towards tether rest
        p_diffuse_ii_fr_self_to_teth_[n_neighbs][x_dub][x] =
            p_diffuse_ii_fr_rest_[n_neighbs][x] * weight_to_teth;
      }
    }
  }
  // KMC STATISTICS BELOW
  double k_on = parameters_->xlinks.k_on;
  double c_xlink = parameters_->xlinks.c_bulk;
  double c_eff_teth = parameters_->motors.c_eff_tether;
  p_avg_bind_i_teth_ = k_on * c_eff_teth * delta_t;
  double c_eff_bind = parameters_->xlinks.c_eff_bind;
  p_avg_bind_ii_ = k_on * c_eff_bind * delta_t;
  double k_off_i = parameters_->xlinks.k_off_i;
  double k_off_ii = parameters_->xlinks.k_off_ii;
  p_bind_i_.resize(max_neighbs_ + 1);
  p_unbind_i_.resize(max_neighbs_ + 1);
  p_unbind_ii_.resize(max_neighbs_ + 1);
  p_unbind_i_teth_.resize(max_neighbs_ + 1);
  p_unbind_ii_to_teth_.resize(max_neighbs_ + 1);
  p_unbind_ii_fr_teth_.resize(max_neighbs_ + 1);
  weight_bind_ii_.resize(max_neighbs_ + 1);
  weight_bind_i_teth_.resize(max_neighbs_ + 1);
  weight_bind_ii_to_teth_.resize(max_neighbs_ + 1);
  weight_bind_ii_fr_teth_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    // dE = E_f - E_i;
    double dE = 0.0 - neighb_energy[n_neighbs];
    // We only consider the energy penalty for unbinding since lambda = 1.0
    double weight_neighb_bind{1.0};
    // dE has units of kbT, so dividing by its numerical value isn't necessary
    double weight_neighb_unbind{exp(-lambda_neighb * dE)};
    p_bind_i_[n_neighbs] = weight_neighb_bind * k_on * c_xlink * delta_t;
    p_unbind_i_[n_neighbs] = weight_neighb_unbind * k_off_i * delta_t;
    p_unbind_ii_[n_neighbs].resize(dist_cutoff_ + 1);
    weight_bind_ii_[n_neighbs].resize(dist_cutoff_ + 1);
    for (int x{0}; x <= dist_cutoff_; x++) {
      // dU = U_f - U_i
      double dU_bind{spring_energy[x] - 0.0};
      double dU_unbind{0.0 - spring_energy[x]};
      double weight_bind{exp(-(1.0 - lambda_spring) * dU_bind / kbT)};
      weight_bind_ii_[n_neighbs][x] = weight_neighb_bind * weight_bind;
      printf("weight_bind_ii_[%i][%i] = %g\n", n_neighbs, x,
             weight_bind_ii_[n_neighbs][x]);
      double weight_unbind{exp(-lambda_spring * dU_unbind / kbT)};
      p_unbind_ii_[n_neighbs][x] =
          weight_neighb_unbind * weight_unbind * k_off_ii * delta_t;
    }
    // Rates involving both xlink & tether spring
    p_unbind_i_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    p_unbind_ii_to_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    p_unbind_ii_fr_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    weight_bind_i_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    weight_bind_ii_to_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    weight_bind_ii_fr_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      // dU = U_f - U_i
      double dU_bind{teth_energy[x_dub] - 0.0};
      double dU_unbind{0.0 - teth_energy[x_dub]};
      double weight_bind{exp(-(1.0 - lambda_teth) * dU_bind / kbT)};
      weight_bind_i_teth_[n_neighbs][x_dub] = weight_neighb_bind * weight_bind;
      double weight_unbind{exp(-lambda_teth * dU_unbind / kbT)};
      p_unbind_i_teth_[n_neighbs][x_dub] =
          weight_unbind * p_unbind_i_[n_neighbs];
      p_unbind_ii_to_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      p_unbind_ii_fr_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      weight_bind_ii_to_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      weight_bind_ii_fr_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      for (int x{0}; x <= dist_cutoff_; x++) {
        // Change in x_dub if 2nd xlink head were to (un)bind at current x_dist
        // (anchor position changes by x/2; then doubled to convert to x_dub)
        int dx_teth_dub = x;
        // Change in x_dub that brings teth towards rest; -1 if extended
        int dir_rest{-1};
        // However, if tether is compressed, increased x_dub is towards rest
        if (x_dub < 2 * rest_dist_teth) {
          dir_rest = 1;
        }
        // Final x_dub if xlink unbinds towards (to) teth rest
        int x_dub_to{x_dub + dir_rest * dx_teth_dub};
        // Final x_dub if xlink unbinds away from (fr) teth rest
        int x_dub_fr{x_dub - dir_rest * dx_teth_dub};
        if (x_dub_fr < 0 or x_dub_fr > 2 * teth_cutoff_) {
          x_dub_fr = 0;
        }
        // dU = U_f - U_i
        double dU_to{teth_energy[x_dub_to] - teth_energy[x_dub]};
        double dU_fr{teth_energy[x_dub_fr] - teth_energy[x_dub]};
        // (Un)binding towards rest is considered an unbinding-type event in
        // regards to Boltzmann factors, since both events let the spring relax
        double weight_to_teth{exp(-lambda_teth * dU_to / kbT)};
        // (Un)binding away from rest is considered a binding-type event in
        // regards to Boltzmann factors, since both events stretch the spring
        double weight_fr_teth{exp(-(1.0 - lambda_teth) * dU_fr / kbT)};
        // If tethers aren't active, all weights are automatically zero
        if (!parameters_->motors.tethers_active) {
          weight_to_teth = 0.0;
          weight_fr_teth = 0.0;
        }
        // To ensure we don't double-count sites for x = 0,
        // only use unbind_ii_fr_teth and ignore to_teth
        if (x == 0) {
          weight_to_teth = 0.0;
        }
        if (x_dub_fr < 2 * comp_cutoff_ or x_dub_fr > 2 * teth_cutoff_) {
          weight_fr_teth = 0.0;
        }
        weight_bind_ii_to_teth_[n_neighbs][x_dub][x] =
            weight_bind_ii_[n_neighbs][x] * weight_to_teth;
        weight_bind_ii_fr_teth_[n_neighbs][x_dub][x] =
            weight_bind_ii_[n_neighbs][x] * weight_fr_teth;
        p_unbind_ii_to_teth_[n_neighbs][x_dub][x] =
            p_unbind_ii_[n_neighbs][x] * weight_to_teth;
        p_unbind_ii_fr_teth_[n_neighbs][x_dub][x] =
            p_unbind_ii_[n_neighbs][x] * weight_fr_teth;
        if (p_unbind_ii_to_teth_[n_neighbs][x_dub][x] > 1.0) {
          printf("WARNING: p_unbind_ii_to_teth_[%i][%i][%i] = %g for xlinks\n",
                 n_neighbs, x_dub, x,
                 p_unbind_ii_to_teth_[n_neighbs][x_dub][x]);
        }
        if (p_unbind_ii_fr_teth_[n_neighbs][x_dub][x] > 1) {
          printf("WARNING: p_unbind_ii_fr_teth_[%i][%i][%i] = %g for xlinks\n",
                 n_neighbs, x_dub, x,
                 p_unbind_ii_fr_teth_[n_neighbs][x_dub][x]);
        }
      }
    }
  }
  double k_tether = parameters_->motors.k_tether;
  p_tether_free_ = k_tether * c_xlink * delta_t;
  double k_untether = parameters_->motors.k_untether;
  p_untether_free_ = k_untether * delta_t;
  // Add these johnnys to our p_theory_ map; need to pad so all are 3-D
  p_theory_["diffuse_i_fwd"] = {{p_diffuse_i_fwd_}};
  p_theory_["diffuse_i_bck"] = {{p_diffuse_i_bck_}};
  if (crosslinking_active_) {
    p_theory_["diffuse_ii_to_rest"] = {p_diffuse_ii_to_rest_};
    p_theory_["diffuse_ii_fr_rest"] = {p_diffuse_ii_fr_rest_};
  }
  if (tethering_active_) {
    p_theory_["diffuse_i_to_teth_rest"] = {p_diffuse_i_to_teth_rest_};
    p_theory_["diffuse_i_fr_teth_rest"] = {p_diffuse_i_fr_teth_rest_};
    if (crosslinking_active_) {
      p_theory_["diffuse_ii_to_both"] = p_diffuse_ii_to_both_;
      p_theory_["diffuse_ii_fr_both"] = p_diffuse_ii_fr_both_;
      p_theory_["diffuse_ii_to_self_fr_teth"] = p_diffuse_ii_to_self_fr_teth_;
      p_theory_["diffuse_ii_fr_self_to_teth"] = p_diffuse_ii_fr_self_to_teth_;
    }
  }
  p_theory_["bind_i"] = {{p_bind_i_}};
  p_theory_["unbind_i"] = {{p_unbind_i_}};
  if (crosslinking_active_) {
    p_theory_["bind_ii"] = {{{p_avg_bind_ii_}}};
    p_theory_["unbind_ii"] = {{p_unbind_ii_}};
  }
  if (tethering_active_) {
    p_theory_["bind_i_teth"] = {{{p_avg_bind_i_teth_}}};
    p_theory_["unbind_i_teth"] = {p_unbind_i_teth_};
    if (crosslinking_active_) {
      p_theory_["bind_ii_teth"] = {{{p_avg_bind_ii_}}};
      p_theory_["unbind_ii_to_teth"] = p_unbind_ii_to_teth_;
      p_theory_["unbind_ii_fr_teth"] = p_unbind_ii_fr_teth_;
    }
    p_theory_["tether_free"] = {{{p_tether_free_}}};
    p_theory_["untether_free"] = {{{p_untether_free_}}};
  }
  /*
  // Manually construct probabilities used in Poisson processes using weights
  Vec<Vec<Vec<double>>> p_bind_ii = {weight_bind_ii_};
  Vec<Vec<Vec<double>>> p_bind_i_teth = {weight_bind_i_teth_};
  Vec<Vec<Vec<double>>> p_bind_ii_to_teth = weight_bind_ii_to_teth_;
  Vec<Vec<Vec<double>>> p_bind_ii_fr_teth = weight_bind_ii_fr_teth_;
  for (int n_neighbs{0}; n_neighbs <= 2; n_neighbs++) {
    for (int x{0}; x <= dist_cutoff_; x++) {
      p_bind_ii[0][n_neighbs][x] *= p_avg_bind_ii_;
    }
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      p_bind_i_teth[0][n_neighbs][x_dub] *= p_avg_bind_i_teth_;
      for (int x{0}; x <= dist_cutoff_; x++) {
        p_bind_ii_to_teth[n_neighbs][x_dub][x] *= p_avg_bind_ii_;
        p_bind_ii_fr_teth[n_neighbs][x_dub][x] *= p_avg_bind_ii_;
      }
    }
  }
  // p_theory_["bind_ii"] = p_bind_ii;
  if (tethering_active_) {
    p_theory_["bind_i_teth"] = p_bind_i_teth;
    p_theory_["bind_ii_to_teth"] = p_bind_ii_to_teth;
    p_theory_["bind_ii_fr_teth"] = p_bind_ii_fr_teth;
  }
  */
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
          if (value[i][j][k] > 1.0) {
            wally_->Log("Error! %s = %g\n", name.c_str(), value[i][j][k]);
            wally_->ErrorExit("AP_MGMT:SetParameters()");
          }
        }
      }
    }
  }
}

void AssociatedProteinManagement::GenerateXLinks() {

  // Since only one head has to be bound, the sim will at most
  // as many xlinks as sites in the bulk (all single-bound)
  int n_mts = parameters_->microtubules.count;
  for (int i_mt{0}; i_mt < n_mts; i_mt++) {
    n_xlinks_ += parameters_->microtubules.length[i_mt];
  }
  xlinks_.resize(n_xlinks_);
  for (int id{0}; id < n_xlinks_; id++) {
    xlinks_[id].Initialize(parameters_, properties_, id);
  }
}

void AssociatedProteinManagement::InitializeLists() {

  // Stats (not a list ok bite me)
  n_bound_i_.resize(max_neighbs_ + 1);
  n_bound_ii_.resize(max_neighbs_ + 1);
  n_bound_i_teth_.resize(max_neighbs_ + 1);
  n_bound_ii_teth_same_.resize(max_neighbs_ + 1);
  n_bound_ii_teth_oppo_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    n_bound_i_[n_neighbs] = 0;
    n_bound_ii_[n_neighbs].resize(dist_cutoff_ + 1);
    for (int x_dist{0}; x_dist <= dist_cutoff_; x_dist++) {
      n_bound_ii_[n_neighbs][x_dist] = 0;
    }
    n_bound_i_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    n_bound_ii_teth_same_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    n_bound_ii_teth_oppo_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    for (int x_dub{0}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      n_bound_i_teth_[n_neighbs][x_dub] = 0;
      n_bound_ii_teth_same_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      n_bound_ii_teth_oppo_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      for (int x{0}; x <= dist_cutoff_; x++) {
        n_bound_ii_teth_same_[n_neighbs][x_dub][x] = 0;
        n_bound_ii_teth_oppo_[n_neighbs][x_dub][x] = 0;
      }
    }
  }
  // Lists
  active_.resize(n_xlinks_);
  bound_unteth_.resize(n_xlinks_);
  free_teth_.resize(n_xlinks_);
  bind_ii_candidates_.resize(n_xlinks_);
  bind_i_teth_candidates_.resize(n_xlinks_);
  bind_ii_teth_candidates_.resize(n_xlinks_);
  bound_i_.resize(max_neighbs_ + 1);
  bound_ii_.resize(max_neighbs_ + 1);
  bound_i_teth_.resize(max_neighbs_ + 1);
  bound_ii_teth_oppo_.resize(max_neighbs_ + 1);
  bound_ii_teth_same_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    bound_i_[n_neighbs].resize(n_xlinks_);
    bound_ii_[n_neighbs].resize(dist_cutoff_ + 1);
    for (int x{0}; x <= dist_cutoff_; x++) {
      bound_ii_[n_neighbs][x].resize(n_xlinks_);
    }
    bound_i_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    bound_ii_teth_oppo_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    bound_ii_teth_same_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    for (int x_dub{0}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      bound_i_teth_[n_neighbs][x_dub].resize(n_xlinks_);
      bound_ii_teth_oppo_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      bound_ii_teth_same_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      for (int x{0}; x <= dist_cutoff_; x++) {
        bound_ii_teth_oppo_[n_neighbs][x_dub][x].resize(n_xlinks_);
        bound_ii_teth_same_[n_neighbs][x_dub][x].resize(n_xlinks_);
      }
    }
  }
}

void AssociatedProteinManagement::InitializeEvents() {

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
    if (n == 0) {
      return 0;
    }
    return properties_->gsl.SampleBinomialDist(p, n);
  };
  /* * Event entries * */
  std::string event_name; // scratch space to construct each event name
  // Diffuse: steps head one site to the left or right
  auto exe_diffuse_fwd = [&](ENTRY_T target) {
    Diffuse(std::get<POP_T *>(target), 1);
  };
  auto exe_diffuse_bck = [&](ENTRY_T target) {
    Diffuse(std::get<POP_T *>(target), -1);
  };
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    // Diffuse_i_fwd
    event_name = "diffuse_i_fwd" + N_NEIGHBS;
    events_.emplace_back(event_name, p_diffuse_i_fwd_[n_neighbs],
                         &n_bound_i_[n_neighbs], &bound_i_[n_neighbs],
                         exe_diffuse_fwd, binomial, set_ran_indices);
    // Diffuse_i_bck
    event_name = "diffuse_i_bck" + N_NEIGHBS;
    events_.emplace_back(event_name, p_diffuse_i_bck_[n_neighbs],
                         &n_bound_i_[n_neighbs], &bound_i_[n_neighbs],
                         exe_diffuse_bck, binomial, set_ran_indices);
    if (!crosslinking_active_) {
      continue;
    }
    for (int x{0}; x <= dist_cutoff_; x++) {
      std::string X_DIST{"_" + std::to_string(x)};
      // Diffuse_ii_to_rest
      event_name = "diffuse_ii_to_rest" + N_NEIGHBS + X_DIST;
      events_.emplace_back(event_name, p_diffuse_ii_to_rest_[n_neighbs][x],
                           &n_bound_ii_[n_neighbs][x], &bound_ii_[n_neighbs][x],
                           exe_diffuse_fwd, binomial, set_ran_indices);
      // Diffuse_ii_fr_rest
      event_name = "diffuse_ii_fr_rest" + N_NEIGHBS + X_DIST;
      events_.emplace_back(event_name, p_diffuse_ii_fr_rest_[n_neighbs][x],
                           &n_bound_ii_[n_neighbs][x], &bound_ii_[n_neighbs][x],
                           exe_diffuse_bck, binomial, set_ran_indices);
    }
  }
  // Diffuse_teth: same as above but for tethered populations
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    if (!tethering_active_) {
      continue;
    }
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      std::string X_DUB{"_" + std::to_string(x_dub)};
      // Diffuse_i_to_teth rest
      event_name = "diffuse_i_to_teth_rest" + N_NEIGHBS + X_DUB;
      events_.emplace_back(
          event_name, p_diffuse_i_to_teth_rest_[n_neighbs][x_dub],
          &n_bound_i_teth_[n_neighbs][x_dub], &bound_i_teth_[n_neighbs][x_dub],
          exe_diffuse_fwd, binomial, set_ran_indices);
      // Diffuse_i_fr_teth_rest
      event_name = "diffuse_i_fr_teth_rest" + N_NEIGHBS + X_DUB;
      events_.emplace_back(
          event_name, p_diffuse_i_fr_teth_rest_[n_neighbs][x_dub],
          &n_bound_i_teth_[n_neighbs][x_dub], &bound_i_teth_[n_neighbs][x_dub],
          exe_diffuse_bck, binomial, set_ran_indices);
      if (!crosslinking_active_) {
        continue;
      }
      for (int x{0}; x <= dist_cutoff_; x++) {
        std::string X_DIST{"_" + std::to_string(x)};
        // Diffuse_ii_to_both_rests
        event_name = "diffuse_ii_to_both" + N_NEIGHBS + X_DUB + X_DIST;
        events_.emplace_back(event_name,
                             p_diffuse_ii_to_both_[n_neighbs][x_dub][x],
                             &n_bound_ii_teth_same_[n_neighbs][x_dub][x],
                             &bound_ii_teth_same_[n_neighbs][x_dub][x],
                             exe_diffuse_fwd, binomial, set_ran_indices);
        // Diffuse_ii_fr_both_rests
        event_name = "diffuse_ii_fr_both" + N_NEIGHBS + X_DUB + X_DIST;
        events_.emplace_back(event_name,
                             p_diffuse_ii_fr_both_[n_neighbs][x_dub][x],
                             &n_bound_ii_teth_same_[n_neighbs][x_dub][x],
                             &bound_ii_teth_same_[n_neighbs][x_dub][x],
                             exe_diffuse_bck, binomial, set_ran_indices);
        // Diffuse_ii_to_self_fr_teth_rest
        event_name = "diffuse_ii_to_self_fr_teth" + N_NEIGHBS + X_DUB + X_DIST;
        events_.emplace_back(event_name,
                             p_diffuse_ii_to_self_fr_teth_[n_neighbs][x_dub][x],
                             &n_bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                             &bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                             exe_diffuse_fwd, binomial, set_ran_indices);
        // Diffuse_ii_fr_self_to_teth_rest
        event_name = "diffuse_ii_fr_self_to_teth" + N_NEIGHBS + X_DUB + X_DIST;
        events_.emplace_back(event_name,
                             p_diffuse_ii_fr_self_to_teth_[n_neighbs][x_dub][x],
                             &n_bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                             &bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                             exe_diffuse_bck, binomial, set_ran_indices);
      }
    }
  }
  // Bind_I: binds first crosslinker head to the microtubule
  auto exe_bind_i = [&](ENTRY_T target) { Bind_I(std::get<SITE_T *>(target)); };
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    event_name = "bind_i" + N_NEIGHBS;
    events_.emplace_back(event_name, p_bind_i_[n_neighbs],
                         &properties_->microtubules.n_unocc_xlink_[n_neighbs],
                         &properties_->microtubules.unocc_xlink_[n_neighbs],
                         exe_bind_i, binomial, set_ran_indices);
  }
  // Bind_I_Teth: same as above but for tethered unbound xlinks (satellites)
  if (tethering_active_) {
    event_name = "bind_i_teth";
    auto exe_bind_i_teth = [&](ENTRY_T target) {
      Bind_I_Teth(std::get<POP_T *>(target));
    };
    // Poisson dist. is used w/ partition function for E-dependent binding
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
    events_.emplace_back(event_name, p_avg_bind_i_teth_,
                         &n_bind_i_teth_candidates_, &bind_i_teth_candidates_,
                         exe_bind_i_teth, poisson_i_teth, set_ran_indices);
  }
  // Bind II: binds the second crosslinker head to adjacent microtubule
  if (crosslinking_active_) {
    event_name = "bind_ii";
    auto exe_bind_ii = [&](ENTRY_T target) {
      try {
        Bind_II(std::get<POP_T *>(target));
      } catch (...) {
        wally_->ErrorExit("AP_MGMT::Bind_II()");
      }
    };
    // Poisson dist. is used w/ partition function for E-dependent binding
    auto poisson_ii = [&](double p, int n) {
      double weight{GetWeight_Bind_II()};
      if (weight == 0.0) {
        return 0;
      }
      int n_expected{properties_->gsl.SamplePoissonDist(p * weight)};
      if (n_expected > 0) {
        if (n_expected > n) {
          wally_->Log("Rescaled n_bind_ii from %i to %i\n", n_expected, n);
          n_expected = n;
        }
        int n_removed{Set_Bind_II_Candidates(n_expected)};
        n_expected -= n_removed;
      }
      return n_expected;
    };
    events_.emplace_back(event_name, p_avg_bind_ii_, &n_bind_ii_candidates_,
                         &bind_ii_candidates_, exe_bind_ii, poisson_ii,
                         set_ran_indices);
  }
  // Bind_II_Teth: same as above but for tethered population
  if (tethering_active_ and crosslinking_active_) {
    event_name = "bind_ii_teth";
    auto exe_bind_ii_teth = [&](ENTRY_T target) {
      Bind_II(std::get<POP_T *>(target));
    };
    // Poisson dist. is used w/ partition function for E-dependent binding
    auto poisson_ii_teth = [&](double p, int n) {
      double weight{GetWeight_Bind_II_Teth()};
      if (weight == 0.0) {
        return 0;
      }
      int n_expected{properties_->gsl.SamplePoissonDist(p * weight)};
      if (n_expected > 0) {
        if (n_expected > n) {
          n_expected = n;
        }
        Set_Bind_II_Teth_Candidates(n_expected);
      }
      return n_expected;
    };
    events_.emplace_back(event_name, p_avg_bind_ii_,
                         &n_bind_ii_teth_candidates_, &bind_ii_teth_candidates_,
                         exe_bind_ii_teth, poisson_ii_teth, set_ran_indices);
  }
  // Unbind II: unbinds a head of doubly-bound crosslinkers
  auto exe_unbind_ii = [&](ENTRY_T target) {
    Unbind_II(std::get<POP_T *>(target));
  };
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    if (!crosslinking_active_) {
      continue;
    }
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    for (int x{0}; x <= dist_cutoff_; x++) {
      std::string X_DIST{"_" + std::to_string(x)};
      event_name = "unbind_ii" + N_NEIGHBS + X_DIST;
      events_.emplace_back(event_name, p_unbind_ii_[n_neighbs][x],
                           &n_bound_ii_[n_neighbs][x], &bound_ii_[n_neighbs][x],
                           exe_unbind_ii, binomial, set_ran_indices);
    }
  }
  // Unbind_II_Teth: same as above but for tethered populations
  auto exe_unbind_ii_teth = [&](ENTRY_T target) {
    Unbind_II(std::get<POP_T *>(target));
  };
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    if (!crosslinking_active_ or !tethering_active_) {
      continue;
    }
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      std::string X_DUB{"_" + std::to_string(x_dub)};
      for (int x{0}; x <= dist_cutoff_; x++) {
        std::string X_DIST{"_" + std::to_string(x)};
        event_name = "unbind_ii_to_teth" + N_NEIGHBS + X_DUB + X_DIST;
        events_.emplace_back(event_name,
                             p_unbind_ii_to_teth_[n_neighbs][x_dub][x],
                             &n_bound_ii_teth_same_[n_neighbs][x_dub][x],
                             &bound_ii_teth_same_[n_neighbs][x_dub][x],
                             exe_unbind_ii_teth, binomial, set_ran_indices);
        event_name = "unbind_ii_fr_teth" + N_NEIGHBS + X_DUB + X_DIST;
        events_.emplace_back(event_name,
                             p_unbind_ii_fr_teth_[n_neighbs][x_dub][x],
                             &n_bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                             &bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                             exe_unbind_ii_teth, binomial, set_ran_indices);
      }
    }
  }
  // Unbind_I: unbinds a head of singly-bound crosslinkres
  auto exe_unbind_i = [&](ENTRY_T target) {
    Unbind_I(std::get<POP_T *>(target));
  };
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    event_name = "unbind_i" + N_NEIGHBS;
    events_.emplace_back(event_name, p_unbind_i_[n_neighbs],
                         &n_bound_i_[n_neighbs], &bound_i_[n_neighbs],
                         exe_unbind_i, binomial, set_ran_indices);
  }
  // Unbind_I_Teth: same as above but for tethered populations
  auto exe_unbind_i_teth = [&](ENTRY_T target) {
    Unbind_I(std::get<POP_T *>(target));
  };
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    if (!tethering_active_) {
      continue;
    }
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      std::string X_DUB{"_" + std::to_string(x_dub)};
      event_name = "unbind_i_teth" + N_NEIGHBS + X_DUB;
      events_.emplace_back(event_name, p_unbind_i_teth_[n_neighbs][x_dub],
                           &n_bound_i_teth_[n_neighbs][x_dub],
                           &bound_i_teth_[n_neighbs][x_dub], exe_unbind_i_teth,
                           binomial, set_ran_indices);
    }
  }
  if (tethering_active_) {
    // Tether_free: tethers free xlinks to bound motors
    event_name = "tether_free";
    auto exe_tether_free = [&](ENTRY_T target) {
      Tether_Free(std::get<ALT_T *>(target));
    };
    events_.emplace_back(event_name, p_tether_free_,
                         &properties_->kinesin4.n_bound_unteth_,
                         &properties_->kinesin4.bound_unteth_, exe_tether_free,
                         binomial, set_ran_indices);
    // Untether_satellite: untethers satellite xlinks from bound motors
    event_name = "untether_satellite";
    auto exe_untether_sat = [&](ENTRY_T target) {
      Untether_Free(std::get<POP_T *>(target));
    };
    events_.emplace_back(event_name, p_untether_free_, &n_free_teth_,
                         &free_teth_, exe_untether_sat, binomial,
                         set_ran_indices);
  }
}

void AssociatedProteinManagement::InitializeTestEvents() {

  printf("Initializing test_%s\n", wally_->test_mode_);
  std::string event_name;
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
    if (n == 0) {
      return 0;
    }
    return properties_->gsl.SampleBinomialDist(p, n);
  };
  if (strcmp(wally_->test_mode_, "bind_ii") == 0) {
    // Bind II: binds the second crosslinker head to adjacent microtubule
    auto exe_bind_ii = [&](ENTRY_T target) {
      POP_T *head{nullptr};
      try {
        head = std::get<POP_T *>(target);
      } catch (...) {
        wally_->ErrorExit("AP_MGMT::Bind_II()");
      }
      Bind_II(head);
      head->xlink_->UpdateExtension();
      bind_ii_stats_[head->xlink_->x_dist_].first++;
    };
    // Poisson dist. is used w/ partition function for E-dependent binding
    auto poisson_ii = [&](double p, int n) {
      for (int x{0}; x <= dist_cutoff_; x++) {
        if (x == 0) {
          // Only one site (directly above) can be bound to for this config
          bind_ii_stats_[x].second += n_bind_ii_candidates_;
        } else {
          // Two sites (to the left or right of above) are eligible for these
          // Include each as a separate opportunity
          bind_ii_stats_[x].second += 2 * n_bind_ii_candidates_;
        }
      }
      double weight{GetWeight_Bind_II()};
      if (weight == 0.0) {
        return 0;
      }
      int n_expected{properties_->gsl.SamplePoissonDist(p * weight)};
      if (n_expected > 0) {
        if (n_expected > n) {
          wally_->Log("Rescaled Bind_II exp from %i to %i\n", n_expected, n);
          n_expected = n;
        }
        int n_removed{Set_Bind_II_Candidates(n_expected)};
        n_expected -= n_removed;
      }
      return n_expected;
    };
    event_name = "bind_ii";
    events_.emplace_back(event_name, p_avg_bind_ii_, &n_bind_ii_candidates_,
                         &bind_ii_candidates_, exe_bind_ii, poisson_ii,
                         set_ran_indices);

    // Unbind II: unbinds a head of doubly-bound crosslinkers
    auto exe_unbind_ii = [&](ENTRY_T target) {
      Unbind_II(std::get<POP_T *>(target));
    };
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
      for (int x{0}; x <= dist_cutoff_; x++) {
        std::string X_DIST{"_" + std::to_string(x)};
        event_name = "unbind_ii" + N_NEIGHBS + X_DIST;
        events_.emplace_back(
            event_name, p_unbind_ii_[n_neighbs][x], &n_bound_ii_[n_neighbs][x],
            &bound_ii_[n_neighbs][x], exe_unbind_ii, binomial, set_ran_indices);
      }
    }
  }
  if (strcmp(wally_->test_mode_, "bind_i_teth") == 0) {
    // Bind_I_Teth: same as above but for tethered unbound xlinks (satellites)
    event_name = "bind_i_teth";
    auto exe_bind_i_teth = [&](ENTRY_T target) {
      Bind_I_Teth(std::get<POP_T *>(target));
    };
    // Poisson dist. is used w/ partition function for E-dependent binding
    auto poisson_i_teth = [&](double p, int n) {
      double weight{GetWeight_Bind_I_Teth()};
      if (weight == 0.0) {
        return 0;
      }
      int n_expected{properties_->gsl.SamplePoissonDist(p * weight)};
      if (n_expected > 0) {
        if (n_expected > n) {
          wally_->Log("Rescaled n_bind_i_teth from %i to %i\n", n_expected, n);
          n_expected = n;
        }
        int n_removed{Set_Bind_I_Teth_Candidates(n_expected)};
        n_expected -= n_removed;
      }
      return n_expected;
    };
    events_.emplace_back(event_name, p_avg_bind_i_teth_,
                         &n_bind_i_teth_candidates_, &bind_i_teth_candidates_,
                         exe_bind_i_teth, poisson_i_teth, set_ran_indices);
    // Unbind_I_Teth: same as above but for tethered populations
    auto exe_unbind_i_teth = [&](ENTRY_T target) {
      Unbind_I(std::get<POP_T *>(target));
    };
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
      for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
        std::string X_DUB{"_" + std::to_string(x_dub)};
        event_name = "unbind_i_teth" + N_NEIGHBS + X_DUB;
        events_.emplace_back(event_name, p_unbind_i_teth_[n_neighbs][x_dub],
                             &n_bound_i_teth_[n_neighbs][x_dub],
                             &bound_i_teth_[n_neighbs][x_dub],
                             exe_unbind_i_teth, binomial, set_ran_indices);
      }
    }
  }
}

void AssociatedProteinManagement::InitializeTestEnvironment() {

  if (strcmp(wally_->test_mode_, "bind_ii") == 0) {
    bind_ii_stats_.resize(dist_cutoff_ + 1);
    for (int x{0}; x <= dist_cutoff_; x++) {
      bind_ii_stats_[x].first = 0;
      bind_ii_stats_[x].second = 0;
    }
    AssociatedProtein *xlink{GetFreeXlink()};
    Tubulin *site{&properties_->microtubules.mt_list_[0].lattice_[50]};
    site->xlink_head_ = &xlink->head_one_;
    site->occupied_ = true;
    xlink->head_one_.site_ = site;
    xlink->heads_active_++;
    AddToActive(xlink);
    properties_->microtubules.FlagForUpdate();
  }
}

void AssociatedProteinManagement::ReportProbabilities() {

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
          wally_->Log("For AP event %s:\n", event->name_.c_str());
          wally_->Log("   p_theory = %g\n", value[i][j][k]);
          wally_->Log("   p_actual = %g", n_exe_tot / n_opp_tot);
          wally_->Log(" (n_exe = %i)\n", event->n_executed_tot_);
        }
      }
    }
  }
  if (wally_->test_mode_ != nullptr) {
    if (strcmp(wally_->test_mode_, "bind_ii") == 0) {
      printf("For bind_ii:\n");
      for (int x{0}; x <= dist_cutoff_; x++) {
        double p{(double)bind_ii_stats_[x].first / bind_ii_stats_[x].second};
        printf("p_theory_[%i] = %g\n", x,
               p_avg_bind_ii_ * weight_bind_ii_[0][x]);
        printf("p_actual_[%i] = %g (n_exe = %i)\n", x, p,
               bind_ii_stats_[x].first);
      }
    }
    if (strcmp(wally_->test_mode_, "bind_i_teth") == 0) {
      printf("For bind_i_teth:\n");
      for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
        double p{(double)bind_i_teth_stats_[x_dub].first /
                 bind_i_teth_stats_[x_dub].second};
        printf("p_theory_[%i] = %g\n", x_dub,
               p_avg_bind_i_teth_ * weight_bind_i_teth_[0][x_dub]);
        printf("p_actual_[%i] = %g (n_exe = %i)\n", x_dub, p,
               bind_i_teth_stats_[x_dub].first);
      }
    }
    if (strcmp(wally_->test_mode_, "bind_ii_teth") == 0) {
      printf("For bind_ii_to_teth:\n");
      for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
        for (int x{0}; x <= dist_cutoff_; x++) {
          double p_to{(double)bind_ii_to_teth_stats_[x_dub][x].first /
                      bind_ii_to_teth_stats_[x_dub][x].second};
          printf("p_theory_[%i][%i] = %g\n", x_dub, x,
                 p_avg_bind_ii_ * weight_bind_ii_to_teth_[0][x_dub][x]);
          printf("p_actual_[%i][%i] = %g (n_exe = %i)\n", x_dub, x, p_to,
                 bind_ii_to_teth_stats_[x_dub][x].first);
        }
      }
      printf("For bind_ii_fr_teth:\n");
      for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
        for (int x{0}; x <= dist_cutoff_; x++) {
          double p_fr{(double)bind_ii_fr_teth_stats_[x_dub][x].first /
                      bind_ii_fr_teth_stats_[x_dub][x].second};
          printf("p_theory_[%i][%i] = %g\n", x_dub, x,
                 p_avg_bind_ii_ * weight_bind_ii_fr_teth_[0][x_dub][x]);
          printf("p_actual_[%i][%i] = %g (n_exe = %i)\n", x_dub, x, p_fr,
                 bind_ii_fr_teth_stats_[x_dub][x].first);
        }
      }
    }
  }
}
AssociatedProtein *AssociatedProteinManagement::GetFreeXlink() {

  // Randomly choose an unbound xlink
  int i_xlink = properties_->gsl.GetRanInt(n_xlinks_);
  AssociatedProtein *xlink = &xlinks_[i_xlink];
  int attempts = 0;
  while (xlink->heads_active_ > 0 || xlink->tethered_) {
    i_xlink++;
    if (i_xlink == n_xlinks_)
      i_xlink = 0;
    xlink = &xlinks_[i_xlink];
    attempts++;
    if (attempts > n_xlinks_) {
      printf("error in get_free_xlink\n");
      exit(1);
    }
  }
  return xlink;
}

void AssociatedProteinManagement::AddToActive(AssociatedProtein *xlink) {

  active_[n_active_] = xlink;
  xlink->active_index_ = n_active_;
  n_active_++;
}

void AssociatedProteinManagement::RemoveFromActive(AssociatedProtein *xlink) {

  AssociatedProtein *last_entry = active_[n_active_ - 1];
  int this_index = xlink->active_index_;
  active_[this_index] = last_entry;
  last_entry->active_index_ = this_index;
  n_active_--;
}

void AssociatedProteinManagement::FlagForUpdate() { lists_up_to_date_ = false; }

void AssociatedProteinManagement::Update_Extensions() {

  if (!crosslinking_active_ and !tethering_active_) {
    return;
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    active_[i_entry]->UpdateExtension();
  }
}

void AssociatedProteinManagement::Update_Bound_I() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting AP_MGMT::Update_Bound_I()\n");
  }
  n_bind_ii_candidates_ = 0;
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    n_bound_i_[n_neighbs] = 0;
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    AssociatedProtein *xlink = active_[i_entry];
    if (xlink->heads_active_ == 1) {
      if (!xlink->tethered_ or xlink->HasSatellite()) {
        AssociatedProtein::Monomer *head = xlink->GetActiveHead();
        int n_neighbs = head->GetPRC1NeighborCount();
        bound_i_[n_neighbs][n_bound_i_[n_neighbs]++] = head;
        bind_ii_candidates_[n_bind_ii_candidates_++] = head;
      }
    }
  }
  if (wally_->test_mode_ != nullptr) {
    /*
    if (strcmp(wally_->test_mode_, "bind_ii") == 0) {
      // for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      //   for (int i_entry{0}; i_entry < n_bound_i_[n_neighbs]; i_entry++) {
      //     POP_T *head{std::get<POP_T *>(bound_i_[n_neighbs][i_entry])};
      //     if (head->site_->mt_ == 0) {
      //       int i_last{--n_bound_i_[n_neighbs]};
      //       bound_i_[n_neighbs][i_entry] = bound_i_[n_neighbs][i_last];
      //     }
      //   }
      // }
      for (int x{0}; x <= dist_cutoff_; x++) {
        if (x == 0) {
          // Only one site (directly above) can be bound to for this config
          bind_ii_stats_[x].second += n_bind_ii_candidates_;
        } else {
          // Two sites (to the left or right of above) are eligible for these
          // Include each as a separate opportunity
          bind_ii_stats_[x].second += 2 * n_bind_ii_candidates_;
        }
      }
    }
    */
  }
  if (verbosity_ >= 1) {
    wally_->Log(" n_bind_ii_candidates = %i\n", n_bind_ii_candidates_);
    for (int i_entry{0}; i_entry < n_bind_ii_candidates_; i_entry++) {
      POP_T *head{std::get<POP_T *>(bind_ii_candidates_[i_entry])};
      wally_->Log("   - xlink %i\n", head->xlink_->id_);
    }
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      wally_->Log(" n_bound_i_[%i] = %i\n", n_neighbs, n_bound_i_[n_neighbs]);
      for (int i_entry{0}; i_entry < n_bound_i_[n_neighbs]; i_entry++) {
        POP_T *head{std::get<POP_T *>(bound_i_[n_neighbs][i_entry])};
        wally_->Log("   - xlink %i\n", head->xlink_->id_);
      }
    }
  }
}

void AssociatedProteinManagement::Update_Bound_I_Teth() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting AP_MGMT::Update_Bound_I_Teth()\n");
  }
  n_bind_ii_teth_candidates_ = 0;
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      n_bound_i_teth_[n_neighbs][x_dub] = 0;
    }
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    AssociatedProtein *xlink = active_[i_entry];
    if (xlink->heads_active_ != 1) {
      continue;
    }
    if (xlink->tethered_ and !xlink->HasSatellite()) {
      int x_dub = xlink->motor_->x_dist_doubled_;
      AssociatedProtein::Monomer *head = xlink->GetActiveHead();
      int n_neighbs = head->GetPRC1NeighborCount();
      int index{n_bound_i_teth_[n_neighbs][x_dub]++};
      bound_i_teth_[n_neighbs][x_dub][index] = head;
      bind_ii_teth_candidates_[n_bind_ii_teth_candidates_++] = head;
    }
  }
}

void AssociatedProteinManagement::Update_Bound_II() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting AP_MGMT::Update_Bound_II()\n");
  }
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    for (int x{0}; x <= dist_cutoff_; x++) {
      n_bound_ii_[n_neighbs][x] = 0;
    }
  }
  for (int i_xlink{0}; i_xlink < n_active_; i_xlink++) {
    AssociatedProtein *xlink = active_[i_xlink];
    if (xlink->heads_active_ != 2) {
      continue;
    }
    if (!xlink->tethered_ or xlink->HasSatellite()) {
      int x = xlink->x_dist_;
      // Head one
      if (wally_->test_mode_ == nullptr) {
        int neighbs_one = xlink->head_one_.GetPRC1NeighborCount();
        bound_ii_[neighbs_one][x][n_bound_ii_[neighbs_one][x]++] =
            &xlink->head_one_;
      }
      // Head two
      int neighbs_two = xlink->head_two_.GetPRC1NeighborCount();
      bound_ii_[neighbs_two][x][n_bound_ii_[neighbs_two][x]++] =
          &xlink->head_two_;
    }
  }
}

void AssociatedProteinManagement::Update_Bound_II_Teth() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting AP_MGMT::Update_Bound_II_Teth()\n");
  }
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      for (int x{0}; x <= dist_cutoff_; x++) {
        n_bound_ii_teth_same_[n_neighbs][x_dub][x] = 0;
        n_bound_ii_teth_oppo_[n_neighbs][x_dub][x] = 0;
      }
    }
  }
  for (int i_xlink{0}; i_xlink < n_active_; i_xlink++) {
    AssociatedProtein *xlink = active_[i_xlink];
    if (xlink->heads_active_ == 2 and xlink->tethered_) {
      // Dont count xlinks tethered to satellite motors
      if (xlink->motor_->heads_active_ > 0) {
        int x = xlink->x_dist_;
        int x_dub = xlink->motor_->x_dist_doubled_;
        // Site one
        Tubulin *site_one = xlink->head_one_.site_;
        int n_one = site_one->GetPRC1NeighborCount();
        if (x == rest_dist_) {
          int i_one = n_bound_ii_teth_same_[n_one][x_dub][x];
          bound_ii_teth_same_[n_one][x_dub][x][i_one] = &xlink->head_one_;
          n_bound_ii_teth_same_[n_one][x_dub][x]++;
          int i_two = n_bound_ii_teth_oppo_[n_one][x_dub][x];
          bound_ii_teth_oppo_[n_one][x_dub][x][i_two] = &xlink->head_one_;
          n_bound_ii_teth_oppo_[n_one][x_dub][x]++;
        } else if (site_one->EquilibriumInSameDirection()) {
          int i = n_bound_ii_teth_same_[n_one][x_dub][x];
          bound_ii_teth_same_[n_one][x_dub][x][i] = &xlink->head_one_;
          n_bound_ii_teth_same_[n_one][x_dub][x]++;
        } else {
          int i = n_bound_ii_teth_oppo_[n_one][x_dub][x];
          bound_ii_teth_oppo_[n_one][x_dub][x][i] = &xlink->head_one_;
          n_bound_ii_teth_oppo_[n_one][x_dub][x]++;
        }
        // Site two
        Tubulin *site_two = xlink->head_two_.site_;
        int n_two = site_two->GetPRC1NeighborCount();
        if (x == rest_dist_) {
          int i_one = n_bound_ii_teth_same_[n_two][x_dub][x];
          bound_ii_teth_same_[n_two][x_dub][x][i_one] = &xlink->head_two_;
          n_bound_ii_teth_same_[n_two][x_dub][x]++;
          int i_two = n_bound_ii_teth_oppo_[n_two][x_dub][x];
          bound_ii_teth_oppo_[n_two][x_dub][x][i_two] = &xlink->head_two_;
          n_bound_ii_teth_oppo_[n_two][x_dub][x]++;
        } else if (site_two->EquilibriumInSameDirection()) {
          int i = n_bound_ii_teth_same_[n_two][x_dub][x];
          bound_ii_teth_same_[n_two][x_dub][x][i] = &xlink->head_two_;
          n_bound_ii_teth_same_[n_two][x_dub][x]++;
        } else {
          int i = n_bound_ii_teth_oppo_[n_two][x_dub][x];
          bound_ii_teth_oppo_[n_two][x_dub][x][i] = &xlink->head_two_;
          n_bound_ii_teth_oppo_[n_two][x_dub][x]++;
        }
      }
    }
  }
}

void AssociatedProteinManagement::Update_Bound_Unteth() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting AP_MGMT::Update_Bound_Unteth()\n");
  }
  n_bound_unteth_ = 0;
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    AssociatedProtein *xlink = active_[i_entry];
    if (xlink->tethered_) {
      continue;
    }
    if (xlink->heads_active_ == 1) {
      bound_unteth_[n_bound_unteth_++] = xlink->GetActiveHead();
    } else {
      bound_unteth_[n_bound_unteth_++] = &xlink->head_one_;
    }
  }
}

void AssociatedProteinManagement::Update_Free_Teth() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting AP_MGMT::Update_Free_Teth()\n");
  }
  n_free_teth_ = 0;
  n_bind_i_teth_candidates_ = 0;
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    AssociatedProtein *xlink = active_[i_entry];
    if (xlink->heads_active_ == 0 and xlink->tethered_) {
      if (xlink->motor_->heads_active_ == 0) {
        wally_->ErrorExit("AP_MGMT::Update_Free_Teth()\n");
      }
      free_teth_[n_free_teth_++] = &xlink->head_one_;
      bind_i_teth_candidates_[n_bind_i_teth_candidates_++] = &xlink->head_one_;
    }
  }
}

double AssociatedProteinManagement::GetWeight_Bind_II() {

  if (verbosity_ >= 1) {
    wally_->Log(" Starting AP_MGMT::GetWeight_Bind_II()\n");
  }
  double weight{0.0};
  for (int i_entry{0}; i_entry < n_bind_ii_candidates_; i_entry++) {
    POP_T *head{std::get<POP_T *>(bind_ii_candidates_[i_entry])};
    weight += head->xlink_->GetTotalWeight_Bind_II();
  }
  return weight;
}

double AssociatedProteinManagement::GetWeight_Bind_I_Teth() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting AP_MGMT::GetWeight_Bind_I_Teth()\n");
  }
  double weight{0.0};
  for (int i_entry{0}; i_entry < n_bind_i_teth_candidates_; i_entry++) {
    POP_T *head{std::get<POP_T *>(bind_i_teth_candidates_[i_entry])};
    weight += head->xlink_->GetTotalWeight_Bind_I_Teth();
  }
  return weight;
}

double AssociatedProteinManagement::GetWeight_Bind_II_Teth() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting AP_MGMT::GetWeight_Bind_II_Teth()\n");
  }
  if (verbosity_ >= 2) {
    wally_->Log("   %i candidates\n", n_bind_ii_teth_candidates_);
  }
  double weight{0.0};
  for (int i_entry{0}; i_entry < n_bind_ii_teth_candidates_; i_entry++) {
    if (verbosity_ >= 2) {
      wally_->Log("   accessing entry #%i\n", i_entry);
    }
    POP_T *head{std::get<POP_T *>(bind_ii_teth_candidates_[i_entry])};
    weight += head->xlink_->GetTotalWeight_Bind_II_Teth();
  }
  return weight;
}

int AssociatedProteinManagement::Set_Bind_II_Candidates(int n_to_set) {

  if (verbosity_ >= 1) {
    wally_->Log(" Starting AP_MGMT::Set_Bind_II_Candidates(%i)\n", n_to_set);
  }
  double weight[n_bind_ii_candidates_];
  double weight_total{0.0};
  for (int i_entry{0}; i_entry < n_bind_ii_candidates_; i_entry++) {
    POP_T *head = std::get<POP_T *>(bind_ii_candidates_[i_entry]);
    double weight_xlink{head->xlink_->GetTotalWeight_Bind_II()};
    // If entry has a weight of 0.0, remove it from candidates
    if (weight_xlink == 0.0) {
      int i_last{--n_bind_ii_candidates_};
      bind_ii_candidates_[i_entry] = bind_ii_candidates_[i_last];
      i_entry--;
      // If we have zero valid candidates, set n_expected to 0 in SampleStats
      if (n_bind_ii_candidates_ == 0) {
        return n_to_set;
      }
      continue;
    }
    weight[i_entry] = weight_xlink;
    weight_total += weight[i_entry];
    if (verbosity_ >= 2) {
      wally_->Log("   - xlink %i has weight %g\n", head->xlink_->id_,
                  weight[i_entry]);
    }
  }
  if (verbosity_ >= 2) {
    wally_->Log("   Total weight is %g\n", weight_total);
  }
  if (weight_total == 0.0) {
    wally_->ErrorExit("AP_MGMT::Set_Bind_II_Candidates()");
  }
  int n_removed{0};
  if (n_to_set > n_bind_ii_candidates_) {
    n_removed = n_to_set - n_bind_ii_candidates_;
    n_to_set = n_bind_ii_candidates_;
  }
  ENTRY_T selected_candidates[n_to_set];
  for (int i_set{0}; i_set < n_to_set; i_set++) {
    double p_cum{0.0};
    double ran{properties_->gsl.GetRanProb()};
    if (verbosity_ >= 2) {
      wally_->Log("   Rolled %g\n", ran);
    }
    for (int i_entry{0}; i_entry < n_bind_ii_candidates_; i_entry++) {
      p_cum += (weight[i_entry] / weight_total);
      if (verbosity_ >= 2) {
        wally_->Log("   p_cum[%i] = %g\n", i_entry, p_cum);
      }
      if (ran < p_cum) {
        selected_candidates[i_set] = bind_ii_candidates_[i_entry];
        weight_total -= weight[i_entry];
        if (verbosity_ >= 2) {
          wally_->Log("   selected entry %i\n", i_entry);
          wally_->Log("   total weight is now %g\n", weight_total);
        }
        bind_ii_candidates_[i_entry] =
            bind_ii_candidates_[n_bind_ii_candidates_ - 1];
        weight[i_entry] = weight[n_bind_ii_candidates_ - 1];
        n_bind_ii_candidates_--;
        break;
      }
    }
  }
  n_bind_ii_candidates_ = n_to_set;
  for (int i_entry{0}; i_entry < n_bind_ii_candidates_; i_entry++) {
    bind_ii_candidates_[i_entry] = selected_candidates[i_entry];
  }
  return n_removed;
}

int AssociatedProteinManagement::Set_Bind_I_Teth_Candidates(int n_to_set) {

  if (verbosity_ >= 1) {
    wally_->Log(" Starting AP_MGMT::Set_Bind_I_Teth_Candidates()\n");
  }
  double weight[n_bind_i_teth_candidates_];
  double weight_total{0.0};
  for (int i_entry{0}; i_entry < n_bind_i_teth_candidates_; i_entry++) {
    POP_T *head = std::get<POP_T *>(bind_i_teth_candidates_[i_entry]);
    double weight_xlink{head->xlink_->GetTotalWeight_Bind_I_Teth()};
    // If entry has a weight of 0.0, remove it from candidates
    if (weight_xlink == 0.0) {
      int i_last{--n_bind_i_teth_candidates_};
      bind_i_teth_candidates_[i_entry] = bind_i_teth_candidates_[i_last];
      i_entry--;
      // If we have zero valid candidates, set n_expected to 0 in SampleStats
      if (n_bind_i_teth_candidates_ == 0) {
        return n_to_set;
      }
      continue;
    }
    weight[i_entry] = weight_xlink;
    weight_total += weight[i_entry];
  }
  if (weight_total == 0.0) {
    wally_->ErrorExit(" AP_MGMT::Set_Bind_I_Teth_Candidates()");
  }
  int n_removed{0};
  if (n_to_set > n_bind_i_teth_candidates_) {
    n_removed = n_to_set - n_bind_i_teth_candidates_;
    n_to_set = n_bind_i_teth_candidates_;
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
  return n_removed;
}

int AssociatedProteinManagement::Set_Bind_II_Teth_Candidates(int n_to_set) {

  if (verbosity_ >= 1) {
    wally_->Log(" Starting AP_MGMT::Set_Bind_II_Teth_Candidates()\n");
  }
  double weight[n_bind_ii_teth_candidates_];
  double weight_total{0.0};
  for (int i_entry{0}; i_entry < n_bind_ii_teth_candidates_; i_entry++) {
    POP_T *head = std::get<POP_T *>(bind_ii_teth_candidates_[i_entry]);
    double weight_xlink{head->xlink_->GetTotalWeight_Bind_II_Teth()};
    // If entry has a weight of 0.0, remove it from candidates
    if (weight_xlink == 0.0) {
      int i_last{--n_bind_ii_teth_candidates_};
      bind_ii_teth_candidates_[i_entry] = bind_ii_teth_candidates_[i_last];
      i_entry--;
      // If we have zero valid candidates, set n_expected to 0 in SampleStats
      if (n_bind_ii_teth_candidates_ == 0) {
        return n_to_set;
      }
      continue;
    }
    weight[i_entry] = weight_xlink;
    weight_total += weight[i_entry];
  }
  if (weight_total == 0.0) {
    wally_->ErrorExit("AP_MGMT::Set_Bind_II_Teth_Candidates()\n");
  }
  int n_removed{0};
  if (n_to_set > n_bind_ii_teth_candidates_) {
    n_removed = n_to_set - n_bind_ii_teth_candidates_;
    n_to_set = n_bind_ii_teth_candidates_;
  }
  ENTRY_T selected_candidates[n_to_set];
  for (int i_set{0}; i_set < n_to_set; i_set++) {
    double p_cum{0.0};
    double ran{properties_->gsl.GetRanProb()};
    for (int i_entry{0}; i_entry < n_bind_ii_teth_candidates_; i_entry++) {
      p_cum += (weight[i_entry] / weight_total);
      if (ran < p_cum) {
        selected_candidates[i_set] = bind_ii_teth_candidates_[i_entry];
        weight_total -= weight[i_entry];
        weight[i_entry] = weight[n_bind_ii_teth_candidates_ - 1];
        n_bind_ii_teth_candidates_--;
        break;
      }
    }
  }
  n_bind_ii_teth_candidates_ = n_to_set;
  for (int i_entry{0}; i_entry < n_bind_ii_teth_candidates_; i_entry++) {
    bind_ii_teth_candidates_[i_entry] = selected_candidates[i_entry];
  }
  return n_removed;
}

void AssociatedProteinManagement::RunKMC() {

  // if (properties_->current_step_ > 4997810) {
  //   verbosity_ = 2;
  // }
  if (!population_active_) {
    return;
  }
  UpdateLists();
  SampleEventStatistics();
  GenerateExecutionSequence();
  ExecuteEvents();
}

void AssociatedProteinManagement::UpdateLists() {

  /*
  if (lists_up_to_date_) {
    return;
  }
  */
  if (verbosity_ >= 1) {
    wally_->Log("\nStarting AP_MGMT::UpdateLists()\n");
  }
  // lists_up_to_date_ = true;
  Update_Extensions();
  properties_->microtubules.UpdateUnoccupied();
  Update_Bound_I();
  if (crosslinking_active_) {
    Update_Bound_II();
  }
  if (tethering_active_) {
    Update_Free_Teth();
    Update_Bound_Unteth();
    Update_Bound_I_Teth();
    if (crosslinking_active_) {
      Update_Bound_II_Teth();
    }
    properties_->kinesin4.Update_Bound_Unteth();
  }
}

void AssociatedProteinManagement::SampleEventStatistics() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting AP_MGMT::SampleEventStatistics()\n");
  }
  n_events_to_exe_ = 0;
  for (auto &&event : events_) {
    if (verbosity_ >= 2) {
      // wally_->Log("Sampling statistics for %s\n", event.name_.c_str());
    }
    n_events_to_exe_ += event.SampleStatistics();
    // Ensure no funny business occurs with statistics
    if (event.n_expected_ > *event.n_avail_) {
      printf("Error; %i events expected but only %i available for %s\n",
             event.n_expected_, *event.n_avail_, event.name_.c_str());
      wally_->ErrorExit("AP_MGMT::SampleEventStatistics()");
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
        wally_->Log("    -target %i: xlink %i\n", i_entry, head->xlink_->id_);
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
  // xlink
  for (int i_entry{0}; i_entry < n_events_to_exe_; i_entry++) {
    EVENT_T *event_i = active_events[i_entry].first;
    AssociatedProtein *xlink_i{nullptr};
    try {
      xlink_i = std::get<POP_T *>(active_events[i_entry].second)->xlink_;
    } catch (...) {
      continue;
    }
    for (int j_entry{i_entry + 1}; j_entry < n_events_to_exe_; j_entry++) {
      EVENT_T *event_j = active_events[j_entry].first;
      AssociatedProtein *xlink_j{nullptr};
      try {
        xlink_j = std::get<POP_T *>(active_events[j_entry].second)->xlink_;
      } catch (...) {
        continue;
      }
      // If event_i and event_j target different xlinks, continue
      if (xlink_i != xlink_j) {
        continue;
      }
      double p_one{event_i->p_occur_};
      double p_two{event_j->p_occur_};
      double ran{properties_->gsl.GetRanProb()};
      if (ran < p_one / (p_one + p_two)) {
        if (verbosity_ >= 2) {
          wally_->Log("Removed xlink #%i from %s targets\n", xlink_i->id_,
                      event_i->name_.c_str());
        }
        event_i->RemoveTarget(active_events[i_entry].second);
        active_events[i_entry] = active_events[n_events_to_exe_ - 1];
        i_entry--;
        n_events_to_exe_--;
        break;
      } else {
        if (verbosity_ >= 2) {
          wally_->Log("Removed xlink #%i from %s targets\n", xlink_j->id_,
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

void AssociatedProteinManagement::GenerateExecutionSequence() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting AP_MGMT::GenerateExecutionSequence()\n");
  }
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
    wally_->ErrorExit("AP_MGMT::GenerateExecutionSequence()");
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

void AssociatedProteinManagement::ExecuteEvents() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting AP_MGMT::ExecuteEvents()\n");
  }
  for (int i_event{0}; i_event < n_events_to_exe_; i_event++) {
    events_to_exe_[i_event]->Execute();
    FlagForUpdate();
  }
}

void AssociatedProteinManagement::Diffuse(POP_T *head, int dir) {

  if (verbosity_ >= 1) {
    wally_->Log("Executing AP_MGMT::Diffuse()\n");
  }
  if (head->site_ == nullptr) {
    printf("woah buddy; u cant diffuse on a NULL site!!!\n");
    return;
  }
  Tubulin *old_site = head->site_;
  int dx{dir * head->GetDirectionTowardRest()};
  // Cannot step off them MTs
  if ((old_site->index_ == 0 and dx == -1) or
      (old_site->index_ == old_site->mt_->n_sites_ - 1 and dx == 1)) {
    return;
  }
  Tubulin *new_site = &old_site->mt_->lattice_[old_site->index_ + dx];
  if (!new_site->occupied_) {
    old_site->xlink_head_ = nullptr;
    old_site->occupied_ = false;
    new_site->xlink_head_ = head;
    new_site->occupied_ = true;
    head->site_ = new_site;
    properties_->microtubules.FlagForUpdate();
  }
}

void AssociatedProteinManagement::Bind_I(SITE_T *site) {

  if (verbosity_ >= 1) {
    wally_->Log("Executing AP_MGMT::Bind_I()\n");
  }
  AssociatedProtein *xlink{GetFreeXlink()};
  site->xlink_head_ = &xlink->head_one_;
  site->occupied_ = true;
  xlink->heads_active_++;
  xlink->head_one_.site_ = site;
  AddToActive(xlink);
  properties_->microtubules.FlagForUpdate();
}

void AssociatedProteinManagement::Bind_I_Teth(POP_T *satellite_head) {

  if (verbosity_ >= 1) {
    wally_->Log("Executing AP_MGMT::Bind_I_Teth()\n");
  }
  AssociatedProtein *xlink{satellite_head->xlink_};
  Tubulin *site{xlink->GetWeightedSite_Bind_I_Teth()};
  if (site == nullptr) {
    printf("failed toBind_I_Free_Tethered (XLINK)\n");
    return;
  }
  site->xlink_head_ = satellite_head;
  site->occupied_ = true;
  satellite_head->site_ = site;
  xlink->heads_active_++;
  properties_->microtubules.FlagForUpdate();
}

void AssociatedProteinManagement::Bind_II(POP_T *bound_head) {

  if (verbosity_ >= 1) {
    wally_->Log("Executing AP_MGMT::Bind_II()\n");
  }
  POP_T *unbound_head{bound_head->GetOtherHead()};
  AssociatedProtein *xlink{bound_head->xlink_};
  Tubulin *site{nullptr};
  bool has_satellite{false};
  if (xlink->tethered_) {
    if (xlink->motor_->heads_active_ == 0) {
      has_satellite = true;
    }
  }
  if (!xlink->tethered_ or has_satellite) {
    site = bound_head->xlink_->GetWeightedSite_Bind_II();
  } else {
    site = bound_head->xlink_->GetWeightedSite_Bind_II_Teth();
  }
  if (site == nullptr) {
    wally_->ErrorExit("AP_MGMT::Bind_II()");
  }
  site->xlink_head_ = unbound_head;
  site->occupied_ = true;
  unbound_head->site_ = site;
  xlink->heads_active_++;
  properties_->microtubules.FlagForUpdate();
}

void AssociatedProteinManagement::Unbind_II(POP_T *head) {

  if (verbosity_ >= 1) {
    wally_->Log("Executing AP_MGMT::Unbind_II()\n");
  }
  AssociatedProtein *xlink{head->xlink_};
  Tubulin *site{head->site_};
  site->xlink_head_ = nullptr;
  site->occupied_ = false;
  head->site_ = nullptr;
  xlink->heads_active_--;
  properties_->microtubules.FlagForUpdate();
}

void AssociatedProteinManagement::Unbind_I(POP_T *head) {

  if (verbosity_ >= 1) {
    wally_->Log("Executing AP_MGMT::Unbind_I()\n");
  }
  AssociatedProtein *xlink{head->xlink_};
  Tubulin *site{head->site_};
  site->xlink_head_ = nullptr;
  site->occupied_ = false;
  head->site_ = nullptr;
  xlink->heads_active_--;
  xlink->UntetherSatellite();
  if (!xlink->tethered_) {
    RemoveFromActive(xlink);
  }
  properties_->microtubules.FlagForUpdate();
}

void AssociatedProteinManagement::Tether_Free(ALT_T *untethered_head) {

  if (verbosity_ >= 1) {
    wally_->Log("Executing AP_MGMT::Tether_Free()\n");
  }
  AssociatedProtein *xlink{GetFreeXlink()};
  Kinesin *motor{untethered_head->motor_};
  xlink->motor_ = motor;
  xlink->tethered_ = true;
  motor->xlink_ = xlink;
  motor->tethered_ = true;
  AddToActive(xlink);
  properties_->kinesin4.FlagForUpdate();
}

void AssociatedProteinManagement::Untether_Free(POP_T *satellite_head) {

  if (verbosity_ >= 1) {
    wally_->Log("Executing AP_MGMT::Untether_Free()\n");
  }
  AssociatedProtein *xlink{satellite_head->xlink_};
  Kinesin *motor{xlink->motor_};
  xlink->motor_ = nullptr;
  xlink->tethered_ = false;
  motor->xlink_ = nullptr;
  motor->tethered_ = false;
  RemoveFromActive(xlink);
  properties_->kinesin4.FlagForUpdate();
}