#include "associated_protein_management.h"
#include "master_header.h"

AssociatedProteinManagement::AssociatedProteinManagement() {}

void AssociatedProteinManagement::Initialize(system_parameters *parameters,
                                             system_properties *properties) {

  verbose_ = false;
  wally_ = &properties->wallace;
  parameters_ = parameters;
  properties_ = properties;
  CalculateCutoffs();
  SetParameters();
  GenerateXLinks();
  InitializeLists();
  InitializeEvents();
}

void AssociatedProteinManagement::CalculateCutoffs() {

  double site_size{parameters_->microtubules.site_size};
  double k_spring{parameters_->xlinks.k_spring};
  double r_y{parameters_->microtubules.y_dist};
  double r_0{parameters_->xlinks.r_0};
  double kbT{parameters_->kbT};
  /* First, calculate rest_dist_ in number of sites */
  // For now, it should always come out to zero
  int approx_rest{(int)(sqrt(r_0 * r_0 - r_y * r_y) / site_size)};
  if (approx_rest == 0) {
    rest_dist_ = approx_rest;
  } else {
    wally_->ErrorExit("AssociatedProteinManagement::CalculateCutoffs()");
  }
  /* next, calculate extension distance cutoff */
  for (int x_dist{rest_dist_}; x_dist < 1000; x_dist++) {
    int r_x = x_dist * site_size;
    double r = sqrt(r_y * r_y + r_x * r_x);
    double dr = r - r_0;
    double U = (k_spring / 2) * dr * dr;
    double boltzmann_weight = exp(U / (2 * kbT));
    if (boltzmann_weight > 1000) {
      dist_cutoff_ = x_dist;
      break;
    }
  }
  if (dist_cutoff_ > 9) {
    printf("As of now, x_dist must be less than 10");
    printf(" (currenly %i with k_spring=%g)\n", dist_cutoff_, k_spring);
    wally_->ErrorExit("AssociatedProteinManagement::CalculateCutoffs() [2]");
  }
  if (parameters_->microtubules.count > 1) {
    wally_->Log("\nFor crosslinkers:\n");
    wally_->Log("  rest_dist is %i\n", rest_dist_);
    wally_->Log("  dist_cutoff is %i\n\n", dist_cutoff_);
  } else {
    wally_->Log("\nCrosslinker double-binding is inactive.\n\n");
  }
}

void AssociatedProteinManagement::SetParameters() {

  /*
    For events that result in a change in energy, we use Boltzmann factors to
    scale rates appropriately. Detailed balance is satisfied with the factors:
                  exp{-(1 - lambda)*(delta_E)/(kB*T)}, and
                  exp{-(lambda)*(delta_E)/(kB*T)}
    for an event and its complement (e.g., binding and unbinding), where lambda
    is a constant that ranges from 0 to 1, delta_E is the change in energy that
    results from the event, kB is Boltzmann's constant, and T is the temperature
  */
  // Lambda = 1.0 means all of the weight goes into unbinding
  double lambda_neighb{1.0}; // Lambda for neighbor interaction energies
  // Lambda = 0.0 means all of the weight goes into binding
  double lambda_spring{0.0}; // Lambda for spring energies
  // Lambda = 0.5 means the weight is equally split between binding & unbinding
  double lambda_teth{0.5}; // lambda for tether spring energies
  // Calculate neighbor interaction energies
  max_neighbs_ = 2;
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
  }
  // Calculate tether extension energies
  double k_teth{parameters_->motors.k_spring};
  double k_slack{parameters_->motors.k_slack};
  double r_0_teth{parameters_->motors.r_0};
  double r_y_teth{parameters_->microtubules.y_dist / 2};
  double rest_dist_teth{properties_->kinesin4.motors_[0].rest_dist_};
  teth_cutoff_ = properties_->kinesin4.motors_[0].dist_cutoff_;
  comp_cutoff_ = properties_->kinesin4.motors_[0].comp_cutoff_;
  std::vector<double> teth_energy(2 * teth_cutoff_ + 1, 0.0);
  for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
    double r_x_teth{((double)x_dub / 2) * site_size};
    double r_teth{sqrt(r_x_teth * r_x_teth + r_y_teth * r_y_teth)};
    double dr_teth{r_teth - r_0_teth};
    // Tether spring energies are positive by definition
    if (dr_teth > 0) {
      // If tether is extended past rest dist, use k_teth as spring constant
      teth_energy[x_dub] = 0.5 * k_teth * dr_teth * dr_teth;
    } else {
      // Otherwise, if tether is compressed past rest dist, use k_slack
      teth_energy[x_dub] = 0.5 * k_slack * dr_teth * dr_teth;
    }
  }
  // Vector that holds the name of each probability and its value
  std::vector<std::pair<std::string, double>> probabilities;
  std::string NEIGHBS; // std::to_string(n_neighbs) for convenience
  std::string X_DIST;  // std::to_string(x) for convenience
  std::string X_DUB;   // std::to_string(x_dub) for convenience
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
  // We consider diff_fwd and diff_bck as two separate events, which effectively
  // doubles the probability to diffuse. To counteract this, divide p_diff by 2.
  p_diffu_i /= 2;
  p_diffu_ii /= 2;
  p_diffuse_i_fwd_.resize(max_neighbs_ + 1);
  p_diffuse_i_bck_.resize(max_neighbs_ + 1);
  p_diffuse_ii_to_rest_.resize(max_neighbs_ + 1);
  p_diffuse_ii_fr_rest_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    NEIGHBS = std::to_string(n_neighbs);
    // For neighb. interactions, we only consider the energy penalty for
    // diffusing away from neighbors (same as unbinding) since lambda = 1.0
    // dE = E_f - E_i
    double dE = 0.0 - neighb_energy[n_neighbs];
    // dE has units of kbT, so dividing by its numerical value isn't necessary
    double weight_neighb{exp(-lambda_neighb * dE)};
    // With two neighbors, the binding head is jammed and cannot diffuse
    if (n_neighbs == 2) {
      weight_neighb = 0.0;
    }
    p_diffuse_i_fwd_[n_neighbs] = weight_neighb * p_diffu_i;
    p_diffuse_i_bck_[n_neighbs] = weight_neighb * p_diffu_i;
    probabilities.emplace_back("diff_i_" + NEIGHBS,
                               p_diffuse_i_fwd_[n_neighbs]);
    p_diffuse_ii_to_rest_[n_neighbs].resize(dist_cutoff_ + 1);
    p_diffuse_ii_fr_rest_[n_neighbs].resize(dist_cutoff_ + 1);
    for (int x{0}; x <= dist_cutoff_; x++) {
      X_DIST = std::to_string(x);
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
      // Diffusing towards rest is considered an unbinding-type event in regards
      // to Boltzmann factors, since both events let the spring relax
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
      probabilities.emplace_back("diff_ii_to_" + NEIGHBS + "_" + X_DIST,
                                 p_diffuse_ii_to_rest_[n_neighbs][x]);
      p_diffuse_ii_fr_rest_[n_neighbs][x] =
          weight_neighb * weight_fr * p_diffu_ii;
      probabilities.emplace_back("diff_ii_fr_" + NEIGHBS + "_" + X_DIST,
                                 p_diffuse_ii_fr_rest_[n_neighbs][x]);
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
    NEIGHBS = std::to_string(n_neighbs);
    p_diffuse_i_to_teth_rest_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    p_diffuse_i_fr_teth_rest_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    p_diffuse_ii_to_both_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    p_diffuse_ii_to_self_fr_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    p_diffuse_ii_fr_self_to_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    p_diffuse_ii_fr_both_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      X_DUB = std::to_string(x_dub);
      // Change in x_dub that brings teth towards rest; -1 if extended (default)
      int dx_rest{-1};
      // However, if tether is compressed, increased x_dub is towards rest
      if (x_dub < 2 * rest_dist_teth) {
        dx_rest = 1;
      }
      // dU (U_f - U_i) for tether if head steps towards (to) tether rest
      double dU_to{teth_energy[x_dub + dx_rest] - teth_energy[x_dub]};
      // dU (U_f - U_i) for tether if head steps away from (fr) tether rest
      double dU_fr{teth_energy[x_dub - dx_rest] - teth_energy[x_dub]};
      // Diffusing towards rest is considered an unbinding-type event in regards
      // to Boltzmann factors, since both events let the spring relax
      double weight_to_teth{exp(-lambda_teth * dU_to / kbT)};
      // Diffusing away from rest is considered a binding-type event in regards
      // to Boltzmann factors, since both events stretch the spring out
      double weight_fr_teth{exp(-(1.0 - lambda_teth) * dU_fr / kbT)};
      if (x_dub == 2 * comp_cutoff_ or x_dub == 2 * teth_cutoff_) {
        weight_fr_teth = 0.0;
      }
      p_diffuse_i_to_teth_rest_[n_neighbs][x_dub] =
          p_diffuse_i_fwd_[n_neighbs] * weight_to_teth;
      probabilities.emplace_back("diff_i_to_teth_" + NEIGHBS + "_" + X_DUB,
                                 p_diffuse_i_to_teth_rest_[n_neighbs][x_dub]);
      // For singly-bound xlinks, diffusing one site changes x_dub by 2, so keep
      // probability 0.0 for 2*cutoff +/- 1
      if (x_dub > 2 * comp_cutoff_ + 1 and x_dub < 2 * teth_cutoff_ - 1) {
        p_diffuse_i_fr_teth_rest_[n_neighbs][x_dub] =
            p_diffuse_i_bck_[n_neighbs] * weight_fr_teth;
        probabilities.emplace_back("diff_i_fr_teth_" + NEIGHBS + "_" + X_DUB,
                                   p_diffuse_i_fr_teth_rest_[n_neighbs][x_dub]);
      }
      p_diffuse_ii_to_both_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      p_diffuse_ii_fr_both_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      p_diffuse_ii_to_self_fr_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      p_diffuse_ii_fr_self_to_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      for (int x{0}; x <= dist_cutoff_; x++) {
        X_DIST = std::to_string(x);
        // Diffuse towards both own and tether rest
        p_diffuse_ii_to_both_[n_neighbs][x_dub][x] =
            p_diffuse_ii_to_rest_[n_neighbs][x] * weight_to_teth;
        probabilities.emplace_back("diff_ii_to_both_" + NEIGHBS + "_" + X_DUB +
                                       "_" + X_DIST,
                                   p_diffuse_ii_to_both_[n_neighbs][x_dub][x]);
        // Diffuse away from both own and tether rest
        p_diffuse_ii_fr_both_[n_neighbs][x_dub][x] =
            p_diffuse_ii_fr_rest_[n_neighbs][x] * weight_fr_teth;
        probabilities.emplace_back("diff_ii_fr_both_" + NEIGHBS + "_" + X_DUB +
                                       "_" + X_DIST,
                                   p_diffuse_ii_fr_both_[n_neighbs][x_dub][x]);
        // Diffuse towards own rest, but away from tether rest
        p_diffuse_ii_to_self_fr_teth_[n_neighbs][x_dub][x] =
            p_diffuse_ii_to_rest_[n_neighbs][x] * weight_fr_teth;
        probabilities.emplace_back(
            "diff_ii_to_self_fr_teth_" + NEIGHBS + "_" + X_DUB + "_" + X_DIST,
            p_diffuse_ii_to_self_fr_teth_[n_neighbs][x_dub][x]);
        // Diffuse away from own rest, but towards tether rest
        p_diffuse_ii_fr_self_to_teth_[n_neighbs][x_dub][x] =
            p_diffuse_ii_fr_rest_[n_neighbs][x] * weight_to_teth;
        probabilities.emplace_back(
            "diff_ii_fr_self_to_teth_" + NEIGHBS + "_" + X_DUB + "_" + X_DIST,
            p_diffuse_ii_fr_self_to_teth_[n_neighbs][x_dub][x]);
      }
    }
  }
  // KMC STATISTICS BELOW
  double k_on = parameters_->xlinks.k_on;
  double c_xlink = parameters_->xlinks.c_bulk;
  double c_eff_teth = parameters_->motors.c_eff_tether;
  if (!parameters_->motors.tethers_active) {
    c_eff_teth = 0;
  }
  p_bind_i_teth_base_ = k_on * c_eff_teth * delta_t;
  double c_eff_bind = parameters_->xlinks.c_eff_bind;
  p_bind_ii_base_ = k_on * c_eff_bind * delta_t;
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
    NEIGHBS = std::to_string(n_neighbs);
    // dE = E_f - E_i;
    double dE = 0.0 - neighb_energy[n_neighbs];
    // We only consider the energy penalty for unbinding since lambda = 1.0
    double weight_neighb_bind{1.0};
    // dE has units of kbT, so dividing by its numerical value isn't necessary
    double weight_neighb_unbind{exp(-lambda_neighb * dE)};
    p_bind_i_[n_neighbs] = k_on * c_xlink * delta_t;
    probabilities.emplace_back("bind_i_" + NEIGHBS, p_bind_i_[n_neighbs]);
    p_unbind_i_[n_neighbs] = weight_neighb_unbind * k_off_i * delta_t;
    probabilities.emplace_back("unbind_i_" + NEIGHBS, p_unbind_i_[n_neighbs]);
    p_unbind_ii_[n_neighbs].resize(dist_cutoff_ + 1);
    weight_bind_ii_[n_neighbs].resize(dist_cutoff_ + 1);
    for (int x{0}; x <= dist_cutoff_; x++) {
      X_DIST = std::to_string(x);
      // dU = U_f - U_i
      double dU_bind{spring_energy[x] - 0.0};
      double dU_unbind{0.0 - spring_energy[x]};
      double weight_bind{exp(-(1.0 - lambda_spring) * dU_bind / kbT)};
      weight_bind_ii_[n_neighbs][x] = weight_neighb_bind * weight_bind;
      double weight_unbind{exp(-lambda_spring * dU_unbind / kbT)};
      p_unbind_ii_[n_neighbs][x] =
          weight_neighb_unbind * weight_unbind * k_off_ii * delta_t;
      probabilities.emplace_back("unbind_ii_" + NEIGHBS + "_" + X_DIST,
                                 p_unbind_ii_[n_neighbs][x]);
    }
    // Rates involving both xlink & tether spring
    p_unbind_i_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    p_unbind_ii_to_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    p_unbind_ii_fr_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    weight_bind_i_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    weight_bind_ii_to_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    weight_bind_ii_fr_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      X_DUB = std::to_string(x_dub);
      // dU = U_f - U_i
      double dU_bind{teth_energy[x_dub] - 0.0};
      double dU_unbind{0.0 - teth_energy[x_dub]};
      double weight_bind{exp(-(1.0 - lambda_teth) * dU_bind / kbT)};
      weight_bind_i_teth_[n_neighbs][x_dub] = weight_neighb_bind * weight_bind;
      double weight_unbind{exp(-lambda_teth * dU_unbind / kbT)};
      p_unbind_i_teth_[n_neighbs][x_dub] =
          weight_unbind * p_unbind_i_[n_neighbs];
      probabilities.emplace_back("unbind_i_teth_" + NEIGHBS + "_" + X_DUB,
                                 p_unbind_i_teth_[n_neighbs][x_dub]);
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
  int align_pos{40};
  // std::sort(probabilities.begin(), probabilities.end(),
  //           [](const std::pair<std::string, double> &a,
  //              const std::pair<std::string, double> &b) {
  //             return a.first.size() < b.first.size();
  //           });

  // for (int i_prob{0}; i_prob < probabilities.size(); i_prob++) {
  //   char probability_report[256];
  //   int whitespace{align_pos - (int)probabilities[i_prob].first.length() +
  //   2}; sprintf(probability_report, "XLINKS: p_%s = %g\n",
  //           probabilities[i_prob].first.c_str(), // whitespace,
  //           probabilities[i_prob].second);
  //   if (probabilities[i_prob].second > 0.0) {
  //     wally_->Log(probability_report);
  //   }
  // }
}

void AssociatedProteinManagement::GenerateXLinks() {

  // Since only one head has to be bound, the sim will at most
  // as many xlinks as sites in the bulk (all single-bound)
  int n_mts = parameters_->microtubules.count;
  for (int i_mt{0}; i_mt < n_mts; i_mt++) {
    n_xlinks_ += parameters_->microtubules.length[i_mt];
  }
  xlinks_.resize(n_xlinks_);
  for (int idx{0}; idx < n_xlinks_; idx++) {
    xlinks_[idx].Initialize(parameters_, properties_, idx);
  }
}

void AssociatedProteinManagement::InitializeLists() {

  scratch_.resize(n_xlinks_);
  // Stats (not a list ok bite me)
  n_bound_i_.resize(max_neighbs_ + 2);
  n_bound_ii_.resize(max_neighbs_ + 2);
  n_bound_i_teth_.resize(max_neighbs_ + 2);
  n_bound_ii_teth_same_.resize(max_neighbs_ + 2);
  n_bound_ii_teth_oppo_.resize(max_neighbs_ + 2);
  // Final entry, [max_neighbs_+1], holds ALL stats regardless of neighbs
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_ + 1; n_neighbs++) {
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
  free_teth_.resize(n_xlinks_);
  bound_unteth_.resize(n_xlinks_);
  bound_i_.resize(max_neighbs_ + 2);
  bound_ii_.resize(max_neighbs_ + 2);
  bound_i_teth_.resize(max_neighbs_ + 2);
  bound_ii_teth_oppo_.resize(max_neighbs_ + 2);
  bound_ii_teth_same_.resize(max_neighbs_ + 2);
  // Final entry, [max_neighbs_+1], holds ALL stats regardless of neighbs
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_ + 1; n_neighbs++) {
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

  /* *** Serialized & unique index of each KMC event *** */
  int ID(0);
  /* *** Basic random integer generator ** */
  auto ran_int = [&](int n) {
    if (n > 0)
      return properties_->gsl.GetRanInt(n);
    else
      return 0;
  };
  /* *** Probability distributions *** */
  // Baseline probability dist. is Binomial since sim is discretized
  auto binomial = [&](double p, int n) {
    if (n > 0)
      return properties_->gsl.SampleBinomialDist(p, n);
    else
      return 0;
  };
  // Poisson dist. is used w/ partition function for E-dependent binding
  auto poisson_ii = [&](double p, int n) {
    if (n > 0) {
      double n_wt = GetWeight_Bind_II();
      if (n_wt > 0) {
        int n_predict = properties_->gsl.SamplePoissonDist(p * n_wt);
        // printf("rolled to bind_ii #%i entries\n", n_predict);
        if (n_predict > n) {
          printf("whoaaa buddy this poisson dist rowdy af\n");
          return n;
        } else
          return n_predict;
      } else
        return 0;
    } else
      return 0;
  };
  auto poisson_i_teth = [&](double p, int n) {
    if (n > 0) {
      double n_wt = GetWeight_Bind_I_Teth();
      if (n_wt > 0)
        return properties_->gsl.SamplePoissonDist(p * n_wt);
      else
        return 0;
    } else
      return 0;
  };
  auto poisson_ii_teth = [&](double p, int n) {
    if (n > 0) {
      double n_wt = GetWeight_Bind_II_Teth();
      if (n_wt > 0)
        return properties_->gsl.SamplePoissonDist(p * n_wt);
      else
        return 0;
    } else
      return 0;
  };
  /* ***    Actual instantiation of event objects; format is:    *** */
  //		(*management, ID, kmc_code, 'event_name', 'target_pop',
  //		  p_occur, *n_avail, *pop_pool, ran_int, p_dist)
  // Bind_ii does not use n_neighbors; append w/ "_ALL" instead
  if (parameters_->microtubules.count > 1) {
    events_.emplace_back(this, ID++, 20, "bind_ii", "bound_i_ALL",
                         p_bind_ii_base_, &n_bound_i_[max_neighbs_ + 1],
                         &bound_i_[max_neighbs_ + 1], ran_int, poisson_ii);
  }
  for (int n_neighbs(0); n_neighbs <= max_neighbs_; n_neighbs++) {
    // Number of neighbors, N, always appends target pop. type
    std::string N = std::to_string(n_neighbs);
    events_.emplace_back(this, ID++, -10, "diff_i_fwd", "bound_i_" + N,
                         p_diffuse_i_fwd_[n_neighbs], &n_bound_i_[n_neighbs],
                         &bound_i_[n_neighbs], ran_int, binomial);
    events_.emplace_back(this, ID++, -11, "diff_i_bck", "bound_i_" + N,
                         p_diffuse_i_bck_[n_neighbs], &n_bound_i_[n_neighbs],
                         &bound_i_[n_neighbs], ran_int, binomial);
    events_.emplace_back(
        this, ID++, 10, "bind_i", "unocc_" + N, p_bind_i_[n_neighbs],
        &properties_->microtubules.n_unoccupied_xl_[n_neighbs],
        &properties_->microtubules.unoccupied_list_xl_[n_neighbs], ran_int,
        binomial);
    events_.emplace_back(this, ID++, 30, "unbind_i", "bound_i_" + N,
                         p_unbind_i_[n_neighbs], &n_bound_i_[n_neighbs],
                         &bound_i_[n_neighbs], ran_int, binomial);
    // Only create doubly-bound events if n_MTs > 1
    if (parameters_->microtubules.count > 1) {
      for (int x(0); x <= dist_cutoff_; x++) {
        events_.emplace_back(this, ID++, -20, "diff_ii_to",
                             "bound_ii_" + std::to_string(x) + "_" + N,
                             p_diffuse_ii_to_rest_[n_neighbs][x],
                             &n_bound_ii_[n_neighbs][x],
                             &bound_ii_[n_neighbs][x], ran_int, binomial);
        events_.emplace_back(this, ID++, -21, "diff_ii_fr",
                             "bound_ii_" + std::to_string(x) + "_" + N,
                             p_diffuse_ii_fr_rest_[n_neighbs][x],
                             &n_bound_ii_[n_neighbs][x],
                             &bound_ii_[n_neighbs][x], ran_int, binomial);
        events_.emplace_back(this, ID++, 40, "*unbind_ii",
                             "bound_ii_" + std::to_string(x) + "_" + N,
                             p_unbind_ii_[n_neighbs][x],
                             &n_bound_ii_[n_neighbs][x],
                             &bound_ii_[n_neighbs][x], ran_int, binomial);
      }
    }
  }
  // If tethering is enabled, add those event populations as well
  if (parameters_->motors.c_bulk > 0.0 and parameters_->motors.tethers_active) {
    events_.emplace_back(this, ID++, 11, "bind_I", "free_teth",
                         p_bind_i_teth_base_, &n_free_teth_, &free_teth_,
                         ran_int, poisson_i_teth);
    // Only create doubly-bound events if n_MTs > 1
    if (parameters_->microtubules.count > 1) {
      events_.emplace_back(this, ID++, 21, "bind_II", "bound_I_ALL",
                           p_bind_ii_base_, &n_bound_i_teth_tot_,
                           &bound_i_teth_[max_neighbs_ + 1][2 * teth_cutoff_],
                           ran_int, poisson_ii_teth);
    }
    for (int n_neighbs(0); n_neighbs <= max_neighbs_; n_neighbs++) {
      //			printf("n neighbs is %i\n", n_neighbs);
      for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
        std::string bound_I = "bound_I_" + std::to_string(x_dub) + "_" +
                              std::to_string(n_neighbs);
        //				std::cout << bound_I << std::endl;
        events_.emplace_back(this, ID++, -30, "diff_I_to", bound_I,
                             p_diffuse_i_to_teth_rest_[n_neighbs][x_dub],
                             &n_bound_i_teth_[n_neighbs][x_dub],
                             &bound_i_teth_[n_neighbs][x_dub], ran_int,
                             binomial);
        events_.emplace_back(this, ID++, -31, "diff_I_fr", bound_I,
                             p_diffuse_i_fr_teth_rest_[n_neighbs][x_dub],
                             &n_bound_i_teth_[n_neighbs][x_dub],
                             &bound_i_teth_[n_neighbs][x_dub], ran_int,
                             binomial);
        events_.emplace_back(this, ID++, 31, "unbind_I", bound_I,
                             p_unbind_i_teth_[n_neighbs][x_dub],
                             &n_bound_i_teth_[n_neighbs][x_dub],
                             &bound_i_teth_[n_neighbs][x_dub], ran_int,
                             binomial);
        // Only create doubly-bound events if n_MTs > 1
        if (parameters_->microtubules.count > 1) {
          for (int x(0); x <= dist_cutoff_; x++) {
            std::string bound_II_same = "bound_II_same_" + std::to_string(x) +
                                        "_" + std::to_string(x_dub) + "_" +
                                        std::to_string(n_neighbs);
            std::string bound_II_oppo = "bound_II_oppo_" + std::to_string(x) +
                                        "_" + std::to_string(x_dub) + "_" +
                                        std::to_string(n_neighbs);
            events_.emplace_back(this, ID++, -40, "diff_II_to", bound_II_same,
                                 p_diffuse_ii_to_both_[n_neighbs][x_dub][x],
                                 &n_bound_ii_teth_same_[n_neighbs][x_dub][x],
                                 &bound_ii_teth_same_[n_neighbs][x_dub][x],
                                 ran_int, binomial);
            events_.emplace_back(this, ID++, -41, "diff_II_fr", bound_II_same,
                                 p_diffuse_ii_fr_both_[n_neighbs][x_dub][x],
                                 &n_bound_ii_teth_same_[n_neighbs][x_dub][x],
                                 &bound_ii_teth_same_[n_neighbs][x_dub][x],
                                 ran_int, binomial);
            events_.emplace_back(
                this, ID++, -50, "diff_II_to_fr", bound_II_oppo,
                p_diffuse_ii_to_self_fr_teth_[n_neighbs][x_dub][x],
                &n_bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                &bound_ii_teth_oppo_[n_neighbs][x_dub][x], ran_int, binomial);
            events_.emplace_back(
                this, ID++, -51, "diff_II_fr_to", bound_II_oppo,
                p_diffuse_ii_fr_self_to_teth_[n_neighbs][x_dub][x],
                &n_bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                &bound_ii_teth_oppo_[n_neighbs][x_dub][x], ran_int, binomial);
            events_.emplace_back(this, ID++, 41, "unbind_II_to", bound_II_same,
                                 p_unbind_ii_to_teth_[n_neighbs][x_dub][x],
                                 &n_bound_ii_teth_same_[n_neighbs][x_dub][x],
                                 &bound_ii_teth_same_[n_neighbs][x_dub][x],
                                 ran_int, binomial);
            events_.emplace_back(this, ID++, 42, "unbind_II_fr", bound_II_oppo,
                                 p_unbind_ii_fr_teth_[n_neighbs][x_dub][x],
                                 &n_bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                                 &bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                                 ran_int, binomial);
          }
        }
      }
    }
    events_.emplace_back(
        this, ID++, 50, "tether_free", "unteth_mots", p_tether_free_,
        &properties_->kinesin4.n_bound_untethered_,
        &properties_->kinesin4.bound_untethered_, ran_int, binomial);
    events_.emplace_back(this, ID++, 60, "untether_free", "free_teth",
                         p_untether_free_, &n_free_teth_, &free_teth_, ran_int,
                         binomial);
  }
  /* ** Segregate events_ into IDs_by_pop_ based on target pop. ** */
  int n_pops = 0;                // Total no. of distinct pops.
  int n_entries[events_.size()]; // No. of entries for each pop.
  // Index of each entry for each population:
  int entries[events_.size()][2 * (2 * teth_cutoff_ + 1) * (dist_cutoff_ + 1)];
  // Certain events have their stats corrected SECONDARY to others,
  // e.g., events that affect all extensions will be corrected after
  // events that affect a specific extension; 'root' refers to the
  // root pop., e.g., bound_i_x w/o neighb info (format can vary)
  int n_roots = 0;
  std::string roots[events_.size()];
  std::string bound_ii("bound_ii_");
  for (int i_entry(0); i_entry < events_.size(); i_entry++) {
    std::string target_pop = events_[i_entry].target_pop_;
    // bound_ii entries go thru a primary AND secondary correction
    // (first over x & n_neighbs, then over x for all n_neighbs)
    if (target_pop.length() >= bound_ii.length()) {
      if (target_pop.substr(0, bound_ii.length()) == bound_ii) {
        int root_length = bound_ii.length() + 1;
        std::string root = target_pop.substr(0, root_length);
        // Check to see if we have recorded this root before
        bool new_root(true);
        for (int i_root(0); i_root < n_roots; i_root++) {
          std::string recorded = roots[i_root];
          if (root == recorded)
            new_root = false;
        }
        if (new_root) {
          roots[n_roots] = root;
          n_roots++;
        }
      }
    }
    // 'ALL' flag at end -> secondary correction only
    if (target_pop.substr(target_pop.length() - 3) == "ALL") {
      std::string root = target_pop.substr(0, target_pop.length() - 3);
      roots[n_roots] = root;
      n_roots++;
    }
    // No 'ALL' flag -> primary scan only
    else {
      bool new_pop(true); // Assume population type is new
      int pop_index(0);   // 1st index of this pop. in entries arry
      for (int i_pop(0); i_pop < n_pops; i_pop++) {
        std::string record = events_[entries[i_pop][0]].target_pop_;
        // If type matches a recorded type; pop. isn't new
        if (target_pop == record) {
          new_pop = false;
          pop_index = i_pop;
        }
      }
      // If indeed a new pop., record it in a new row
      if (new_pop) {
        entries[n_pops][0] = events_[i_entry].ID_;
        n_entries[n_pops] = 1;
        n_pops++;
      }
      // Otherwise, add entry to row of already-recorded pop.
      else {
        int entry_no = n_entries[pop_index];
        entries[pop_index][entry_no] = events_[i_entry].ID_;
        n_entries[pop_index]++;
      }
    }
  }
  // Scan through primary pops. & place them into IDs_by_pop_
  // As we transfer primary pops, tally up all secondary events
  int n_sec_entries[n_roots];
  for (int i_root(0); i_root < n_roots; i_root++)
    n_sec_entries[i_root] = 0;
  int sec_entries[n_roots][2 * (2 * teth_cutoff_ + 1) * (dist_cutoff_ + 1)];
  // 1st index of IDs_by_pop_ corresponds to population type
  IDs_by_pop_.resize(n_pops);
  for (int i_pop(0); i_pop < n_pops; i_pop++) {
    std::string target_pop = events_[entries[i_pop][0]].target_pop_;
    // 2nd index corresponds to entry of events that target this pop.
    IDs_by_pop_[i_pop].resize(n_entries[i_pop]);
    for (int i_entry(0); i_entry < n_entries[i_pop]; i_entry++) {
      IDs_by_pop_[i_pop][i_entry] = entries[i_pop][i_entry];
    }
    // Check to see if this population matches any secondary-scan roots
    for (int i_root(0); i_root < n_roots; i_root++) {
      std::string target_root;
      std::string recorded_root = roots[i_root];
      if (target_pop.length() >= recorded_root.length())
        target_root = target_pop.substr(0, recorded_root.length());
      else
        target_root = std::string("nope");
      // If it does match, add all entries of this pop. to sec_entries
      if (target_root == recorded_root) {
        for (int i_entry(0); i_entry < n_entries[i_pop]; i_entry++) {
          sec_entries[i_root][n_sec_entries[i_root]] = entries[i_pop][i_entry];
          n_sec_entries[i_root]++;
        }
      }
    }
  }
  // Finally, scan through secondary pops. & place them into IDs_by_root_
  IDs_by_root_.resize(n_roots);
  n_avail_by_root_.resize(n_roots);
  for (int i_root(0); i_root < n_roots; i_root++) {
    std::string root = roots[i_root];
    if (root.substr(0, 8) == "bound_i_") {
      n_avail_by_root_[i_root] = &n_bound_i_[max_neighbs_ + 1];
    } else if (root.substr(0, 9) == "bound_ii_") {
      int x = std::stoi(root.substr(9, 1));
      n_avail_by_root_[i_root] = &n_bound_ii_[max_neighbs_ + 1][x];
    } else if (root.substr(0, 8) == "bound_I_") {
      n_avail_by_root_[i_root] = &n_bound_i_teth_tot_;
    } else {
      printf("???? in ROOT portion of GEN XLINK EVENTS");
    }
    IDs_by_root_[i_root].resize(n_sec_entries[i_root]);
    for (int i_entry(0); i_entry < n_sec_entries[i_root]; i_entry++) {
      int entry_ID = sec_entries[i_root][i_entry];
      // If target event has prefix of "*", make IDs negative to
      // flag that entries are coupled & n_avail should be halved
      bool coupled(false);
      if (events_[entry_ID].name_.substr(0, 1) == "*") {
        coupled = true;
      }
      if (coupled)
        IDs_by_root_[i_root][i_entry] = -1 * entry_ID;
      else
        IDs_by_root_[i_root][i_entry] = entry_ID;
    }
  }
  if (verbose_) {
    for (int i_pop{0}; i_pop < IDs_by_pop_.size(); i_pop++) {
      printf("Pop partition #%i: ", i_pop);
      int n_events = IDs_by_pop_[i_pop].size();
      for (int i_event{0}; i_event < n_events; i_event++) {
        int event_ID = IDs_by_pop_[i_pop][i_event];
        std::cout << events_[event_ID].name_;
        std::cout << " (targets " << events_[event_ID].target_pop_ << "), ";
      }
      printf("\n");
    }
    for (int i_root{0}; i_root < IDs_by_root_.size(); i_root++) {
      printf("Root partition #%i: ", i_root);
      int n_events = IDs_by_root_[i_root].size();
      for (int i_event{0}; i_event < n_events; i_event++) {
        if (i_event != 0) {
          printf("                   ");
        }
        int event_ID = IDs_by_root_[i_root][i_event];
        if (event_ID < 0) {
          event_ID = abs(event_ID);
          printf("[COUPLED] ");
        }
        std::cout << events_[event_ID].name_;
        std::cout << " (targets " << events_[event_ID].target_pop_ << "), ";
        printf("\n");
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

AssociatedProtein *AssociatedProteinManagement::GetBoundUntetheredXlink() {

  Update_Bound_Unteth();
  if (n_bound_unteth_ > 0) {
    int i_entry = properties_->gsl.GetRanInt(n_bound_unteth_);
    AssociatedProtein *xlink = bound_unteth_[i_entry];
    return xlink;
  } else {
    printf("Error in GetUnTetheredXlink: no untethered xlinks!\n");
    exit(1);
  }
}

double AssociatedProteinManagement::GetWeight_Bind_II() {

  double weights_summed = 0;
  // Sum over all single-bound xlinks
  for (int i_xlink = 0; i_xlink < n_bound_i_[max_neighbs_ + 1]; i_xlink++) {
    POP_T *head = std::get<POP_T *>(bound_i_[max_neighbs_ + 1][i_xlink]);
    head->xlink_->UpdateNeighborSites_II();
    // Get weight of every possible orientation w/ neighbors
    int n_neighbors = head->xlink_->n_neighbor_sites_;
    // printf("%i neighbs for xlink #%i\n", n_neighbors, head->xlink_->ID_);
    for (int i_neighb = 0; i_neighb < n_neighbors; i_neighb++) {
      Tubulin *site = head->xlink_->neighbor_sites_[i_neighb];
      double weight = head->xlink_->GetBindingWeight_II(site);
      // printf("returned wt is %g\n", weight);
      weights_summed += weight;
    }
  }
  // printf("total weight is %g\n", weights_summed);
  return weights_summed;
}

double AssociatedProteinManagement::GetWeight_Bind_I_Teth() {

  double weights_summed = 0;
  // Sum over all free_tethered xlinks
  for (int i_xlink = 0; i_xlink < n_free_teth_; i_xlink++) {
    POP_T *head = std::get<POP_T *>(free_teth_[i_xlink]);
    head->xlink_->UpdateNeighborSites_I_Teth();
    // Get weight of every possible orientation with neighbors
    int n_neighbs = head->xlink_->n_teth_neighbor_sites_;
    for (int i_neighb = 0; i_neighb < n_neighbs; i_neighb++) {
      Tubulin *site = head->xlink_->teth_neighbor_sites_[i_neighb];
      double weight = head->xlink_->GetBindingWeight_I_Teth(site);
      weights_summed += weight;
    }
  }
  return weights_summed;
}

double AssociatedProteinManagement::GetWeight_Bind_II_Teth() {

  double weights_summed = 0;
  // Sum over all single-bound tethered xlink extensions
  for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
    int n_bound = n_bound_i_teth_[max_neighbs_ + 1][x_dub];
    // Sum over xlinks at this specific extension
    for (int i_xlink(0); i_xlink < n_bound; i_xlink++) {
      AssociatedProtein *xlink =
          std::get<POP_T *>(bound_i_teth_[max_neighbs_ + 1][x_dub][i_xlink])
              ->xlink_;
      xlink->UpdateNeighborSites_II_Teth();
      int n_neighbs = xlink->n_teth_neighbor_sites_ii_;
      // Get weight of every possible orientation with neighbors
      for (int i_neighb(0); i_neighb < n_neighbs; i_neighb++) {
        Tubulin *site = xlink->teth_neighbor_sites_ii_[i_neighb];
        double weight = xlink->GetBindingWeight_II_Teth(site);
        weights_summed += weight;
      }
    }
  }
  if (verbose_ and weights_summed > 0.0) {
    printf("weights for bind_ii_teth is %g\n", weights_summed);
  }
  return weights_summed;
}

void AssociatedProteinManagement::Update_All_Lists() {

  // Update all extensions to ensure lists are accurate
  sys_timepoint start = sys_clock::now();
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    active_[i_entry]->UpdateExtension();
    if (active_[i_entry]->tethered_) {
      active_[i_entry]->motor_->UpdateExtension();
    }
  }
  sys_timepoint finish_ext = sys_clock::now();
  properties_->wallace.t_xlinks_[5] += (finish_ext - start).count();
  // Run through lists & update them
  properties_->microtubules.UpdateUnoccupied();
  sys_timepoint finish_unocc = sys_clock::now();
  properties_->wallace.t_xlinks_[6] += (finish_unocc - finish_ext).count();
  Update_Bound_I();
  sys_timepoint finish_bound = sys_clock::now();
  properties_->wallace.t_xlinks_[7] += (finish_bound - finish_unocc).count();
  if (parameters_->microtubules.count > 1) {
    Update_Bound_II();
  }
  if (parameters_->motors.c_bulk > 0.0 and parameters_->motors.tethers_active) {
    Update_Free_Teth();
    Update_Bound_Unteth();
    Update_Bound_I_Teth();
    if (parameters_->microtubules.count > 1) {
      Update_Bound_II_Teth();
    }
    properties_->kinesin4.UpdateBoundUntethered();
  }
}

void AssociatedProteinManagement::Update_List_Relay(std::string event,
                                                    std::string target_pop) {

  if (target_pop.substr(0, 5) == "unocc") {
    properties_->microtubules.UpdateUnoccupied();
  } else if (target_pop.substr(0, 8) == "bound_i_") {
    Update_Bound_I();
    if (event.substr(0, 7) == "bind_ii") {
      properties_->microtubules.UpdateUnoccupied();
      Update_Bind_II_Candidate();
    }
  } else if (target_pop.substr(0, 9) == "bound_ii_") {
    Update_Bound_II();
  } else if (parameters_->motors.tethers_active) {
    if (target_pop.substr(0, 5) == "free_") {
      Update_Free_Teth();
      if (event.substr(0, 6) == "bind_I") {
        properties_->microtubules.UpdateUnoccupied();
        Update_Bind_I_Teth_Candidate();
      }
    } else if (target_pop.substr(0, 8) == "bound_I_") {
      Update_Bound_I_Teth();
      if (event.substr(0, 7) == "bind_II") {
        properties_->microtubules.UpdateUnoccupied();
        Update_Bind_II_Teth_Candidate();
      }
    } else if (target_pop.substr(0, 9) == "bound_II_") {
      Update_Bound_II_Teth();
    } else if (target_pop.substr(0, 6) == "unteth") {
      properties_->kinesin4.UpdateBoundUntethered();
    } else {
      printf("somethin SPOOKY TWO in update_relay (XLINKS) - EXITING\n");
      std::cout << target_pop << std::endl;
      exit(1);
    }
  } else {
    printf("somethin SPOOKY in update_relay (XLINKS) - EXITING\n");
    std::cout << target_pop << std::endl;
    exit(1);
  }
}

void AssociatedProteinManagement::Update_Free_Teth() {

  n_free_teth_ = 0;
  for (int i_xlink = 0; i_xlink < n_active_; i_xlink++) {
    AssociatedProtein *xlink = active_[i_xlink];
    if (xlink->heads_active_ == 0 && xlink->tethered_ == true) {
      if (xlink->motor_->heads_active_ > 0) {
        free_teth_[n_free_teth_] = &xlink->head_one_;
        n_free_teth_++;
      } else {
        printf("woah. error in update_free_teth_list (XLINKS)\n");
        exit(1);
      }
    }
  }
}

void AssociatedProteinManagement::Update_Bound_Unteth() {

  n_bound_unteth_ = 0;
  for (int i_xlink = 0; i_xlink < n_active_; i_xlink++) {
    AssociatedProtein *xlink = active_[i_xlink];
    if (xlink->heads_active_ > 0 && xlink->tethered_ == false) {
      bound_unteth_[n_bound_unteth_] = xlink;
      n_bound_unteth_++;
    }
  }
}

void AssociatedProteinManagement::Update_Bound_I() {

  for (int n_neighbs = 0; n_neighbs <= max_neighbs_ + 1; n_neighbs++) {
    n_bound_i_[n_neighbs] = 0;
  }
  for (int i_xlink = 0; i_xlink < n_active_; i_xlink++) {
    AssociatedProtein *xlink = active_[i_xlink];
    bool relevant_entry(false);
    if (xlink->heads_active_ == 1) {
      if (!xlink->tethered_)
        relevant_entry = true;
      else if (xlink->motor_->heads_active_ == 0)
        relevant_entry = true;
    }
    if (relevant_entry) {
      AssociatedProtein::Monomer *head = xlink->GetActiveHead();
      int n_neighbs = head->site_->GetPRC1NeighborCount();
      bound_i_[n_neighbs][n_bound_i_[n_neighbs]] = head;
      n_bound_i_[n_neighbs]++;
      // Last entry in n_neighbs index holds ALL entries
      bound_i_[max_neighbs_ + 1][n_bound_i_[max_neighbs_ + 1]] = head;
      n_bound_i_[max_neighbs_ + 1]++;
    }
  }
}

void AssociatedProteinManagement::Update_Bound_I_Teth() {

  n_bound_i_teth_tot_ = 0;
  for (int n_neighbs(0); n_neighbs <= max_neighbs_ + 1; n_neighbs++) {
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
      n_bound_i_teth_[n_neighbs][x_dub] = 0;
    }
  }
  for (int i_xlink = 0; i_xlink < n_active_; i_xlink++) {
    AssociatedProtein *xlink = active_[i_xlink];
    if (xlink->heads_active_ == 1 && xlink->tethered_) {
      if (xlink->motor_->heads_active_ > 0) {
        xlink->motor_->UpdateExtension();
        // Ensure we didn't force an untether event
        if (xlink->tethered_) {
          int x_dub = xlink->motor_->x_dist_doubled_;
          AssociatedProtein::Monomer *head = xlink->GetActiveHead();
          int n_neighbs = head->site_->GetPRC1NeighborCount();
          //				printf("FOUND %i NEBS\n", n_neighbs);
          int index = n_bound_i_teth_[n_neighbs][x_dub];
          bound_i_teth_[n_neighbs][x_dub][index] = head;
          n_bound_i_teth_[n_neighbs][x_dub]++;
          index = n_bound_i_teth_[max_neighbs_ + 1][x_dub];
          bound_i_teth_[max_neighbs_ + 1][x_dub][index] = head;
          n_bound_i_teth_[max_neighbs_ + 1][x_dub]++;
          n_bound_i_teth_tot_++;
        }
      }
    }
  }
}

void AssociatedProteinManagement::Update_Bound_II() {

  for (int n_neighbs{0}; n_neighbs <= max_neighbs_ + 1; n_neighbs++) {
    for (int x{0}; x <= dist_cutoff_; x++) {
      n_bound_ii_[n_neighbs][x] = 0;
    }
  }
  for (int i_xlink = 0; i_xlink < n_active_; i_xlink++) {
    AssociatedProtein *xlink = active_[i_xlink];
    bool relevant_entry(false);
    if (xlink->heads_active_ == 2) {
      if (!xlink->tethered_) {
        relevant_entry = true;
      } else if (xlink->motor_->heads_active_ == 0) {
        relevant_entry = true;
      }
    }
    if (relevant_entry) {
      xlink->UpdateExtension();
      // Ensure we didn't force an unbind event
      if (xlink->heads_active_ == 2) {
        int x = xlink->x_dist_;
        // Head one
        int n_one = xlink->head_one_.GetPRC1NeighbCount();
        int i_one = n_bound_ii_[n_one][x];
        bound_ii_[n_one][x][i_one] = &xlink->head_one_;
        n_bound_ii_[n_one][x]++;
        // Last entry in n_neighbs index holds ALL entries
        i_one = n_bound_ii_[max_neighbs_ + 1][x]++;
        bound_ii_[max_neighbs_ + 1][x][i_one] = &xlink->head_one_;
        // Head two
        int n_two = xlink->head_two_.GetPRC1NeighbCount();
        int i_two = n_bound_ii_[n_two][x];
        bound_ii_[n_two][x][i_two] = &xlink->head_two_;
        n_bound_ii_[n_two][x]++;
        // Last entry in n_neighbs index holds ALL entries
        i_two = n_bound_ii_[max_neighbs_ + 1][x]++;
        bound_ii_[max_neighbs_ + 1][x][i_two] = &xlink->head_two_;
      }
    }
  }
}

void AssociatedProteinManagement::Update_Bound_II_Teth() {

  for (int n_neighbs{0}; n_neighbs <= max_neighbs_ + 1; n_neighbs++) {
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      for (int x{0}; x <= dist_cutoff_; x++) {
        n_bound_ii_teth_same_[n_neighbs][x_dub][x] = 0;
        n_bound_ii_teth_oppo_[n_neighbs][x_dub][x] = 0;
      }
    }
  }
  for (int i_xlink = 0; i_xlink < n_active_; i_xlink++) {
    AssociatedProtein *xlink = active_[i_xlink];
    if (xlink->heads_active_ == 2 and xlink->tethered_) {
      // Dont count xlinks tethered to satellite motors
      if (xlink->motor_->heads_active_ > 0) {
        xlink->UpdateExtension();
        xlink->motor_->UpdateExtension();
        // Ensure we didn't force an untether or unbind event
        if (xlink->heads_active_ == 2 and xlink->tethered_) {
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
}

void AssociatedProteinManagement::Update_Bind_II_Candidate() {

  int all_neighbs = max_neighbs_ + 1;
  int n_bound = n_bound_i_[all_neighbs];
  double weight_total = 0;
  double weight_local[n_bound];
  for (int i_entry(0); i_entry < n_bound; i_entry++) {
    weight_local[i_entry] = 0;
    POP_T *head = std::get<POP_T *>(bound_i_[all_neighbs][i_entry]);
    head->xlink_->UpdateNeighborSites_II();
    int n_sites = head->xlink_->n_neighbor_sites_;
    for (int i_site(0); i_site < n_sites; i_site++) {
      Tubulin *site = head->xlink_->neighbor_sites_[i_site];
      double weight = head->xlink_->GetBindingWeight_II(site);
      weight_local[i_entry] += weight;
    }
    weight_total += weight_local[i_entry];
  }
  if (weight_total > 0) {
    // Normalize local weights to get relative probabilities
    // Using these relative probs, randomly pick an entry
    bool failed(true);
    double p_cum(0);
    double ran = properties_->gsl.GetRanProb();
    for (int i_entry(0); i_entry < n_bound; i_entry++) {
      p_cum += weight_local[i_entry] / weight_total;
      if (p_cum > ran) {
        failed = false;
        bound_i_[all_neighbs][0] = bound_i_[all_neighbs][i_entry];
        n_bound_i_[all_neighbs] = 1;
        return;
      }
    }
    if (failed) {
      printf("nope in update_bind_ii_candidate; EXIT \n");
      exit(1);
    }
  } else
    n_bound_i_[all_neighbs] = 0;
}

void AssociatedProteinManagement::Update_Bind_I_Teth_Candidate() {

  int all_neighbs = max_neighbs_ + 1;
  int n_bound = n_free_teth_;
  double weight_total = 0;
  double weight_local[n_bound];
  for (int i_entry(0); i_entry < n_bound; i_entry++) {
    weight_local[i_entry] = 0;
    POP_T *head = std::get<POP_T *>(free_teth_[i_entry]);
    head->xlink_->UpdateNeighborSites_I_Teth();
    int n_sites = head->xlink_->n_teth_neighbor_sites_;
    for (int i_site(0); i_site < n_sites; i_site++) {
      Tubulin *site = head->xlink_->teth_neighbor_sites_[i_site];
      double weight = head->xlink_->GetBindingWeight_I_Teth(site);
      weight_local[i_entry] += weight;
    }
    weight_total += weight_local[i_entry];
  }
  if (weight_total > 0) {
    // Normalize local weights to get relative probabilities
    // Using these relative probs, randomly pick an entry
    bool failed(true);
    double p_cum(0);
    double ran = properties_->gsl.GetRanProb();
    for (int i_entry(0); i_entry < n_bound; i_entry++) {
      p_cum += weight_local[i_entry] / weight_total;
      if (p_cum > ran) {
        failed = false;
        free_teth_[0] = free_teth_[i_entry];
        n_free_teth_ = 1;
        return;
      }
    }
    if (failed) {
      printf("nope in update_bind_i_teth_candidate; EXIT \n");
      exit(1);
    }
  } else
    n_free_teth_ = 0;
}

void AssociatedProteinManagement::Update_Bind_II_Teth_Candidate() {

  int all_neighbs = max_neighbs_ + 1;
  int n_bound = n_bound_i_teth_tot_;
  int i_entry = 0;
  double weight_total = 0.0;
  double weight_local[n_bound];
  if (verbose_) {
    printf("For Bind_II_Teth_Candidate: %i xlinks to choose from\n", n_bound);
  }
  for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
    for (int i(0); i < n_bound_i_teth_[all_neighbs][x_dub]; i++) {
      weight_local[i_entry] = 0.0;
      POP_T *head = std::get<POP_T *>(bound_i_teth_[all_neighbs][x_dub][i]);
      head->xlink_->UpdateNeighborSites_II_Teth();
      int n_sites = head->xlink_->n_teth_neighbor_sites_ii_;
      for (int i_site(0); i_site < n_sites; i_site++) {
        Tubulin *site = head->xlink_->teth_neighbor_sites_ii_[i_site];
        double weight = head->xlink_->GetBindingWeight_II_Teth(site);
        weight_local[i_entry] += weight;
      }
      if (verbose_ and weight_local[i_entry] > 0.0) {
        printf("  Entry #%i has a weight of %g\n", i_entry,
               weight_local[i_entry]);
      }
      weight_total += weight_local[i_entry];
      i_entry++;
    }
  }
  if (verbose_) {
    printf("Total weight is %g\n", weight_total);
  }
  // Make sure we didn't bamboozle ourselves here
  if (i_entry != n_bound) {
    printf("check Get_Bind_II_Teth_Candidate()");
  }
  if (weight_total > 0.0) {
    // Normalize local weights to get relative probabilities
    // Using these relative probs, randomly pick an entry
    i_entry = 0;
    bool failed(true);
    double p_cum(0);
    double ran = properties_->gsl.GetRanProb();
    if (verbose_) {
      printf("ran is %g\n", ran);
    }
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
      for (int i(0); i < n_bound_i_teth_[all_neighbs][x_dub]; i++) {
        p_cum += weight_local[i_entry] / weight_total;
        if (verbose_) {
          printf("  p_cum is %g after entry %i\n", p_cum, i_entry);
        }
        if (p_cum > ran) {
          if (verbose_) {
            printf("  Picked entry #%i as candidate ", i_entry);
            printf("  (x_dub = %i)\n", x_dub);
          }
          failed = false;
          // Store candidate (arbitrary indices; needs to match w/ event)
          bound_i_teth_[all_neighbs][2 * teth_cutoff_][0] =
              bound_i_teth_[all_neighbs][x_dub][i];
          n_bound_i_teth_tot_ = 1;
          return;
        }
        i_entry++;
      }
    }
    if (failed) {
      printf("nope in update_bind_ii_teth_candidate; EXIT \n");
      exit(1);
    }
  } else {
    n_bound_i_teth_tot_ = 0;
  }
}

AssociatedProtein::Monomer *
AssociatedProteinManagement::CheckScratchFor(std::string pop) {

  if (verbose_) {
    printf("Scratch was checked for ");
    std::cout << pop << std::endl;
  }
  POP_T *head(nullptr);
  // Tethered xlinks at self-rest are counted as both same and oppo sites
  // However, only one label can be held, so use this hack for now
  std::string default_root{"bound_II_same_0_"};
  std::string alternate_root{"bound_II_oppo_0_"};
  std::string corrected_pop{"NUNYA"};
  bool pop_was_corrected{false};
  if (pop.substr(0, alternate_root.length()) == alternate_root) {
    pop_was_corrected = true;
    std::string suffix{pop.substr(default_root.length(), pop.length())};
    corrected_pop = default_root + suffix;
    if (verbose_) {
      std::cout << "corrected_pop: " << corrected_pop << std::endl;
    }
  }
  // Scan through scratch; check if any entries match desired pop
  for (int i_entry(0); i_entry < n_scratched_; i_entry++) {
    if (verbose_) {
      printf("   entry #%i: ", i_entry);
      std::cout << scratch_[i_entry]->state_ << std::endl;
    }
    if (pop_was_corrected and scratch_[i_entry]->state_ == corrected_pop) {
      head = scratch_[i_entry];
      if (verbose_) {
        printf("      FOUND! at index %i (W/ CORRECTION)\n", i_entry);
      }
      break;
    } else if (scratch_[i_entry]->state_ == pop) {
      head = scratch_[i_entry];
      if (verbose_) {
        printf("      FOUND! at index %i\n", i_entry);
      }
      break;
    }
  }
  return head;
}

void AssociatedProteinManagement::SaveToScratch(POP_T *head) {

  // head->xlink_->is_outdated_ = true;
  if (!head->in_scratch_) {
    if (verbose_) {
      std::cout << head->state_;
    }
    head->UpdateState();
    if (verbose_) {
      std::cout << " was updated to " << head->state_ << std::endl;
    }
    scratch_[n_scratched_] = head;
    n_scratched_++;
    head->in_scratch_ = true;
  } else if (verbose_) {
    std::cout << head->state_ << " is already in scratch!\n";
  }
  int n_neighbs = head->GetPRC1NeighbCount();
  if (n_neighbs == 1) {
    POP_T *neighb;
    int i_site = head->site_->index_;
    if (i_site == 0)
      neighb = head->site_->mt_->lattice_[i_site + 1].xlink_head_;
    else if (i_site == head->site_->mt_->n_sites_ - 1)
      neighb = head->site_->mt_->lattice_[i_site - 1].xlink_head_;
    else if (head->site_->mt_->lattice_[i_site + 1].xlink_head_ != nullptr)
      neighb = head->site_->mt_->lattice_[i_site + 1].xlink_head_;
    else if (head->site_->mt_->lattice_[i_site - 1].xlink_head_ != nullptr)
      neighb = head->site_->mt_->lattice_[i_site - 1].xlink_head_;
    if (!neighb->in_scratch_) {
      if (verbose_) {
        std::cout << neighb->state_;
      }
      neighb->UpdateState();
      if (verbose_) {
        std::cout << " (NB) was updated to " << neighb->state_ << std::endl;
      }
      scratch_[n_scratched_] = neighb;
      n_scratched_++;
      neighb->in_scratch_ = true;
    } else if (verbose_) {
      std::cout << neighb->state_ << " (NB) is already in scratch!\n";
    }
  } else if (n_neighbs == 2) {
    POP_T *neighb_one, *neighb_two;
    int i_site = head->site_->index_;
    neighb_one = head->site_->mt_->lattice_[i_site + 1].xlink_head_;
    if (!neighb_one->in_scratch_) {
      if (verbose_) {
        std::cout << neighb_one->state_;
      }
      neighb_one->UpdateState();
      if (verbose_) {
        std::cout << " (NB1) was updated to" << neighb_one->state_ << std::endl;
      }
      scratch_[n_scratched_] = neighb_one;
      n_scratched_++;
      neighb_one->in_scratch_ = true;
    } else if (verbose_) {
      std::cout << neighb_one->state_ << " (NB1) is already in scratch!\n";
    }
    neighb_two = head->site_->mt_->lattice_[i_site - 1].xlink_head_;
    if (!neighb_two->in_scratch_) {
      if (verbose_) {
        std::cout << neighb_two->state_;
      }
      neighb_two->UpdateState();
      if (verbose_) {
        std::cout << " (NB2) was updated to" << neighb_two->state_ << std::endl;
      }
      scratch_[n_scratched_] = neighb_two;
      n_scratched_++;
      neighb_two->in_scratch_ = true;
    } else if (verbose_) {
      std::cout << neighb_two->state_ << " (NB2) is already in scratch!\n";
    }
  }
}

void AssociatedProteinManagement::SaveNeighbsToScratch(SITE_T *site) {

  int n_neighbs = site->GetPRC1NeighborCount();
  if (n_neighbs == 1) {
    POP_T *neighb;
    int i_site = site->index_;
    if (i_site == 0)
      neighb = site->mt_->lattice_[i_site + 1].xlink_head_;
    else if (i_site == site->mt_->n_sites_ - 1)
      neighb = site->mt_->lattice_[i_site - 1].xlink_head_;
    else if (site->mt_->lattice_[i_site + 1].xlink_head_ != nullptr)
      neighb = site->mt_->lattice_[i_site + 1].xlink_head_;
    else if (site->mt_->lattice_[i_site - 1].xlink_head_ != nullptr)
      neighb = site->mt_->lattice_[i_site - 1].xlink_head_;
    if (!neighb->in_scratch_) {
      neighb->UpdateState();
      scratch_[n_scratched_] = neighb;
      n_scratched_++;
      neighb->in_scratch_ = true;
    }
  } else if (n_neighbs == 2) {
    POP_T *neighb_one, *neighb_two;
    int i_site = site->index_;
    neighb_one = site->mt_->lattice_[i_site + 1].xlink_head_;
    if (!neighb_one->in_scratch_) {
      neighb_one->UpdateState();
      scratch_[n_scratched_] = neighb_one;
      n_scratched_++;
      neighb_one->in_scratch_ = true;
    }
    neighb_two = site->mt_->lattice_[i_site - 1].xlink_head_;
    if (!neighb_two->in_scratch_) {
      neighb_two->UpdateState();
      scratch_[n_scratched_] = neighb_two;
      n_scratched_++;
      neighb_two->in_scratch_ = true;
    }
  }
}

void AssociatedProteinManagement::Run_KMC() {

  // if (verbose_) {
  //   printf("\n Start of new xlink KMC cycle\n");
  // }
  if (parameters_->xlinks.c_bulk == 0.0) {
    return;
  }
  sys_timepoint start = sys_clock::now();
  Update_All_Lists();
  sys_timepoint finish_list = sys_clock::now();
  properties_->wallace.t_xlinks_[1] += (finish_list - start).count();
  Refresh_Populations();
  sys_timepoint finish_pops = sys_clock::now();
  properties_->wallace.t_xlinks_[2] += (finish_pops - finish_list).count();
  Generate_Execution_Sequence();
  sys_timepoint finish_seq = sys_clock::now();
  properties_->wallace.t_xlinks_[3] += (finish_seq - finish_pops).count();
  if (verbose_ and !IDs_to_exe_.empty()) {
    printf("Executing %lu events\n", IDs_to_exe_.size());
  }
  for (int i_event = 0; i_event < IDs_to_exe_.size(); i_event++) {
    events_[IDs_to_exe_[i_event]].Execute();
  }
  sys_timepoint finish_all = sys_clock::now();
  properties_->wallace.t_xlinks_[4] += (finish_all - finish_seq).count();
  properties_->wallace.t_xlinks_[0] += (finish_all - start).count();
}

void AssociatedProteinManagement::Refresh_Populations() {

  // Clear scratch list
  for (int i_entry{0}; i_entry < n_scratched_; i_entry++) {
    scratch_[i_entry]->in_scratch_ = false;
  }
  n_scratched_ = 0;

  /*
  // Run through all active entries and update labels if necessary
  for (int i_entry(0); i_entry < n_active_; i_entry++) {
    AssociatedProtein *entry = active_[i_entry];
    entry->is_outdated_ = false;
    std::string state;
    std::string root, x, x_dub, neighbs;
    int n_neighbs;
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
      POP_T *active_head = entry->GetActiveHead();
      n_neighbs = active_head->GetPRC1NeighbCount();
      neighbs = std::to_string(n_neighbs);
      if (entry->tethered_) {
        if (entry->motor_->heads_active_ > 0) {
          x_dub = std::to_string(entry->motor_->x_dist_doubled_);
          root = std::string("bound_I_") + x_dub;
        } else {
          root = std::string("bound_i");
        }
      } else {
        root = std::string("bound_i");
      }
      state = root + std::string("_") + neighbs;
      if (active_head == &entry->head_one_) {
        entry->head_one_.state_ = state;
        entry->head_two_.state_ = std::string("docked");
      } else {
        entry->head_one_.state_ = std::string("docked");
        entry->head_two_.state_ = state;
      }
    } else if (entry->heads_active_ == 2) {
      x = std::to_string(entry->x_dist_);
      // Head one
      if (entry->tethered_) {
        if (entry->motor_->heads_active_ > 0) {
          x_dub = std::to_string(entry->motor_->x_dist_doubled_);
          // For x = rest dist, store as "same" regardless
          // Note: this accounted for in CheckScratchFor()
          if (entry->x_dist_ == rest_dist_) {
            root = "bound_II_same_" + x + "_" + x_dub;
          } else if (entry->head_one_.site_->EquilibriumInSameDirection()) {
            root = "bound_II_same_" + x + "_" + x_dub;
          } else {
            root = "bound_II_oppo_" + x + "_" + x_dub;
          }
        } else {
          root = std::string("bound_ii_") + x;
        }
      } else {
        root = std::string("bound_ii_") + x;
      }
      n_neighbs = entry->head_one_.GetPRC1NeighbCount();
      neighbs = std::to_string(n_neighbs);
      state = root + std::string("_") + neighbs;
      entry->head_one_.state_ = state;
      // Head two
      if (entry->tethered_) {
        if (entry->motor_->heads_active_ > 0) {
          x_dub = std::to_string(entry->motor_->x_dist_doubled_);
          // For x = rest dist, store as "same" regardless
          // Note: this accounted for in CheckScratchFor()
          if (entry->x_dist_ == rest_dist_) {
            root = "bound_II_same_" + x + "_" + x_dub;
          } else if (entry->head_two_.site_->EquilibriumInSameDirection()) {
            root = "bound_II_same_" + x + "_" + x_dub;
          } else {
            root = "bound_II_oppo_" + x + "_" + x_dub;
          }
        } else {
          root = std::string("bound_ii_") + x;
        }
      } else {
        root = std::string("bound_ii_") + x;
      }
      n_neighbs = entry->head_two_.GetPRC1NeighbCount();
      neighbs = std::to_string(n_neighbs);
      state = root + std::string("_") + neighbs;
      entry->head_two_.state_ = state;
    } else {
      printf("Check XLINK Refresh_pop; EXITING\n");
      exit(1);
    }
    //			printf("STATE IS ");
    //			std::cout << state << std::endl;
  }
  */
}

void AssociatedProteinManagement::Generate_Execution_Sequence() {

  int n_events = Sample_Event_Statistics();
  // Construct KMC list if any events are predicted, otherwise clear it
  if (n_events > 0) {
    int pre_array[n_events];
    int i_array(0);
    for (int i_event(0); i_event < events_.size(); i_event++) {
      int n_expected = events_[i_event].n_expected_;
      for (int i_entry(0); i_entry < n_expected; i_entry++) {
        if (i_array >= n_events) {
          printf("WOAH BUDDY. check XLINK genKMC\n");
          exit(1);
        }
        // Store & pass unique index of this event
        pre_array[i_array] = events_[i_event].ID_;
        i_array++;
      }
    }
    if (i_array != n_events) {
      printf("NOT SURE in xlink GEN KMC \n");
      exit(1);
    }
    // Shuffle for some guuuud random exe order
    if (n_events > 1)
      gsl_ran_shuffle(properties_->gsl.rng_, pre_array, n_events, sizeof(int));
    // Transfer info to permanent vector owned by MGMT
    IDs_to_exe_.resize(n_events);
    for (int i_entry(0); i_entry < n_events; i_entry++) {
      IDs_to_exe_[i_entry] = pre_array[i_entry];
    }
  } else
    IDs_to_exe_.clear();
}

int AssociatedProteinManagement::Sample_Event_Statistics() {

  // Scan through all events & get expected number of occurrences
  int n_events_tot = 0;
  for (int i_event(0); i_event < events_.size(); i_event++) {
    events_[i_event].SampleStatistics();
    n_events_tot += events_[i_event].n_expected_;
    if (events_[i_event].n_expected_ > 0) {
      if (verbose_) {
        printf("%i expected for ", events_[i_event].n_expected_);
        std::cout << events_[i_event].name_;
        printf(" (code: %i, target: ", events_[i_event].code_);
        std::cout << events_[i_event].target_pop_ << ")" << std::endl;
      }
    }
  }
  /* Scan through all target pops. & ensure none will become negative */
  // Primary scan: compare pops. segregated by ext. & n_neighbs
  for (int i_pop(0); i_pop < IDs_by_pop_.size(); i_pop++) {
    int n_competitors = IDs_by_pop_[i_pop].size();
    if (n_competitors > 1) {
      double p_tot(0);
      int n_events_loc(0);
      for (int i_entry(0); i_entry < n_competitors; i_entry++) {
        int i_competitor = IDs_by_pop_[i_pop][i_entry];
        n_events_loc += events_[i_competitor].n_expected_;
        p_tot += (events_[i_competitor].n_expected_ *
                  events_[i_competitor].p_occur_);
      }
      // By convention, always use first entry for total avail pop
      // (Ensure emplace_back pattern in Init_Events() matches this)
      int n_avail_loc = *events_[IDs_by_pop_[i_pop][0]].n_avail_;
      if (n_avail_loc == 0 && n_events_loc > 0) {
        printf("uhhhh ?? - %i\n", n_events_loc);
        std::cout << events_[IDs_by_pop_[i_pop][0]].name_ << std::endl;
        std::cout << events_[IDs_by_pop_[i_pop][0]].target_pop_ << std::endl;
      }
      while (n_events_loc > n_avail_loc) {
        double p_cum = 0;
        double ran = properties_->gsl.GetRanProb();
        for (int i_entry(0); i_entry < n_competitors; i_entry++) {
          int i_competitor = IDs_by_pop_[i_pop][i_entry];
          p_cum += (events_[i_competitor].n_expected_ *
                    events_[i_competitor].p_occur_) /
                   p_tot;
          if (ran < p_cum and events_[i_competitor].n_expected_ > 0) {
            events_[i_competitor].n_expected_--;
            n_events_tot--;
            n_events_loc--;
            p_tot -= events_[i_competitor].p_occur_;
            break;
          }
        }
      }
    }
  }
  // Secondary scan: compares exts. over ALL neighbors
  for (int i_root(0); i_root < IDs_by_root_.size(); i_root++) {
    int n_competitors = IDs_by_root_[i_root].size();
    if (n_competitors > 1) {
      double p_tot(0);
      bool coupled(false);
      int n_events_loc(0);
      for (int i_entry(0); i_entry < n_competitors; i_entry++) {
        int i_competitor = IDs_by_root_[i_root][i_entry];
        if (i_competitor < 0) {
          coupled = true;
          i_competitor = abs(i_competitor);
        }
        n_events_loc += events_[i_competitor].n_expected_;
        p_tot += (events_[i_competitor].n_expected_ *
                  events_[i_competitor].p_occur_);
      }
      int n_avail_loc = *n_avail_by_root_[i_root];
      if (coupled and n_avail_loc > 0 and n_events_loc > 0) {
        if (verbose_) {
          printf("For ");
          std::cout << events_[abs(IDs_by_root_[i_root][0])].target_pop_;
          printf(": %i avail; %i events expected\n", n_avail_loc, n_events_loc);
        }
        n_avail_loc /= 2;
        if (verbose_) {
          printf("n avail scaled down to %i\n", n_avail_loc);
        }
      }
      if (n_avail_loc == 0 && n_events_loc > 0) {
        printf("uhhhh TWO ?? - %i\n", n_events_loc);
        std::cout << events_[abs(IDs_by_root_[i_root][0])].name_ << std::endl;
        std::cout << events_[abs(IDs_by_root_[i_root][0])].target_pop_
                  << std::endl;
        exit(1);
      }
      while (n_events_loc > n_avail_loc) {
        double p_cum = 0;
        double ran = properties_->gsl.GetRanProb();
        for (int i_entry(0); i_entry < n_competitors; i_entry++) {
          int i_competitor = abs(IDs_by_root_[i_root][i_entry]);
          p_cum += (events_[i_competitor].n_expected_ *
                    events_[i_competitor].p_occur_) /
                   p_tot;
          if (ran < p_cum and events_[i_competitor].n_expected_ > 0) {
            events_[i_competitor].n_expected_--;
            if (verbose_) {
              printf("removed 1 expected from ");
              std::cout << events_[i_competitor].name_ << std::endl;
            }
            n_events_tot--;
            n_events_loc--;
            p_tot -= events_[i_competitor].p_occur_;
            break;
          }
        }
      }
    }
  }
  return n_events_tot;
}

void AssociatedProteinManagement::Execute_Function_Relay(ENTRY_T target,
                                                         int code) {

  bool verbose = verbose_;
  POP_T *head(nullptr);
  SITE_T *site(nullptr);
  ALT_T *motor_head(nullptr);
  try {
    head = std::get<POP_T *>(target);
  } catch (...) {
    try {
      site = std::get<SITE_T *>(target);
    } catch (...) {
      motor_head = std::get<ALT_T *>(target);
    }
  }
  if (head == nullptr && site == nullptr && motor_head == nullptr) {
    printf("What the FUCK!! RUN TO THE CENTER!!! (...kmc_relay)\n");
    printf("code: %i\n", code);
    if (verbose_) {
      properties_->wallace.PauseSim(10);
    }
    return;
    //		exit(1);
  }
  bool diffusion_event(false);
  // Diffusion events are encoded as negative values
  if (code < 0) {
    diffusion_event = true;
    code = abs(code);
  }
  if (diffusion_event) {
    switch (code) {
    case 10:
      if (verbose) {
        printf("diffuse_i_fwd on ");
        std::cout << head->state_ << std::endl;
      }
      Diffuse(head, 1);
      break;
    case 11:
      if (verbose) {
        printf("diffuse_i_bck on ");
        std::cout << head->state_ << std::endl;
      }
      Diffuse(head, -1);
      break;
    case 20:
      if (verbose) {
        printf("diffuse_ii_to_rest on ");
        std::cout << head->state_ << std::endl;
      }
      Diffuse(head, 1);
      break;
    case 21:
      if (verbose) {
        printf("diffuse_ii_fr_rest on ");
        std::cout << head->state_ << std::endl;
      }
      Diffuse(head, -1);
      break;
    case 30:
      if (verbose) {
        printf("diffuse_i_to_teth on ");
        std::cout << head->state_ << std::endl;
      }
      Diffuse(head, 1);
      break;
    case 31:
      if (verbose) {
        printf("diffuse_i_fr_rest on ");
        std::cout << head->state_ << std::endl;
      }
      Diffuse(head, -1);
      break;
    case 40:
      if (verbose) {
        printf("diffuse_ii_to_both on ");
        std::cout << head->state_ << std::endl;
      }
      Diffuse(head, 1);
      break;
    case 41:
      if (verbose) {
        printf("diffuse_ii_fr_both on ");
        std::cout << head->state_ << std::endl;
      }
      Diffuse(head, -1);
      break;
    case 50:
      if (verbose) {
        printf("diffuse_ii_to_self_fr_teth on ");
        std::cout << head->state_ << std::endl;
      }
      Diffuse(head, -1);
      break;
    case 51:
      if (verbose) {
        printf("diffuse_ii_fr_self_to_teth on ");
        std::cout << head->state_ << std::endl;
      }
      Diffuse(head, 1);
      break;
    }
  } else {
    switch (code) {
    case 10:
      if (verbose)
        printf("bind_i");
      Bind_I(site);
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
        printf("bind_ii on ");
        std::cout << head->state_ << std::endl;
      }
      Bind_II(head);
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
        printf("unbind_i on ");
        std::cout << head->state_ << std::endl;
      }
      Unbind_I(head);
      break;
    case 31:
      if (verbose) {
        printf("unbind_i_teth on ");
        std::cout << head->state_ << std::endl;
      }
      Unbind_I(head);
      break;
    case 40:
      if (verbose) {
        printf("unbind_ii on ");
        std::cout << head->state_ << std::endl;
      }
      Unbind_II(head);
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
        printf("tether_free on motor ID #%i\n", motor_head->motor_->ID_);
      }
      Tether_Free(motor_head);
      break;
    case 60:
      if (verbose) {
        printf("untether_free on ");
        std::cout << head->state_ << std::endl;
      }
      Untether_Free(head);
      break;
    }
  }
}

void AssociatedProteinManagement::Diffuse(POP_T *head, int dir) {

  if (head->site_ == nullptr) {
    printf("woah buddy; u cant diffuse on a NULL site!!!\n");
    return;
  }
  int dx{0};
  if (head->xlink_->tethered_) {
    // xlink->GetDirToRest is ill-defined for x = 0, so use motor's
    // function instead (multiply by -1 since xlink is the one moving)
    if (head->xlink_->motor_->heads_active_ > 0) {
      dx = dir * -1 * head->xlink_->motor_->GetDirectionTowardRest();
    }
    // if tethered to a satellite motor, the above is not a concern
    else {
      dx = dir * head->GetDirectionToRest();
    }
  }
  // if not tethered, the above is not a concern
  else {
    dx = dir * head->GetDirectionToRest();
  }
  if (dx == 0) {
    printf("Error in xlink diffusion brobro\n");
    return;
  }
  Tubulin *old_site = head->site_;
  // Cannot step off them MTs
  if (!(old_site->index_ == 0 and dx == -1) and
      !(old_site->index_ == old_site->mt_->n_sites_ - 1 and dx == 1)) {
    Tubulin *new_site = &old_site->mt_->lattice_[old_site->index_ + dx];
    if (!new_site->occupied_) {
      // Save xlink head(s) to scratch_ before modifying it (them)
      SaveToScratch(head);
      if (head->xlink_->heads_active_ == 2) {
        SaveToScratch(head->GetOtherHead());
      }
      // Save neighbors of new site to scratch_ before modifying them
      SaveNeighbsToScratch(new_site);
      old_site->xlink_head_ = nullptr;
      old_site->occupied_ = false;
      new_site->xlink_head_ = head;
      new_site->occupied_ = true;
      head->site_ = new_site;
      properties_->microtubules.FlagForUpdate();
      if (verbose_) {
        head->xlink_->UpdateExtension();
        printf("extension is %i after diffusion\n", head->xlink_->x_dist_);
      }
      // properties_->microtubules.DeleteFromUnoccupied(new_site);
      // properties_->microtubules.PushToUnoccupied(old_site);
    }
  }
}

void AssociatedProteinManagement::Bind_I(SITE_T *site) {

  // Save neighbors of unbound site to scratch_ before modifying them
  SaveNeighbsToScratch(site);
  // Randomly choose a free xlink
  AssociatedProtein *xlink = GetFreeXlink();
  // Place xlink onto site
  site->xlink_head_ = &xlink->head_one_;
  site->occupied_ = true;
  // Update xlink details
  xlink->heads_active_++;
  xlink->head_one_.site_ = site;
  // Add this xlink to active_ list
  active_[n_active_] = xlink;
  xlink->active_index_ = n_active_;
  n_active_++;
  properties_->microtubules.FlagForUpdate();
  // properties_->microtubules.DeleteFromUnoccupied(site);
}

void AssociatedProteinManagement::Bind_II(POP_T *bound_head) {

  // Save bound head to scratch_ before modifying it (by making it doubly-bound)
  SaveToScratch(bound_head);
  POP_T *unbound_head{bound_head->GetOtherHead()};
  Tubulin *site{nullptr};
  if (bound_head->xlink_->tethered_) {
    if (bound_head->xlink_->motor_->heads_active_ > 0) {
      site = bound_head->xlink_->GetWeightedSite_Bind_II_Teth();
    } else {
      site = bound_head->xlink_->GetWeightedSite_Bind_II();
    }
  } else {
    site = bound_head->xlink_->GetWeightedSite_Bind_II();
  }
  if (site != nullptr) {
    // printf("bind_ii went thru\n");
    // Save neighbors of site to scratch_ before modifying them
    SaveNeighbsToScratch(site);
    // Place xlink's second head onto site
    site->xlink_head_ = unbound_head;
    site->occupied_ = true;
    // Update xlink details
    unbound_head->site_ = site;
    unbound_head->xlink_->heads_active_++;
    properties_->microtubules.FlagForUpdate();
    // properties_->microtubules.DeleteFromUnoccupied(site);
  } else {
    printf("failed to bind_ii in XLINK MGMT\n");
    if (bound_head->xlink_->tethered_) {
      printf("    (tethered)\n");
    }
    // properties_->wallace.PauseSim(10);
    exit(1);
  }
}

void AssociatedProteinManagement::Unbind_I(POP_T *head) {

  // Get site
  Tubulin *site = head->site_;
  // Save neighbors of bound site to scratch_ before modifying them
  SaveNeighbsToScratch(site);
  // Remove xlink from site
  site->xlink_head_ = nullptr;
  site->occupied_ = false;
  // Update xlink details
  head->site_ = nullptr;
  head->xlink_->heads_active_--;
  properties_->microtubules.FlagForUpdate();
  // properties_->microtubules.PushToUnoccupied(site);
  // Check to see if xlink has become inactive
  bool now_inactive{false};
  if (head->xlink_->heads_active_ == 0) {
    // Fully unbound but tethered xlinks remain active ...
    if (head->xlink_->tethered_) {
      // ... unless they were tethered to a satellite motor
      if (head->xlink_->motor_->heads_active_ == 0) {
        now_inactive = true;
        head->xlink_->UntetherSatellite();
      }
    }
    // If untethered & unbound, automatically goes inactive
    else {
      now_inactive = true;
    }
  }
  // If xlink is now inactive, remove it from the active list
  // and replace it with the last entry to keep list updated
  if (now_inactive) {
    AssociatedProtein *last_entry = active_[n_active_ - 1];
    int this_index = head->xlink_->active_index_;
    active_[this_index] = last_entry;
    last_entry->active_index_ = this_index;
    n_active_--;
  }
}

void AssociatedProteinManagement::Unbind_II(POP_T *head) {

  // Save other head (that will remain bound) to scratch_ before modifying it
  SaveToScratch(head->GetOtherHead());
  // Get site
  Tubulin *site = head->site_;
  // Save neighbors of site to scratch_ before modifying them
  SaveNeighbsToScratch(site);
  // Remove xlink from site
  site->xlink_head_ = nullptr;
  site->occupied_ = false;
  // Update xlink details
  head->site_ = nullptr;
  head->xlink_->heads_active_--;
  properties_->microtubules.FlagForUpdate();
  // properties_->microtubules.PushToUnoccupied(site);
}

void AssociatedProteinManagement::Bind_I_Teth(POP_T *satellite_head) {

  // Save this xlink head to scratch_ before modifying it
  SaveToScratch(satellite_head);
  // Randomly choose a weighted neighbor site
  Tubulin *site = satellite_head->xlink_->GetWeightedSite_Bind_I_Teth();
  // Save neighbors of site to scratch_ before modifying thme
  SaveNeighbsToScratch(site);
  if (site != nullptr) {
    // Update site details
    site->xlink_head_ = satellite_head;
    site->occupied_ = true;
    // Update xlink details
    satellite_head->xlink_->heads_active_++;
    satellite_head->site_ = site;
    properties_->microtubules.FlagForUpdate();
    // properties_->microtubules.DeleteFromUnoccupied(site);
  } else {
    printf("failed to XLINK Bind_I_Free_Tethered\n");
  }
}

void AssociatedProteinManagement::Tether_Free(ALT_T *untethered_head) {

  Kinesin *motor = untethered_head->motor_;
  AssociatedProtein *xlink = GetFreeXlink();
  // Update xlink and motor details
  xlink->motor_ = motor;
  xlink->tethered_ = true;
  motor->xlink_ = xlink;
  motor->tethered_ = true;
  // Add xlink to active_ list
  active_[n_active_] = xlink;
  xlink->active_index_ = n_active_;
  n_active_++;
}

void AssociatedProteinManagement::Untether_Free(POP_T *satellite_head) {

  SaveToScratch(satellite_head);
  AssociatedProtein *xlink = satellite_head->xlink_;
  Kinesin *motor = xlink->motor_;
  // Update xlink and motor details
  xlink->motor_ = nullptr;
  xlink->tethered_ = false;
  motor->xlink_ = nullptr;
  motor->tethered_ = false;
  // Remove this xlink from active_, replace with last entry
  AssociatedProtein *last_entry = active_[n_active_ - 1];
  int this_index = xlink->active_index_;
  if (this_index != n_active_ - 1) {
    active_[this_index] = last_entry;
    last_entry->active_index_ = this_index;
  }
  n_active_--;
}