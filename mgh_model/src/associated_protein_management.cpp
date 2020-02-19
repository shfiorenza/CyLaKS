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
    wally_->ErrorExit("AP_MGMT::CalculateCutoffs() [2]");
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

  if (parameters_->xlinks.c_bulk > 0.0) {
    population_active_ = true;
  }
  if (parameters_->microtubules.count > 1) {
    crosslinking_active_ = true;
  }
  if (parameters_->motors.tethers_active and parameters_->motors.c_bulk > 0.0) {
    tethering_active_ = true;
  }
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
  teth_cutoff_ = properties_->kinesin4.motors_[0].teth_cutoff_;
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
      p_diffuse_ii_fr_rest_[n_neighbs][x] =
          weight_neighb * weight_fr * p_diffu_ii;
    }
  }
  // p_theory_.emplace("diffuse_i_fwd", p_diffuse_i_fwd_);
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
      // For singly-bound xlinks, diffusing one site changes x_dub by 2, so keep
      // probability 0.0 for 2*cutoff +/- 1
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
    // dE = E_f - E_i;
    double dE = 0.0 - neighb_energy[n_neighbs];
    // We only consider the energy penalty for unbinding since lambda = 1.0
    double weight_neighb_bind{1.0};
    // dE has units of kbT, so dividing by its numerical value isn't necessary
    double weight_neighb_unbind{exp(-lambda_neighb * dE)};
    p_bind_i_[n_neighbs] = k_on * c_xlink * delta_t;
    p_unbind_i_[n_neighbs] = weight_neighb_unbind * k_off_i * delta_t;
    p_unbind_ii_[n_neighbs].resize(dist_cutoff_ + 1);
    weight_bind_ii_[n_neighbs].resize(dist_cutoff_ + 1);
    for (int x{0}; x <= dist_cutoff_; x++) {
      // dU = U_f - U_i
      double dU_bind{spring_energy[x] - 0.0};
      double dU_unbind{0.0 - spring_energy[x]};
      double weight_bind{exp(-(1.0 - lambda_spring) * dU_bind / kbT)};
      weight_bind_ii_[n_neighbs][x] = weight_neighb_bind * weight_bind;
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
    POP_T *head = std::get<POP_T *>(target);
    Diffuse(head, 1);
  };
  auto exe_diffuse_bck = [&](ENTRY_T target) {
    POP_T *head = std::get<POP_T *>(target);
    Diffuse(head, -1);
  };
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    event_name = "diffuse_i_fwd" + N_NEIGHBS;
    events_.emplace_back(event_name, p_diffuse_i_fwd_[n_neighbs],
                         &n_bound_i_[n_neighbs], &bound_i_[n_neighbs],
                         exe_diffuse_fwd, binomial, set_ran_indices);
    event_name = "diffuse_i_bck" + N_NEIGHBS;
    events_.emplace_back(event_name, p_diffuse_i_bck_[n_neighbs],
                         &n_bound_i_[n_neighbs], &bound_i_[n_neighbs],
                         exe_diffuse_bck, binomial, set_ran_indices);
    if (!crosslinking_active_) {
      continue;
    }
    for (int x{0}; x <= dist_cutoff_; x++) {
      std::string X_DIST{"_" + std::to_string(x)};
      event_name = "diffuse_ii_to_rest" + N_NEIGHBS + X_DIST;
      events_.emplace_back(event_name, p_diffuse_ii_to_rest_[n_neighbs][x],
                           &n_bound_ii_[n_neighbs][x], &bound_ii_[n_neighbs][x],
                           exe_diffuse_fwd, binomial, set_ran_indices);
      event_name = "diffuse_ii_fr_rest" + N_NEIGHBS + X_DIST;
      events_.emplace_back(event_name, p_diffuse_ii_fr_rest_[n_neighbs][x],
                           &n_bound_ii_[n_neighbs][x], &bound_ii_[n_neighbs][x],
                           exe_diffuse_bck, binomial, set_ran_indices);
    }
  }
  // Diffuse_teth: same as above but for tethered populations
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    if (!tethering_active_) {
      continue;
    }
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      std::string X_DUB{"_" + std::to_string(x_dub)};
      event_name = "diffuse_i_to_teth_rest" + N_NEIGHBS + X_DUB;
      events_.emplace_back(
          event_name, p_diffuse_i_to_teth_rest_[n_neighbs][x_dub],
          &n_bound_i_teth_[n_neighbs][x_dub], &bound_i_teth_[n_neighbs][x_dub],
          exe_diffuse_fwd, binomial, set_ran_indices);
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
        event_name = "diffuse_ii_to_both" + N_NEIGHBS + X_DUB + X_DIST;
        events_.emplace_back(event_name,
                             p_diffuse_ii_to_both_[n_neighbs][x_dub][x],
                             &n_bound_ii_teth_same_[n_neighbs][x_dub][x],
                             &bound_ii_teth_same_[n_neighbs][x_dub][x],
                             exe_diffuse_fwd, binomial, set_ran_indices);
        event_name = "diffuse_ii_fr_both" + N_NEIGHBS + X_DUB + X_DIST;
        events_.emplace_back(event_name,
                             p_diffuse_ii_fr_both_[n_neighbs][x_dub][x],
                             &n_bound_ii_teth_same_[n_neighbs][x_dub][x],
                             &bound_ii_teth_same_[n_neighbs][x_dub][x],
                             exe_diffuse_bck, binomial, set_ran_indices);
        event_name = "diffuse_ii_to_self_fr_teth" + N_NEIGHBS + X_DUB + X_DIST;
        events_.emplace_back(event_name,
                             p_diffuse_ii_to_self_fr_teth_[n_neighbs][x_dub][x],
                             &n_bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                             &bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                             exe_diffuse_fwd, binomial, set_ran_indices);
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
    events_.emplace_back(
        event_name, p_bind_i_[n_neighbs],
        &properties_->microtubules.n_unoccupied_xl_[n_neighbs],
        &properties_->microtubules.unoccupied_list_xl_[n_neighbs], exe_bind_i,
        binomial, set_ran_indices);
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
      Set_Bind_I_Teth_Candidates(n_expected);
      return n_expected;
    };
    events_.emplace_back(event_name, p_bind_i_teth_base_,
                         &n_bind_i_teth_candidates_, &bind_i_teth_candidates_,
                         exe_bind_i_teth, poisson_i_teth, set_ran_indices);
  }
  // Bind II: binds the second crosslinker head to adjacent microtubule
  if (crosslinking_active_) {
    event_name = "bind_ii";
    auto exe_bind_ii = [&](ENTRY_T target) {
      Bind_II(std::get<POP_T *>(target));
    };
    // Poisson dist. is used w/ partition function for E-dependent binding
    auto poisson_ii = [&](double p, int n) {
      double weight{GetWeight_Bind_II()};
      if (weight == 0.0) {
        return 0;
      }
      int n_expected{properties_->gsl.SamplePoissonDist(p * weight)};
      Set_Bind_II_Candidates(n_expected);
      return n_expected;
    };
    events_.emplace_back(event_name, p_bind_ii_base_, &n_bind_ii_candidates_,
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
      Set_Bind_II_Teth_Candidates(n_expected);
      return n_expected;
    };
    events_.emplace_back(event_name, p_bind_ii_base_,
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
                           &n_bound_ii_[n_neighbs][x], bound_ii_[n_neighbs][x],
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
                         &properties_->kinesin4.n_bound_untethered_,
                         &properties_->kinesin4.bound_untethered_,
                         exe_tether_free, binomial, set_ran_indices);
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

void AssociatedProteinManagement::AddToActive(AssociatedProtein *xlink) {

  active_[n_active_] = xlink;
  xlink->active_index_ = n_active_;
  n_active_++;
}

void AssociatedProteinManagement::RemoveFromActive(AssociatedProtein *xlink) {

  AssociatedProtein *last_entry = active_[n_active_ - 1];
  int this_index = xlink->active_index_;
  if (this_index != n_active_ - 1) {
    active_[this_index] = last_entry;
    last_entry->active_index_ = this_index;
  }
  n_active_--;
}

void AssociatedProteinManagement::Update_Free_Teth() {

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

void AssociatedProteinManagement::Update_Bound_Unteth() {

  n_bound_unteth_ = 0;
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    AssociatedProtein *xlink = active_[i_entry];
    if (xlink->heads_active_ > 0 and !xlink->tethered_) {
      bound_unteth_[n_bound_unteth_++] = xlink;
    }
  }
}

void AssociatedProteinManagement::Update_Bound_I() {

  n_bind_ii_candidates_ = 0;
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    n_bound_i_[n_neighbs] = 0;
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    AssociatedProtein *xlink = active_[i_entry];
    if (!xlink->tethered_ or xlink->HasSatellite()) {
      AssociatedProtein::Monomer *head = xlink->GetActiveHead();
      int n_neighbs = head->GetPRC1NeighborCount();
      bound_i_[n_neighbs][n_bound_i_[n_neighbs]++] = head;
      bind_ii_candidates_[n_bind_ii_candidates_++] = head;
    }
  }
}

void AssociatedProteinManagement::Update_Bound_I_Teth() {

  n_bind_ii_teth_candidates_ = 0;
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      n_bound_i_teth_[n_neighbs][x_dub] = 0;
    }
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    AssociatedProtein *xlink = active_[i_entry];
    if (xlink->heads_active_ == 1 and xlink->tethered_) {
      if (xlink->motor_->heads_active_ > 0) {
        int x_dub = xlink->motor_->x_dist_doubled_;
        AssociatedProtein::Monomer *head = xlink->GetActiveHead();
        int n_neighbs = head->GetPRC1NeighborCount();
        bound_i_teth_[n_neighbs][x_dub][n_bound_i_teth_[n_neighbs][x_dub]++] =
            head;
        bind_ii_teth_candidates_[n_bind_ii_teth_candidates_++] = head;
      }
    }
  }
}

void AssociatedProteinManagement::Update_Bound_II() {

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
      int neighbs_one = xlink->head_one_.GetPRC1NeighborCount();
      bound_ii_[neighbs_one][x][n_bound_ii_[neighbs_one][x]++] =
          &xlink->head_one_;
      // Head two
      int neighbs_two = xlink->head_two_.GetPRC1NeighborCount();
      bound_ii_[neighbs_two][x][n_bound_ii_[neighbs_two][x]++] =
          &xlink->head_two_;
    }
  }
}

void AssociatedProteinManagement::Update_Bound_II_Teth() {

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

double AssociatedProteinManagement::GetWeight_Bind_II() {

  double weight{0.0};
  for (int i_entry{0}; i_entry < n_bind_ii_candidates_; i_entry++) {
    POP_T *head{std::get<POP_T *>(bind_ii_candidates_[i_entry])};
    weight += head->xlink_->GetTotalWeight_Bind_II();
  }
  return weight;
}

double AssociatedProteinManagement::GetWeight_Bind_I_Teth() {

  double weight{0.0};
  for (int i_entry{0}; i_entry < n_bind_i_teth_candidates_; i_entry++) {
    POP_T *head{std::get<POP_T *>(bind_i_teth_candidates_[i_entry])};
    weight += head->xlink_->GetTotalWeight_Bind_I_Teth();
  }
  return weight;
}

double AssociatedProteinManagement::GetWeight_Bind_II_Teth() {

  double weight{0.0};
  for (int i_entry{0}; i_entry < n_bind_ii_teth_candidates_; i_entry++) {
    POP_T *head{std::get<POP_T *>(bind_ii_teth_candidates_[i_entry])};
    weight += head->xlink_->GetTotalWeight_Bind_II_Teth();
  }
  return weight;
}

void AssociatedProteinManagement::Set_Bind_II_Candidates(int n_to_set) {

  double weight[n_bind_ii_candidates_];
  double weight_total{0.0};
  for (int i_entry{0}; i_entry < n_bind_ii_candidates_; i_entry++) {
    POP_T *head = std::get<POP_T *>(bind_ii_candidates_[i_entry]);
    weight[i_entry] = head->xlink_->GetTotalWeight_Bind_II();
    weight_total += weight[i_entry];
  }
  if (weight_total == 0.0) {
    wally_->ErrorExit("AP_MGMT::Set_Bind_II_Candidate()\n");
  }
  ENTRY_T selected_candidates[n_to_set];
  for (int i_set{0}; i_set < n_to_set; i_set++) {
    double p_cum{0.0};
    double ran{properties_->gsl.GetRanProb()};
    for (int i_entry{0}; i_entry < n_bind_ii_candidates_; i_entry++) {
      p_cum += (weight[i_entry] / weight_total);
      if (ran < p_cum) {
        selected_candidates[i_set] = bind_ii_candidates_[i_entry];
        weight_total -= weight[i_entry];
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
}

void AssociatedProteinManagement::Set_Bind_I_Teth_Candidates(int n_to_set) {

  double weight[n_bind_i_teth_candidates_];
  double weight_total{0.0};
  for (int i_entry{0}; i_entry < n_bind_i_teth_candidates_; i_entry++) {
    POP_T *head = std::get<POP_T *>(bind_i_teth_candidates_[i_entry]);
    weight[i_entry] = head->xlink_->GetTotalWeight_Bind_I_Teth();
    weight_total += weight[i_entry];
  }
  if (weight_total == 0.0) {
    wally_->ErrorExit("AP_MGMT::Set_Bind_I_Teth_Candidate()\n");
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

void AssociatedProteinManagement::Set_Bind_II_Teth_Candidates(int n_to_set) {

  double weight[n_bind_ii_teth_candidates_];
  double weight_total{0.0};
  for (int i_entry{0}; i_entry < n_bind_ii_teth_candidates_; i_entry++) {
    POP_T *head = std::get<POP_T *>(bind_ii_teth_candidates_[i_entry]);
    weight[i_entry] = head->xlink_->GetTotalWeight_Bind_II_Teth();
    weight_total += weight[i_entry];
  }
  if (weight_total == 0.0) {
    wally_->ErrorExit("AP_MGMT::Set_Bind_II_Teth_Candidate()\n");
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
}

void AssociatedProteinManagement::RunKMC() {

  if (!population_active_) {
    return;
  }
  UpdateLists();
  SampleEventStatistics();
  GenerateExecutionSequence();
  ExecuteEvents();
}

void AssociatedProteinManagement::UpdateLists() {

  properties_->microtubules.UpdateUnoccupied();
  Update_Bound_I();
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
    properties_->kinesin4.Update_Bound_Unteth();
  }
}

void AssociatedProteinManagement::SampleEventStatistics() {

  n_events_to_exe_ = 0;
  for (auto &&event : events_) {
    n_events_to_exe_ += event.SampleStatistics();
    // Ensure no funny business occurs with statistics
    if (event.n_expected_ > *event.n_avail_) {
      printf("Error; %i events expected but only %i available for %s\n",
             event.n_expected_, *event.n_avail_, event.name_.c_str());
      wally_->ErrorExit("AP_MGMT::SampleEventStatistics");
    }
  }
  if (n_events_to_exe_ <= 1) {
    return;
  }
  // Put all events with > 1 target into an active_ array
  EVENT_T *active_events[n_events_to_exe_];
  int i_active{0};
  for (auto &&event : events_) {
    for (int i_tar{0}; i_tar < event.n_expected_; i_tar++) {
      active_events[i_active++] = &event;
    }
  }
  // Scan through all active events to ensure that no two target the same
  // xlink
  int n_removed_tot{0};
  for (int i_event{0}; i_event < n_events_to_exe_; i_event++) {
    EVENT_T *event_i = active_events[i_event];
    // For each event, compare each target
    int n_removed_i{0};
    for (int i_tar{0}; i_tar < event_i->n_expected_; i_tar++) {
      i_tar -= n_removed_i;
      AssociatedProtein *xlink_i{nullptr};
      try {
        xlink_i = std::get<POP_T *>(event_i->targets_[i_tar])->xlink_;
      } catch (...) {
        continue;
      }
      // Compare w/ events starting at i + 1 to avoid double counting
      for (int j_event{i_event + 1}; j_event < n_events_to_exe_; j_event++) {
        EVENT_T *event_j = active_events[j_event];
        // If events are identical (i.e., n_expected > 1 for an event), skip
        if (event_i == event_j) {
          continue;
        }
        // Otherwise, scan over targets of event_j
        int n_removed_j{0};
        for (int j_tar{0}; j_tar < event_j->n_expected_; j_tar++) {
          j_tar -= n_removed_j;
          AssociatedProtein *xlink_j{nullptr};
          try {
            xlink_j = std::get<POP_T *>(event_j->targets_[j_tar])->xlink_;
          } catch (...) {
            continue;
          }
          // If event_i and event_j have the same target, roll to remove one
          if (xlink_i == xlink_j) {
            double p_one{event_i->p_occur_};
            double p_two{event_j->p_occur_};
            double ran{properties_->gsl.GetRanProb()};
            if (ran < p_one / (p_one + p_two)) {
              event_i->RemoveTarget(i_tar);
              n_removed_i++;
            } else {
              event_j->RemoveTarget(j_tar);
              n_removed_j++;
            }
            n_removed_tot++;
          }
        }
      }
    }
  }
  n_events_to_exe_ -= n_removed_tot;
}

void AssociatedProteinManagement::GenerateExecutionSequence() {

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
    wally_->ErrorExit("AssociatedProteinMGMT::GenerateExecutionSequence()");
  }
  if (n_events_to_exe_ > 1) {
    properties_->gsl.Shuffle(pre_array, n_events_to_exe_, sizeof(EVENT_T *));
  }
  for (int i_event{0}; i_event < n_events_to_exe_; i_event++) {
    events_to_exe_[i_event] = pre_array[i_event];
  }
}

void AssociatedProteinManagement::ExecuteEvents() {

  for (int i_event{0}; i_event < n_events_to_exe_; i_event++) {
    events_to_exe_[i_event]->Execute();
  }
}

void AssociatedProteinManagement::Diffuse(POP_T *head, int dir) {

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
    head->xlink_->UpdateExtension();
    properties_->microtubules.FlagForUpdate();
  }
}

void AssociatedProteinManagement::Bind_I(SITE_T *site) {

  AssociatedProtein *xlink{GetFreeXlink()};
  site->xlink_head_ = &xlink->head_one_;
  site->occupied_ = true;
  xlink->heads_active_++;
  xlink->head_one_.site_ = site;
  xlink->UpdateExtension();
  AddToActive(xlink);
  properties_->microtubules.FlagForUpdate();
}

void AssociatedProteinManagement::Bind_I_Teth(POP_T *satellite_head) {

  AssociatedProtein *xlink{satellite_head->xlink_};
  Tubulin *site{xlink->GetWeightedSite_Bind_I_Teth()};
  if (site == nullptr) {
    printf("failed to XLINK Bind_I_Free_Tethered\n");
    return;
  }
  site->xlink_head_ = satellite_head;
  site->occupied_ = true;
  satellite_head->site_ = site;
  xlink->heads_active_++;
  xlink->UpdateExtension();
  properties_->microtubules.FlagForUpdate();
}

void AssociatedProteinManagement::Bind_II(POP_T *bound_head) {

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
  xlink->UpdateExtension();
  properties_->microtubules.FlagForUpdate();
}

void AssociatedProteinManagement::Unbind_II(POP_T *head) {

  AssociatedProtein *xlink{head->xlink_};
  Tubulin *site{head->site_};
  site->xlink_head_ = nullptr;
  site->occupied_ = false;
  head->site_ = nullptr;
  xlink->heads_active_--;
  xlink->UpdateExtension();
  properties_->microtubules.FlagForUpdate();
}

void AssociatedProteinManagement::Unbind_I(POP_T *head) {

  AssociatedProtein *xlink{head->xlink_};
  Tubulin *site{head->site_};
  site->xlink_head_ = nullptr;
  site->occupied_ = false;
  head->site_ = nullptr;
  xlink->heads_active_--;
  xlink->UpdateExtension();
  xlink->UntetherSatellite();
  if (!xlink->tethered_) {
    RemoveFromActive(xlink);
  }
  properties_->microtubules.FlagForUpdate();
}

void AssociatedProteinManagement::Tether_Free(ALT_T *untethered_head) {

  AssociatedProtein *xlink{GetFreeXlink()};
  Kinesin *motor{untethered_head->motor_};
  xlink->motor_ = motor;
  xlink->tethered_ = true;
  motor->xlink_ = xlink;
  motor->tethered_ = true;
  AddToActive(xlink);
}

void AssociatedProteinManagement::Untether_Free(POP_T *satellite_head) {

  AssociatedProtein *xlink{satellite_head->xlink_};
  Kinesin *motor{xlink->motor_};
  xlink->motor_ = nullptr;
  xlink->tethered_ = false;
  motor->xlink_ = nullptr;
  motor->tethered_ = false;
  RemoveFromActive(xlink);
}