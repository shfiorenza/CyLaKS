#include "associated_protein_management.h"
#include "master_header.h"

AssociatedProteinManagement::AssociatedProteinManagement() {}

void AssociatedProteinManagement::Initialize(system_parameters *parameters,
                                             system_properties *properties) {

  wally_ = &properties->wallace;
  gsl_ = &properties->gsl;
  parameters_ = parameters;
  properties_ = properties;
  CalculateCutoffs();
  SetParameters();
  GenerateXLinks();
  InitializeLists();
  if (wally_->test_mode_ == nullptr) {
    InitializeEvents();
  } else {
    InitializeTestEnvironment();
    InitializeTestEvents();
  }
}

void AssociatedProteinManagement::CalculateCutoffs() {

  double kbT{parameters_->kbT};                           // in pN * nm
  double k_spring{parameters_->xlinks.k_spring};          // in pN / nm
  double site_size{parameters_->microtubules.site_size};  // in nm
  r_0_ = parameters_->xlinks.r_0 / site_size;             // in n_sites
  r_y_ = parameters_->microtubules.y_dist / site_size;    // in n_sites
  k_spring_eff_ = k_spring * site_size * site_size / kbT; // in 1 / (n_sites)^2
  /*
      For events that result in a change in energy, we use Boltzmann factors to
      scale rates appropriately. Detailed balance is satisfied with the factors:
                    exp{-(1 - lambda)*(delta_E)/(kB*T)}, and
                    exp{(lambda)*(delta_E)/(kB*T)}
      for an event and its complement (e.g., binding and unbinding), where
     lambda is a constant that ranges from 0 to 1, delta_E is the change in
     energy that results from the event, kB is Boltzmann's constant, and T is
     the temperature
    */
  // Lambda = 1.0 means all of the weight goes into unbinding
  // Lambda = 0.5 means the weight is equally split between binding & unbinding
  // Lambda = 0.0 means all of the weight goes into binding
  lambda_neighb_ = 1.0; // Lambda for neighbor interaction energies
  lambda_spring_ = 0.5; // Lambda for spring energies
  lambda_teth_ = 0.5;   // lambda for tether spring energies
  if (lambda_neighb_ != 1.0) {
    wally_->Log("Lambda != 1.0 for neighb interactions not implemented yet!\n");
    wally_->ErrorExit("AssociatedProteinManagement::SetParameters()\n");
  }
  double f_cutoff{0.0};
  for (int x_dist{0}; x_dist < 1000; x_dist++) {
    double r_x{x_dist};
    double r{sqrt(r_x * r_x + r_y_ * r_y_)};
    double dr{r - r_0_};                                 // in n_sites
    double energy_factor{0.5 * k_spring_eff_ * dr * dr}; // unitless energy
    double boltzmann_weight = exp(lambda_spring_ * energy_factor);
    if (boltzmann_weight > 1e2) {
      f_cutoff = fabs(dr * site_size * k_spring);
      dist_cutoff_ = x_dist;
      break;
    }
  }
  if (parameters_->microtubules.count > 1) {
    crosslinking_active_ = true;
    wally_->Log("\nFor crosslinkers:\n");
    wally_->Log("  rest_dist is %i\n", 0);
    wally_->Log("  dist_cutoff is %i (f = %g pN)\n\n", dist_cutoff_, f_cutoff);
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
  // Construct array of neighbor weights
  weight_neighb_bind_.resize(max_neighbs_ + 1);
  weight_neighb_unbind_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    // ! energy is in units of kbT
    // value in params is aboslute value; need to multiply by -1
    double energy{-1 * n_neighbs * parameters_->xlinks.interaction_energy};
    weight_neighb_bind_[n_neighbs] = exp(-(1.0 - lambda_neighb_) * energy);
    weight_neighb_unbind_[n_neighbs] = exp(lambda_neighb_ * energy);
  }
  auto xlink_test{strstr(wally_->test_mode_, "xlink")};
  if (xlink_test != nullptr) {
    properties_->microtubules.InitializeTestEnvironment();
  }
  // Construct array of spring extensions and cosines for each x_dist value
  possible_cosines_.resize(2 * dist_cutoff_ + 1);
  possible_extensions_.resize(2 * dist_cutoff_ + 1);
  // if (wally_->test_mode_ == nullptr) {
  Update_Extensions();
  // } else if (strcmp(wally_->test_mode_, "xlink_diffuse_ii") == 0) {
  //   Update_Extensions();
  // }
  double kbT{parameters_->kbT};
  double delta_t{parameters_->delta_t};
  double site_size{parameters_->microtubules.site_size};
  double x_squared{(site_size / 1000) * (site_size / 1000)}; //! in um^2
  // Characteristic timescale for doubly-bound xlink to diffuse 1 site
  double tau_ii{x_squared / (2 * parameters_->xlinks.diffu_coeff_ii)};
  double p_diffu_ii{delta_t / tau_ii};
  // We consider diff_fwd and diff_bck as two separate events, which
  // effectively doubles the probability to diffuse. To counteract this,
  // divide p_diff by 2.
  p_diffu_ii /= 2;
  p_avg_diffuse_ii_ = p_diffu_ii;
  double k_off_ii{parameters_->xlinks.k_off_ii};
  p_avg_unbind_ii_ = k_off_ii * delta_t;
  p_diffuse_ii_to_rest_.resize(max_neighbs_ + 1);
  p_diffuse_ii_fr_rest_.resize(max_neighbs_ + 1);
  p_unbind_ii_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    p_diffuse_ii_to_rest_[n_neighbs].resize(2 * dist_cutoff_ + 1);
    p_diffuse_ii_fr_rest_[n_neighbs].resize(2 * dist_cutoff_ + 1);
    p_unbind_ii_[n_neighbs].resize(2 * dist_cutoff_ + 1);
  }
  weight_spring_bind_.resize(2 * dist_cutoff_ + 1);
  weight_spring_unbind_.resize(2 * dist_cutoff_ + 1);
  // if (xlink_test != nullptr) {
  Update_Weights();
  // }

  /*
  // Calculate neighbor interaction energies
  std::vector<double> neighb_energy(max_neighbs_ + 1, 0.0);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    // Neighbor energies are negative since its an attractive potential
    neighb_energy[n_neighbs] = n_neighbs * interaction_energy; //! in kBT
  }
  // Calculate spring extension energies
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
    // wally_->Log("spring_energy[%i] = %g\n", x, spring_energy[x]);
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
  */
  // [DIFFUSION STATISTICS FOR CROSSLINKER W/O TETH BELOW] //  //
  // Characteristic timescale for singly-bound xlink to diffuse 1 site
  double tau_i{x_squared / (2 * parameters_->xlinks.diffu_coeff_i)};
  double p_diffu_i{delta_t / tau_i};
  p_diffu_i /= 2;
  p_diffuse_i_fwd_.resize(max_neighbs_ + 1);
  p_diffuse_i_bck_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    // currently, diffusion is thought of as an unbinding event
    double weight_neighb{weight_neighb_unbind_[n_neighbs]};
    // With two neighbors, the binding head is jammed and cannot diffuse
    if (n_neighbs == 2) {
      weight_neighb = 0.0;
    }
    p_diffuse_i_fwd_[n_neighbs] = weight_neighb * p_diffu_i;
    p_diffuse_i_bck_[n_neighbs] = weight_neighb * p_diffu_i;
    /*
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
      // Diffusing away from rest is considered a binding-type event in
    regards
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
    */
  }
  /*
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
      // Diffusing away from rest is considered a binding-type event in
  regards
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
      p_diffuse_ii_to_self_fr_teth_[n_neighbs][x_dub].resize(dist_cutoff_ +
  1); p_diffuse_ii_fr_self_to_teth_[n_neighbs][x_dub].resize(dist_cutoff_ +
  1); for (int x{0}; x <= dist_cutoff_; x++) {
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
  */
  // KMC STATISTICS BELOW
  double k_on = parameters_->xlinks.k_on;
  double c_xlink = parameters_->xlinks.c_bulk;
  double p_bind_i{k_on * c_xlink * delta_t};
  double c_eff_teth = parameters_->motors.c_eff_tether;
  p_avg_bind_i_teth_ = k_on * c_eff_teth * delta_t;
  double c_eff_bind = parameters_->xlinks.c_eff_bind;
  p_avg_bind_ii_ = k_on * c_eff_bind * delta_t;
  double k_off_i = parameters_->xlinks.k_off_i;
  double p_unbind_i{k_off_i * delta_t};
  p_bind_i_.resize(max_neighbs_ + 1);
  p_unbind_i_.resize(max_neighbs_ + 1);
  p_unbind_i_teth_.resize(max_neighbs_ + 1);
  /*
    p_unbind_ii_to_teth_.resize(max_neighbs_ + 1);
    p_unbind_ii_fr_teth_.resize(max_neighbs_ + 1);
    weight_bind_ii_.resize(max_neighbs_ + 1);
    weight_bind_i_teth_.resize(max_neighbs_ + 1);
    weight_bind_ii_to_teth_.resize(max_neighbs_ + 1);
    weight_bind_ii_fr_teth_.resize(max_neighbs_ + 1);
    */
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    p_bind_i_[n_neighbs] = weight_neighb_bind_[n_neighbs] * p_bind_i;
    p_unbind_i_[n_neighbs] = weight_neighb_unbind_[n_neighbs] * p_unbind_i;
    /*
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
      weight_bind_i_teth_[n_neighbs][x_dub] = weight_neighb_bind *
    weight_bind; double weight_unbind{exp(-lambda_teth * dU_unbind / kbT)};
      p_unbind_i_teth_[n_neighbs][x_dub] =
          weight_unbind * p_unbind_i_[n_neighbs];
      p_unbind_ii_to_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      p_unbind_ii_fr_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      weight_bind_ii_to_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      weight_bind_ii_fr_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      for (int x{0}; x <= dist_cutoff_; x++) {
        // Change in x_dub if 2nd xlink head were to (un)bind at current
    x_dist
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
        // regards to Boltzmann factors, since both events let the spring
    relax double weight_to_teth{exp(-lambda_teth * dU_to / kbT)};
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
      }
    }
    */
  }
  /*
  double k_tether = parameters_->motors.k_tether;
  p_tether_free_ = k_tether * c_xlink * delta_t;
  double k_untether = parameters_->motors.k_untether;
  p_untether_free_ = k_untether * delta_t;
  */
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
    /*
    if (crosslinking_active_) {
      p_theory_["diffuse_ii_to_both"] = p_diffuse_ii_to_both_;
      p_theory_["diffuse_ii_fr_both"] = p_diffuse_ii_fr_both_;
      p_theory_["diffuse_ii_to_self_fr_teth"] = p_diffuse_ii_to_self_fr_teth_;
      p_theory_["diffuse_ii_fr_self_to_teth"] = p_diffuse_ii_fr_self_to_teth_;
    }
    */
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
      /*
      p_theory_["unbind_ii_to_teth"] = p_unbind_ii_to_teth_;
      p_theory_["unbind_ii_fr_teth"] = p_unbind_ii_fr_teth_;
      */
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
          // if (verbosity_ >= 3) {
          wally_->Log("%s = %g\n", name.c_str(), value[i][j][k]);
          // }
          if (value[i][j][k] > 1.0) {
            wally_->Log("Error! %s = %g\n", name.c_str(), value[i][j][k]);
            // wally_->ErrorExit("AP_MGMT:SetParameters()");
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
  for (int i_entry{0}; i_entry < n_xlinks_; i_entry++) {
    xlinks_[i_entry].Initialize(parameters_, properties_, i_entry);
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
    n_bound_ii_[n_neighbs].resize(2 * dist_cutoff_ + 1);
    for (int x{0}; x <= 2 * dist_cutoff_; x++) {
      n_bound_ii_[n_neighbs][x] = 0;
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
    bound_ii_[n_neighbs].resize(2 * dist_cutoff_ + 1);
    for (int x{0}; x <= 2 * dist_cutoff_; x++) {
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

  // Various GSL functions packaged in lambda expressions so that they
  // are stand-alone and easily passable to event class by reference
  auto binomial = [&](double p, int n) {
    if (n == 0) {
      return 0;
    }
    return gsl_->SampleBinomialDist(p, n);
  };
  auto poisson = [&](double p, int n) {
    if (p == 0.0) {
      return 0;
    }
    return gsl_->SamplePoissonDist(p);
  };
  auto get_ran_prob = [&]() { return gsl_->GetRanProb(); };
  auto set_ran_indices = [&](int *indices, int n, int m) {
    if (m > 1) {
      gsl_->SetRanIndices(indices, n, m);
    } else {
      indices[0] = 0;
    }
  };
  /* * Event entries * */
  std::string event_name; // scratch space to construct each event name
  // Diffuse: steps head one site to the left or right
  auto exe_diffuse_fwd = [&](ENTRY_T target) {
    try {
      Diffuse(std::get<POP_T *>(target), 1);
    } catch (...) {
      wally_->ErrorExit("AP_MGMT::Diffuse_FWD");
    }
  };
  auto exe_diffuse_bck = [&](ENTRY_T target) {
    try {
      Diffuse(std::get<POP_T *>(target), -1);
    } catch (...) {
      wally_->ErrorExit("AP_MGMT::Diffuse_BCK");
    }
  };
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    // Diffuse_i_fwd
    event_name = "diffuse_i_fwd" + N_NEIGHBS;
    events_.emplace_back(event_name, &p_diffuse_i_fwd_[n_neighbs],
                         &n_bound_i_[n_neighbs], &bound_i_[n_neighbs],
                         exe_diffuse_fwd, binomial, set_ran_indices);
    // Diffuse_i_bck
    event_name = "diffuse_i_bck" + N_NEIGHBS;
    events_.emplace_back(event_name, &p_diffuse_i_bck_[n_neighbs],
                         &n_bound_i_[n_neighbs], &bound_i_[n_neighbs],
                         exe_diffuse_bck, binomial, set_ran_indices);
    if (!crosslinking_active_) {
      continue;
    }
    for (int x{0}; x <= 2 * dist_cutoff_; x++) {
      std::string X_DIST{"_" + std::to_string(x)};
      event_name = "diffuse_ii_to_rest" + N_NEIGHBS + X_DIST;
      events_.emplace_back(event_name, &p_diffuse_ii_to_rest_[n_neighbs][x],
                           &n_bound_ii_[n_neighbs][x], &bound_ii_[n_neighbs][x],
                           exe_diffuse_fwd, binomial, set_ran_indices);
      event_name = "diffuse_ii_fr_rest" + N_NEIGHBS + X_DIST;
      events_.emplace_back(event_name, &p_diffuse_ii_fr_rest_[n_neighbs][x],
                           &n_bound_ii_[n_neighbs][x], &bound_ii_[n_neighbs][x],
                           exe_diffuse_bck, binomial, set_ran_indices);
    }
  }
  /*
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
  */
  // Bind_I: binds first crosslinker head to the microtubule
  auto exe_bind_i = [&](ENTRY_T target) { Bind_I(std::get<SITE_T *>(target)); };
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    event_name = "bind_i" + N_NEIGHBS;
    events_.emplace_back(event_name, &p_bind_i_[n_neighbs],
                         &properties_->microtubules.n_unocc_xlink_[n_neighbs],
                         &properties_->microtubules.unocc_xlink_[n_neighbs],
                         exe_bind_i, binomial, set_ran_indices);
  }
  // Bind_I_Teth: same as above but for tethered unbound xlinks (satellites)
  if (tethering_active_) {
    /*
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
    */
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
    auto weight_bind_ii = [&](ENTRY_T target) {
      try {
        return std::get<POP_T *>(target)->xlink_->GetTotalWeight_Bind_II();
      } catch (...) {
        wally_->ErrorExit("AP_MGMT::Weight_Bind_II()");
      }
    };
    events_.emplace_back(event_name, p_avg_bind_ii_, &n_bind_ii_candidates_,
                         &bind_ii_candidates_, exe_bind_ii, weight_bind_ii,
                         poisson, get_ran_prob, set_ran_indices);
  }
  // Bind_II_Teth: same as above but for tethered population
  if (tethering_active_ and crosslinking_active_) {
    /*
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
    */
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
    for (int x{0}; x <= 2 * dist_cutoff_; x++) {
      std::string X_DIST{"_" + std::to_string(x)};
      event_name = "unbind_ii" + N_NEIGHBS + X_DIST;
      events_.emplace_back(event_name, &p_unbind_ii_[n_neighbs][x],
                           &n_bound_ii_[n_neighbs][x], &bound_ii_[n_neighbs][x],
                           exe_unbind_ii, binomial, set_ran_indices);
    }
  }
  /*
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
  */
  // Unbind_I: unbinds a head of singly-bound crosslinkres
  auto exe_unbind_i = [&](ENTRY_T target) {
    Unbind_I(std::get<POP_T *>(target));
  };
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    event_name = "unbind_i" + N_NEIGHBS;
    events_.emplace_back(event_name, &p_unbind_i_[n_neighbs],
                         &n_bound_i_[n_neighbs], &bound_i_[n_neighbs],
                         exe_unbind_i, binomial, set_ran_indices);
  }
  /*
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
  */
  /*
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
  */
}

void AssociatedProteinManagement::InitializeTestEnvironment() {

  if (strcmp(wally_->test_mode_, "xlink_bind_ii") == 0) {
    wally_->Log("Initializing xlink_bind_ii test in AP_MGMT\n");
    for (int x{0}; x <= 2 * dist_cutoff_; x++) {
      bind_ii_stats_.emplace_back(0, 0);
      unbind_ii_stats_.emplace_back(0, 0);
    }
    int i_site{dist_cutoff_ + 1};
    Tubulin *site{&properties_->microtubules.mt_list_[0].lattice_[i_site]};
    Bind_I(site);
  }
  if (strcmp(wally_->test_mode_, "xlink_diffuse_ii") == 0) {
    wally_->Log("Initializing xlink_diffuse_ii test in AP_MGMT\n");
    for (int x{0}; x <= 2 * dist_cutoff_; x++) {
      diff_ii_to_stats_.emplace_back(0, 0);
      diff_ii_fr_stats_.emplace_back(0, 0);
    }
    int i_site{properties_->microtubules.mt_list_[0].n_sites_ / 2};
    Tubulin *site{&properties_->microtubules.mt_list_[0].lattice_[i_site]};
    Bind_I(site);
    Update_Extensions();
    Update_Weights();
    Bind_II(site->xlink_head_);
  }
}

void AssociatedProteinManagement::InitializeTestEvents() {

  std::string event_name;
  auto binomial = [&](double p, int n) {
    if (n == 0) {
      return 0;
    }
    return properties_->gsl.SampleBinomialDist(p, n);
  };
  auto set_ran_indices = [&](int *indices, int n, int m) {
    if (m > 1) {
      properties_->gsl.SetRanIndices(indices, n, m);
    } else {
      indices[0] = 0;
    }
  };
  auto get_ran_prob = [&]() { return properties_->gsl.GetRanProb(); };
  if (strcmp(wally_->test_mode_, "xlink_bind_ii") == 0) {
    // Bind II: binds the second crosslinker head to adjacent microtubule
    event_name = "bind_ii";
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
    auto weight_bind_ii = [&](ENTRY_T target) {
      try {
        return std::get<POP_T *>(target)->xlink_->GetTotalWeight_Bind_II();
      } catch (...) {
        wally_->ErrorExit("AP_MGMT::Weight_Bind_II()");
      }
    };
    // Poisson dist. is used w/ partition function for E-dependent binding
    auto poisson_bind_ii = [&](double p, int n) {
      double offset{properties_->microtubules.GetSiteOffset()};
      for (int x{0}; x <= 2 * dist_cutoff_; x++) {
        double factor{1.0};
        // With offset = 0.0 or 0.5, there is a symmetry in available sites
        if (x != 0 and (offset == 0.0 or offset == 0.5)) {
          factor = 2.0;
        }
        bind_ii_stats_[x].second += factor * n_bind_ii_candidates_;
      }
      if (p == 0.0) {
        return 0;
      }
      return properties_->gsl.SamplePoissonDist(p);
    };
    events_.emplace_back(event_name, p_avg_bind_ii_, &n_bind_ii_candidates_,
                         &bind_ii_candidates_, exe_bind_ii, weight_bind_ii,
                         poisson_bind_ii, get_ran_prob, set_ran_indices);

    // Unbind II: unbinds a head of doubly-bound crosslinkers
    auto exe_unbind_ii = [&](ENTRY_T target) {
      POP_T *head{nullptr};
      try {
        head = std::get<POP_T *>(target);
      } catch (...) {
        wally_->ErrorExit("AP_MGMT::Unbind_II()");
      }
      unbind_ii_stats_[head->xlink_->x_dist_].first++;
      Unbind_II(head);
    };
    auto binomial_unbind_ii = [&](double p, int n) {
      // There will always be one active xlink
      if (active_[0]->heads_active_ == 2) {
        unbind_ii_stats_[active_[0]->x_dist_].second += n;
      }
      if (p == 0.0 or n == 0) {
        return 0;
      }
      return properties_->gsl.SampleBinomialDist(p, n);
    };
    for (int x{0}; x <= 2 * dist_cutoff_; x++) {
      std::string N_NEIGHBS{"_" + std::to_string(0)};
      std::string X_DIST{"_" + std::to_string(x)};
      event_name = "unbind_ii" + N_NEIGHBS + X_DIST;
      events_.emplace_back(event_name, &p_unbind_ii_[0][x], &n_bound_ii_[0][x],
                           &bound_ii_[0][x], exe_unbind_ii, binomial_unbind_ii,
                           set_ran_indices);
    }
  }
  if (strcmp(wally_->test_mode_, "xlink_diffuse_ii") == 0) {
    auto exe_diff_to = [&](ENTRY_T target) {
      POP_T *head{std::get<POP_T *>(target)};
      diff_ii_to_stats_[head->xlink_->x_dist_].first++;
      Diffuse(head, 1);
      head->xlink_->UpdateExtension();
    };
    auto exe_diff_fr = [&](ENTRY_T target) {
      POP_T *head{std::get<POP_T *>(target)};
      diff_ii_fr_stats_[head->xlink_->x_dist_].first++;
      Diffuse(head, -1);
      head->xlink_->UpdateExtension();
    };
    auto binomial_to = [&](double p, int n) {
      diff_ii_to_stats_[active_[0]->x_dist_].second += n;
      if (p == 0.0 or n == 0) {
        return 0;
      }
      return properties_->gsl.SampleBinomialDist(p, n);
    };
    auto binomial_fr = [&](double p, int n) {
      diff_ii_fr_stats_[active_[0]->x_dist_].second += n;
      if (p == 0.0 or n == 0) {
        return 0;
      }
      return properties_->gsl.SampleBinomialDist(p, n);
    };
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
      for (int x{0}; x <= 2 * dist_cutoff_; x++) {
        std::string X_DIST{"_" + std::to_string(x)};
        event_name = "diffuse_ii_to_rest" + N_NEIGHBS + X_DIST;
        events_.emplace_back(event_name, &p_diffuse_ii_to_rest_[n_neighbs][x],
                             &n_bound_ii_[n_neighbs][x],
                             &bound_ii_[n_neighbs][x], exe_diff_to, binomial_to,
                             set_ran_indices);
        event_name = "diffuse_ii_fr_rest" + N_NEIGHBS + X_DIST;
        events_.emplace_back(event_name, &p_diffuse_ii_fr_rest_[n_neighbs][x],
                             &n_bound_ii_[n_neighbs][x],
                             &bound_ii_[n_neighbs][x], exe_diff_fr, binomial_fr,
                             set_ran_indices);
      }
    }
  }
}

void AssociatedProteinManagement::ReportProbabilities() {

  if (!population_active_) {
    return;
  }

  // if (wally_->test_mode_ == nullptr) {
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
          double n_opp_tot{(double)event->n_opportunities_tot_};
          if (n_opp_tot < 0.0) {
            wally_->Log("WHAT?? n_opp = %g\n", n_opp_tot);
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
  // } else
  if (strcmp(wally_->test_mode_, "xlink_bind_ii") == 0) {
    wally_->Log("For xlink_bind_ii test:\n");
    for (int x{0}; x <= 2 * dist_cutoff_; x++) {
      wally_->Log("p_bind_theory_[%i] = %g\n", x,
                  p_avg_bind_ii_ * weight_spring_bind_[x]);
      wally_->Log("p_bind_actual_[%i] = %g (n_exe = %i)\n", x,
                  double(bind_ii_stats_[x].first) / bind_ii_stats_[x].second,
                  bind_ii_stats_[x].first);
    }
    for (int x{0}; x <= 2 * dist_cutoff_; x++) {
      wally_->Log("p_unbind_theory_[%i] = %g\n", x,
                  p_avg_unbind_ii_ * weight_spring_unbind_[x]);
      wally_->Log("p_unbind_actual_[%i] = %g (n_exe = %i)\n", x,
                  double(unbind_ii_stats_[x].first) /
                      unbind_ii_stats_[x].second,
                  unbind_ii_stats_[x].first);
    }
  } else if (strcmp(wally_->test_mode_, "xlink_diffuse_ii") == 0) {
    wally_->Log("For xlink_diffuse_ii test:\n");
    for (int x{0}; x <= 2 * dist_cutoff_; x++) {
      int x_to{x - 1};
      // x_dist == 1 is handled by diffuse_fr for x == 0
      if (x == 0) {
        x_to = dist_cutoff_ + 1;
      }
      if (x == dist_cutoff_ + 1) {
        x_to = 0;
      }
      wally_->Log("p_to_theory_[%i] = %g\n", x,
                  p_avg_diffuse_ii_ * weight_spring_unbind_[x] /
                      weight_spring_unbind_[x_to]);
      wally_->Log("p_to_actual_[%i] = %g (n_exe = %i)\n", x,
                  double(diff_ii_to_stats_[x].first) /
                      diff_ii_to_stats_[x].second,
                  diff_ii_to_stats_[x].first);
    }
    for (int x{0}; x <= 2 * dist_cutoff_; x++) {
      double p_theory{0.0};
      if (x != dist_cutoff_ and x != 2 * dist_cutoff_) {
        p_theory = p_avg_diffuse_ii_ * weight_spring_bind_[x + 1] /
                   weight_spring_bind_[x];
      }
      wally_->Log("p_fr_theory_[%i] = %g\n", x, p_theory);
      wally_->Log("p_fr_actual_[%i] = %g (n_exe = %i)\n", x,
                  double(diff_ii_fr_stats_[x].first) /
                      diff_ii_fr_stats_[x].second,
                  diff_ii_fr_stats_[x].first);
    }
  }
  /*
  if (strcmp(wally_->test_mode_, "xlink_bind_i_teth") == 0) {
    wally_->Log("For bind_i_teth:\n");
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      double p{(double)bind_i_teth_stats_[x_dub].first /
               bind_i_teth_stats_[x_dub].second};
      wally_->Log("p_theory_[%i] = %g\n", x_dub,
                  p_avg_bind_i_teth_ * weight_bind_i_teth_[0][x_dub]);
      wally_->Log("p_actual_[%i] = %g (n_exe = %i)\n", x_dub, p,
                  bind_i_teth_stats_[x_dub].first);
    }
  }
  if (strcmp(wally_->test_mode_, "xlink_bind_ii_teth") == 0) {
    wally_->Log("For bind_ii_to_teth:\n");
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      for (int x{0}; x <= dist_cutoff_; x++) {
        double p_to{(double)bind_ii_to_teth_stats_[x_dub][x].first /
                    bind_ii_to_teth_stats_[x_dub][x].second};
        wally_->Log("p_theory_[%i][%i] = %g\n", x_dub, x,
                    p_avg_bind_ii_ * weight_bind_ii_to_teth_[0][x_dub][x]);
        wally_->Log("p_actual_[%i][%i] = %g (n_exe = %i)\n", x_dub, x, p_to,
                    bind_ii_to_teth_stats_[x_dub][x].first);
      }
    }
    wally_->Log("For bind_ii_fr_teth:\n");
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      for (int x{0}; x <= dist_cutoff_; x++) {
        double p_fr{(double)bind_ii_fr_teth_stats_[x_dub][x].first /
                    bind_ii_fr_teth_stats_[x_dub][x].second};
        wally_->Log("p_theory_[%i][%i] = %g\n", x_dub, x,
                    p_avg_bind_ii_ * weight_bind_ii_fr_teth_[0][x_dub][x]);
        wally_->Log("p_actual_[%i][%i] = %g (n_exe = %i)\n", x_dub, x, p_fr,
                    bind_ii_fr_teth_stats_[x_dub][x].first);
      }
    }
  }
  */
}

AssociatedProtein *AssociatedProteinManagement::GetFreeXlink() {

  // Randomly choose an unbound xlink
  int i_xlink = gsl_->GetRanInt(n_xlinks_);
  AssociatedProtein *xlink = &xlinks_[i_xlink];
  int attempts{0};
  while (xlink->heads_active_ > 0 or xlink->tethered_) {
    attempts++;
    if (attempts > n_xlinks_) {
      wally_->ErrorExit("AP_MGMT::GetFreeXlink()");
    }
    i_xlink++;
    if (i_xlink == n_xlinks_) {
      i_xlink = 0;
    }
    xlink = &xlinks_[i_xlink];
  }
  return xlink;
}

void AssociatedProteinManagement::AddToActive(AssociatedProtein *xlink) {

  active_[n_active_] = xlink;
  xlink->active_index_ = n_active_;
  n_active_++;
}

void AssociatedProteinManagement::RemoveFromActive(AssociatedProtein *xlink) {

  AssociatedProtein *last_entry{active_[n_active_ - 1]};
  int this_index{xlink->active_index_};
  active_[this_index] = last_entry;
  last_entry->active_index_ = this_index;
  n_active_--;
}

void AssociatedProteinManagement::FlagForUpdate() { lists_up_to_date_ = false; }

void AssociatedProteinManagement::Update_Extensions() {

  /*
  if (!crosslinking_active_ and !tethering_active_) {
    return;
  }
  */
  // Offset ranges from 0 to 0.5; corresponds to the fractional misalignment
  // of MTs in n_sites, e.g., offset = 0 means they are perfectly aligned and
  // offset = 0.25 means all sites are misaligned by a 0.25*site_size
  double offset{properties_->microtubules.GetSiteOffset()};
  for (int x_dist{0}; x_dist <= 2 * dist_cutoff_; x_dist++) {
    // printf("\nx_dist is %i\n", x_dist);
    double r_x{x_dist + offset};
    if (x_dist > dist_cutoff_) {
      r_x = (x_dist - dist_cutoff_) - offset;
    }
    // printf("r_x is %g\n", r_x);
    double r{sqrt(r_x * r_x + r_y_ * r_y_)};
    possible_extensions_[x_dist] = r - r_0_;
    // printf("extension is %g\n", possible_extensions_[x_dist]);
    possible_cosines_[x_dist] = r_x / r;
    // printf("cosine is %g\n", possible_cosines_[x_dist]);
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    active_[i_entry]->UpdateExtension();
  }
}

void AssociatedProteinManagement::Update_Weights() {

  /*
    if (!crosslinking_active_ or !population_active_) {
      printf("boink\n");
      exit(1);
      return;
    }
    */
  for (int x_dist{0}; x_dist <= 2 * dist_cutoff_; x_dist++) {
    double dr{possible_extensions_[x_dist]};             // in n_sites
    double spring_energy{0.5 * k_spring_eff_ * dr * dr}; // unitless
    weight_spring_bind_[x_dist] = exp(-(1.0 - lambda_spring_) * spring_energy);
    weight_spring_unbind_[x_dist] = exp(lambda_spring_ * spring_energy);
  }
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    for (int x_dist{0}; x_dist <= 2 * dist_cutoff_; x_dist++) {
      p_unbind_ii_[n_neighbs][x_dist] = p_avg_unbind_ii_ *
                                        weight_spring_unbind_[x_dist] *
                                        weight_neighb_unbind_[n_neighbs];
      if (n_neighbs == max_neighbs_) {
        // printf("[%i][%i] -> 0\n", n_neighbs, x_dist);
        p_diffuse_ii_to_rest_[n_neighbs][x_dist] = 0.0;
        p_diffuse_ii_fr_rest_[n_neighbs][x_dist] = 0.0;
        continue;
      }
      int x_fr{x_dist + 1};
      int x_to{x_dist - 1};
      if (x_dist == 0) {
        x_to = dist_cutoff_ + 1;
      } else if (x_dist == dist_cutoff_ + 1) {
        x_to = 0;
      }
      // printf("x_to = %i | x_fr = %i\n", x_to, x_fr);
      // Dividing in the following way results in:
      // p_diff_ii * exp{lambda_spring_ * -(E_f - E_i)} * weight_neighb
      // Which will get larger with x_dist since E_f < E_i generally
      p_diffuse_ii_to_rest_[n_neighbs][x_dist] =
          p_avg_diffuse_ii_ * weight_spring_unbind_[x_dist] /
          weight_spring_unbind_[x_to] * weight_neighb_unbind_[n_neighbs];
      // printf("p_to = %g * %g / %g * %g\n", p_avg_diffuse_ii_,
      //        weight_spring_unbind_[x_dist], weight_spring_unbind_[x_to],
      //        weight_neighb_unbind_[n_neighbs]);
      // printf("     = %g\n", p_diffuse_ii_to_rest_[n_neighbs][x_dist]);
      if (x_dist % dist_cutoff_ == 0 and x_dist != 0) {
        p_diffuse_ii_fr_rest_[n_neighbs][x_dist] = 0.0;
        // printf("x_fr = %i --> 0.0\n", x_fr);
        continue;
      }
      // Dividing in the following way results in:
      // p_diff_ii * exp{-(1.0 - lambda_spring_) * (E_f - E_i)} *
      // weight_neighb Which will get smaller with x_dist since E_f > E_i
      p_diffuse_ii_fr_rest_[n_neighbs][x_dist] =
          p_avg_diffuse_ii_ * weight_spring_bind_[x_fr] /
          weight_spring_bind_[x_dist] * weight_neighb_unbind_[n_neighbs];
      // printf("p_fr = %g * %g / %g * %g\n", p_avg_diffuse_ii_,
      //        weight_spring_bind_[x_fr], weight_spring_bind_[x_dist],
      //        weight_neighb_unbind_[n_neighbs]);
      // printf("     = %g\n", p_diffuse_ii_fr_rest_[n_neighbs][x_dist]);
    }
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
  if (verbosity_ < 2) {
    return;
  }
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
    for (int x{0}; x <= 2 * dist_cutoff_; x++) {
      n_bound_ii_[n_neighbs][x] = 0;
    }
  }

  for (int i_xlink{0}; i_xlink < n_active_; i_xlink++) {
    AssociatedProtein *xlink = active_[i_xlink];
    if (xlink->heads_active_ != 2) {
      continue;
    }
    if (!xlink->tethered_ or xlink->HasSatellite()) {
      int x{xlink->x_dist_};
      // If we're testing bind_ii or bind_ii_teth, head one remains fixed
      if (wally_->test_mode_ == nullptr) {
        POP_T *head_one{&xlink->head_one_};
        int neighbs_one{head_one->GetPRC1NeighborCount()};
        bound_ii_[neighbs_one][x][n_bound_ii_[neighbs_one][x]++] = head_one;
      }
      // Always add head two
      POP_T *head_two{&xlink->head_two_};
      int neighbs_two{head_two->GetPRC1NeighborCount()};
      bound_ii_[neighbs_two][x][n_bound_ii_[neighbs_two][x]++] = head_two;
    }
  }
  if (verbosity_ < 2) {
    return;
  }
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    for (int x{0}; x <= 2 * dist_cutoff_; x++) {
      wally_->Log(" n_bound_ii_[%i][%i]  = %i\n", n_neighbs, x,
                  n_bound_ii_[n_neighbs][x]);
      for (int i_entry{0}; i_entry < n_bound_ii_[n_neighbs][x]; i_entry++) {
        POP_T *head{std::get<POP_T *>(bound_ii_[n_neighbs][x][i_entry])};
        wally_->Log("   - xlink %i\n", head->xlink_->id_);
      }
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
        if (x == 0) {
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
        if (x == 0) {
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

void AssociatedProteinManagement::RunKMC() {

  // if (properties_->current_step_ > 94445) {
  //   verbosity_ = 3;
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
  Update_Weights();
  properties_->microtubules.UpdateUnoccupied();
  Update_Bound_I();
  if (crosslinking_active_) {
    Update_Bound_II();
  }
  /*
  if (tethering_active_) {
    Update_Free_Teth();
    Update_Bound_Unteth();
    Update_Bound_I_Teth();
    if (crosslinking_active_) {
      Update_Bound_II_Teth();
    }
    properties_->kinesin4.Update_Bound_Unteth();
  }
  */
}

void AssociatedProteinManagement::SampleEventStatistics() {

  if (verbosity_ >= 1) {
    wally_->Log("Starting AP_MGMT::SampleEventStatistics()\n");
  }
  n_events_to_exe_ = 0;
  for (auto &&event : events_) {
    if (verbosity_ >= 3) {
      wally_->Log("Sampling statistics for %s\n", event.name_.c_str());
    }
    n_events_to_exe_ += event.SampleStatistics();
    // Ensure no funny business occurs with statistics
    if (event.n_expected_ > *event.n_avail_) {
      wally_->Log("Error; %i events expected but only %i available for %s\n",
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
          wally_->Log("    -target %i: NULL (i_var = %lu)\n", i_entry, index);
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
  // Put all events with >0 targets into an active_ array
  std::pair<EVENT_T *, ENTRY_T> active_events[n_events_to_exe_];
  int i_active{0};
  for (auto &&event : events_) {
    // Add a ptr to the event for each target it has
    for (int i_tar{0}; i_tar < event.n_expected_; i_tar++) {
      active_events[i_active++] = std::make_pair(&event, event.targets_[i_tar]);
      if (verbosity_ >= 2) {
        wally_->Log(" %s is an active_event\n", event.name_.c_str());
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
    gsl_->Shuffle(pre_array, n_events_to_exe_, sizeof(EVENT_T *));
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
    wally_->Log("woah buddy; u cant diffuse on a NULL site!!!\n");
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
    wally_->Log("failed toBind_I_Free_Tethered (XLINK)\n");
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
  if (!xlink->tethered_ or xlink->HasSatellite()) {
    site = bound_head->xlink_->GetWeightedSite_Bind_II();
  } else {
    site = bound_head->xlink_->GetWeightedSite_Bind_II_Teth();
  }
  // In case another crosslinker binds and takes up the only available
  // neighbor
  if (site == nullptr) {
    return;
    wally_->Log("woah XL.\n");
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