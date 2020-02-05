#include "associated_protein_management.h"
#include "master_header.h"

AssociatedProteinManagement::AssociatedProteinManagement() {}

void AssociatedProteinManagement::Initialize(system_parameters *parameters,
                                             system_properties *properties) {
  parameters_ = parameters;
  properties_ = properties;
  GenerateXLinks();
  SetParameters();
  InitializeLists();
  // InitializeEvents();
}

void AssociatedProteinManagement::GenerateXLinks() {

  int n_mts = parameters_->microtubules.count;
  n_xlinks_ = 0;
  for (int i_mt = 0; i_mt < n_mts; i_mt++) {
    n_xlinks_ += parameters_->microtubules.length[i_mt];
  }
  // Since only one head has to be bound, the sim will at most
  // as many xlinks as sites in the bulk (all single-bound)
  xlinks_.resize(n_xlinks_);
  for (int id = 0; id < n_xlinks_; id++) {
    xlinks_[id].Initialize(parameters_, properties_, id);
  }
}

void AssociatedProteinManagement::SetParameters() {

  /*    Assume lambda = 0.5 & 0.5 for binding & unbinding,		  *
   *  	whereas lambda = 1 & 0 for diffusing away/towards	  */
  max_neighbs_ = 2; // can have one in front & one behind
  interaction_energy_ = -1 * parameters_->xlinks.interaction_energy;
  double delta_t = parameters_->delta_t;
  double site_size = parameters_->microtubules.site_size;
  // DIFFUSION STATISTICS FOR SELF BELOW
  double D_coeff = parameters_->xlinks.diffusion_coeff;
  double x_squared = (site_size / 1000) * (site_size / 1000); // in um^2
  double tau = x_squared / (2 * D_coeff);
  p_diffuse_i_fwd_.resize(max_neighbs_ + 1);
  p_diffuse_i_bck_.resize(max_neighbs_ + 1);
  if (parameters_->xlinks.c_bulk > 0.0) {
    printf("add kbt to bolzmann factor in xlinks YA DOOFUS\n");
    exit(1);
  }
  for (int n_neighbs(0); n_neighbs <= max_neighbs_; n_neighbs++) {
    double tot_E = n_neighbs * interaction_energy_;
    // Lambda = 1 & 0 for diffusion
    double weight = exp(tot_E);
    if (n_neighbs == 2)
      weight = 0;
    p_diffuse_i_fwd_[n_neighbs] = (delta_t / tau) * weight;
    p_diffuse_i_bck_[n_neighbs] = (delta_t / tau) * weight;
  }
  // Generate different stepping rates based on changes in
  // potential energy (dU) associated with that step
  dist_cutoff_ = xlinks_[0].dist_cutoff_;
  rest_dist_ = xlinks_[0].rest_dist_;
  teth_cutoff_ = properties_->kinesin4.teth_cutoff_;
  comp_cutoff_ = properties_->kinesin4.comp_cutoff_;
  if (parameters_->motors.tethers_active) {
    properties_->wallace.Log("\nFor crosslinkers:\n");
    properties_->wallace.Log("  rest_dist is %i\n", rest_dist_);
    properties_->wallace.Log("  dist_cutoff is %i\n\n", dist_cutoff_);
  }
  double kbT = parameters_->kbT;
  double r_0 = xlinks_[0].r_0_;
  double k_spring = xlinks_[0].k_spring_;
  double r_y = parameters_->microtubules.y_dist;
  p_diffuse_ii_to_rest_.resize(max_neighbs_ + 1);
  p_diffuse_ii_fr_rest_.resize(max_neighbs_ + 1);
  for (int n_neighbs(0); n_neighbs <= max_neighbs_; n_neighbs++) {
    p_diffuse_ii_to_rest_[n_neighbs].resize(dist_cutoff_ + 1);
    p_diffuse_ii_fr_rest_[n_neighbs].resize(dist_cutoff_ + 1);
    double tot_E = n_neighbs * interaction_energy_;
    double weight_neighb = exp(tot_E);
    if (n_neighbs == max_neighbs_)
      weight_neighb = 0;
    for (int x_dist = 0; x_dist <= dist_cutoff_; x_dist++) {
      double r_x = x_dist * site_size;
      double r_x_to = (x_dist - 1) * site_size;
      double r_x_fr = (x_dist + 1) * site_size;
      double r = sqrt(r_x * r_x + r_y * r_y);
      double r_to = sqrt(r_x_to * r_x_to + r_y * r_y);
      double r_fr = sqrt(r_x_fr * r_x_fr + r_y * r_y);
      // Get extension for current dist and steps to/from spring rest
      double dr = r - r_0;
      double dr_to = r_to - r_0;
      double dr_fr = r_fr - r_0;
      if (dr >= 0) {
        // Get corresponding changes in potential energy
        double dU_to_rest = (k_spring / 2) * (dr_to * dr_to - dr * dr);
        double dU_fr_rest = (k_spring / 2) * (dr_fr * dr_fr - dr * dr);
        // Weights according to Lanksy et al.
        double weight_to = exp(-dU_to_rest / (2 * kbT));
        double weight_fr = exp(-dU_fr_rest / (2 * kbT));
        double p_to = weight_neighb * weight_to * delta_t / tau;
        double p_fr = weight_neighb * weight_fr * delta_t / tau;
        if (x_dist == 0) {
          p_diffuse_ii_to_rest_[n_neighbs][x_dist] = 0;
          p_diffuse_ii_fr_rest_[n_neighbs][x_dist] = 2 * p_fr;
        } else if (x_dist == dist_cutoff_) {
          p_diffuse_ii_to_rest_[n_neighbs][x_dist] = p_to;
          p_diffuse_ii_fr_rest_[n_neighbs][x_dist] = 0;
        } else {
          p_diffuse_ii_to_rest_[n_neighbs][x_dist] = p_to;
          p_diffuse_ii_fr_rest_[n_neighbs][x_dist] = p_fr;
        }
        if (p_to > 1)
          printf("WARNING: p_diffuse_to_rest=%g for x=%i\n", p_to, x_dist);
        if (2 * p_fr > 1)
          printf("WARNING: 2*p_diffuse_fr_rest=%g for x=%i\n", 2 * p_fr,
                 x_dist);
      } else {
        printf("woah mayne. xlink set parameters \n");
        exit(1);
      }
    }
  }
  // DIFFUSION STATISTICS INVOLVING TETHER BELOW
  int teth_cutoff = properties_->kinesin4.motors_[0].teth_cutoff_;
  double k_teth_spring = properties_->kinesin4.motors_[0].k_spring_;
  double k_teth_slack = properties_->kinesin4.motors_[0].k_slack_;
  double r_0_teth = properties_->kinesin4.motors_[0].r_0_;
  double r_y_teth = parameters_->microtubules.y_dist / 2;
  double rest_teth = properties_->kinesin4.motors_[0].rest_dist_;
  double r_rest_teth =
      sqrt(site_size * rest_teth * site_size * rest_teth + r_y_teth * r_y_teth);
  p_diffuse_i_to_teth_rest_.resize(max_neighbs_ + 1);
  p_diffuse_i_fr_teth_rest_.resize(max_neighbs_ + 1);
  p_diffuse_ii_to_both_.resize(max_neighbs_ + 1);
  p_diffuse_ii_to_self_fr_teth_.resize(max_neighbs_ + 1);
  p_diffuse_ii_fr_self_to_teth_.resize(max_neighbs_ + 1);
  p_diffuse_ii_fr_both_.resize(max_neighbs_ + 1);
  for (int n_neighbs = 0; n_neighbs <= max_neighbs_; n_neighbs++) {
    double tot_E = n_neighbs * interaction_energy_;
    // Lambda = 1 & 0 for diffusion
    double weight_neighb = exp(tot_E);
    if (n_neighbs == 2)
      weight_neighb = 0;
    p_diffuse_i_to_teth_rest_[n_neighbs].resize(2 * teth_cutoff + 1);
    p_diffuse_i_fr_teth_rest_[n_neighbs].resize(2 * teth_cutoff + 1);
    p_diffuse_ii_to_both_[n_neighbs].resize(2 * teth_cutoff + 1);
    p_diffuse_ii_to_self_fr_teth_[n_neighbs].resize(2 * teth_cutoff + 1);
    p_diffuse_ii_fr_self_to_teth_[n_neighbs].resize(2 * teth_cutoff + 1);
    p_diffuse_ii_fr_both_[n_neighbs].resize(2 * teth_cutoff + 1);
    for (int x_dub = 0; x_dub <= 2 * teth_cutoff; x_dub++) {
      // Calc x-distances (in nm) for tether
      double r_x_teth = x_dub * site_size / 2;
      double r_x_teth_bck = (x_dub - 1) * site_size / 2;
      double r_x_teth_fwd = (x_dub + 1) * site_size / 2;
      // Calc total r values
      double r_teth = sqrt(r_x_teth * r_x_teth + r_y_teth * r_y_teth);
      double r_teth_bck =
          sqrt(r_x_teth_bck * r_x_teth_bck + r_y_teth * r_y_teth);
      double r_teth_fwd =
          sqrt(r_x_teth_fwd * r_x_teth_fwd + r_y_teth * r_y_teth);
      // Calc tether exts for current dist and stepping to/from rest
      double dr_teth = r_teth - r_0_teth, dr_teth_to, dr_teth_fr;
      if (dr_teth >= 0) {
        dr_teth_to = r_teth_bck - r_0_teth;
        dr_teth_fr = r_teth_fwd - r_0_teth;
      } else {
        dr_teth_to = r_teth_fwd - r_0_teth;
        dr_teth_fr = r_teth_bck - r_0_teth;
      }
      double dU_fr_teth, dU_to_teth;
      if (x_dub == 2 * rest_teth) {
        if (r_0_teth > r_rest_teth) {
          dU_fr_teth = (k_teth_slack / 2) *
                       (dr_teth_fr * dr_teth_fr - dr_teth * dr_teth);
          dU_to_teth = (0.5) * (k_teth_spring * dr_teth_to * dr_teth_to -
                                k_teth_slack * dr_teth * dr_teth);
        } else {
          dU_fr_teth = (k_teth_spring / 2) *
                       (dr_teth_fr * dr_teth_fr - dr_teth * dr_teth);
          dU_to_teth = (0.5) * (k_teth_slack * dr_teth_to * dr_teth_to -
                                k_teth_spring * dr_teth * dr_teth);
        }
      } else if (dr_teth > 0) {
        dU_fr_teth =
            (k_teth_spring / 2) * (dr_teth_fr * dr_teth_fr - dr_teth * dr_teth);
        dU_to_teth =
            (k_teth_spring / 2) * (dr_teth_to * dr_teth_to - dr_teth * dr_teth);
      } else {
        dU_fr_teth =
            (k_teth_slack / 2) * (dr_teth_fr * dr_teth_fr - dr_teth * dr_teth);
        dU_to_teth =
            (k_teth_slack / 2) * (dr_teth_to * dr_teth_to - dr_teth * dr_teth);
      }
      double weight_to_teth = exp(-dU_to_teth / (2 * kbT));
      double weight_fr_teth = exp(-dU_fr_teth / (2 * kbT));
      if (x_dub < 2 * comp_cutoff_ || !parameters_->motors.tethers_active) {
        weight_to_teth = 0;
        weight_fr_teth = 0;
      }
      if (x_dub < 2 * (comp_cutoff_ + 1) || x_dub > 2 * (teth_cutoff - 1)) {
        weight_fr_teth = 0;
      }
      double p_to_teth_i = weight_neighb * weight_to_teth * delta_t / tau;
      double p_fr_teth_i = weight_neighb * weight_fr_teth * delta_t / tau;
      // Input probabilities for stage_i / tethered xlinks
      p_diffuse_i_to_teth_rest_[n_neighbs][x_dub] = p_to_teth_i;
      p_diffuse_i_fr_teth_rest_[n_neighbs][x_dub] = p_fr_teth_i;
      if (p_to_teth_i > 1)
        printf("WARNING: p_diffuse_to_teth_i=%g for 2x=%i\n", p_to_teth_i,
               x_dub);
      if (p_fr_teth_i > 1)
        printf("WARNING: p_diffuse_fr_teth_i=%g for 2x=%i\n", p_fr_teth_i,
               x_dub);
      // Run through x to get probs for stage_ii tethered xlinks
      p_diffuse_ii_to_both_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      p_diffuse_ii_to_self_fr_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      p_diffuse_ii_fr_self_to_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      p_diffuse_ii_fr_both_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      for (int x_dist = 0; x_dist <= dist_cutoff_; x_dist++) {
        double r_x = x_dist * site_size;
        double r_x_to = (x_dist - 1) * site_size;
        double r_x_fr = (x_dist + 1) * site_size;
        double r = sqrt(r_x * r_x + r_y * r_y);
        double r_to = sqrt(r_x_to * r_x_to + r_y * r_y);
        double r_fr = sqrt(r_x_fr * r_x_fr + r_y * r_y);
        // Get extension for current dist and steps to/from rest
        double dr = r - r_0;
        double dr_to = r_to - r_0;
        double dr_fr = r_fr - r_0;
        if (dr >= 0) {
          // Get corresponding changes in potential energy
          double dU_to_rest = (k_spring / 2) * (dr_to * dr_to - dr * dr);
          double dU_fr_rest = (k_spring / 2) * (dr_fr * dr_fr - dr * dr);
          // Weights according to Lanksy et al.
          double weight_to;
          double weight_fr;
          if (x_dist == rest_dist_) {
            weight_to = 0;
            weight_fr = 2 * exp(-dU_fr_rest / (2 * kbT));
          } else if (x_dist == dist_cutoff_) {
            weight_to = exp(-dU_to_rest / (2 * kbT));
            weight_fr = 0;
          } else {
            weight_to = exp(-dU_to_rest / (2 * kbT));
            weight_fr = exp(-dU_fr_rest / (2 * kbT));
          }
          // Convolve these bitches
          double p_to_both =
              weight_neighb * weight_to * weight_to_teth * delta_t / tau;
          double p_to_self_fr_teth =
              weight_neighb * weight_to * weight_fr_teth * delta_t / tau;
          double p_fr_self_to_teth =
              weight_neighb * weight_fr * weight_to_teth * delta_t / tau;
          double p_fr_both =
              weight_neighb * weight_fr * weight_fr_teth * delta_t / tau;
          p_diffuse_ii_to_both_[n_neighbs][x_dub][x_dist] = p_to_both;
          p_diffuse_ii_to_self_fr_teth_[n_neighbs][x_dub][x_dist] =
              p_to_self_fr_teth;
          p_diffuse_ii_fr_self_to_teth_[n_neighbs][x_dub][x_dist] =
              p_fr_self_to_teth;
          p_diffuse_ii_fr_both_[n_neighbs][x_dub][x_dist] = p_fr_both;

          if (p_to_both > 1) {
            printf("WARNING: p_diff_to_both=%g", p_to_both);
            printf(" for 2x=%i, x=%i\n", x_dub, x_dist);
          }
          if (p_to_self_fr_teth > 1) {
            printf("WARNING: p_diff_to_self_fr_teth=%g", p_to_self_fr_teth);
            printf(" for 2x=%i, x=%i\n", x_dub, x_dist);
          }
          if (p_fr_self_to_teth > 1) {
            printf("WARNING: p_diff_fr_self_to_teth=%g", p_fr_self_to_teth);
            printf(" for 2x=%i, x=%i\n", x_dub, x_dist);
          }
          if (p_fr_both > 1) {
            printf("WARNING: p_diff_fr_both=%g", p_fr_both);
            printf(" for 2x=%i, x=%i\n", x_dub, x_dist);
          }
        } else {
          printf("woah mayne. xlink set parameters TWOO \n");
          exit(1);
        }
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
  double k_off = parameters_->xlinks.k_off;
  p_bind_i_.resize(max_neighbs_ + 1);
  p_unbind_i_.resize(max_neighbs_ + 1);
  p_unbind_ii_.resize(max_neighbs_ + 1);
  p_unbind_i_teth_.resize(max_neighbs_ + 1);
  p_unbind_ii_to_teth_.resize(max_neighbs_ + 1);
  p_unbind_ii_fr_teth_.resize(max_neighbs_ + 1);
  for (int n_neighbs(0); n_neighbs <= max_neighbs_; n_neighbs++) {
    double tot_E = n_neighbs * interaction_energy_;
    double weight_neighb = exp(tot_E / 2);
    p_bind_i_[n_neighbs] = (k_on * c_xlink * delta_t) / weight_neighb;
    p_unbind_i_[n_neighbs] = k_off * delta_t * weight_neighb;
    // Rates involving xlink spring only
    p_unbind_ii_[n_neighbs].resize(dist_cutoff_ + 1);
    for (int x = 0; x <= dist_cutoff_; x++) {
      double r_x = x * site_size;
      double r = sqrt(r_y * r_y + r_x * r_x);
      double dr = r - r_0;
      double U_xlink = (k_spring / 2) * dr * dr;
      double unbind_weight = exp(U_xlink / (2 * kbT));
      if (x == rest_dist_)
        unbind_weight = 1;
      p_unbind_ii_[n_neighbs][x] =
          weight_neighb * unbind_weight * k_off * delta_t;
    }
    // Rates involving both xlink & tether spring
    p_unbind_i_teth_[n_neighbs].resize(2 * teth_cutoff + 1);
    p_unbind_ii_to_teth_[n_neighbs].resize(2 * teth_cutoff + 1);
    p_unbind_ii_fr_teth_[n_neighbs].resize(2 * teth_cutoff + 1);
    for (int x_dub = 0; x_dub <= 2 * teth_cutoff; x_dub++) {
      // Calc x-distances (in nm) for tether
      double r_x_teth = x_dub * site_size / 2;
      // Calc total r values
      double r_teth = sqrt(r_x_teth * r_x_teth + r_y_teth * r_y_teth);
      // Calc tether exts for current dist and stepping to/from rest
      double dr_teth = r_teth - r_0_teth;
      double U_at_teth;
      if (dr_teth > 0)
        U_at_teth = (k_teth_spring / 2) * dr_teth * dr_teth;
      else
        U_at_teth = (k_teth_slack / 2) * dr_teth * dr_teth;
      double weight_at_teth = exp(U_at_teth / (2 * kbT));
      if (x_dub < 2 * comp_cutoff_ || !parameters_->motors.tethers_active) {
        weight_at_teth = 0;
      }
      p_unbind_i_teth_[n_neighbs][x_dub] =
          weight_neighb * weight_at_teth * k_off * delta_t;
      if (p_unbind_i_teth_[n_neighbs][x_dub] > 1)
        printf("WARNING: p_unbind_teth (XLINK)=%g for 2x=%i\n",
               p_unbind_i_teth_[n_neighbs][x_dub], x_dub);
      // Run through x to get probs for stage_ii / tethered xlinks
      p_unbind_ii_to_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      p_unbind_ii_fr_teth_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      for (int x = 0; x <= dist_cutoff_; x++) {
        // change in teth extension if 2nd xlink head were to unbind
        double dx_teth = (double)x / 2;
        double dr_x_teth = dx_teth * site_size;
        double r_x_teth_to, r_x_teth_fr;
        if (dr_teth >= 0) {
          r_x_teth_to = r_x_teth - dr_x_teth;
          r_x_teth_fr = r_x_teth + dr_x_teth;
        } else {
          r_x_teth_to = r_x_teth + dr_x_teth;
          r_x_teth_fr = r_x_teth - dr_x_teth;
        }
        double r_teth_to =
            sqrt(r_x_teth_to * r_x_teth_to + r_y_teth * r_y_teth);
        double r_teth_fr =
            sqrt(r_x_teth_fr * r_x_teth_fr + r_y_teth * r_y_teth);
        double dr_teth_to = r_teth_to - r_0_teth;
        double dr_teth_fr = r_teth_fr - r_0_teth;
        double dU_to_teth, dU_fr_teth;
        int x_fr_rest_dub = abs(2 * rest_teth - x_dub);
        int dx_teth_dub = 2 * dx_teth;
        // Check if we're crossing over equil. point of tether
        if (dx_teth_dub > x_fr_rest_dub) {
          if (r_teth < r_0_teth) {
            dU_fr_teth = (k_teth_slack / 2) *
                         (dr_teth_fr * dr_teth_fr - dr_teth * dr_teth);
            dU_to_teth = (0.5) * (k_teth_spring * dr_teth_to * dr_teth_to -
                                  k_teth_slack * dr_teth * dr_teth);
          } else {
            dU_fr_teth = (k_teth_spring / 2) *
                         (dr_teth_fr * dr_teth_fr - dr_teth * dr_teth);
            dU_to_teth = (0.5) * (k_teth_slack * dr_teth_to * dr_teth_to -
                                  k_teth_spring * dr_teth * dr_teth);
          }
        } else if (dr_teth > 0) {
          dU_fr_teth = (k_teth_spring / 2) *
                       (dr_teth_fr * dr_teth_fr - dr_teth * dr_teth);
          dU_to_teth = (k_teth_spring / 2) *
                       (dr_teth_to * dr_teth_to - dr_teth * dr_teth);
        } else {
          dU_fr_teth = (k_teth_slack / 2) *
                       (dr_teth_fr * dr_teth_fr - dr_teth * dr_teth);
          dU_to_teth = (k_teth_slack / 2) *
                       (dr_teth_to * dr_teth_to - dr_teth * dr_teth);
        }
        double weight_to_teth = exp(-dU_to_teth / (2 * kbT));
        double weight_fr_teth = exp(-dU_fr_teth / (2 * kbT));
        // To ensure we don't double-count sites for x = 0,
        // only use unbind_ii_fr_teth and ignore to_teth
        if (x == rest_dist_) {
          weight_to_teth = 0;
          weight_fr_teth = 1;
        }
        if (x_dub < 2 * comp_cutoff_ || !parameters_->motors.tethers_active) {
          weight_to_teth = 0;
          weight_fr_teth = 0;
        }
        if (x_dub < 2 * (comp_cutoff_ + 1) || x_dub > 2 * (teth_cutoff - 1)) {
          weight_fr_teth = 0;
        }
        p_unbind_ii_to_teth_[n_neighbs][x_dub][x] =
            weight_neighb * weight_to_teth * k_off * delta_t;
        p_unbind_ii_fr_teth_[n_neighbs][x_dub][x] =
            weight_neighb * weight_fr_teth * k_off * delta_t;
        if (p_unbind_ii_to_teth_[n_neighbs][x_dub][x] > 1)
          printf("WARNING: p_unbind_to = %g for 2x=%ix, x=%i\n",
                 p_unbind_ii_to_teth_[n_neighbs][x_dub][x], x_dub, x);
        if (p_unbind_ii_fr_teth_[n_neighbs][x_dub][x] > 1)
          printf("WARNING: p_unbind_fr = %g for 2x=%ix, x=%i\n",
                 p_unbind_ii_fr_teth_[n_neighbs][x_dub][x], x_dub, x);
      }
    }
  }
  double k_teth = parameters_->motors.k_tether;
  double k_unteth = parameters_->motors.k_untether;
  p_tether_free_ = k_teth * c_xlink * delta_t;
  p_untether_free_ = k_unteth * delta_t;
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
    for (int x_dist = 0; x_dist <= dist_cutoff_; x_dist++)
      n_bound_ii_[n_neighbs][x_dist] = 0;
    n_bound_i_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    n_bound_ii_teth_same_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    n_bound_ii_teth_oppo_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    for (int x_dub = 0; x_dub <= 2 * teth_cutoff_; x_dub++) {
      n_bound_i_teth_[n_neighbs][x_dub] = 0;
      n_bound_ii_teth_same_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      n_bound_ii_teth_oppo_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      for (int x = 0; x <= dist_cutoff_; x++) {
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
  for (int n_neighbs(0); n_neighbs <= max_neighbs_ + 1; n_neighbs++) {
    bound_i_[n_neighbs].resize(n_xlinks_);
    bound_ii_[n_neighbs].resize(dist_cutoff_ + 1);
    for (int x = 0; x <= dist_cutoff_; x++)
      bound_ii_[n_neighbs][x].resize(n_xlinks_);
    bound_i_teth_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    bound_ii_teth_oppo_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    bound_ii_teth_same_[n_neighbs].resize(2 * teth_cutoff_ + 1);
    for (int x_dub = 0; x_dub <= 2 * teth_cutoff_; x_dub++) {
      bound_i_teth_[n_neighbs][x_dub].resize(n_xlinks_);
      bound_ii_teth_oppo_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      bound_ii_teth_same_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      for (int x = 0; x <= dist_cutoff_; x++) {
        bound_ii_teth_oppo_[n_neighbs][x_dub][x].resize(n_xlinks_);
        bound_ii_teth_same_[n_neighbs][x_dub][x].resize(n_xlinks_);
      }
    }
  }
}

void AssociatedProteinManagement::InitializeEvents() {

  /* *** Serialized & unique index of each KMC event *** */
  int id(0);
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
      if (n_wt > 0)
        return properties_->gsl.SamplePoissonDist(p * n_wt);
      else
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
  //		(*management, id, kmc_code, 'event_name', 'target_pop',
  //		  p_occur, *n_avail, *pop_pool, ran_int, p_dist)
  // Bind_ii does not use n_neighbors; append w/ "_ALL" instead
  /*
  if (parameters_->microtubules.count > 1) {
    events_.emplace_back(this, id++, 20, "bind_ii", "bound_i_ALL",
                         p_bind_ii_base_, &n_bound_i_[max_neighbs_ + 1],
                         &bound_i_[max_neighbs_ + 1], ran_int, poisson_ii);
  }
  for (int n_neighbs(0); n_neighbs <= max_neighbs_; n_neighbs++) {
    // Number of neighbors, N, always appends target pop. type
    std::string N = std::to_string(n_neighbs);
    events_.emplace_back(this, id++, -10, "diff_i_fwd", "bound_i_" + N,
                         p_diffuse_i_fwd_[n_neighbs], &n_bound_i_[n_neighbs],
                         &bound_i_[n_neighbs], ran_int, binomial);
    events_.emplace_back(this, id++, -11, "diff_i_bck", "bound_i_" + N,
                         p_diffuse_i_bck_[n_neighbs], &n_bound_i_[n_neighbs],
                         &bound_i_[n_neighbs], ran_int, binomial);
    events_.emplace_back(
        this, id++, 10, "bind_i", "unocc_" + N, p_bind_i_[n_neighbs],
        &properties_->microtubules.n_unoccupied_xl_[n_neighbs],
        &properties_->microtubules.unoccupied_list_xl_[n_neighbs], ran_int,
        binomial);
    events_.emplace_back(this, id++, 30, "unbind_i", "bound_i_" + N,
                         p_unbind_i_[n_neighbs], &n_bound_i_[n_neighbs],
                         &bound_i_[n_neighbs], ran_int, binomial);
    // Only create doubly-bound events if n_MTs > 1
    if (parameters_->microtubules.count > 1) {
      for (int x(0); x <= dist_cutoff_; x++) {
        events_.emplace_back(this, id++, -20, "diff_ii_to",
                             "bound_ii_" + std::to_string(x) + "_" + N,
                             p_diffuse_ii_to_rest_[n_neighbs][x],
                             &n_bound_ii_[n_neighbs][x],
                             &bound_ii_[n_neighbs][x], ran_int, binomial);
        events_.emplace_back(this, id++, -21, "diff_ii_fr",
                             "bound_ii_" + std::to_string(x) + "_" + N,
                             p_diffuse_ii_fr_rest_[n_neighbs][x],
                             &n_bound_ii_[n_neighbs][x],
                             &bound_ii_[n_neighbs][x], ran_int, binomial);
        events_.emplace_back(this, id++, 40, "*unbind_ii",
                             "bound_ii_" + std::to_string(x) + "_" + N,
                             p_unbind_ii_[n_neighbs][x],
                             &n_bound_ii_[n_neighbs][x],
                             &bound_ii_[n_neighbs][x], ran_int, binomial);
      }
    }
  }
  // If tethering is enabled, add those event populations as well
  if (parameters_->motors.tethers_active) {
    events_.emplace_back(this, id++, 11, "bind_I", "free_teth",
                         p_bind_i_teth_base_, &n_free_teth_, &free_teth_,
                         ran_int, poisson_i_teth);
    events_.emplace_back(this, id++, 21, "bind_II", "bound_I_ALL",
                         p_bind_ii_base_, &n_bound_i_teth_tot_,
                         &bound_i_teth_[max_neighbs_ + 1][0], ran_int,
                         poisson_ii_teth);
    for (int n_neighbs(0); n_neighbs <= max_neighbs_; n_neighbs++) {
      //			printf("n neighbs is %i\n", n_neighbs);
      for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
        std::string bound_I = "bound_I_" + std::to_string(x_dub) + "_" +
                              std::to_string(n_neighbs);
        //				std::cout << bound_I << std::endl;
        events_.emplace_back(this, id++, -30, "diff_I_to", bound_I,
                             p_diffuse_i_to_teth_rest_[n_neighbs][x_dub],
                             &n_bound_i_teth_[n_neighbs][x_dub],
                             &bound_i_teth_[n_neighbs][x_dub], ran_int,
                             binomial);
        events_.emplace_back(this, id++, -31, "diff_I_fr", bound_I,
                             p_diffuse_i_fr_teth_rest_[n_neighbs][x_dub],
                             &n_bound_i_teth_[n_neighbs][x_dub],
                             &bound_i_teth_[n_neighbs][x_dub], ran_int,
                             binomial);
        events_.emplace_back(this, id++, 31, "unbind_I", bound_I,
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
            events_.emplace_back(this, id++, -40, "diff_II_to", bound_II_same,
                                 p_diffuse_ii_to_both_[n_neighbs][x_dub][x],
                                 &n_bound_ii_teth_same_[n_neighbs][x_dub][x],
                                 &bound_ii_teth_same_[n_neighbs][x_dub][x],
                                 ran_int, binomial);
            events_.emplace_back(this, id++, -41, "diff_II_fr", bound_II_same,
                                 p_diffuse_ii_fr_both_[n_neighbs][x_dub][x],
                                 &n_bound_ii_teth_same_[n_neighbs][x_dub][x],
                                 &bound_ii_teth_same_[n_neighbs][x_dub][x],
                                 ran_int, binomial);
            events_.emplace_back(
                this, id++, -50, "diff_II_to_fr", bound_II_oppo,
                p_diffuse_ii_to_self_fr_teth_[n_neighbs][x_dub][x],
                &n_bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                &bound_ii_teth_oppo_[n_neighbs][x_dub][x], ran_int, binomial);
            events_.emplace_back(
                this, id++, -50, "diff_II_fr_to", bound_II_oppo,
                p_diffuse_ii_fr_self_to_teth_[n_neighbs][x_dub][x],
                &n_bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                &bound_ii_teth_oppo_[n_neighbs][x_dub][x], ran_int, binomial);
            events_.emplace_back(this, id++, 41, "unbind_II_to", bound_II_same,
                                 p_unbind_ii_to_teth_[n_neighbs][x_dub][x],
                                 &n_bound_ii_teth_same_[n_neighbs][x_dub][x],
                                 &bound_ii_teth_same_[n_neighbs][x_dub][x],
                                 ran_int, binomial);
            events_.emplace_back(this, id++, 42, "unbind_II_fr", bound_II_oppo,
                                 p_unbind_ii_fr_teth_[n_neighbs][x_dub][x],
                                 &n_bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                                 &bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                                 ran_int, binomial);
          }
        }
      }
    }
    events_.emplace_back(
        this, id++, 50, "tether_free", "unteth_mots", p_tether_free_,
        &properties_->kinesin4.n_bound_untethered_,
        &properties_->kinesin4.bound_untethered_, ran_int, binomial);
    events_.emplace_back(this, id++, 60, "untether_free", "free_teth",
                         p_untether_free_, &n_free_teth_, &free_teth_, ran_int,
                         binomial);
  }
  */
  /* ** Segregate events_ into IDs_by_pop_ based on target pop. ** */
  /*
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
    std::string target_pop = events_[i_entry].targets_[0];
    // bound_ii entries go thru a primary AND secondary correction
    // (first over x & n_neighbs, then over x for all n_neighbs)
    if (target_pop.length() >= bound_ii.length()) {
      if (target_pop.substr(0, bound_ii.length()) == bound_ii) {
        // std::cout << target_pop;
        // printf(" will be corrected in both scans; root is ");
        int root_length = bound_ii.length() + 1;
        std::string root = target_pop.substr(0, root_length);
        // std::cout << root;
        // Check to see if we have recorded this root before
        bool new_root(true);
        for (int i_root(0); i_root < n_roots; i_root++) {
          std::string recorded = roots[i_root];
          if (root == recorded)
            new_root = false;
        }
        if (new_root) {
          //					printf(" - NEW!\n");
          roots[n_roots] = root;
          n_roots++;
        }
        //				else printf("\n");
      }
    }
    // 'ALL' flag at end -> secondary correction only
    if (target_pop.substr(target_pop.length() - 3) == "ALL") {
      //			std::cout << target_pop;
      //			printf(" will be corrected in secondary scan
      // only; root is ");
      std::string root = target_pop.substr(0, target_pop.length() - 3);
      //			std::cout << root << std::endl;
      roots[n_roots] = root;
      n_roots++;
    }
    // No 'ALL' flag -> primary scan only
    else {
      bool new_pop(true); // Assume population type is new
      int pop_index(0);   // 1st index of this pop. in entries arry
      for (int i_pop(0); i_pop < n_pops; i_pop++) {
        std::string record = events_[entries[i_pop][0]].targets_[0];
        // If type matches a recorded type; pop. isn't new
        if (target_pop == record) {
          new_pop = false;
          pop_index = i_pop;
        }
      }
      // If indeed a new pop., record it in a new row
      if (new_pop) {
        entries[n_pops][0] = events_[i_entry].id_;
        n_entries[n_pops] = 1;
        n_pops++;
      }
      // Otherwise, add entry to row of already-recorded pop.
      else {
        int entry_no = n_entries[pop_index];
        entries[pop_index][entry_no] = events_[i_entry].id_;
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
    std::string target_pop = events_[entries[i_pop][0]].targets_[0];
    printf("KMC Pop #%i - ", i_pop);
    std::cout << target_pop;
    printf(":");
    // 2nd index corresponds to entry of events that target this pop.
    IDs_by_pop_[i_pop].resize(n_entries[i_pop]);
    for (int i_entry(0); i_entry < n_entries[i_pop]; i_entry++) {
      IDs_by_pop_[i_pop][i_entry] = entries[i_pop][i_entry];
      //			printf(" %i,", entries[i_pop][i_entry]);
    }
    //		printf("\n");
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
  // Finally, scan through secondary pops. & place them into eventIDs_by_root_
  eventIDs_by_root_.resize(n_roots);
  n_avail_by_root_.resize(n_roots);
  for (int i_root(0); i_root < n_roots; i_root++) {
    std::string root = roots[i_root];
    printf("KMC ROOT #%i - ", i_root);
    std::cout << root;
    printf(":");
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
    eventIDs_by_root_[i_root].resize(n_sec_entries[i_root]);
    for (int i_entry(0); i_entry < n_sec_entries[i_root]; i_entry++) {
      int entry_ID = sec_entries[i_root][i_entry];
      // If target event has prefix of "*", make IDs negative to
      // flag that entries are coupled & n_avail should be halved
      bool coupled(false);
      if (events_[entry_ID].name_.substr(0, 1) == "*") {
        coupled = true;
      }
      if (coupled)
        eventIDs_by_root_[i_root][i_entry] = -1 * entry_ID;
      else
        eventIDs_by_root_[i_root][i_entry] = entry_ID;
      //			printf(" %i,",
      // eventIDs_by_root_[i_root][i_entry]);
    }
    //		printf("\n");
  }
  */
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
    for (int i_neighb = 0; i_neighb < n_neighbors; i_neighb++) {
      Tubulin *site = head->xlink_->neighbor_sites_[i_neighb];
      double weight = head->xlink_->GetBindingWeight_II(site);
      weights_summed += weight;
    }
  }
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
  return weights_summed;
}

AssociatedProtein::Monomer *
AssociatedProteinManagement::CheckScratchFor(std::string pop) {

  //	printf("Scratch was checked for ");
  //	std::cout << pop << std::endl;
  POP_T *head(nullptr);
  for (int i_entry(0); i_entry < n_scratched_; i_entry++) {
    //		printf("   entry #%i: ", i_entry);
    //		std::cout << scratch_[i_entry]->state_ << std::endl;
    if (scratch_[i_entry]->state_ == pop) {
      head = scratch_[i_entry];
      break;
    }
  }
  return head;
}

void AssociatedProteinManagement::SaveToScratch(POP_T *head) {

  head->xlink_->is_outdated_ = true;
  scratch_[n_scratched_] = head;
  n_scratched_++;
  int n_neighbs = head->GetPRC1NeighbCount();
  if (n_neighbs == 1) {
    POP_T *neighb;
    int i_site = head->site_->index_;
    int n_sites = head->site_->mt_->n_sites_;
    if (i_site == 0)
      neighb = head->site_->mt_->lattice_[i_site + 1].xlink_head_;
    else if (i_site == n_sites - 1)
      neighb = head->site_->mt_->lattice_[i_site - 1].xlink_head_;
    else if (head->site_->mt_->lattice_[i_site + 1].xlink_head_ != nullptr)
      neighb = head->site_->mt_->lattice_[i_site + 1].xlink_head_;
    else if (head->site_->mt_->lattice_[i_site - 1].xlink_head_ != nullptr)
      neighb = head->site_->mt_->lattice_[i_site - 1].xlink_head_;
    neighb->xlink_->is_outdated_ = true;
    scratch_[n_scratched_] = neighb;
    n_scratched_++;
  } else if (n_neighbs == 2) {
    POP_T *neighb_one, *neighb_two;
    int i_site = head->site_->index_;
    neighb_one = head->site_->mt_->lattice_[i_site + 1].xlink_head_;
    scratch_[n_scratched_] = neighb_one;
    n_scratched_++;
    neighb_one->xlink_->is_outdated_ = true;
    neighb_two = head->site_->mt_->lattice_[i_site - 1].xlink_head_;
    scratch_[n_scratched_] = neighb_two;
    n_scratched_++;
    neighb_two->xlink_->is_outdated_ = true;
  }
}

void AssociatedProteinManagement::Update_Relay(std::string event,
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
    } else if (target_pop.substr(0, 12) == "bound_I_") {
      Update_Bound_I_Teth();
      if (event.substr(0, 7) == "bind_II") {
        properties_->microtubules.UpdateUnoccupied();
        Update_Bind_II_Teth_Candidate();
      }
    } else if (target_pop.substr(0, 13) == "bound_II_") {
      Update_Bound_II_Teth();
    } else if (target_pop.substr(0, 6) == "unteth") {
      properties_->kinesin4.Update_Bound_Unteth();
    }
  } else {
    printf("somethin SPOOKY in update_relay (XLINKS) - EXITING\n");
    std::cout << target_pop << std::endl;
    exit(1);
  }
}

void AssociatedProteinManagement::Update_All_Lists() {

  properties_->microtubules.UpdateUnoccupied();
  Update_Bound_I();
  if (parameters_->microtubules.count > 1) {
    Update_Bound_II();
  }
  if (parameters_->motors.tethers_active) {
    Update_Free_Teth();
    Update_Bound_Unteth();
    Update_Bound_I_Teth();
    if (parameters_->microtubules.count > 1) {
      Update_Bound_II_Teth();
    }
    properties_->kinesin4.Update_Bound_Unteth();
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
        int x_dub = xlink->motor_->x_dist_doubled_;
        AssociatedProtein::Monomer *head = xlink->GetActiveHead();
        int n_neighbs = head->site_->GetPRC1NeighborCount();
        // printf("FOUND %i NEBS\n", n_neighbs);
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

void AssociatedProteinManagement::Update_Bound_II() {

  for (int n_neighbs(0); n_neighbs <= max_neighbs_ + 1; n_neighbs++) {
    for (int x = 0; x <= dist_cutoff_; x++) {
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

void AssociatedProteinManagement::Update_Bound_II_Teth() {

  for (int n_neighbs(0); n_neighbs <= max_neighbs_ + 1; n_neighbs++) {
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
      for (int x(0); x <= dist_cutoff_; x++) {
        n_bound_ii_teth_same_[n_neighbs][x_dub][x] = 0;
        n_bound_ii_teth_oppo_[n_neighbs][x_dub][x] = 0;
      }
    }
  }
  for (int i_xlink = 0; i_xlink < n_active_; i_xlink++) {
    AssociatedProtein *xlink = active_[i_xlink];
    if (xlink->heads_active_ == 2 && xlink->tethered_ == true) {
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
        break;
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
        break;
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
      weight_total += weight_local[i_entry];
      i_entry++;
    }
  }
  if (weight_total > 0.0) {
    // Normalize local weights to get relative probabilities
    // Using these relative probs, randomly pick an entry
    i_entry = 0;
    bool failed(true);
    double p_cum(0);
    double ran = properties_->gsl.GetRanProb();
    int chosen_x_dub = -1;
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
      n_bound_i_teth_[all_neighbs][x_dub] = 0;
      for (int i(0); i < n_bound_i_teth_[all_neighbs][x_dub]; i++) {
        p_cum += weight_local[i_entry] / weight_total;
        if (p_cum > ran) {
          failed = false;
          bound_i_teth_[all_neighbs][x_dub][0] =
              bound_i_teth_[all_neighbs][x_dub][i];
          n_bound_i_teth_[all_neighbs][x_dub] = 1;
          chosen_x_dub = x_dub;
          break;
        }
      }
    }
    if (failed) {
      printf("nope in update_bind_ii_teth_candidate; EXIT \n");
      exit(1);
    }
  } else {
    for (int x_dub(2 * comp_cutoff_); x_dub <= 2 * teth_cutoff_; x_dub++) {
      n_bound_i_teth_[all_neighbs][x_dub] = 0;
    }
  }
}

void AssociatedProteinManagement::Run_KMC() {

  if (parameters_->xlinks.c_bulk == 0) {
    return;
  }
  sys_time start = sys_clock::now();
  Update_All_Lists();
  sys_time finish_list = sys_clock::now();
  properties_->wallace.t_xlinks_[1] += (finish_list - start).count();
  Refresh_Populations();
  sys_time finish_pops = sys_clock::now();
  properties_->wallace.t_xlinks_[2] += (finish_pops - finish_list).count();
  Generate_Execution_Sequence();
  sys_time finish_seq = sys_clock::now();
  properties_->wallace.t_xlinks_[3] += (finish_seq - finish_pops).count();
  //	printf("Start of xlink KMC cycle\n");
  //	if(IDs_to_exe_.size() > 0) printf("%lu EVENTS\n", IDs_to_exe_.size());
  for (int i_event = 0; i_event < IDs_to_exe_.size(); i_event++) {
    events_[IDs_to_exe_[i_event]].Execute();
  }
  sys_time finish_all = sys_clock::now();
  properties_->wallace.t_xlinks_[4] += (finish_all - finish_seq).count();
  properties_->wallace.t_xlinks_[0] += (finish_all - start).count();
}

void AssociatedProteinManagement::Refresh_Populations() {

  // Clear scratch list (setting n_scratch = 0 is all that is needed)
  n_scratched_ = 0;
  // If we have more than 1 MT, update all extensions
  if (parameters_->microtubules.count > 1) {
    for (int i_entry(0); i_entry < n_active_; i_entry++) {
      active_[i_entry]->UpdateExtension();
      if (active_[i_entry]->tethered_) {
        active_[i_entry]->motor_->UpdateExtension();
      }
    }
  }
  // Run through all active entries and update labels if necessary
  for (int i_entry(0); i_entry < n_active_; i_entry++) {
    AssociatedProtein *entry = active_[i_entry];
    if (entry->is_outdated_) {
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
            int X_DUB = entry->motor_->x_dist_doubled_;
            x_dub = std::to_string(X_DUB);
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
        if (entry->tethered_) {
          if (entry->motor_->heads_active_ > 0) {
            int X_DUB = entry->motor_->x_dist_doubled_;
            x_dub = std::to_string(X_DUB);
            root = std::string("bound_II_") + x + std::string("_") + x_dub;
          } else {
            root = std::string("bound_ii_") + x;
          }
        } else {
          root = std::string("bound_ii_") + x;
        }
        // Head one
        n_neighbs = entry->head_one_.GetPRC1NeighbCount();
        neighbs = std::to_string(n_neighbs);
        state = root + std::string("_") + neighbs;
        entry->head_one_.state_ = state;
        // Head two
        n_neighbs = entry->head_two_.GetPRC1NeighbCount();
        neighbs = std::to_string(n_neighbs);
        state = root + std::string("_") + neighbs;
        entry->head_two_.state_ = state;
      }
      //			printf("STATE IS ");
      //			std::cout << state << std::endl;
    }
  }
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
        pre_array[i_array] = 0; // events_[i_event].id_;
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
        // printf("uhhhh ?? - %i\n",
        // n_events_loc);
        std::cout << events_[IDs_by_pop_[i_pop][0]].name_ << std::endl;
        // std::cout << events_[IDs_by_pop_[i_pop][0]].targets_[0] << std::endl;
      }
      while (n_events_loc > n_avail_loc) {
        double p_cum = 0;
        double ran = properties_->gsl.GetRanProb();
        for (int i_entry(0); i_entry < n_competitors; i_entry++) {
          int i_competitor = IDs_by_pop_[i_pop][i_entry];
          p_cum += (events_[i_competitor].n_expected_ *
                    events_[i_competitor].p_occur_) /
                   p_tot;
          if (p_cum >= ran && events_[i_competitor].n_expected_ > 0) {
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
  for (int i_root(0); i_root < eventIDs_by_root_.size(); i_root++) {
    int n_competitors = eventIDs_by_root_[i_root].size();
    if (n_competitors > 1) {
      double p_tot(0);
      bool coupled(false);
      int n_events_loc(0);
      for (int i_entry(0); i_entry < n_competitors; i_entry++) {
        int i_competitor = eventIDs_by_root_[i_root][i_entry];
        if (i_competitor < 0) {
          coupled = true;
          i_competitor = abs(i_competitor);
        }
        n_events_loc += events_[i_competitor].n_expected_;
        p_tot += (events_[i_competitor].n_expected_ *
                  events_[i_competitor].p_occur_);
      }
      int n_avail_loc = *n_avail_by_root_[i_root];
      if (coupled) {
        if (n_avail_loc > 0) {
          std::cout << events_[abs(eventIDs_by_root_[i_root][0])].name_
                    << std::endl;
          // std::cout << events_[abs(eventIDs_by_root_[i_root][0])].targets_[0]
          //           << std::endl;
        }
        n_avail_loc /= 2;
      }
      if (n_avail_loc == 0 && n_events_loc > 0) {
        printf("uhhhh TWO ?? - %i\n", n_events_loc);
        std::cout << events_[abs(eventIDs_by_root_[i_root][0])].name_
                  << std::endl;
        // std::cout << events_[abs(eventIDs_by_root_[i_root][0])].targets_[0]
        //           << std::endl;
      }
      while (n_events_loc > n_avail_loc) {
        double p_cum = 0;
        double ran = properties_->gsl.GetRanProb();
        for (int i_entry(0); i_entry < n_competitors; i_entry++) {
          int i_competitor = abs(eventIDs_by_root_[i_root][i_entry]);
          p_cum += (events_[i_competitor].n_expected_ *
                    events_[i_competitor].p_occur_) /
                   p_tot;
          if (p_cum >= ran && events_[i_competitor].n_expected_ > 0) {
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
  return n_events_tot;
}

void AssociatedProteinManagement::Execution_Relay(ENTRY_T target, int code) {

  bool verbose(false);
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
        printf("tether_free\n");
        //					std::cout << head->state_ <<
        // std::endl;
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

  Tubulin *site = head->site_;
  int dx(0);
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
  } else {
    // if not tethered, the above is not a concern
    dx = dir * head->GetDirectionToRest();
  }
  // Cannot step off them MTs
  if (dx != 0 && !(site->index_ == 0 && dx == -1) &&
      !(site->index_ == site->mt_->n_sites_ - 1 && dx == 1)) {
    Tubulin *old_site = site;
    Tubulin *new_site = &site->mt_->lattice_[site->index_ + dx];
    if (!new_site->occupied_) {
      SaveToScratch(head);
      if (head->xlink_->heads_active_ == 2) {
        SaveToScratch(head->GetOtherHead());
      }
      if (!(new_site->index_ == 0 && dx == -1) &&
          !(new_site->index_ == site->mt_->n_sites_ - 1 && dx == 1)) {
        Tubulin *new_neighb = &site->mt_->lattice_[new_site->index_ + dx];
        if (new_neighb->xlink_head_ != nullptr) {
          SaveToScratch(new_neighb->xlink_head_);
        }
      }
      old_site->xlink_head_ = nullptr;
      old_site->occupied_ = false;
      new_site->xlink_head_ = head;
      new_site->occupied_ = true;
      head->site_ = new_site;
    }
  }
}

void AssociatedProteinManagement::Bind_I(SITE_T *site) {

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
  // For binding (and only binding), save to scratch AFTER event
  SaveToScratch(&xlink->head_one_);
}

void AssociatedProteinManagement::Bind_II(POP_T *bound_head) {

  // XXX do i like this? idk think about it
  SaveToScratch(bound_head);
  POP_T *second_head = bound_head->GetOtherHead();
  Tubulin *site(nullptr);
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
    // Place xlink's second head onto site
    site->xlink_head_ = second_head;
    site->occupied_ = true;
    // Update xlink details
    second_head->site_ = site;
    second_head->xlink_->heads_active_++;
    SaveToScratch(second_head);
  } else {
    printf("failed to bind_ii in XLINK MGMT\n");
  }
}

void AssociatedProteinManagement::Unbind_I(POP_T *head) {

  SaveToScratch(head);
  // Get site
  Tubulin *site = head->site_;
  // Remove xlink from site
  site->xlink_head_ = nullptr;
  site->occupied_ = false;
  // Update xlink details
  head->site_ = nullptr;
  head->xlink_->heads_active_--;
  // Check to see if xlink has become inactive
  bool now_inactive(false);
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
    else
      now_inactive = true;
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

  SaveToScratch(head);
  SaveToScratch(head->GetOtherHead());
  // Get site
  Tubulin *site = head->site_;
  // Remove xlink from site
  site->xlink_head_ = nullptr;
  site->occupied_ = false;
  // Update xlink details
  head->site_ = nullptr;
  head->xlink_->heads_active_--;
}

void AssociatedProteinManagement::Bind_I_Teth(POP_T *satellite_head) {

  // Randomly choose a weighted neighbor site
  Tubulin *site = satellite_head->xlink_->GetWeightedSite_Bind_I_Teth();
  if (site != nullptr) {
    // Update site details
    site->xlink_head_ = satellite_head;
    site->occupied_ = true;
    // Update xlink details
    satellite_head->xlink_->heads_active_++;
    satellite_head->site_ = site;
    SaveToScratch(satellite_head);
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