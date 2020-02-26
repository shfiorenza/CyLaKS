#include "microtubule_management.h"
#include "master_header.h"

MicrotubuleManagement::MicrotubuleManagement() {}

void MicrotubuleManagement::Initialize(system_parameters *parameters,
                                       system_properties *properties) {

  parameters_ = parameters;
  properties_ = properties;

  SetParameters();
  GenerateMicrotubules();
  UpdateNeighbors();
  UpdateUnoccupied();
}

void MicrotubuleManagement::SetParameters() {

  for (int i_mt{0}; i_mt < parameters_->microtubules.count; i_mt++) {
    n_sites_tot_ += parameters_->microtubules.length[i_mt];
  }
  n_unocc_xlink_.resize(max_neighbs_xlink_ + 1);
  unocc_xlink_.resize(max_neighbs_xlink_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_xlink_; n_neighbs++) {
    n_unocc_xlink_[n_neighbs] = 0;
    unocc_xlink_[n_neighbs].resize(n_sites_tot_);
  }
  n_unocc_motor_.resize(n_affinities_);
  unocc_motor_.resize(n_affinities_);
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    n_unocc_motor_[i_aff].resize(max_neighbs_motor_ + 1);
    unocc_motor_[i_aff].resize(max_neighbs_motor_ + 1);
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_motor_; n_neighbs++) {
      n_unocc_motor_[i_aff][n_neighbs] = 0;
      unocc_motor_[i_aff][n_neighbs].resize(n_sites_tot_);
    }
  }
}

void MicrotubuleManagement::GenerateMicrotubules() {

  mt_list_.resize(parameters_->microtubules.count);
  for (int i_mt{0}; i_mt < parameters_->microtubules.count; i_mt++) {
    mt_list_[i_mt].Initialize(parameters_, properties_, i_mt);
  }
}

void MicrotubuleManagement::FlagForUpdate() { lists_up_to_date_ = false; }

void MicrotubuleManagement::UpdateNeighbors() {

  int n_mts = parameters_->microtubules.count;
  if (n_mts > 1) {
    for (int i_mt = 0; i_mt < n_mts; i_mt += 2) {
      Microtubule *mt = &mt_list_[i_mt];
      Microtubule *mt_adj = &mt_list_[i_mt + 1];
      mt->neighbor_ = mt_adj;
      mt_adj->neighbor_ = mt;
    }
  } else
    mt_list_[0].neighbor_ = nullptr;
}

void MicrotubuleManagement::UpdateUnoccupied() {

  if (lists_up_to_date_) {
    return;
  }
  lists_up_to_date_ = true;
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_xlink_; n_neighbs++) {
    n_unocc_xlink_[n_neighbs] = 0;
  }
  for (int i_aff{0}; i_aff < n_affinities_; i_aff++) {
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_motor_; n_neighbs++) {
      n_unocc_motor_[i_aff][n_neighbs] = 0;
    }
  }
  for (int i_mt{0}; i_mt < parameters_->microtubules.count; i_mt++) {
    for (int i_site{0}; i_site < mt_list_[i_mt].n_sites_; i_site++) {
      Tubulin *site{&mt_list_[i_mt].lattice_[i_site]};
      site->UpdateAffinity();
      if (site->occupied_) {
        continue;
      }
      int n_neighbs_xl{site->GetPRC1NeighborCount()};
      unocc_xlink_[n_neighbs_xl][n_unocc_xlink_[n_neighbs_xl]++] = site;
      int n_neighbs_mot{site->GetKif4ANeighborCount()};
      int tub_aff{site->affinity_};
      int index{n_unocc_motor_[tub_aff][n_neighbs_mot]++};
      unocc_motor_[tub_aff][n_neighbs_mot][index] = site;
    }
  }
}

void MicrotubuleManagement::RunDiffusion() {

  // If diffusion is disabled for this simulation, immediately return
  if (!parameters_->microtubules.diffusion_on) {
    return;
  }
  int n_mts{parameters_->microtubules.count};
  bool mts_inactive{true};
  int current_step{properties_->current_step_};
  double delta_t{parameters_->delta_t};
  double current_time{current_step * delta_t};
  // Check that at least one microtubule is active (not immobilized)
  for (int i_mt{0}; i_mt < n_mts; i_mt++) {
    if (
        // current_step % mt_list_[i_mt].steps_per_iteration_ == 0 and
        current_time >= parameters_->microtubules.immobile_until[i_mt]) {
      mts_inactive = false;
    }
  }
  // If no microtubules are active, return
  if (mts_inactive) {
    return;
  }
  double kbT{parameters_->kbT};
  double site_size{parameters_->microtubules.site_size};
  // Sum up all forces exerted on each microtubule
  double forces_summed[n_mts];
  properties_->prc1.Update_Extensions();
  properties_->kinesin4.Update_Extensions();
  for (int i_mt = 0; i_mt < n_mts; i_mt++) {
    forces_summed[i_mt] = mt_list_[i_mt].GetNetForce();
    forces_summed[i_mt] += parameters_->microtubules.applied_force;
  }
  // Check for symmetry
  double delta{forces_summed[0] + forces_summed[1]};
  double tolerance{0.0001};
  if (delta > tolerance) {
    printf("Error in RunMTDiffusion\n");
    printf(" *** EXITING *** \n");
    exit(1);
  }
  // Calculate the displacement of each microtubule (in n_sites)
  int displacement[n_mts];
  for (int i_mt{0}; i_mt < n_mts; i_mt++) {
    // If microtubule is still immobilized, set displacement to 0 and continue
    if (current_time < parameters_->microtubules.immobile_until[i_mt]) {
      //  or current_step % mt_list_[i_mt].steps_per_iteration_ != 0) {
      displacement[i_mt] = 0;
      continue;
    }
    // Calculate instanteous velocity due to forces using drag coefficient
    double velocity{forces_summed[i_mt] / mt_list_[i_mt].gamma_};
    double dx_mean{velocity * delta_t};
    // if (current_step % mt_list_[i_mt].steps_per_iteration_ != 0) {
    double n_sites{dx_mean / site_size};
    int n_sites_whole{(int)n_sites};
    displacement[i_mt] = n_sites_whole;
    double leftover{fabs(n_sites - n_sites_whole)};
    double ran{properties_->gsl.GetRanProb()};
    if (ran < leftover) {
      if (n_sites > 0) {
        displacement[i_mt]++;
      } else {
        displacement[i_mt]--;
      }
    }
    continue;
    // }
    /*
    // Add gaussan noise, meant to represent thermal motion
    double delta_t_eff{delta_t * mt_list_[i_mt].steps_per_iteration_};
    double dx_sigma{sqrt(2 * kbT * delta_t_eff / mt_list_[i_mt].gamma_)};
    // Convert mean and sigma from nm to n_sites
    dx_mean /= site_size;
    dx_sigma /= site_size;
    int sigma_cutoff{3};
    // Construct a discrete gaussian cumulative distribution table
    int range{(int)ceil(fabs(dx_mean) + dx_sigma * sigma_cutoff)};
    double discrete_cdf[2 * range + 1];
    double p_cum{0.0};
    for (int i_bin{0}; i_bin < 2 * range + 1; i_bin++) {
      // dx goes from -range to +range
      int dx{i_bin - range};
      p_cum += properties_->gsl.GetGaussianPDF(dx - dx_mean, dx_sigma);
      discrete_cdf[i_bin] = p_cum;
    }
    // Normalize table so that last entry is 1.0
    for (int i_bin{0}; i_bin < 2 * range + 1; i_bin++) {
      discrete_cdf[i_bin] /= discrete_cdf[2 * range];
    }
    // Roll a random number to determine which dx to choose
    double ran{properties_->gsl.GetRanProb()};
    for (int i_bin{0}; i_bin < 2 * range + 1; i_bin++) {
      if (ran < discrete_cdf[i_bin]) {
        displacement[i_mt] = i_bin - range;
        break;
      }
    }
    */
  }
  /*  Run through MT list and update displacementsi */
  for (int i_mt = 0; i_mt < n_mts; i_mt++) {
    if (displacement[i_mt] == 0) {
      continue;
    }
    mt_list_[i_mt].coord_ += displacement[i_mt];
    properties_->prc1.FlagForUpdate();
    properties_->kinesin4.FlagForUpdate();
  }
}