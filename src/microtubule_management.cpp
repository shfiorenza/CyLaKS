#include "microtubule_management.h"
#include "master_header.h"

MicrotubuleManagement::MicrotubuleManagement() {}

void MicrotubuleManagement::Initialize(system_parameters *parameters,
                                       system_properties *properties) {

  wally_ = &properties->wallace;
  parameters_ = parameters;
  properties_ = properties;
  // K_ or AP_MGMT will initialize the MT environment for each specific test
  if (wally_->test_mode_ != nullptr) {
    return;
  }
  SetParameters();
  GenerateMicrotubules();
  InitializeLists();
  UpdateNeighbors();
}

void MicrotubuleManagement::InitializeTestEnvironment() {

  if (strcmp(properties_->wallace.test_mode_, "motor_lattice_bind") == 0) {
    int lattice_cutoff{properties_->kinesin4.lattice_cutoff_};
    parameters_->microtubules.count = 1;
    parameters_->microtubules.length[0] = 2 * lattice_cutoff + 1;
    parameters_->microtubules.diffusion_on = false;
  }
  if (strcmp(properties_->wallace.test_mode_, "motor_lattice_step") == 0) {
    parameters_->microtubules.count = 1;
    parameters_->microtubules.length[0] = 50;
    parameters_->microtubules.diffusion_on = false;
  }
  if (strcmp(properties_->wallace.test_mode_, "xlink_bind_ii") == 0) {
    int xlink_cutoff{properties_->prc1.dist_cutoff_};
    parameters_->microtubules.count = 2;
    parameters_->microtubules.length[0] = 2 * xlink_cutoff + 1;
    parameters_->microtubules.length[1] = 2 * xlink_cutoff + 1;
    parameters_->microtubules.diffusion_on = false;
  }
  wally_->Log("   MT count set to %i\n", parameters_->microtubules.count);
  for (int i_mt{0}; i_mt < parameters_->microtubules.count; i_mt++) {
    wally_->Log("   MT #%i length set to %i\n", i_mt,
                parameters_->microtubules.length[i_mt]);
  }
  SetParameters();
  GenerateMicrotubules();
  InitializeLists();
  UpdateNeighbors();
  UpdateUnoccupied();
}

void MicrotubuleManagement::SetParameters() {

  for (int i_mt{0}; i_mt < parameters_->microtubules.count; i_mt++) {
    n_sites_tot_ += parameters_->microtubules.length[i_mt];
  }
}

void MicrotubuleManagement::GenerateMicrotubules() {

  mt_list_.resize(parameters_->microtubules.count);
  for (int i_mt{0}; i_mt < parameters_->microtubules.count; i_mt++) {
    mt_list_[i_mt].Initialize(parameters_, properties_, i_mt);
  }
}

void MicrotubuleManagement::InitializeLists() {

  unocc_motor_.resize(n_sites_tot_);
  n_unocc_xlink_.resize(max_neighbs_xlink_ + 1);
  unocc_xlink_.resize(max_neighbs_xlink_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_xlink_; n_neighbs++) {
    n_unocc_xlink_[n_neighbs] = 0;
    unocc_xlink_[n_neighbs].resize(n_sites_tot_);
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
  n_unocc_motor_ = 0;
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_xlink_; n_neighbs++) {
    n_unocc_xlink_[n_neighbs] = 0;
  }
  KinesinManagement *motors{&properties_->kinesin4};
  for (int i_mt{0}; i_mt < mt_list_.size(); i_mt++) {
    for (int i_site{0}; i_site < mt_list_[i_mt].n_sites_; i_site++) {
      Tubulin *site{&mt_list_[i_mt].lattice_[i_site]};
      // Uncomment for motor_lattice_step w/ delta > 0
      if (wally_->test_mode_ == nullptr or properties_->current_step_ == 0) {
        int n_neighbs_mot{site->GetKif4ANeighborCount()};
        site->weight_bind_ = motors->weight_neighbs_bind_[n_neighbs_mot];
        site->weight_unbind_ = motors->weight_neighbs_unbind_[n_neighbs_mot];
      }
      if (site->occupied_) {
        continue;
      }
      unocc_motor_[n_unocc_motor_++] = site;
      int n_neighbs_xl{site->GetPRC1NeighborCount()};
      unocc_xlink_[n_neighbs_xl][n_unocc_xlink_[n_neighbs_xl]++] = site;
    }
  }
  motors->UpdateLatticeWeights();
}

void MicrotubuleManagement::RunDiffusion() {

  // If diffusion is disabled for this simulation, immediately return
  if (!parameters_->microtubules.diffusion_on) {
    return;
  }
  bool mts_inactive{true};
  int n_mts{parameters_->microtubules.count};
  unsigned long current_step{properties_->current_step_};
  double delta_t{parameters_->delta_t};
  double current_time{current_step * delta_t};
  // Check that at least one microtubule is active (not immobilized)
  for (int i_mt{0}; i_mt < n_mts; i_mt++) {
    if (current_step % mt_list_[i_mt].steps_per_iteration_ == 0 and
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
  properties_->prc1.UpdateExtensions();
  properties_->kinesin4.UpdateExtensions();
  for (int i_mt = 0; i_mt < n_mts; i_mt++) {
    forces_summed[i_mt] = mt_list_[i_mt].GetNetForce();
    forces_summed[i_mt] += parameters_->microtubules.applied_force;
  }
  // Check for symmetry
  double delta{forces_summed[0] + forces_summed[1]};
  double tolerance{0.0001};
  if (delta > tolerance) {
    properties_->wallace.ErrorExit("MT_MGMT::RunDiffusion() [1]");
  }
  // Calculate the displacement of each microtubule (in n_sites)
  int displacement[n_mts];
  for (int i_mt{0}; i_mt < n_mts; i_mt++) {
    // If microtubule is still immobilized, set displacement to 0 and continue
    if (current_time < parameters_->microtubules.immobile_until[i_mt] or
        current_step % mt_list_[i_mt].steps_per_iteration_ != 0) {
      displacement[i_mt] = 0;
      continue;
    }
    // Effective dt for this microtubule due to steps_per_iteration_
    double delta_t_eff{delta_t * mt_list_[i_mt].steps_per_iteration_};
    // Calculate instanteous velocity due to forces using drag coefficient
    double velocity{forces_summed[i_mt] / mt_list_[i_mt].gamma_};
    // In order to properly model Brownian dynamics, we use the displacement of
    // this velocity as the mean of a gaussian distribution. The width of this
    // gaussian represents thermal motion, i.e., diffusion of the MTs
    double dx_mean{velocity * delta_t_eff};
    // double dx_sigma{sqrt(2 * kbT * delta_t_eff / mt_list_[i_mt].gamma_)};
    double dx_sigma{0.0};
    if (dx_sigma == 0.0) {
      double dx_sites{dx_mean / site_size};
      int dx_sites_whole{(int)dx_sites};
      double dx_sites_leftover{fabs(dx_sites - dx_sites_whole)};
      double ran{properties_->gsl.GetRanProb()};
      if (ran < dx_sites_leftover) {
        if (dx_sites > 0) {
          dx_sites_whole++;
        } else {
          dx_sites_whole++;
        }
      }
      displacement[i_mt] = dx_sites_whole;
      continue;
    }
    // Convert mean and sigma from nm to n_sites
    dx_mean /= site_size;
    dx_sigma /= site_size;
    // Construct a discrete gaussian cumulative distribution table
    int sigma_cutoff{3};
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
  }
  // Run through MT list and update displacements
  for (int i_mt = 0; i_mt < n_mts; i_mt++) {
    // printf("dx = %i for mt #%i\n", displacement[i_mt], i_mt);
    mt_list_[i_mt].coord_ += displacement[i_mt];
  }
  // wally_->PauseSim(2);
}