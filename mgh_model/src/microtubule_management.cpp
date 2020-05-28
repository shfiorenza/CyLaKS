#include "microtubule_management.h"
#include "master_header.h"

MicrotubuleManagement::MicrotubuleManagement() {}

void MicrotubuleManagement::Initialize(system_parameters *parameters,
                                       system_properties *properties) {

  wally_ = &properties->wallace;
  gsl_ = &properties->gsl;
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
  UpdateUnoccupied();
}

void MicrotubuleManagement::InitializeTestEnvironment() {

  if (strcmp(wally_->test_mode_, "motor_lattice_coop") == 0) {
    // Set parameters
    /*
    int lattice_cutoff{properties_->kinesin4.lattice_cutoff_};
    parameters_->microtubules.count = 1;
    parameters_->microtubules.length[0] = 2 * lattice_cutoff + 1;
    parameters_->microtubules.diffusion_on = false;
    */
  }
  if (strcmp(wally_->test_mode_, "xlink_bind_ii") == 0) {
    int xlink_cutoff{properties_->prc1.dist_cutoff_};
    parameters_->microtubules.count = 2;
    parameters_->microtubules.length[0] = 2 * (xlink_cutoff + 1) + 1;
    parameters_->microtubules.length[1] = 2 * (xlink_cutoff + 1) + 1;
    parameters_->microtubules.diffusion_on = false;
  }
  if (strcmp(wally_->test_mode_, "xlink_diffuse_ii") == 0) {
    parameters_->microtubules.count = 2;
    parameters_->microtubules.length[0] = 500;
    parameters_->microtubules.length[1] = 500;
    parameters_->microtubules.diffusion_on = false;
  }
  SetParameters();
  GenerateMicrotubules();
  InitializeLists();
  UpdateNeighbors();
  UpdateUnoccupied();
}

void MicrotubuleManagement::SetParameters() {

  max_neighbs_xlink_ = 2;
  // max_neighbs_motor_ = 0;
  n_diffusion_iterations_ = parameters_->microtubules.n_iterations;
  dt_eff_ = parameters_->delta_t / n_diffusion_iterations_;
  for (int i_mt{0}; i_mt < parameters_->microtubules.count; i_mt++) {
    n_sites_tot_ += parameters_->microtubules.length[i_mt];
  }
}

void MicrotubuleManagement::GenerateMicrotubules() {

  mt_list_.resize(parameters_->microtubules.count);
  for (int i_mt{0}; i_mt < mt_list_.size(); i_mt++) {
    mt_list_[i_mt].Initialize(parameters_, properties_, i_mt);
  }
}

void MicrotubuleManagement::InitializeLists() {

  n_unocc_xlink_.resize(max_neighbs_xlink_ + 1);
  unocc_xlink_.resize(max_neighbs_xlink_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_xlink_; n_neighbs++) {
    n_unocc_xlink_[n_neighbs] = 0;
    unocc_xlink_[n_neighbs].resize(n_sites_tot_);
  }
  /*
  n_unocc_motor_.resize(max_neighbs_motor_ + 1);
  unocc_motor_.resize(max_neighbs_motor_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_motor_; n_neighbs++) {
    n_unocc_motor_[n_neighbs] = 0;
    unocc_motor_[n_neighbs].resize(n_sites_tot_);
  }
  */
}

void MicrotubuleManagement::FlagForUpdate() { lists_up_to_date_ = false; }

void MicrotubuleManagement::UpdateNeighbors() {

  int n_mts{(int)mt_list_.size()};
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
  /*
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_motor_; n_neighbs++) {
    n_unocc_motor_[n_neighbs] = 0;
  }
  */
  for (int i_mt{0}; i_mt < mt_list_.size(); i_mt++) {
    for (int i_site{0}; i_site < mt_list_[i_mt].n_sites_; i_site++) {
      Tubulin *site{&mt_list_[i_mt].lattice_[i_site]};
      if (site->occupied_) {
        continue;
      }
      int n_neighbs_xl{site->GetPRC1NeighborCount()};
      unocc_xlink_[n_neighbs_xl][n_unocc_xlink_[n_neighbs_xl]++] = site;
      // int n_neighbs_mot{site->GetKif4ANeighborCount()};
      // unocc_motor_[n_neighbs_mot][n_unocc_motor_[n_neighbs_mot]++] = site;
    }
  }
}

double MicrotubuleManagement::GetSiteOffset() {

  if (mt_list_.size() != 2) {
    return 0.0;
  }
  double site_diff{mt_list_[0].coord_ - mt_list_[1].coord_};
  // Ranges from 0.0 to 0.5
  double site_offset{fabs(std::round(site_diff) - site_diff)};
  // printf("offset is %g\n", site_offset);
  return site_offset;
}

void MicrotubuleManagement::RunDiffusion() {

  // If diffusion is disabled for this simulation, immediately return
  if (!parameters_->microtubules.diffusion_on) {
    return;
  }
  bool mts_inactive{true};
  int n_mts{(int)mt_list_.size()};
  double current_time{properties_->current_step_ * parameters_->delta_t};
  // Check that at least one microtubule is active (not immobilized)
  for (int i_mt{0}; i_mt < n_mts; i_mt++) {
    if (current_time >= mt_list_[i_mt].immobile_until_) {
      mts_inactive = false;
    }
  }
  // If no microtubules are active, return
  if (mts_inactive) {
    return;
  }
  for (int i_itr{0}; i_itr < n_diffusion_iterations_; i_itr++) {
    // Sum up all forces exerted on each microtubule
    double forces_summed[n_mts];
    properties_->prc1.Update_Extensions();
    properties_->kinesin4.Update_Extensions();
    for (int i_mt{0}; i_mt < n_mts; i_mt++) {
      // If microtubule is still immobilized, simply continue
      if (current_time < mt_list_[i_mt].immobile_until_) {
        continue;
      }
      forces_summed[i_mt] = mt_list_[i_mt].GetNetForce();
      // Calculate instanteous velocity due to forces using drag coefficient
      double velocity{forces_summed[i_mt] / mt_list_[i_mt].gamma_};
      // Add gaussian noise to represent theraml motion, i.e., diffusion
      double noise{gsl_->GetGaussianNoise(mt_list_[i_mt].sigma_)};
      mt_list_[i_mt].coord_ += velocity * dt_eff_ + noise;
      // printf("coord[%i] is %g\n", i_mt, mt_list_[i_mt].coord_);
    }
  }
  // printf("\n  -- \n");
}