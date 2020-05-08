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
    // Set parameters
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
  for (int i_mt{0}; i_mt < mt_list_.size(); i_mt++) {
    for (int i_site{0}; i_site < mt_list_[i_mt].n_sites_; i_site++) {
      Tubulin *site{&mt_list_[i_mt].lattice_[i_site]};
      // Uncomment for motor_lattice_bind and self-coop test to work
      // if (wally_->test_mode_ == nullptr) {
      int n_neighbs{site->GetKif4ANeighborCount()};
      site->weight_bind_ =
          properties_->kinesin4.weight_neighbs_bind_[n_neighbs];
      site->weight_unbind_ =
          properties_->kinesin4.weight_neighbs_unbind_[n_neighbs];
      // }
      if (site->occupied_) {
        continue;
      }
      unocc_motor_[n_unocc_motor_++] = site;
      int n_neighbs_xl{site->GetPRC1NeighborCount()};
      unocc_xlink_[n_neighbs_xl][n_unocc_xlink_[n_neighbs_xl]++] = site;
    }
  }
  properties_->kinesin4.Update_Weights();
}

/*
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
  }
  //  Run through MT list and update displacements
  for (int i_mt = 0; i_mt < n_mts; i_mt++) {
    if (displacement[i_mt] == 0) {
      continue;
    }
    mt_list_[i_mt].coord_ += displacement[i_mt];
  }
}
*/