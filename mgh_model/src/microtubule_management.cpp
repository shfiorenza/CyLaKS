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

  int n_mts = parameters_->microtubules.count;
  n_sites_tot_ = 0;
  for (int i_mt = 0; i_mt < n_mts; i_mt++) {
    n_sites_tot_ += parameters_->microtubules.length[i_mt];
  }
  // int n_sites_bulk = n_sites_tot_ - 2*n_mts;
  unoccupied_list_.resize(n_sites_tot_);
  n_unoccupied_xl_.resize(3);
  unoccupied_list_xl_.resize(3);
  for (int n_neighbs(0); n_neighbs < 3; n_neighbs++) {
    n_unoccupied_xl_[n_neighbs] = 0;
    unoccupied_list_xl_[n_neighbs].resize(n_sites_tot_);
  }
}

void MicrotubuleManagement::GenerateMicrotubules() {

  int n_mts = parameters_->microtubules.count;
  mt_list_.resize(n_mts);
  for (int i_mt = 0; i_mt < n_mts; i_mt++) {
    mt_list_[i_mt].Initialize(parameters_, properties_, i_mt);
    // for (int i_site{0}; i_site < mt_list_[i_mt].n_sites_; i_site++) {
    //   PushToUnoccupied(&mt_list_[i_mt].lattice_[i_site]);
    // }
  }
}

void MicrotubuleManagement::UnoccupiedCheck(Tubulin *site) {

  if (site->motor_head_ != nullptr || site->xlink_head_ != nullptr) {
    printf("Error @ site %i_%i: should be unoccupied\n", site->mt_->index_,
           site->index_);
    exit(1);
  }
}

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
  n_unoccupied_ = 0;
  for (int n_neighbs(0); n_neighbs < 3; n_neighbs++) {
    n_unoccupied_xl_[n_neighbs] = 0;
  }
  for (int i_mt = 0; i_mt < parameters_->microtubules.count; i_mt++) {
    int n_sites = parameters_->microtubules.length[i_mt];
    int i_plus = mt_list_[i_mt].plus_end_;
    int i_minus = mt_list_[i_mt].minus_end_;
    int dx = mt_list_[i_mt].delta_x_;
    for (int i_site = 0; i_site < n_sites; i_site++) {
      Tubulin *site = &mt_list_[i_mt].lattice_[i_site];
      if (site->occupied_ == false) {
        unoccupied_list_[n_unoccupied_] = site;
        n_unoccupied_++;
        int n_neighbs = 0;
        if (i_site == i_plus) {
          if (mt_list_[i_mt].lattice_[i_site - dx].xlink_head_ != nullptr)
            n_neighbs++;
        } else if (i_site == i_minus) {
          if (mt_list_[i_mt].lattice_[i_site + dx].xlink_head_ != nullptr)
            n_neighbs++;
        } else {
          if (mt_list_[i_mt].lattice_[i_site - dx].xlink_head_ != nullptr)
            n_neighbs++;
          if (mt_list_[i_mt].lattice_[i_site + dx].xlink_head_ != nullptr)
            n_neighbs++;
        }
        unoccupied_list_xl_[n_neighbs][n_unoccupied_xl_[n_neighbs]] = site;
        n_unoccupied_xl_[n_neighbs]++;
      }
    }
  }
}

void MicrotubuleManagement::FlagForUpdate() { lists_up_to_date_ = false; }

void MicrotubuleManagement::PushToUnoccupied(Tubulin *site) {

  // Place in generic unocc_ list first
  unoccupied_list_[n_unoccupied_] = site;
  site->unocc_index_[0] = n_unoccupied_++;
  // Next, place in neighbor-specific unocc_ list
  int n_neighbs = site->GetPRC1NeighborCount();
  unoccupied_list_xl_[n_neighbs][n_unoccupied_xl_[n_neighbs]] = site;
  site->unocc_index_[1] = n_unoccupied_xl_[n_neighbs]++;
}

void MicrotubuleManagement::DeleteFromUnoccupied(Tubulin *site) {

  // Delete from generic unocc_ first
  Tubulin *last_entry = unoccupied_list_[n_unoccupied_ - 1];
  int site_index = site->unocc_index_[0];
  unoccupied_list_[site_index] = last_entry;
  last_entry->unocc_index_[0] = site_index;
  site->unocc_index_[0] = -1;
  n_unoccupied_--;
  // Next, delete from neighor_specific unocc_
  int n_neighbs = site->GetPRC1NeighborCount();
  int n_unocc = n_unoccupied_xl_[n_neighbs];
  last_entry = std::get<Tubulin *>(unoccupied_list_xl_[n_neighbs][n_unocc - 1]);
  site_index = site->unocc_index_[1];
  unoccupied_list_xl_[n_neighbs][site_index] = last_entry;
  last_entry->unocc_index_[1] = site_index;
  site->unocc_index_[1] = -1;
  n_unoccupied_xl_[n_neighbs]--;
}

Tubulin *MicrotubuleManagement::GetUnoccupiedSite() {

  UpdateUnoccupied();
  int n_unoccupied = n_unoccupied_;
  // Make sure an unoccupied site exists
  if (n_unoccupied > 0) {
    int i_entry = properties_->gsl.GetRanInt(n_unoccupied);
    Tubulin *site = unoccupied_list_[i_entry];
    UnoccupiedCheck(site);
    return site;
  } else {
    printf("Error: GetUnoccupiedSite called, but no unoccupied sites\n");
    exit(1);
  }
}

Tubulin *MicrotubuleManagement::GetUnoccupiedSite(int n_neighbs) {

  UpdateUnoccupied();
  int n_unoccupied = n_unoccupied_xl_[n_neighbs];
  if (n_unoccupied > 0) {
    int i_entry = properties_->gsl.GetRanInt(n_unoccupied);
    Tubulin *site =
        std::get<Tubulin *>(unoccupied_list_xl_[n_neighbs][i_entry]);
    UnoccupiedCheck(site);
    return site;
  } else {
    printf("Error: GetUnoccupiedSiteNEIGHB called, but no sites\n");
    exit(1);
  }
}

void MicrotubuleManagement::RunDiffusion() {

  if (!parameters_->microtubules.diffusion_on) {
    return;
  }
  double delta_t = parameters_->delta_t;
  int n_mts = parameters_->microtubules.count;
  double current_time = properties_->current_step_ * delta_t;
  bool mts_active{false};
  for (int i_mt{0}; i_mt < n_mts; i_mt++) {
    if (current_time >= parameters_->microtubules.immobile_until[i_mt]) {
      mts_active = true;
    }
  }
  if (!mts_active) {
    return;
  }
  int n_iterations = 10;
  double delta_t_eff = delta_t / n_iterations;
  double site_size = parameters_->microtubules.site_size;
  double kbT = parameters_->kbT;
  for (int i_it{0}; i_it < n_iterations; i_it++) {
    sys_timepoint start1 = sys_clock::now();
    /*	Sum up all forces exerted on the MTs 	*/
    double forces_summed[n_mts];
    for (int i_mt = 0; i_mt < n_mts; i_mt++) {
      mt_list_[i_mt].UpdateExtensions();
    }
    for (int i_mt = 0; i_mt < n_mts; i_mt++) {
      forces_summed[i_mt] = mt_list_[i_mt].GetNetForce();
      // if (mt_list_[i_mt].GetNetForce_Motors() != 0) {
      //   printf("motor force for mt %i is %g\n", i_mt,
      //          mt_list_[i_mt].GetNetForce_Motors());
      // }
      //       printf("force is %g for mt #%i\n", forces_summed[i_mt],
      //       i_mt);
    }
    // Check for symmetry
    double tolerance = 0.000001;
    for (int i_mt{0}; i_mt < n_mts; i_mt += 2) {
      double delta = fabs(forces_summed[i_mt] + forces_summed[i_mt + 1]);
      if (delta > tolerance) {
        printf("aw man in MT diffusion\n");
        printf("for mt %i: %g, for mt %i: %g\n", i_mt, forces_summed[i_mt],
               i_mt + 1, -forces_summed[i_mt + 1]);
        exit(1);
      }
    }
    sys_timepoint finish = sys_clock::now();
    auto elapsed = std::chrono::duration_cast<t_unit>(finish - start1);
    properties_->wallace.t_MTs_[1] += elapsed.count();
    sys_timepoint start2 = sys_clock::now();
    elapsed = std::chrono::duration_cast<t_unit>(finish - start2);
    /*	Calculate MT displacements for this timestep */
    int displacement[n_mts];
    for (int i_mt{0}; i_mt < n_mts; i_mt++) {
      if (current_time >= parameters_->microtubules.immobile_until[i_mt]) {
        double velocity = forces_summed[i_mt] / mt_list_[i_mt].gamma_;
        // gaussian noise is added
        double sigma = sqrt(2 * kbT * delta_t_eff / mt_list_[i_mt].gamma_);
        double noise = properties_->gsl.GetGaussianNoise(sigma);
        double raw_displacement = velocity * delta_t_eff + noise;
        double site_displacement = raw_displacement / site_size;
        // Get number of sites MT is expected to move
        int n_steps = (int)site_displacement;
        // Use leftover as a prob. to roll for another step
        double leftover = fabs(site_displacement - n_steps);
        double ran = properties_->gsl.GetRanProb();
        if (ran < leftover and site_displacement > 0)
          n_steps++;
        else if (ran < leftover and site_displacement < 0)
          n_steps--;
        // Store value in array
        displacement[i_mt] = n_steps;
      } else
        displacement[i_mt] = 0;
    }
    finish = sys_clock::now();
    elapsed = std::chrono::duration_cast<t_unit>(finish - start2);
    properties_->wallace.t_MTs_[2] += elapsed.count();
    start2 = sys_clock::now();
    /*  Run through MT list and update displacementsi */
    for (int i_mt = 0; i_mt < n_mts; i_mt++) {
      Microtubule *mt = &mt_list_[i_mt];
      Microtubule *neighb = &mt_list_[mt->mt_index_adj_];
      int n_steps = displacement[i_mt];
      if (n_steps != 0) {
        int dx = 0;
        if (n_steps > 0)
          dx = 1;
        else if (n_steps < 0) {
          dx = -1;
          n_steps = abs(n_steps);
        }
        for (int i_step = 0; i_step < n_steps; i_step++) {
          mt->coord_ += dx;
        }
        mt->UpdateExtensions();
        neighb->UpdateExtensions();
      }
    }
    finish = sys_clock::now();
    elapsed = std::chrono::duration_cast<t_unit>(finish - start2);
    properties_->wallace.t_MTs_[3] += elapsed.count();
    elapsed = std::chrono::duration_cast<t_unit>(finish - start1);
    properties_->wallace.t_MTs_[0] += elapsed.count();
  }
}
