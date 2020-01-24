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

  // If diffusion is disabled for this simulation, immediately return
  if (!parameters_->microtubules.diffusion_on) {
    return;
  }
  int n_mts{parameters_->microtubules.count};
  bool mts_inactive{true};
  double delta_t{parameters_->delta_t};
  double current_time{properties_->current_step_ * delta_t};
  // Check that at least one microtubule is active (not immobilized)
  for (int i_mt{0}; i_mt < n_mts; i_mt++) {
    if (current_time >= parameters_->microtubules.immobile_until[i_mt]) {
      mts_inactive = false;
    }
  }
  // If no microtubules are active, return
  if (mts_inactive) {
    return;
  }
  double kbT{parameters_->kbT};
  double site_size{parameters_->microtubules.site_size};
  // Number of iterations to split diffusion step into
  int n_iterations = 10;
  // Adjust dt so that total time over all iterations is equal to one timestep
  double delta_t_eff = delta_t / n_iterations;
  // Loop thru iterations
  for (int i_it{0}; i_it < n_iterations; i_it++) {
    // Sum up all forces exerted on each microtubule
    double forces_summed[n_mts];
    for (int i_mt = 0; i_mt < n_mts; i_mt++) {
      forces_summed[i_mt] = 0.0; // mt_list_[i_mt].GetNetForce();
    }
    // Calculate the displacement of each microtubule
    int displacement[n_mts];
    for (int i_mt{0}; i_mt < n_mts; i_mt++) {
      // If microtubule is still immobilized, set displacement to 0 and continue
      if (current_time < parameters_->microtubules.immobile_until[i_mt]) {
        displacement[i_mt] = 0;
        continue;
      }
      double x_sq{site_size * site_size};
      double D{kbT / mt_list_[i_mt].gamma_};
      double tau{x_sq / (2 * D)};
      double p_diffu{delta_t_eff / tau};
      double ran0 = properties_->gsl.GetRanProb();
      if (ran0 < p_diffu) {
        double ran1 = properties_->gsl.GetRanProb();
        if (ran1 < 0.5) {
          displacement[i_mt] = -1;
        } else {
          displacement[i_mt] = 1;
        }
      } else {
        displacement[i_mt] = 0;
      }
      /*
       // Calculate instanteous velocity due to forces using drag coefficient
       double velocity{forces_summed[i_mt] / mt_list_[i_mt].gamma_};
       // Add gaussian noise
       double sigma{sqrt(2 * kbT * delta_t_eff / mt_list_[i_mt].gamma_)};
       double noise{properties_->gsl.GetGaussianNoise(sigma)};
       // Calculate raw displacement in nanometers
       double raw_displacement{velocity * delta_t_eff + noise};
       // Convert nm displacement to no. of sites
       double site_displacement{raw_displacement / site_size};
       int n_steps{(int)site_displacement};
       // Use leftover as a prob. to roll for another step
       double leftover{fabs(site_displacement - n_steps)};
       double ran{properties_->gsl.GetRanProb()};
       if (ran < leftover) {
         if (site_displacement > 0) {
           n_steps++;
         } else if (site_displacement < 0) {
           n_steps--;
         } else {
           printf("WHAT\n");
         }
       }
       displacement[i_mt] = n_steps;
       */
      // printf("Noise is %g: ", noise);
      // printf("leftover is %g | ", leftover);
      // printf("displacement was %i\n", n_steps);
    }
    /*  Run through MT list and update displacementsi */
    for (int i_mt = 0; i_mt < n_mts; i_mt++) {
      if (displacement[i_mt] == 0) {
        continue;
      }
      // Microtubule *mt = &mt_list_[i_mt];
      // Microtubule *neighb = mt->neighbor_;
      mt_list_[i_mt].coord_ += displacement[i_mt];
      // mt->UpdateExtensions();
      // neighb->UpdateExtensions();
    }
  }
}