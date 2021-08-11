#include "cylaks/motor.hpp"
#include "cylaks/protofilament.hpp"

void Motor::ChangeConformation() {

  if (n_heads_active_ != 1) {
    Sys::ErrorExit("Motor::ChangeConformation()");
  }
  BindingSite *site{GetActiveHead()->site_};
  if (site == nullptr) {
    Sys::ErrorExit("Motor::ChangeConformation() [2]");
  }
  if (site == site->filament_->plus_end_ and
      Params::Motors::endpausing_active) {
    return;
    // if (Sys::test_mode_.empty() or site->filament_->index_ == 0) {
    //   return;
    // }
  }
  head_one_.trailing_ = !head_one_.trailing_;
  head_two_.trailing_ = !head_two_.trailing_;
}

BindingSite *Motor::GetDockSite() {

  CatalyticHead *active_head{GetActiveHead()};
  if (active_head->ligand_ != CatalyticHead::Ligand::ADPP) {
    return nullptr;
  }
  BindingSite *site{active_head->site_};
  int dir{active_head->trailing_ ? 1 : -1};
  int i_dock{(int)site->index_ + dir * site->filament_->dx_};
  if (i_dock < 0 or i_dock > site->filament_->sites_.size() - 1) {
    if (Sys::test_mode_.empty()) {
      return nullptr;
    }
    if (Sys::test_mode_ != "filament_ablation") {
      return nullptr;
    }
    if (Sys::i_step_ > Sys::ablation_step_) {
      //   Params::Motors::endpausing_active = false;
      return nullptr;
    }
    if (i_dock == -1 and site->filament_->index_ == 1) {
      BindingSite *minus_end{site->filament_->neighbor_->minus_end_};
      return site->filament_->neighbor_->minus_end_;
    } else {
      return nullptr;
    }
  }
  return &site->filament_->sites_[i_dock];
}

CatalyticHead *Motor::GetDockedHead() {

  if (n_heads_active_ == 1) {
    BindingSite *dock{GetDockSite()};
    if (dock != nullptr) {
      if (dock->occupant_ == nullptr) {
        return GetActiveHead()->GetOtherHead();
      }
    }
  }
  return nullptr;
}

void Motor::ApplyLatticeDeformation() {

  if (n_heads_active_ == 0) {
    return;
  }
  // Epicenter is the location of the motor; deformation propagates from here
  BindingSite *epicenter{nullptr};
  // If singly bound, there is no ambiguity to epicenter location
  if (n_heads_active_ == 1) {
    epicenter = GetActiveHead()->site_;
    // Otherwise, by convention, always use the trailing head when doubly bound
  } else if (head_one_.trailing_) {
    epicenter = head_one_.site_;
  } else {
    epicenter = head_two_.site_;
  }
  // Get index of epicenter in microtubule lattice
  int i_epicenter{(int)epicenter->index_};
  // Starting from +/- 1 site from the epicenter, apply lattice deformations
  for (int delta{1}; delta <= Sys::lattice_cutoff_; delta++) {
    for (int dir{-1}; dir <= 1; dir += 2) {
      // Index of site we're currently applying the deformation to
      int i_scan{i_epicenter + dir * delta};
      int mt_length{(int)epicenter->filament_->sites_.size() - 1};
      // Do not access lattice sites that don't exist in memory
      if (i_scan < 0 or i_scan > mt_length) {
        continue;
        /*
        if (Sys::test_mode_.empty()) {
          continue;
        }
        if (Sys::test_mode_ != "filament_ablation" or
            Sys::i_step_ >= Sys::ablation_step_) {
          continue;
        }
        if (epicenter->filament_->index_ == 0 and i_scan > mt_length) {
          int i_adj{i_scan - mt_length};
          // printf("i_adj = %i\n", i_adj);
          auto other_mt{epicenter->filament_->neighbor_};
          if (i_adj > other_mt->sites_.size() - 1) {
            continue;
          }
          BindingSite *site{&other_mt->sites_[i_adj]};
          if (site == nullptr) {
            Sys::ErrorExit("Motor::ApplyLatticeDeformation [0]");
          }
          site->AddWeight_Bind(Sys::weight_lattice_bind_[delta]);
          site->AddWeight_Unbind(Sys::weight_lattice_unbind_[delta]);
          continue;
        }
        if (epicenter->filament_->index_ == 1 and i_scan < 0) {
          auto other_mt{epicenter->filament_->neighbor_};
          int i_adj{(int)other_mt->sites_.size() + i_scan};
          // printf("i_adj = %i\n", i_adj);
          if (i_adj < 0) {
            continue;
          }
          BindingSite *site{&other_mt->sites_[i_adj]};
          if (site == nullptr) {
            Sys::ErrorExit("Motor::ApplyLatticeDeformation [1]");
          }
          site->AddWeight_Bind(Sys::weight_lattice_bind_[delta]);
          site->AddWeight_Unbind(Sys::weight_lattice_unbind_[delta]);
          continue;
        }
        */
      }
      if (i_scan >= 0 and i_scan < epicenter->filament_->sites_.size()) {
        BindingSite *site{&epicenter->filament_->sites_[i_scan]};
        site->AddWeight_Bind(Sys::weight_lattice_bind_[delta]);
        site->AddWeight_Unbind(Sys::weight_lattice_unbind_[delta]);
      }
    }
  }
}

double Motor::GetWeight_Diffuse(CatalyticHead *head, int dir) { return 0.0; }

double Motor::GetWeight_Bind_II() {

  BindingSite *dock{GetDockSite()};
  if (dock == nullptr) {
    return 0.0;
  }
  if (dock->occupant_ != nullptr) {
    return 0.0;
  }
  // Use n_neighbs from dock, but lattice weight from bound_site
  CatalyticHead *bound_head{GetActiveHead()};
  // By using bound_head for lattice weight, we avoid any self-coop
  double weight_site{bound_head->site_->GetWeight_Bind()};
  // FIXME  _lambda_neighb = 1.0 now, but need to include neighb coop otherwise
  return weight_site;
}

double Motor::GetWeight_BindATP_II(CatalyticHead *head) {

  if (!Params::Motors::gaussian_stepping_coop or
      Params::Motors::gaussian_range == 0) {
    return 0.0;
  }
  if (n_heads_active_ != 2) {
    printf("bruh\n");
  }
  // Ensure we use weight from trailing head to avoid self-coop from lattice
  if (!head->trailing_) {
    head = head->GetOtherHead();
  }
  double weight_site{head->site_->GetWeight_Bind()};
  // int n_neighbs{head->site_->GetNumNeighborsOccupied()};
  // Remove contribution from neighb mechanism
  return weight_site; // / Sys::weight_neighb_bind_[n_neighbs];
}

double Motor::GetWeight_Unbind_II(CatalyticHead *head) {

  if (head->ligand_ != CatalyticHead::Ligand::ADPP) {
    Sys::ErrorExit("Motor::GetWeight_Unbind_II()");
  }

  if (n_heads_active_ != 2) {
    printf("bruh\n");
  }
  // Ensure we use weight from trailing head to avoid self-coop from lattice
  if (!head->trailing_) {
    head = head->GetOtherHead();
  }
  double weight_site{head->site_->GetWeight_Unbind()};
  // Divide out the weight from one neighbor, since it's the motor's other foot
  // if (Sys::test_mode_ != "motor_lattice_step") {
  // weight_site /= Sys::weight_neighb_unbind_[1];
  // }
  // Disregard effects from internal force if it's disabled
  if (!Params::Motors::gaussian_stepping_coop or
      Params::Motors::gaussian_range == 0) {
    return weight_site;
  }
  double weight_sq{Square(weight_site)};
  int n_neighbs{head->site_->GetNumNeighborsOccupied()};
  if (n_neighbs == 2) {
    Sys::ErrorExit("bruh");
  }
  weight_sq /= Sys::weight_neighb_unbind_[n_neighbs];
  // FIXME make function that automatically divides out neighbors to avoid
  // this Divide out sq'd contribution if motor has a real neighb (i.e., not
  // itself) if (head->GetNumNeighborsOccupied() == 2) {
  //   weight_sq /= Sys::weight_neighb_unbind_[1];
  //   // printf("%zu and %zu\n", Sys::i_step_, Sys::i_datapoint_);
  // }
  // printf("wt is %g\n", weight_sq);
  return weight_sq;
}

double Motor::GetWeight_Unbind_I() {
  //   printf("%g\n", GetActiveHead()->site_->GetWeight_Unbind());
  return GetActiveHead()->site_->GetWeight_Unbind();
}

bool Motor::Diffuse(CatalyticHead *head, int dir) {

  BindingSite *old_site = head->site_;
  int i_new{(int)old_site->index_ + dir};
  if (i_new < 0 or i_new > old_site->filament_->sites_.size() - 1) {
    return false;
  }
  BindingSite *new_site{&old_site->filament_->sites_[i_new]};
  if (new_site->occupant_ != nullptr) {
    return false;
  }
  old_site->occupant_ = nullptr;
  new_site->occupant_ = head;
  head->site_ = new_site;
  return true;
}

// FIXME see if we can down-cast to bindinghead & call Protein's Bind() funct
bool Motor::Bind(BindingSite *site, CatalyticHead *head) {

  if (site->occupant_ != nullptr) {
    return false;
  }
  site->occupant_ = head;
  head->site_ = site;
  head->ligand_ = CatalyticHead::Ligand::NONE;
  // If we bound from bulk solution, initialize head direction
  if (n_heads_active_ == 0) {
    head->trailing_ = false;
    head->GetOtherHead()->trailing_ = true;
  }
  n_heads_active_++;
  if (Sys::test_mode_.empty()) {
    return true;
  }
  if (Sys::test_mode_ == "kinesin_mutant") {
    if (head == &head_two_) {
      head->ligand_ = CatalyticHead::Ligand::ADPP;
    }
  }
  return true;
}

bool Motor::Bind_ATP(CatalyticHead *head) {

  if (head->ligand_ != CatalyticHead::Ligand::NONE) {
    Sys::ErrorExit("Motor::Bind_ATP");
  }
  head->ligand_ = CatalyticHead::Ligand::ATP;
  // Do not change conformation of trailing heads (for end-pausing)
  //   printf("ATP bound to head of motor %i on site %i %s\n",
  //          head->parent_->GetID(), head->site_->index_,
  //          head->trailing_ ? "(trailing)" : "");
  if (!head->trailing_) {
    ChangeConformation();
  }
  return true;
}

bool Motor::Hydrolyze(CatalyticHead *head) {

  if (head->ligand_ != CatalyticHead::Ligand::ATP) {
    Sys::ErrorExit("Motor::Hydrolyze()");
  }
  head->ligand_ = CatalyticHead::Ligand::ADPP;
  return true;
}

// FIXME same as Bind();
bool Motor::Unbind(CatalyticHead *head) {

  BindingSite *site{head->site_};
  site->occupant_ = nullptr;
  head->site_ = nullptr;
  head->ligand_ = CatalyticHead::Ligand::ADP;
  // If we are about to completely unbind, record this motor run
  if (n_heads_active_ == 1) {
    if (!Sys::equilibrating_) {
      Sys::n_runs_recorded_++;
      if (Sys::n_runs_recorded_ >= Params::Motors::n_runs_to_exit) {
        Sys::early_exit_triggered_ = true;
      }
    }
  }
  n_heads_active_--;
  if (Sys::test_mode_.empty()) {
    return true;
  }
  if (Sys::test_mode_ == "kinesin_mutant") {
    if (n_heads_active_ == 1 and head == &head_one_) {
      ChangeConformation();
    }
  }
  return true;
}

bool Motor::Tether() { return false; }

bool Motor::Untether() { return false; }