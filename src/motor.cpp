#include "motor.hpp"
#include "protofilament.hpp"

void Motor::ChangeConformation() {

  if (n_heads_active_ != 1) {
    Sys::ErrorExit("Motor::ChangeConformation()");
  }
  BindingSite *site{GetActiveHead()->site_};
  if (site == nullptr) {
    Sys::ErrorExit("Motor::ChangeConformation() [2]");
  }
  //   printf("site %i\n", site->index_);
  if (site == site->filament_->plus_end_ and
      Params::Motors::endpausing_active) {
    if (Sys::test_mode_.empty() or site->filament_->index_ == 0) {
      return;
    }
  }
  //   printf("ChangeConformation for motor %i!\n", GetID());
  head_one_.trailing_ = !head_one_.trailing_;
  head_two_.trailing_ = !head_two_.trailing_;
  //   printf("wut\n");
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
  //   printf("hi\n");
  BindingSite *epicenter{nullptr};
  if (n_heads_active_ == 1) {
    epicenter = GetActiveHead()->site_;
  } else if (head_one_.trailing_) {
    epicenter = head_one_.site_;
  } else {
    epicenter = head_two_.site_;
  }
  int i_epicenter{(int)epicenter->index_};
  for (int delta{1}; delta <= Sys::lattice_cutoff_; delta++) {
    for (int dir{-1}; dir <= 1; dir += 2) {
      int i_scan{i_epicenter + dir * delta};
      int mt_length{(int)epicenter->filament_->sites_.size() - 1};
      if (i_scan < 0 or i_scan > mt_length) {
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
            printf("NO\n");
            exit(1);
          }
          site->AddWeight_Bind(Sys::weight_lattice_bind_[delta]);
          // printf("nop\n");
          site->AddWeight_Unbind(Sys::weight_lattice_unbind_[delta]);
          // printf("yop\n");
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
            printf("NO\n");
            exit(1);
          }
          site->AddWeight_Bind(Sys::weight_lattice_bind_[delta]);
          // printf("yop\n");
          site->AddWeight_Unbind(Sys::weight_lattice_unbind_[delta]);
          // printf("nop\n");
          continue;
        }
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

  if (Params::Motors::internal_force == 0.0 or
      Params::Motors::gaussian_range == 0) {
    return 0.0;
  }
  double weight_site{head->site_->GetWeight_Bind()};
  int n_neighbs{head->site_->GetNumNeighborsOccupied()};
  // Remove contribution from neighb mechanism
  return weight_site / Sys::weight_neighb_bind_[n_neighbs];
}

double Motor::GetWeight_Unbind_II(CatalyticHead *head) {

  if (head->ligand_ != CatalyticHead::Ligand::ADPP) {
    Sys::ErrorExit("Motor::GetWeight_Unbind_II()");
  }
  double weight_site{head->site_->GetWeight_Unbind()};
  // Ensure we use weight from trailing head to avoid self-coop
  if (!head->trailing_) {
    weight_site = head->GetOtherHead()->site_->GetWeight_Unbind();
  }
  // Divide out the weight from one neighbor, since it's the motor's other foot
  // weight_site /= Sys::weight_neighb_unbind_[1];
  // Disregard effects from internal force if it's disabled
  if (Params::Motors::internal_force == 0.0 or
      Params::Motors::gaussian_range == 0) {
    return weight_site;
  }
  double weight_sq{Square(weight_site)};
  // We only want the lattice contribution to be squared; divide out neighb term
  // int n_neighbs{head->site_->GetNumNeighborsOccupied() - 1};
  int n_neighbs{1};
  //   double weight{weight_sq / Sys::weight_neighb_unbind_[n_neighbs]};
  //   printf("weight = %g\n", weight);
  return weight_sq / Sys::weight_neighb_unbind_[n_neighbs];
}

double Motor::GetWeight_Unbind_I() {
  //   printf("%g\n", GetActiveHead()->site_->GetWeight_Unbind());
  return GetActiveHead()->site_->GetWeight_Unbind();
}

bool Motor::Diffuse(CatalyticHead *head, int dir) {

  BindingSite *old_site = head->site_;
  // printf("no\n");
  int i_new{(int)old_site->index_ + dir};
  // printf("i_old: %i | i_new: %i\n", old_site->index_, i_new);
  if (i_new < 0 or i_new > old_site->filament_->sites_.size() - 1) {
    return false;
  }
  // printf("chaching\n");
  BindingSite *new_site{&old_site->filament_->sites_[i_new]};
  if (new_site->occupant_ != nullptr) {
    return false;
  }
  old_site->occupant_ = nullptr;
  new_site->occupant_ = head;
  head->site_ = new_site;
  // printf("frfr\n\n");
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
      // n_runs_rec++ GOES HERE
    }
  }
  n_heads_active_--;
  if (Sys::test_mode_.empty()) {
    return true;
  }
  if (Sys::test_mode_ == "kinesin_mutant") {
    if (n_heads_active_ == 1 and head == &head_one_) {
      //   printf("bang\n");
      //   printf("head_")
      ChangeConformation();
    }
  }
  return true;
}

bool Motor::Tether() { return false; }

bool Motor::Untether() { return false; }