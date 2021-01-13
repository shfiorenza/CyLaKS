#include "motor.hpp"
#include "protofilament.hpp"

void Motor::ChangeConformation() {

  if (n_heads_active_ != 1) {
    Sys::ErrorExit("Motor::ChangeConformation()");
  }
  BindingSite *site{GetActiveHead()->site_};
  if (site == site->filament_->plus_end_ and
      Params::Motors::endpausing_active) {
    return;
  }
  //   printf("ChangeConformation for motor %i!\n", GetID());
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
  int i_dock{site->index_ + dir * site->filament_->dx_};
  if (i_dock < 0 or i_dock > site->filament_->sites_.size() - 1) {
    return nullptr;
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

double Motor::GetWeight_Unbind_II(CatalyticHead *head) {

  if (head->ligand_ != CatalyticHead::Ligand::ADPP) {
    Sys::ErrorExit("Motor::GetWeight_Unbind_II()");
  }
  double weight_site{head->site_->GetWeight_Unbind()};
  // Disregard effects from internal force if it's disabled
  if (Params::Motors::internal_force == 0.0) {
    return weight_site;
  }
  double weight_sq{Square(weight_site)};
  // We only want the lattice contribution to be squared; divide out neighb term
  int n_neighbs{head->site_->GetNumNeighborsOccupied()};
  return weight_sq / Sys::weight_neighb_unbind_[n_neighbs];
}

double Motor::GetWeight_Unbind_I() {
  //   printf("%g\n", GetActiveHead()->site_->GetWeight_Unbind());
  return GetActiveHead()->site_->GetWeight_Unbind();
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
  return true;
}

bool Motor::Tether() {}

bool Motor::Untether() {}