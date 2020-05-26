#include "kinesin.h"
#include "master_header.h"

int Kinesin::Monomer::GetKif4ANeighborCount() {

  if (site_ == nullptr or motor_->properties_->kinesin4.max_neighbs_ == 0) {
    return 0;
  }
  int n_neighbs{site_->GetKif4ANeighborCount()};
  if (motor_->heads_active_ == 2) {
    n_neighbs--;
  }
  return n_neighbs;
}

Kinesin::Kinesin() {}

void Kinesin::Initialize(system_parameters *parameters,
                         system_properties *properties, int id) {

  id_ = id;
  wally_ = &properties->wallace;
  parameters_ = parameters;
  properties_ = properties;
}

void Kinesin::ChangeConformation() {

  if (heads_active_ == 0) {
    wally_->ErrorExit("Kinesin::ChangeConformation()");
  }
  Tubulin *site{GetActiveHead()->site_};
  if (site->index_ == site->mt_->plus_end_ and
      parameters_->motors.endpausing_active) {
    return;
  }
  head_one_.trailing_ = !head_one_.trailing_;
  head_two_.trailing_ = !head_two_.trailing_;
}

double Kinesin::GetWeight_Bind_II() {

  // Use n_neighbs weight from dock site, use lattice weight from bound site
  Tubulin *dock_site{GetDockSite()};
  if (dock_site->occupied_) {
    return 0.0;
  }
  Monomer *bound_head{GetActiveHead()};
  double corrected_weight{bound_head->site_->weight_bind_};
  int neighbs_dock{dock_site->GetKif4ANeighborCount() - 1};
  int neighbs_bound{bound_head->site_->GetKif4ANeighborCount()};
  corrected_weight /= properties_->kinesin4.weight_neighbs_bind_[neighbs_bound];
  corrected_weight *= properties_->kinesin4.weight_neighbs_bind_[neighbs_dock];
  return corrected_weight;
}

double Kinesin::GetWeight_Unbind_II() {

  Monomer *chosen_head{nullptr};
  bool head_found{false};
  if (head_one_.ligand_ == "ADPP") {
    chosen_head = &head_one_;
    head_found = true;
  } else if (head_two_.ligand_ == "ADPP") {
    chosen_head = &head_two_;
    if (head_found) {
      wally_->ErrorExit("Kinesin::GetWeight_Unbind_II() [1]");
    }
  }
  if (chosen_head == nullptr) {
    wally_->ErrorExit("Kinesin::GetWeight_Unbind_II() [2]");
  }
  double corrected_weight{chosen_head->site_->weight_unbind_};
  corrected_weight *= corrected_weight;
  // Divide out one of the squared neighbor weights
  int n_neighbs{chosen_head->site_->GetKif4ANeighborCount()};
  corrected_weight /= properties_->kinesin4.weight_neighbs_unbind_[n_neighbs];
  return corrected_weight;
}

double Kinesin::GetWeight_BindATP_II() {

  Monomer *chosen_head{nullptr};
  bool head_found{false};
  if (head_one_.ligand_ == "NULL") {
    chosen_head = &head_one_;
    head_found = true;
  } else if (head_two_.ligand_ == "NULL") {
    chosen_head = &head_two_;
    if (head_found) {
      wally_->ErrorExit("Kinesin::GetWeight_BindATP_II()");
    }
  }
  if (chosen_head == nullptr) {
    wally_->ErrorExit("Kinesin::GetWeight_BindATP_II()");
  }
  double corrected_weight{chosen_head->site_->weight_bind_};
  // Remove contribution from neighb mechanism
  int n_neighbs{chosen_head->site_->GetKif4ANeighborCount()};
  corrected_weight /= properties_->kinesin4.weight_neighbs_bind_[n_neighbs];
  return corrected_weight;
}

Kinesin::Monomer *Kinesin::GetActiveHead() {

  if (heads_active_ != 1) {
    wally_->ErrorExit("Kinesin::GetActiveHead()");
  }
  if (head_one_.site_ != nullptr) {
    return &head_one_;
  } else {
    return &head_two_;
  }
}

Kinesin::Monomer *Kinesin::GetDockedHead() {

  Kinesin::Monomer *active_head{GetActiveHead()};
  if (active_head->ligand_ != "ADPP") {
    wally_->ErrorExit("Kinesin::GetDockedHead() [1]");
  }
  Kinesin::Monomer *docked_head{active_head->GetOtherHead()};
  if (docked_head->ligand_ != "ADP") {
    wally_->ErrorExit("Kinesin::GetDockedHead() [2]");
  }
  return docked_head;
}

Tubulin *Kinesin::GetDockSite() {

  Kinesin::Monomer *active_head{GetActiveHead()};
  Tubulin *site{active_head->site_};
  if (active_head->ligand_ != "ADPP") {
    wally_->ErrorExit("Kinesin::GetDockSite()");
  }
  int dir{-1};
  if (active_head->trailing_) {
    dir = 1;
  }
  int i_dock{site->index_ + dir * site->mt_->delta_x_};
  if (i_dock >= 0 and i_dock <= site->mt_->n_sites_ - 1) {
    return &site->mt_->lattice_[i_dock];
  } else {
    return nullptr;
  }
}
