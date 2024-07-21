#include "cylaks/binding_site.hpp"
#include "cylaks/binding_head.hpp"
#include "cylaks/protofilament.hpp"

int BindingSite::GetNumNeighborsOccupied_Tot() {

  if (_n_neighbs_max == 0) {
    return 0;
  }
  int n_neighbs{0};
  for (auto const &site : neighbors_) {
    if (site->occupant_ != nullptr) {
      if (occupant_ != nullptr) {
        if (occupant_->parent_ != site->occupant_->parent_) {
          n_neighbs++;
        }
      } else {
        n_neighbs++;
      }
    }
  }
  return n_neighbs;
}

int BindingSite::GetNumNeighborsOccupied_Tot_Side() {

  int n_neighbs{0};
  if (filament_->top_neighb_ != nullptr) {
    if (filament_->top_neighb_->sites_[index_].occupant_ != nullptr) {
      n_neighbs++;
    }
  }
  if (filament_->bot_neighb_ != nullptr) {
    if (filament_->bot_neighb_->sites_[index_].occupant_ != nullptr) {
      n_neighbs++;
    }
  }
  return n_neighbs;
}

int BindingSite::GetNumNeighborsOccupied_Xlink() {

  if (_n_neighbs_max == 0) {
    return 0;
  }
  int n_neighbs{0};
  for (auto const &site : neighbors_) {
    if (site->occupant_ != nullptr) {
      if (occupant_ != nullptr) {
        if (occupant_->parent_ != site->occupant_->parent_) {
          if (site->occupant_->GetSpeciesID() == _id_xlink) {
            n_neighbs++;
          }
        }
      } else if (site->occupant_->GetSpeciesID() == _id_xlink) {
        n_neighbs++;
      }
    }
  }
  return n_neighbs;
}

int BindingSite::GetNumNeighborsOccupied_Xlink_Side() {

  int n_neighbs{0};
  if (filament_->top_neighb_ != nullptr) {
    if (filament_->top_neighb_->sites_[index_].occupant_ != nullptr) {
      if (filament_->top_neighb_->sites_[index_].occupant_->GetSpeciesID() ==
          _id_xlink) {
        n_neighbs++;
      }
    }
  }
  if (filament_->bot_neighb_ != nullptr) {
    if (filament_->bot_neighb_->sites_[index_].occupant_ != nullptr) {
      if (filament_->bot_neighb_->sites_[index_].occupant_->GetSpeciesID() ==
          _id_xlink) {
        n_neighbs++;
      }
    }
  }
  return n_neighbs;
}

int BindingSite::GetNumNeighborsOccupied_Motor() {

  if (_n_neighbs_max == 0) {
    return 0;
  }
  int n_neighbs{0};
  for (auto const &site : neighbors_) {
    if (site->occupant_ != nullptr) {
      if (occupant_ != nullptr) {
        if (occupant_->parent_ != site->occupant_->parent_) {
          if (site->occupant_->GetSpeciesID() == _id_motor) {
            n_neighbs++;
          }
        }
      } else if (site->occupant_->GetSpeciesID() == _id_motor) {
        n_neighbs++;
      }
    }
  }
  return n_neighbs;
}

int BindingSite::GetNumNeighborsOccupied_Motor_Side() {

  int n_neighbs{0};
  if (filament_->top_neighb_ != nullptr) {
    if (filament_->top_neighb_->sites_[index_].occupant_ != nullptr) {
      if (filament_->top_neighb_->sites_[index_].occupant_->GetSpeciesID() ==
          _id_motor) {
        n_neighbs++;
      }
    }
  }
  if (filament_->bot_neighb_ != nullptr) {
    if (filament_->bot_neighb_->sites_[index_].occupant_ != nullptr) {
      if (filament_->top_neighb_->sites_[index_].occupant_->GetSpeciesID() ==
          _id_motor) {
        n_neighbs++;
      }
    }
  }
  return n_neighbs;
}

BindingSite *BindingSite::GetNeighbor(int dir) {

  if (dir != 1 and dir != -1) {
    Sys::ErrorExit("BindingSite::GetNeighb()");
  }
  for (auto const &neighb : neighbors_) {
    if (index_ + dir == neighb->index_) {
      return neighb;
    }
  }
  return nullptr;
}

BindingSite *BindingSite::GetNeighbor_Side(int dir) {

  if (dir == 1) {
    if (filament_->top_neighb_ == nullptr) {
      return nullptr;
    }
    return &filament_->top_neighb_->sites_[index_];
  } else if (dir == -1) {
    if (filament_->bot_neighb_ == nullptr) {
      return nullptr;
    }
    return &filament_->bot_neighb_->sites_[index_];
  } else {
    Sys::ErrorExit("BindingSite::GetNeighb_Side()");
  }
  return nullptr;
}

void BindingSite::AddForce(Vec<double> f) { filament_->AddForce(this, f); }
void BindingSite::AddTorque(double tq) { filament_->AddTorque(tq); }