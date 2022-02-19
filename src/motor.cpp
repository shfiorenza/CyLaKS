#include "cylaks/motor.hpp"
#include "cylaks/protofilament.hpp"

double Motor::GetAnchorCoordinate(int i_dim) {

  // printf("mot called\n");
  if (n_heads_active_ == 0) {
    Sys::ErrorExit("Motor::GetAnchorCoord()");
  }
  if (n_heads_active_ == 1) {
    CatalyticHead *active_head{GetActiveHead()};
    int dir{active_head->trailing_ ? 1 : -1};
    int dx{active_head->site_->filament_->dx_};
    double delta{Params::Filaments::site_size / 2.0};
    return active_head->site_->pos_[i_dim] + (dir * dx * delta);
  }
  return (head_one_.site_->pos_[i_dim] + head_two_.site_->pos_[i_dim]) / 2.0;
}

void Motor::UpdateNeighbors_Bind_I_Teth() {

  if (n_heads_active_ != 0) {
    Sys::ErrorExit("Protein::UpdateNeighbors_Bind_I_Teth() [1]");
  }
  if (!IsTethered()) {
    Sys::ErrorExit("Protein::UpdateNeighbors_Bind_I_Teth() [2]");
  }
  if (teth_partner_->GetNumHeadsActive() == 0) {
    Sys::ErrorExit("Protein::UpdateNeighbors_Bind_I_Teth() [3]");
  }
  n_neighbors_bind_i_teth_ = 0;
  // ! FIXME add n_mt = 2 case
  BindingSite *anchor{teth_partner_->GetHeadOne()->site_};
  if (anchor == nullptr) {
    anchor = teth_partner_->GetHeadTwo()->site_;
  }
  for (int dx{Sys::teth_x_min_}; dx <= Sys::teth_x_max_; dx++) {
    for (int dir{-1}; dir <= 1; dir += 2) {
      int i_neighb{(int)anchor->index_ + dir * dx};
      // printf("%i + (%i)(%i) = %i\n", (int)anchor->index_, dir, dx, i_neighb);
      if (i_neighb < 0 or i_neighb >= anchor->filament_->sites_.size()) {
        continue;
      }
      BindingSite *neighb{&anchor->filament_->sites_[i_neighb]};
      if (!neighb->IsOccupied()) {
        // printf("%zu\n", neighbors_bind_i_teth_.size());
        // printf("%i\n", n_neighbors_bind_i_teth_);
        neighbors_bind_i_teth_[n_neighbors_bind_i_teth_++] = neighb;
      }
    }
  }
}

double Motor::GetSoloWeight_Bind_I_Teth(BindingSite *target) {

  double r_x{teth_partner_->GetAnchorCoordinate(0) - target->pos_[0]};
  double r_y{teth_partner_->GetAnchorCoordinate(1) - target->pos_[1]};
  double r{sqrt(Square(r_x) + Square(r_y))};
  if (r < tether_.r_min_ or r > tether_.r_max_) {
    return 0.0;
  }
  double weight_teth{tether_.GetWeight_Bind(r)};
  double weight_site{target->GetWeight_Bind()};
  return weight_teth * weight_site;
}

BindingSite *Motor::GetNeighbor_Bind_I_Teth() {

  double weight_tot{GetWeight_Bind_I_Teth()};
  double ran{SysRNG::GetRanProb()};
  double p_cum{0.0};
  Sys::Log(2, "%i NEIGHBS\n", n_neighbors_bind_i_teth_);
  Sys::Log(2, "ran = %g\n", ran);
  for (int i_neighb{0}; i_neighb < n_neighbors_bind_i_teth_; i_neighb++) {
    BindingSite *neighb{neighbors_bind_i_teth_[i_neighb]};
    p_cum += GetSoloWeight_Bind_I_Teth(neighb) / weight_tot;
    Sys::Log(2, "p_cum = %g\n", p_cum);
    if (ran < p_cum) {
      return neighb;
    }
  }
  return nullptr;
}

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
    // ! FIXME sort out this test nonsense
    // if (Sys::test_mode_.empty() or site->filament_->index_ == 0) {
    //   return;
    // }
  }
  head_one_.trailing_ = !head_one_.trailing_;
  head_two_.trailing_ = !head_two_.trailing_;
}

BindingSite *Motor::GetDockSite() {

  CatalyticHead *active_head{GetActiveHead()};
  if (active_head->ligand_ != Ligand::ADPP) {
    return nullptr;
  }
  BindingSite *site{active_head->site_};
  int dir{active_head->trailing_ ? 1 : -1};
  int i_dock{(int)site->index_ + dir * site->filament_->dx_};
  if (i_dock < 0 or i_dock > site->filament_->sites_.size() - 1) {
    return nullptr;
    // ! FIXME sort out this test nonsense
    /*
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
    */
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
        // ! FIXME figure out this test nonsense
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
  // ! FIXME  _lambda_neighb = 1.0 now; need to include neighb coop otherwise
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
  // ! FIXME  _lambda_neighb = 1.0 now; need to include neighb coop otherwise
  // int n_neighbs{head->site_->GetNumNeighborsOccupied()};
  // Remove contribution from neighb mechanism
  return weight_site; // / Sys::weight_neighb_bind_[n_neighbs];
}

double Motor::GetWeight_Unbind_II(CatalyticHead *head) {

  if (head->ligand_ != Ligand::ADPP) {
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
  return weight_sq;
}

double Motor::GetWeight_Unbind_I() {
  return GetActiveHead()->site_->GetWeight_Unbind();
}

double Motor::GetWeight_Bind_I_Teth() {

  double tot_weight{0.0};
  UpdateNeighbors_Bind_I_Teth();
  for (int i_neighb{0}; i_neighb < n_neighbors_bind_i_teth_; i_neighb++) {
    tot_weight += GetSoloWeight_Bind_I_Teth(neighbors_bind_i_teth_[i_neighb]);
  }
  // if (tot_weight != 0.0) {
  //   printf("%g\n", tot_weight);
  // }
  return tot_weight;
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

// ! FIXME see if we can down-cast to bindinghead & call Protein's Bind() funct
bool Motor::Bind(BindingSite *site, CatalyticHead *head) {

  if (site->occupant_ != nullptr) {
    return false;
  }
  site->occupant_ = head;
  head->site_ = site;
  head->ligand_ = Ligand::NONE;
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
      head->ligand_ = Ligand::ADPP;
    }
  }
  return true;
}

bool Motor::Bind_ATP(CatalyticHead *head) {

  if (head->ligand_ != Ligand::NONE) {
    Sys::ErrorExit("Motor::Bind_ATP");
  }
  head->ligand_ = Ligand::ATP;
  // Do not change conformation of trailing heads (for end-pausing)
  // Sys::Log()
  // ! FIXME convert this to Log() -- verbosity of 2 maybe?
  //   printf("ATP bound to head of motor %i on site %i %s\n",
  //          head->parent_->GetID(), head->site_->index_,
  //          head->trailing_ ? "(trailing)" : "");
  if (!head->trailing_) {
    ChangeConformation();
  }
  return true;
}

bool Motor::Hydrolyze(CatalyticHead *head) {

  if (head->ligand_ != Ligand::ATP) {
    Sys::ErrorExit("Motor::Hydrolyze()");
  }
  head->ligand_ = Ligand::ADPP;
  return true;
}

// ! FIXME same as Bind();
bool Motor::Unbind(CatalyticHead *head) {

  BindingSite *site{head->site_};
  site->occupant_ = nullptr;
  head->site_ = nullptr;
  head->ligand_ = Ligand::ADP;
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

bool Motor::Tether(Protein *target) {

  if (target->IsTethered()) {
    Sys::ErrorExit("Protein::Tether()");
  }
  teth_partner_ = target;
  target->teth_partner_ = this;
  return true;
}

// bool Motor::Untether() { return false; }