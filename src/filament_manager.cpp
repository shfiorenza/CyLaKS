#include "filament_manager.hpp"
#include "protein_manager.hpp"

void FilamentManager::SetParameters() {

  using namespace Params;
  using namespace Filaments;
  // Filaments are mobile so long as at least 1 dimension is enabled
  for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
    if (Params::Filaments::translation_enabled[i_dim]) {
      mobile_ = true;
    }
  }
  if (Params::Filaments::rotation_enabled) {
    mobile_ = true;
  }
  n_bd_iterations_ = Params::Filaments::n_bd_per_kmc;
  dt_eff_ = Params::dt / n_bd_iterations_;
  /*
  auto is_unoccupied = [](BindingSite *site) {
    if (site->occupant_ == nullptr) {
      return true;
    }
    return false;
  };
  size_t max_list_size{0};
  for (int i_fil{0}; i_fil < count; i_fil++) {
    max_list_size += n_sites[i_fil];
  }
  unoccupied_.emplace("motors", Population<BindingSite>("motors", is_unoccupied,
                                                        max_list_size));
  Vec<int> i_min{0, 0, 0};
  Vec<size_t> max_size{1, 1, _n_neighbs_max + 1, max_list_size};
  auto get_n_neighbs = [](BindingSite *site) {
    Vec<int> indices_vec{(int)site->GetNumNeighborsOccupied()};
    return indices_vec;
  };
  unoccupied_.emplace("xlinks",
                      Population<BindingSite>("xlinks", is_unoccupied, max_size,
                                              i_min, get_n_neighbs));
                                              */
}

void FilamentManager::GenerateFilaments() {

  proto_.resize(Params::Filaments::count);
  for (int i_fil{0}; i_fil < proto_.size(); i_fil++) {
    proto_[i_fil].Initialize(_id_site, Sys::n_unique_objects_++, i_fil);
  }
  for (auto &&pf : proto_) {
    for (auto &&site : pf.sites_) {
      sites_.emplace_back(&site);
    }
  }
  if (proto_.size() == 2) {
    proto_[0].neighbor_ = &proto_[1];
    proto_[1].neighbor_ = &proto_[0];
  }
  /*
  Sys::Log("***\n");
  Sys::Log(" %zu \n", proto_.size());
  for (int i_fil{0}; i_fil < Params::Filaments::count; i_fil++) {
    Sys::Log("mt #%i\n", proto_[i_fil].index_);
    Sys::Log("%zu & %zu\n", proto_[i_fil].plus_end_->pos_.size(),
             proto_[i_fil].minus_end_->pos_.size());
    Sys::Log("plus-end: (%g, %g)\n", proto_[i_fil].plus_end_->pos_[0],
             proto_[i_fil].plus_end_->pos_[1]);
    Sys::Log("minus_end: (%g, %g)\n", proto_[i_fil].minus_end_->pos_[0],
             proto_[i_fil].minus_end_->pos_[1]);
  }
  exit(1);
  */
}

bool FilamentManager::AllFilamentsImmobile() {

  if (!mobile_) {
    return true;
  }
  for (auto const &pf : proto_) {
    if (Sys::i_step_ >= pf.immobile_until_) {
      return false;
    }
  }
  return true;
}

void FilamentManager::UpdateForces() {
  for (auto &&pf : proto_) {
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      pf.force_[i_dim] = 0.0;
    }
    pf.torque_ = 0.0;
  }
  proteins_->UpdateExtensions();
}

void FilamentManager::UpdateLattice() { proteins_->UpdateLatticeDeformation(); }