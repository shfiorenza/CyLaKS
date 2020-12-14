#include "filament_manager.hpp"
#include "protein_manager.hpp"

void FilamentManager::SetParameters() {

  // Filaments are mobile so long as at least 1 dimension is enabled
  for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
    if (Params::Filaments::dimension_enabled[i_dim]) {
      mobile_ = true;
    }
  }
}

void FilamentManager::GenerateFilaments() {

  proto_.resize(Params::Filaments::count);
  for (int i_fil{0}; i_fil < proto_.size(); i_fil++) {
    proto_[i_fil].Initialize(_id_site, Sys::n_unique_objects_++, i_fil);
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

void FilamentManager::InitializeTestEnvironment() {}

bool FilamentManager::NoMobileFilamentsYet() {

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

void FilamentManager::UpdateProteins() { proteins_->UpdateExtensions(); }

void FilamentManager::UpdateLattice() { proteins_->UpdateLatticeDeformation(); }