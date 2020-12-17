#include "protein_manager.hpp"
#include "binding_site.hpp"
#include "curator.hpp"
#include "filament_manager.hpp"
#include "system_namespace.hpp"
#include "system_rng.hpp"

void ProteinManager::GenerateReservoirs() {

  size_t reservoir_size{0};
  for (int i_mt{0}; i_mt < Params::Filaments::count; i_mt++) {
    reservoir_size += Params::Filaments::n_sites[i_mt];
  }
  size_t motor_step_active{size_t(Params::Motors::t_active / Params::dt)};
  size_t xlink_step_active(size_t(Params::Xlinks::t_active / Params::dt));
  motors_.Initialize(_id_motor, reservoir_size, motor_step_active);
  xlinks_.Initialize(_id_xlink, reservoir_size, xlink_step_active);
}

void ProteinManager::InitializeWeights() {

  /*
    For events that result in a change in energy dE, we use Boltzmann factors to
    scale rates appropriately. Detailed balance is satisfied by the factors:
                  exp{-(1.0 - lambda) * dE / kBT}, and
                  exp{lambda * dE / kBT}
    for forward and reverse pathways, respectively, (e.g., binding and
    unbinding), where lambda is a constant that ranges from 0 to 1, and kBT is
    thermal energy of the system. 3 values of lambda demonstrate its function:
        Lambda = 0.0 means all energy dependences is in binding
        Lambda = 0.5 means energy dependence is equal for binding and unbinding
        Lambda = 1.0 means all energy dependence is in unbinding
   */
  // Neighbor stuff
  Str name{"neighbs"};
  motors_.AddWeight(name, _n_neighbs_max + 1);
  xlinks_.AddWeight(name, _n_neighbs_max + 1);
  double lambda_neighb{1.0};
  double dE{0.0};
  for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
    dE = -1 * Params::Motors::neighb_neighb_energy * n_neighbs;
    motors_.weights_[name].bind_[n_neighbs] = exp(-(1.0 - lambda_neighb) * dE);
    motors_.weights_[name].unbind_[n_neighbs] = exp(lambda_neighb * dE);
    dE = -1 * Params::Xlinks::neighb_neighb_energy * n_neighbs;
    xlinks_.weights_[name].bind_[n_neighbs] = exp(-(1.0 - lambda_neighb) * dE);
    xlinks_.weights_[name].unbind_[n_neighbs] = exp(lambda_neighb * dE);
  }
}

void ProteinManager::SetParameters() {

  using namespace Params;

  // Bind I -- motors
  auto is_unoccupied = [&](Object *base) { return !base->IsOccupied(); };
  // filaments_->AddPop("motors", is_unoccupied);
  // motors_.AddProb("bind_i", Motors::k_on * Motors::c_bulk * dt);
  // Bind I -- xlinks
  Vec<int> i_min{0, 0, 0};
  Vec<size_t> dim_size{1, 1, _n_neighbs_max + 1};
  auto get_n_neighbs = [](Object *site) {
    Vec<int> indices_vec{site->GetNeighborCount()};
    return indices_vec;
  };
  filaments_->AddPop("xlinks", is_unoccupied, dim_size, i_min, get_n_neighbs);
  xlinks_.AddProb("bind_i", Xlinks::k_on * Xlinks::c_bulk * dt, "neighbs", 0);

  // Bind_I_Teth
  // Bind_II
  // Unbind_II
  // Unbind_I -- motors
  auto is_docked = [](Object *base) {
    if (base->GetNumHeadsActive() == 1) {
      return true;
    }
    return false;
  };
  // motors_.AddPop("bound_ADPP_i", is_docked);
  // motors_.AddProb("unbind_i", Motors::k_off_i * dt);
  // Unbind I -- xlinks
  xlinks_.AddPop("bound_i", is_docked, dim_size, i_min, get_n_neighbs);
  xlinks_.AddProb("unbind_i", Xlinks::k_off_i * dt, "neighbs", 1);

  // motors_.sorted_.emplace(pop_name,
  // Population<Object>(pop_name, is_docked, n_max));

  /*
  for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
    printf("bind[%i]: %g\n", n_neighbs,
           xlinks_.p_event_.at("bind_i").GetVal(n_neighbs));
  }
  for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
    printf("unbind[%i]: %g\n", n_neighbs,
           xlinks_.p_event_.at("unbind_i").GetVal(n_neighbs));
  }
  */
  // Unbind_I_Teth
  // Tether_Free
  // Untether_Free

  // vv motors only vv
  // Bind_ATP
  // Hydrolyze_ATP
  // Tether_Bound
  // Untether_Bound

  // vv xlinks only vv
  // Diffusion
  // Bind_II_Teth
  // Unbind_II_Teth
}

void ProteinManager::InitializeTestEnvironment() {}

void ProteinManager::InitializeTestEvents() {}

void ProteinManager::InitializeEvents() {

  Str name{"bruh"};
  //  Binomial probabilitiy distribution; sampled to predict most events
  auto binomial = [&](double p, int n) {
    if (n > 0) {
      return SysRNG::SampleBinomial(p, n);
    } else {
      return 0;
    }
  };
  // Poisson distribution; sampled to predict events w/ variable probabilities
  auto poisson = [&](double p, int n) {
    if (p > 0.0) {
      return SysRNG::SamplePoisson(p);
    } else {
      return 0;
    }
  };

  // Bind_I: Bind first head of a protein
  name = "bind_i";
  /*
  auto exe_bind_i = [&](auto *site, auto *population) {
    auto entry{population->GetFreeEntry()};
    entry->Bind(site, &entry->head_one_);
    population->AddToActive(entry);
    filaments_->FlagForUpdate();
  };
  */
  /*
   auto weight_bind_i = [&](Object *base) {
     return dynamic_cast<BindingSite *>(base)->GetWeight_Bind();
   };
   */
  auto exe_bind_i = [&](Object *base) {
    if (Sys::i_step_ < xlinks_.step_active_) {
      return;
    }
    auto site{dynamic_cast<BindingSite *>(base)};
    auto xlink{xlinks_.GetFreeEntry()};
    // printf("xlink id is %i\n", xlink->GetID());
    // printf("site id is %i\n", site->GetID());
    // printf("%li active entries\n", xlinks_.n_active_entries_);
    bool executed{xlink->Bind(site, &xlink->head_one_)};
    if (executed) {
      xlinks_.AddToActive(xlink);
      filaments_->FlagForUpdate();
    }
  };
  /*
  kmc_.events_.emplace_back(name, motors_.p_event_.at(name).GetVal(),
                            &filaments_->unoccupied_.at("motors").size_,
                            &filaments_->unoccupied_.at("motors").entries_,
                            poisson, weight_bind_i, exe_bind_i);
  */
  for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
    kmc_.events_.emplace_back(
        name, xlinks_.p_event_.at("bind_i").GetVal(n_neighbs),
        &filaments_->unoccupied_.at("xlinks").bin_size_[0][0][n_neighbs],
        &filaments_->unoccupied_.at("xlinks").bin_entries_[0][0][n_neighbs],
        binomial, exe_bind_i);
  }
  // Bind_I_Teth
  // Bind_II
  // Unbind_II
  // Unbind_I: Unbind first (singly bound) head of a protein
  /*
  auto exe_unbind_i = [&](auto *head, auto *population) {
    auto entry{head->GetParent()};
    entry->Unbind(head);
    entry->UntetherSatellite();
    population->RemoveFromActive(entry);
    filaments_->FlagForUpdate();
  };
  */
  name = "unbind_i";
  /*
  auto weight_unbind_i = [&](Object *head) {
    return dynamic_cast<CatalyticHead *>(head)->site_->GetWeight_Unbind();
  };
  */

  // ACTIVE HOLDS PROTS/MOTS NOT HEADS ATM!!
  auto exe_unbind_i = [&](Object *base) {
    // auto head{dynamic_cast<BindingHead *>(base)};
    auto xlink{dynamic_cast<Protein *>(base)};
    bool executed{xlink->Unbind(xlink->GetActiveHead())};
    if (executed) {
      xlink->UntetherSatellite();
      xlinks_.RemoveFromActive(xlink);
      filaments_->FlagForUpdate();
    }
  };
  /*
  kmc_.events_.emplace_back(name, motors_.p_event_.at(name).GetVal(),
                            &motors_.sorted_.at(name).size_,
                            &motors_.sorted_.at(name).entries_, poisson,
                            weight_unbind_i, exe_unbind_i);
  */
  for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
    kmc_.events_.emplace_back(
        name, xlinks_.p_event_.at("unbind_i").GetVal(n_neighbs),
        &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
        &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs], binomial,
        exe_unbind_i);
  }
  // Unbind_I_Teth
  // Tether_Free
  // Untether_Free

  // vv motors only vv
  // Bind_ATP
  // Hydrolyze_ATP
  // Tether_Bound
  // Untether_Bound

  // vv xlinks only vv
  // Diffusion
  // Bind_II_Teth
  // Unbind_II_Teth

  for (auto const &event : kmc_.events_) {
    printf("%s: %g\n", event.name_.c_str(), event.p_occur_);
  }
}

void ProteinManager::UpdateFilaments() { filaments_->UpdateUnoccupied(); }
