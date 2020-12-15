#include "protein_manager.hpp"
#include "binding_site.hpp"
#include "curator.hpp"
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
  motors_.AddWeight(name, _n_neighbs_max);
  xlinks_.AddWeight(name, _n_neighbs_max);
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

  Str name{"bruh"};
  double dt{Params::dt};
  double kbT{Params::kbT};

  // Bind I
  name = "bind_i";
  motors_.p_event_[name].val_ =
      Params::Motors::k_on * Params::Motors::c_bulk * dt;
  xlinks_.p_event_[name].InitVals(_n_neighbs_max + 1);
  for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
    xlinks_.p_event_[name].vals_[0][0][n_neighbs] =
        xlinks_.weights_["neighbs"].bind_[n_neighbs] * Params::Xlinks::k_on *
        Params::Xlinks::c_bulk * dt;
  }
  // Bind_I_Teth
  // Bind_II
  // Unbind_II
  // Unbind_I
  name = "unbind_i";
  motors_.p_event_[name].val_ = Params::Motors::k_off_i * dt;
  xlinks_.p_event_[name].InitVals(_n_neighbs_max + 1);
  for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
    xlinks_.p_event_[name].vals_[0][0][n_neighbs] =
        xlinks_.weights_["neighbs"].unbind_[n_neighbs] *
        Params::Xlinks::k_off_i * dt;
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

  /*
  // Bind_I: Bind first head of a protein
  name = "bind_i";
  auto exe_bind_i = [&](auto *site, auto *population) {
    auto entry{population->GetFreeEntry()};
    entry->Bind(site, &entry->head_one_);
    population->AddToActive(entry);
    filaments_->FlagForUpdate();
  };
  auto weight_bind_i = [&](BindingSite *site) {
    return site->GetWeight_Bind();
  };
  kmc_.events_.emplace_back(
      name, motors_.p_event_[name].GetVal(),
      &filaments_->unocc_["motors"].size_,
      &filaments_->unocc_["motors"].entries_, poisson, weight_bind_i,
      [&](BindingSite *site) { exe_bind_i(site, &motors_); }, gsl_);
  for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
    kmc_.events_.emplace_back(
        name, xlinks_.p_event_[name].GetVal(n_neighbs),
        &filaments_->unocc_["xlinks"].bin_size_[0][0][n_neighbs],
        &filaments_->unocc_["xlinks"].bin_entries_[0][0][n_neighbs], binomial,
        [&](BindingSite *site) { exe_bind_i(site, &xlinks_); }, gsl_);
  }
  */
  // Bind_I_Teth
  // Bind_II
  // Unbind_II
  // Unbind_I: Unbind first (singly bound) head of a protein
  /*
  name = "unbind_i";
  auto exe_unbind_i = [&](auto *head, auto *population) {
    auto entry{head->GetParent()};
    entry->Unbind(head);
    entry->UntetherSatellite();
    population->RemoveFromActive(entry);
    filaments_->FlagForUpdate();
  };
  auto weight_unbind_i = [&](auto head) {
    return head->site_->GetWeight_Unbind();
  };
  kmc_.events_.emplace_back(
      name, motors_.p_event_[name].GetVal(), &motors_.sorted_[name].size_,
      &motors_.sorted_[name].entries_, poisson, weight_unbind_i,
      [&](CatalyticHead *head) { exe_unbind_i(head, &motors_); }, gsl_);
      */
  /*
  for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
    kmc_.events_.emplace_back(
        name, xlinks_.p_event_[name].GetVal(n_neighbs),
        &xlinks_.sorted_[name].bin_size_[0][0][n_neighbs],
        &xlinks_.sorted_[name].bin_entries_[0][0][n_neighbs], binomial,
        [&](BindingHead *head) { exe_unbind_i(head, &xlinks_); }, gsl_);
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
