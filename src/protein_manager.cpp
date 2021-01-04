#include "protein_manager.hpp"
#include "binding_site.hpp"
#include "curator.hpp"
#include "filament_manager.hpp"
#include "system_namespace.hpp"
#include "system_rng.hpp"
#include <iostream>
#include <string>

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
  auto is_unocc = [&](Object *site) -> Vec<Object *> {
    if (!site->IsOccupied()) {
      return {site};
    }
    return {};
  };
  // filaments_->AddPop("motors", get_unocc);
  // motors_.AddProb("bind_i", Motors::k_on * Motors::c_bulk * dt);
  // Bind I -- xlinks
  Vec<int> i_min{0, 0, 0};
  Vec<size_t> dim_size{1, 1, _n_neighbs_max + 1};
  auto get_n_neighbs = [](Object *entry) {
    Vec<int> indices_vec{entry->GetNumNeighborsOccupied()};
    return indices_vec;
  };
  filaments_->AddPop("xlinks", is_unocc, dim_size, i_min, get_n_neighbs);

  xlinks_.AddProb("bind_i", Xlinks::k_on * Xlinks::c_bulk * dt, "neighbs", 0);

  // Bind_I_Teth
  // Bind_II
  auto is_singly_bound = [](Object *protein) -> Vec<Object *> {
    if (protein->GetNumHeadsActive() == 1) {
      return {protein->GetActiveHead()};
    }
    return {};
  };
  xlinks_.AddPop("bind_ii_candidates", is_singly_bound);
  xlinks_.AddProb("bind_ii", Xlinks::k_on * Xlinks::c_eff_bind * dt);
  // Unbind_II
  auto is_doubly_bound = [&](Object *protein) -> Vec<Object *> {
    if (protein->GetNumHeadsActive() == 2) {
      return {protein->GetHeadOne(), protein->GetHeadTwo()};
    }
    return {};
  };
  xlinks_.AddPop("bound_ii", is_doubly_bound);
  xlinks_.AddProb("unbind_ii", Xlinks::k_off_ii * dt);

  // Unbind_I -- motors
  // motors_.AddPop("bound_ADPP_i", is_docked);
  // motors_.AddProb("unbind_i", Motors::k_off_i * dt);
  // Unbind I -- xlinks
  xlinks_.AddPop("bound_i", is_singly_bound, dim_size, i_min, get_n_neighbs);
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
  double x_sq{Square(Filaments::site_size / 1000)}; // in um^2
  double tau_i{x_sq / (2 * Xlinks::d_i)};
  double tau_ii{x_sq / (2 * Xlinks::d_ii)};
  double p_diffuse_i{dt / tau_i};
  double p_diffuse_ii{dt / tau_ii};
  // diff_fwd and diff_bck are two separate events, which effectively
  // doubles the probability to diffuse. Thus we divide p_diff by 2.
  p_diffuse_i /= 2.0;
  p_diffuse_ii /= 2.0;
  Vec3D<double> vec_diff_i(
      1, Vec2D<double>(1, Vec<double>(_n_neighbs_max + 1, p_diffuse_i)));
  vec_diff_i[0][0][1] *= xlinks_.weights_.at("neighbs").unbind_[1];
  vec_diff_i[0][0][2] = 0.0;
  xlinks_.AddProb("diffuse_i_fwd", vec_diff_i);
  xlinks_.AddProb("diffuse_i_bck", vec_diff_i);
  xlinks_.AddPop("diffuse_ii_to_rest_candidates", is_doubly_bound);
  xlinks_.AddPop("diffuse_ii_fr_rest_candidates", is_doubly_bound);
  xlinks_.AddProb("diffuse_ii_to_rest", p_diffuse_ii);
  xlinks_.AddProb("diffuse_ii_fr_rest", p_diffuse_ii);
  // Bind_II_Teth
  // Unbind_II_Teth
}

void ProteinManager::InitializeTestEnvironment() {

  Sys::Log(" ** Initializing test '%s'\n", Sys::test_mode_.c_str());
  using namespace Params;
  if (Sys::test_mode_ == "xlink_bind_ii") {
    Sys::Log("Enter xlink spring cutoff in nm: ");
    Str input;
    std::getline(std::cin, input);
    double dist_cutoff{std::stod(input)};
    printf("input = %g\n", dist_cutoff);
    // Initialize statistics
    double r_y{std::fabs(Filaments::y_initial[0] - Filaments::y_initial[1])};
    if (r_y > dist_cutoff) {
      Sys::Log("Error; filament distance greater than spring cutoff\n");
      exit(1);
    }
    double r_x_cutoff{sqrt(Square(dist_cutoff) - Square(r_y))};
    printf("r_x_max = %g\n", r_x_cutoff);
    size_t x_max((size_t)std::ceil(r_x_cutoff / Filaments::site_size));
    printf("x_max = %zu\n", x_max);
    // Calculate event probability (sans energy effects) and store it
    xlinks_.AddProb("bind_ii", Xlinks::k_on * Xlinks::c_eff_bind * dt);
    xlinks_.AddProb("unbind_ii", Xlinks::k_off_ii * dt);
    Vec<double> p_bind(x_max + 1, xlinks_.p_event_.at("bind_ii").GetVal());
    Vec<double> p_unbind(x_max + 1, xlinks_.p_event_.at("unbind_ii").GetVal());
    for (int x{0}; x <= x_max; x++) {
      double lambda{0.5};
      double r_x{x * Filaments::site_size};
      double r{sqrt(Square(r_x) + Square(r_y))};
      double dr{r - Params::Xlinks::r_0};
      double dE{0.5 * Params::Xlinks::k_spring * Square(dr)};
      // printf("dE[%i] = %.3g\n", x, dE);
      p_bind[x] *= exp(-(1.0 - lambda) * dE / Params::kbT);
      p_unbind[x] *= exp(lambda * dE / Params::kbT);
      printf("p[%i] = %.3g (%.3g)\n", x, p_bind[x], p_unbind[x]);
    }
    test_ref_.emplace("bind_ii", p_bind);
    test_ref_.emplace("unbind_ii", p_unbind);
    Vec<Pair<size_t, size_t>> zeros(x_max + 1, {0, 0});
    test_stats_.emplace("bind_ii", zeros);
    test_stats_.emplace("unbind_ii", zeros);
    // Initialize filament environment
    Filaments::count = 2;
    Filaments::n_sites[0] = Filaments::n_sites[1] = 2 * x_max + 1;
    Filaments::translation_enabled[0] = false;
    Filaments::translation_enabled[1] = false;
    Filaments::rotation_enabled = false;
    filaments_->Initialize(this);
    // FIXME -- shouldn't matter, but we cant have this in SetParams
    // filaments_->AddPop("xlinks", get_unocc, dim_size, i_min, get_n_neighbs);

    // Place first xlink head on lower MT; remains static for entire sim
    int i_site{x_max};
    BindingSite *site{&filaments_->proto_[0].sites_[i_site]};
    Protein *xlink{xlinks_.GetFreeEntry()};
    bool executed{xlink->Bind(site, &xlink->head_one_)};
    if (executed) {
      xlinks_.AddToActive(xlink);
      filaments_->FlagForUpdate();
    } else {
      Sys::ErrorExit("ProteinManager::InitializeTestEnvironment()");
    }
  }
}

void ProteinManager::InitializeTestEvents() {

  using namespace Params;
  // Str name{};
  if (Sys::test_mode_ == "xlink_bind_ii") {
    // KMC Event -- Bind_II
    // name = "bind_ii";
    // Add population tracker for potential targets of Bind_II event
    auto is_singly_bound = [](Object *protein) -> Vec<Object *> {
      if (protein->GetNumHeadsActive() == 1) {
        return {protein->GetActiveHead()};
      }
      return {};
    };
    xlinks_.AddPop("bind_ii", is_singly_bound);
    // Helper functions for Bind_II event structure
    auto poisson_bind_ii = [&](double p, int n) {
      for (int x{0}; x < test_stats_.at("bind_ii").size(); x++) {
        int c{x == 0 ? 1 : 2}; // 1 site at x = 0. Otherwise, 2 due to
        int n_entries{c * xlinks_.sorted_.at("bind_ii").size_};
        test_stats_.at("bind_ii")[x].second += n_entries;
      }
      if (p == 0.0) {
        return 0;
      }
      return SysRNG::SamplePoisson(p);
    };
    auto get_weight_bind_ii = [](Object *base) {
      auto head{dynamic_cast<BindingHead *>(base)};
      return head->parent_->GetTotalWeight_Bind_II();
    };
    auto exe_bind_ii = [&](Object *base_head) {
      auto bound_head{dynamic_cast<BindingHead *>(base_head)};
      auto head{bound_head->GetOtherHead()};
      auto site{head->parent_->GetNeighbor_Bind_II()};
      auto executed{head->parent_->Bind(site, head)};
      if (executed) {
        xlinks_.FlagForUpdate();
        filaments_->FlagForUpdate();
        head->parent_->UpdateExtension();
        double r_x{head->pos_[0] - bound_head->pos_[0]};
        size_t x{std::abs(std::round(r_x / Params::Filaments::site_size))};
        test_stats_.at("bind_ii")[x].first++;
      } else {
        Sys::ErrorExit("Bind_II (TEST)");
      }
    };
    // Construct KMC event fr Bind_II
    kmc_.events_.emplace_back("bind_ii",
                              xlinks_.p_event_.at("bind_ii").GetVal(),
                              &xlinks_.sorted_.at("bind_ii").size_,
                              &xlinks_.sorted_.at("bind_ii").entries_,
                              poisson_bind_ii, get_weight_bind_ii, exe_bind_ii);
    // KMC event -- Unbind_II
    // name = "unbind_ii";
    auto is_doubly_bound = [&](Object *protein) -> Vec<Object *> {
      // Only ever unbind second head
      if (protein->GetNumHeadsActive() == 2) {
        return {protein->GetHeadTwo()};
      }
      return {};
    };
    xlinks_.AddPop("unbind_ii", is_doubly_bound);
    auto poisson_unbind_ii = [&](double p, int n) {
      if (xlinks_.sorted_.at("unbind_ii").size_ > 0) {
        auto head{xlinks_.sorted_.at("unbind_ii").entries_[0]};
        double r_x{head->pos_[0] - head->GetOtherHead()->pos_[0]};
        size_t x{std::abs(std::round(r_x / Params::Filaments::site_size))};
        test_stats_.at("unbind_ii")[x].second += 1;
      }
      if (p == 0.0) {
        return 0;
      }
      return SysRNG::SamplePoisson(p);
    };
    auto get_weight_unbind_ii = [](Object *base) {
      auto head{dynamic_cast<BindingHead *>(base)};
      return head->GetWeight_Unbind_II();
    };
    auto exe_unbind_ii = [&](Object *base) {
      auto head{dynamic_cast<BindingHead *>(base)};
      bool executed{head->Unbind()};
      if (executed) {
        xlinks_.FlagForUpdate();
        filaments_->FlagForUpdate();
        double r_x{head->pos_[0] - head->GetOtherHead()->pos_[0]};
        size_t x{std::abs(std::round(r_x / Params::Filaments::site_size))};
        test_stats_.at("unbind_ii")[x].first++;
        head->parent_->UpdateExtension();
      } else {
        Sys::ErrorExit("Unbind_II (TEST)");
      }
    };
    kmc_.events_.emplace_back(
        "unbind_ii", xlinks_.p_event_.at("unbind_ii").GetVal(),
        &xlinks_.sorted_.at("unbind_ii").size_,
        &xlinks_.sorted_.at("unbind_ii").entries_, poisson_unbind_ii,
        get_weight_unbind_ii, exe_unbind_ii);
  }
}

void ProteinManager::InitializePopulations() {}

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
  auto exe_bind_ii = [&](Object *base_head) {
    auto bound_head{dynamic_cast<BindingHead *>(base_head)};
    auto head{bound_head->GetOtherHead()};
    auto site{head->parent_->GetNeighbor_Bind_II()};
    auto executed{head->parent_->Bind(site, head)};
    if (executed) {
      head->parent_->UpdateExtension();
      xlinks_.FlagForUpdate();
      filaments_->FlagForUpdate();
    }
  };
  auto get_weight_bind_ii = [](Object *base) {
    auto head{dynamic_cast<BindingHead *>(base)};
    return head->parent_->GetTotalWeight_Bind_II();
  };
  if (xlinks_.crosslinking_active_) {
    kmc_.events_.emplace_back(
        "bind_ii", xlinks_.p_event_.at("bind_ii").GetVal(),
        &xlinks_.sorted_.at("bind_ii_candidates").size_,
        &xlinks_.sorted_.at("bind_ii_candidates").entries_, poisson,
        get_weight_bind_ii, exe_bind_ii);
  }
  // Unbind_II
  auto exe_unbind_ii = [&](Object *base) {
    auto head{dynamic_cast<BindingHead *>(base)};
    bool executed{head->Unbind()};
    if (executed) {
      xlinks_.FlagForUpdate();
      filaments_->FlagForUpdate();
    }
  };
  auto get_weight_unbind_ii = [](Object *base) {
    auto head{dynamic_cast<BindingHead *>(base)};
    return head->GetWeight_Unbind_II();
  };
  if (xlinks_.crosslinking_active_) {
    kmc_.events_.emplace_back("unbind_ii",
                              xlinks_.p_event_.at("unbind_ii").GetVal(),
                              &xlinks_.sorted_.at("bound_ii").size_,
                              &xlinks_.sorted_.at("bound_ii").entries_, poisson,
                              get_weight_unbind_ii, exe_unbind_ii);
  }

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

  auto exe_unbind_i = [&](Object *base) {
    auto head{dynamic_cast<BindingHead *>(base)};
    bool executed{head->Unbind()};
    if (executed) {
      head->UntetherSatellite();
      xlinks_.RemoveFromActive(head->parent_);
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
  auto exe_diffuse_fwd = [&](Object *base) {
    auto head{dynamic_cast<BindingHead *>(base)};
    bool executed{head->Diffuse(1)};
    if (executed) {
      head->parent_->UpdateExtension();
      filaments_->FlagForUpdate();
    }
  };
  auto exe_diffuse_bck = [&](Object *base) {
    auto head{dynamic_cast<BindingHead *>(base)};
    bool executed{head->Diffuse(-1)};
    if (executed) {
      head->parent_->UpdateExtension();
      filaments_->FlagForUpdate();
    }
  };
  for (int n_neighbs{0}; n_neighbs < _n_neighbs_max; n_neighbs++) {
    kmc_.events_.emplace_back(
        "diffuse_i_fwd", xlinks_.p_event_.at("diffuse_i_fwd").GetVal(n_neighbs),
        &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
        &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs], binomial,
        exe_diffuse_fwd);
    kmc_.events_.emplace_back(
        "diffuse_i_bck", xlinks_.p_event_.at("diffuse_i_bck").GetVal(n_neighbs),
        &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
        &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs], binomial,
        exe_diffuse_bck);
  }
  auto get_weight_diff_ii_to = [](Object *base) {
    auto head{dynamic_cast<BindingHead *>(base)};
    return head->GetWeight_Diffuse(1);
  };
  auto get_weight_diff_ii_fr = [](Object *base) {
    auto head{dynamic_cast<BindingHead *>(base)};
    return head->GetWeight_Diffuse(-1);
  };
  if (xlinks_.crosslinking_active_) {
    kmc_.events_.emplace_back(
        "diffuse_ii_to_rest",
        xlinks_.p_event_.at("diffuse_ii_to_rest").GetVal(),
        &xlinks_.sorted_.at("diffuse_ii_to_rest_candidates").size_,
        &xlinks_.sorted_.at("diffuse_ii_to_rest_candidates").entries_, poisson,
        get_weight_diff_ii_to, exe_diffuse_fwd);
    kmc_.events_.emplace_back(
        "diffuse_ii_fr_rest",
        xlinks_.p_event_.at("diffuse_ii_fr_rest").GetVal(),
        &xlinks_.sorted_.at("diffuse_ii_fr_rest_candidates").size_,
        &xlinks_.sorted_.at("diffuse_ii_fr_rest_candidates").entries_, poisson,
        get_weight_diff_ii_fr, exe_diffuse_bck);
  }
  // Bind_II_Teth
  // Unbind_II_Teth

  for (auto const &event : kmc_.events_) {
    printf("%s: %g\n", event.name_.c_str(), event.p_occur_);
  }
}

void ProteinManager::UpdateFilaments() { filaments_->UpdateUnoccupied(); }
