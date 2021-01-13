#include "protein_manager.hpp"
#include "binding_site.hpp"
#include "curator.hpp"
#include "filament_manager.hpp"
#include "system_namespace.hpp"
#include "system_rng.hpp"
#include <iostream>
#include <limits>
#include <string>

void ProteinManager::GenerateReservoirs() {

  size_t reservoir_size{0};
  for (int i_mt{0}; i_mt < Params::Filaments::count; i_mt++) {
    reservoir_size += Params::Filaments::n_sites[i_mt];
  }
  size_t motor_step_active{size_t(Params::Motors::t_active / Params::dt)};
  if (Params::Motors::c_bulk == 0.0) {
    motor_step_active = std::numeric_limits<size_t>::max();
  }
  motors_.Initialize(_id_motor, reservoir_size, motor_step_active);
  size_t xlink_step_active{size_t(Params::Xlinks::t_active / Params::dt)};
  if (Params::Xlinks::c_bulk == 0.0) {
    xlink_step_active = std::numeric_limits<size_t>::max();
  }
  xlinks_.Initialize(_id_xlink, reservoir_size, xlink_step_active);
}

void ProteinManager::InitializeWeights() {

  /*
    For events that result in a change in energy dE, we use Boltzmann factors
    to scale rates appropriately. Detailed balance is satisfied by the
    factors: exp{-(1.0 - lambda) * dE / kBT}, and exp{lambda * dE / kBT} for
    forward and reverse pathways, respectively, (e.g., binding and unbinding),
    where lambda is a constant that ranges from 0 to 1, and kBT is thermal
    energy of the system. 3 values of lambda demonstrate its function: Lambda
    = 0.0 means all energy dependences is in binding Lambda = 0.5 means energy
    dependence is equal for binding and unbinding Lambda = 1.0 means all
    energy dependence is in unbinding
   */
  // Neighbor stuff
  Str name{"neighbs"};
  Sys::weight_neighb_bind_.resize(_n_neighbs_max + 1);
  Sys::weight_neighb_unbind_.resize(_n_neighbs_max + 1);
  motors_.AddWeight(name, _n_neighbs_max + 1);
  xlinks_.AddWeight(name, _n_neighbs_max + 1);
  double lambda_neighb{1.0};
  double dE{0.0};
  for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
    dE = -1 * Params::Motors::neighb_neighb_energy * n_neighbs;
    motors_.weights_[name].bind_[n_neighbs] = exp(-(1.0 - _lambda_neighb) * dE);
    motors_.weights_[name].unbind_[n_neighbs] = exp(_lambda_neighb * dE);
    Sys::weight_neighb_bind_[n_neighbs] = exp(-(1.0 - _lambda_neighb) * dE);
    Sys::weight_neighb_unbind_[n_neighbs] = exp(_lambda_neighb * dE);
    dE = -1 * Params::Xlinks::neighb_neighb_energy * n_neighbs;
    xlinks_.weights_[name].bind_[n_neighbs] = exp(-(1.0 - _lambda_neighb) * dE);
    xlinks_.weights_[name].unbind_[n_neighbs] = exp(_lambda_neighb * dE);
  }
}

void ProteinManager::SetParameters() {

  using namespace Params;
  // Bind I
  xlinks_.AddProb("bind_i", Xlinks::k_on * Xlinks::c_bulk * dt, "neighbs", 0);
  motors_.AddProb("bind_i", Motors::k_on * Motors::c_bulk * dt);
  // Bind_I_Teth

  // Bind_ATP -- motors only
  double p_bind_ATP{Motors::k_on_ATP * Motors::c_ATP * dt};
  motors_.AddProb("bind_ATP_i", p_bind_ATP);
  double wt_ATP_ii{exp(-Motors::internal_force * Motors::sigma_ATP / kbT)};
  if (Motors::internal_force == 0.0) {
    wt_ATP_ii = 0.0;
  }
  motors_.AddProb("bind_ATP_ii", p_bind_ATP * wt_ATP_ii);
  // Hydrolyze_ATP -- motors only
  motors_.AddProb("hydrolyze", Motors::k_hydrolyze * dt);
  // Bind_II
  xlinks_.AddProb("bind_ii", Xlinks::k_on * Xlinks::c_eff_bind * dt);
  motors_.AddProb("bind_ii", Motors::k_on * Motors::c_eff_bind * dt);
  // Bind_II_Teth -- xlinks only

  // Unbind_II
  xlinks_.AddProb("unbind_ii", Xlinks::k_off_ii * dt);
  motors_.AddProb("unbind_ii", Motors::k_off_ii * dt);
  // Unbind_II_Teth -- xlinks only

  // Unbind I
  xlinks_.AddProb("unbind_i", Xlinks::k_off_i * dt, "neighbs", 1);
  motors_.AddProb("unbind_i", Motors::k_off_i * dt);
  // Unbind_I_Teth

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
  Vec3D<double> weight_diff{{Vec<double>(_n_neighbs_max + 1, p_diffuse_i)}};
  weight_diff[0][0][1] *= xlinks_.weights_.at("neighbs").unbind_[1];
  weight_diff[0][0][2] *= 0.0;
  xlinks_.AddProb("diffuse_i_fwd", weight_diff);
  xlinks_.AddProb("diffuse_i_bck", weight_diff);
  xlinks_.AddProb("diffuse_ii_to_rest", p_diffuse_ii);
  xlinks_.AddProb("diffuse_ii_fr_rest", p_diffuse_ii);
  // Tether_Free

  // Tether_Bound -- motors only

  // Untether_Free

  // Untether_Bound -- motors only
}

void ProteinManager::InitializeTestEnvironment() {

  Sys::Log(" ** Initializing test '%s'\n", Sys::test_mode_.c_str());
  using namespace Params;
  if (Sys::test_mode_ == "xlink_bind_ii") {
    double r_y{std::fabs(Filaments::y_initial[0] - Filaments::y_initial[1])};
    double r_max{xlinks_.r_max_};
    printf("r_max = %g\n", r_max);
    double r_x_max{sqrt(Square(r_max) - Square(r_y))};
    printf("r_x_max = %g\n", r_x_max);
    size_t x_max((size_t)std::ceil(r_x_max / Filaments::site_size));
    printf("x_max = %zu\n", x_max);
    Vec<double> p_bind(x_max + 1, xlinks_.p_event_.at("bind_ii").GetVal());
    Vec<double> p_unbind(x_max + 1, xlinks_.p_event_.at("unbind_ii").GetVal());
    for (int x{0}; x <= x_max; x++) {
      double lambda{0.5};
      double r_x{x * Filaments::site_size};
      double r{sqrt(Square(r_x) + Square(r_y))};
      if (r < xlinks_.r_min_ or r > xlinks_.r_max_) {
        p_bind[x] *= 0.0;
        p_unbind[x] *= 0.0;
        continue;
      }
      double dr{r - Params::Xlinks::r_0};
      double dE{0.5 * Params::Xlinks::k_spring * Square(dr)};
      p_bind[x] *= exp(-(1.0 - lambda) * dE / Params::kbT);
      p_unbind[x] *= 100 * exp(lambda * dE / Params::kbT);
    }
    test_ref_.emplace("bind_ii", p_bind);
    test_ref_.emplace("unbind_ii", p_unbind);
    Vec<Pair<size_t, size_t>> zeros(x_max + 1, {0, 0});
    test_stats_.emplace("bind_ii", zeros);
    test_stats_.emplace("unbind_ii", zeros);
    // Initialize filament environment
    Filaments::count = 2;
    Filaments::n_sites[0] = Filaments::n_sites[1] = 2 * x_max + 1;
    Sys::Log("  N_SITES[0] = %i\n", Filaments::n_sites[0]);
    Sys::Log("  N_SITES[1] = %i\n", Filaments::n_sites[1]);
    Filaments::translation_enabled[0] = false;
    Filaments::translation_enabled[1] = false;
    Filaments::rotation_enabled = false;
    filaments_->Initialize(this);
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
  } else if (Sys::test_mode_ == "filament_separation") {
    // Initialize filament environment
    if (Filaments::n_sites[0] != Filaments::n_sites[1]) {
      printf("\nError! Filaments must be the same length.\n");
      exit(1);
    }
    Filaments::immobile_until[0] = 0.0;
    Filaments::immobile_until[1] = 0.0;
    Filaments::translation_enabled[0] = true;
    Filaments::translation_enabled[1] = true;
    Filaments::rotation_enabled = false; // true;
    filaments_->Initialize(this);
    Str response;
    printf("\nEnter number of crosslinkers: ");
    std::getline(std::cin, response);
    int n_xlinks{(int)std::stoi(response)};
    Sys::Log("%i crosslinkers initialized.\n", n_xlinks);
    int n_places{(int)filaments_->sites_.size() / 2};
    if (n_xlinks > n_places) {
      printf("\nError! Too many crosslinkers for filament length used.\n");
      exit(1);
    }
    // Randomly place crosslinkers on filaments w/ x = 0
    int site_indices[n_places];
    for (int index{0}; index < n_places; index++) {
      site_indices[index] = index;
    }
    SysRNG::Shuffle(site_indices, n_places, sizeof(int));
    for (int i_xlink{0}; i_xlink < n_xlinks; i_xlink++) {
      Protein *xlink{xlinks_.GetFreeEntry()};
      int i_site{site_indices[i_xlink]};
      BindingSite *site_one{&filaments_->proto_[0].sites_[i_site]};
      BindingSite *site_two{&filaments_->proto_[1].sites_[i_site]};
      bool exe_one{xlink->Bind(site_one, &xlink->head_one_)};
      bool exe_two{xlink->Bind(site_two, &xlink->head_two_)};
      if (exe_one and exe_two) {
        bool still_attached{xlink->UpdateExtension()};
        if (still_attached) {
          xlinks_.AddToActive(xlink);
          filaments_->FlagForUpdate();
        } else {
          Sys::ErrorExit("ProteinManager::InitializeTestEnvironment() [2]");
        }
      } else {
        Sys::ErrorExit("ProteinManager::InitializeTestEnvironment() [1]");
      }
    }
  }
}

void ProteinManager::InitializeTestEvents() {

  using namespace Params;
  if (Sys::test_mode_ == "xlink_bind_ii") {
    // KMC Event -- Bind_II
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
      return head->parent_->GetWeight_Bind_II();
    };
    auto exe_bind_ii = [&](Object *base_head) {
      auto bound_head{dynamic_cast<BindingHead *>(base_head)};
      auto head{bound_head->GetOtherHead()};
      auto site{head->parent_->GetNeighbor_Bind_II()};
      auto executed{head->parent_->Bind(site, head)};
      if (executed) {
        xlinks_.FlagForUpdate();
        filaments_->FlagForUpdate();
        bool still_attached{head->parent_->UpdateExtension()};
        if (!still_attached) {
          return;
        }
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
        // bool head->parent_->UpdateExtension();
      } else {
        Sys::ErrorExit("Unbind_II (TEST)");
      }
    };
    kmc_.events_.emplace_back(
        "unbind_ii", 100 * xlinks_.p_event_.at("unbind_ii").GetVal(),
        &xlinks_.sorted_.at("unbind_ii").size_,
        &xlinks_.sorted_.at("unbind_ii").entries_, poisson_unbind_ii,
        get_weight_unbind_ii, exe_unbind_ii);
  } else if (Sys::test_mode_ == "filament_separation") {
    // Poisson distribution; sampled to predict events w/ variable probabilities
    auto poisson = [&](double p, int n) {
      if (p > 0.0) {
        return SysRNG::SamplePoisson(p);
      } else {
        return 0;
      }
    };
    auto is_doubly_bound = [&](Object *protein) -> Vec<Object *> {
      if (protein->GetNumHeadsActive() == 2) {
        return {protein->GetHeadOne(), protein->GetHeadTwo()};
      }
      return {};
    };
    xlinks_.AddPop("diffuse_ii_to_rest", is_doubly_bound);
    xlinks_.AddPop("diffuse_ii_fr_rest", is_doubly_bound);
    auto exe_diffuse_fwd = [&](Object *base) {
      auto head{dynamic_cast<BindingHead *>(base)};
      bool executed{head->Diffuse(1)};
      if (executed) {
        bool still_attached{head->parent_->UpdateExtension()};
        if (!still_attached) {
        }
        // FIXME had to move this from if statement above -- why ?
        xlinks_.FlagForUpdate();
        filaments_->FlagForUpdate();
      }
    };
    auto exe_diffuse_bck = [&](Object *base) {
      auto head{dynamic_cast<BindingHead *>(base)};
      bool executed{head->Diffuse(-1)};
      if (executed) {
        bool still_attached{head->parent_->UpdateExtension()};
        if (!still_attached) {
        }
        // FIXME had to move this from if statement above -- why ?
        xlinks_.FlagForUpdate();
        filaments_->FlagForUpdate();
      }
    };
    auto get_weight_diff_ii_to = [](Object *base) {
      auto head{dynamic_cast<BindingHead *>(base)};
      return head->GetWeight_Diffuse(1);
    };
    auto get_weight_diff_ii_fr = [](Object *base) {
      auto head{dynamic_cast<BindingHead *>(base)};
      return head->GetWeight_Diffuse(-1);
    };
    kmc_.events_.emplace_back(
        "diffuse_ii_to_rest",
        xlinks_.p_event_.at("diffuse_ii_to_rest").GetVal(),
        &xlinks_.sorted_.at("diffuse_ii_to_rest").size_,
        &xlinks_.sorted_.at("diffuse_ii_to_rest").entries_, poisson,
        get_weight_diff_ii_to, exe_diffuse_fwd);
    kmc_.events_.emplace_back(
        "diffuse_ii_fr_rest",
        xlinks_.p_event_.at("diffuse_ii_fr_rest").GetVal(),
        &xlinks_.sorted_.at("diffuse_ii_fr_rest").size_,
        &xlinks_.sorted_.at("diffuse_ii_fr_rest").entries_, poisson,
        get_weight_diff_ii_fr, exe_diffuse_bck);
  }
  /*
  for (auto const &event : kmc_.events_) {
    printf("%s: %g\n", event.name_.c_str(), event.p_occur_);
  }
  */
}

void ProteinManager::InitializeEvents() {

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
  auto exe_bind_i = [&](auto *site, auto *pop, auto *fil) {
    if (Sys::i_step_ < pop->step_active_) {
      return;
    }
    auto entry{pop->GetFreeEntry()};
    bool executed{entry->Bind(site, &entry->head_one_)};
    if (executed) {
      pop->AddToActive(entry);
      fil->FlagForUpdate();
    }
  };
  auto weight_bind_i = [](auto *site) { return site->GetWeight_Bind(); };
  auto is_unocc = [](Object *site) -> Vec<Object *> {
    if (!site->IsOccupied()) {
      return {site};
    }
    return {};
  };
  Vec<int> i_min{0, 0, 0};
  Vec<size_t> dim_size{1, 1, _n_neighbs_max + 1};
  auto get_n_neighbs = [](Object *entry) {
    Vec<int> indices_vec{entry->GetNumNeighborsOccupied()};
    return indices_vec;
  };
  if (xlinks_.active_) {
    filaments_->AddPop("xlinks", is_unocc, dim_size, i_min, get_n_neighbs);
    for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
      kmc_.events_.emplace_back(
          "bind_i", xlinks_.p_event_.at("bind_i").GetVal(n_neighbs),
          &filaments_->unoccupied_.at("xlinks").bin_size_[0][0][n_neighbs],
          &filaments_->unoccupied_.at("xlinks").bin_entries_[0][0][n_neighbs],
          binomial, [&](Object *base) {
            exe_bind_i(dynamic_cast<BindingSite *>(base), &xlinks_, filaments_);
          });
    }
  }
  if (motors_.active_) {
    filaments_->AddPop("motors", is_unocc);
    kmc_.events_.emplace_back(
        "bind_i", motors_.p_event_.at("bind_i").GetVal(),
        &filaments_->unoccupied_.at("motors").size_,
        &filaments_->unoccupied_.at("motors").entries_, poisson,
        [&](Object *base) {
          return weight_bind_i(dynamic_cast<BindingSite *>(base));
        },
        [&](Object *base) {
          exe_bind_i(dynamic_cast<BindingSite *>(base), &motors_, filaments_);
        });
  }
  // Bind_I_Teth

  // Bind_II
  auto exe_bind_ii = [](auto *bound_head, auto *pop, auto *fil) {
    // auto bound_head{dynamic_cast<BindingHead *>(base_head)};
    auto head{bound_head->GetOtherHead()};
    auto site{head->parent_->GetNeighbor_Bind_II()};
    auto executed{head->parent_->Bind(site, head)};
    if (executed) {
      bool still_attached{head->parent_->UpdateExtension()};
      pop->FlagForUpdate();
      fil->FlagForUpdate();
    }
  };
  auto weight_bind_ii = [](auto *head) {
    return head->parent_->GetWeight_Bind_II();
  };
  auto is_singly_bound = [](Object *protein) -> Vec<Object *> {
    if (protein->GetNumHeadsActive() == 1) {
      return {protein->GetActiveHead()};
    }
    return {};
  };
  if (xlinks_.crosslinking_active_) {
    xlinks_.AddPop("bind_ii", is_singly_bound);
    kmc_.events_.emplace_back(
        "bind_ii", xlinks_.p_event_.at("bind_ii").GetVal(),
        &xlinks_.sorted_.at("bind_ii").size_,
        &xlinks_.sorted_.at("bind_ii").entries_, poisson,
        [&](Object *base) {
          return weight_bind_ii(dynamic_cast<BindingHead *>(base));
        },
        [&](Object *base) {
          exe_bind_ii(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_);
        });
  }
  if (motors_.active_) {
    auto is_docked = [](auto *motor) -> Vec<Object *> {
      auto *docked_head{motor->GetDockedHead()};
      if (docked_head != nullptr) {
        return {docked_head->GetOtherHead()};
      }
      return {};
    };
    motors_.AddPop("bind_ii", [&](Object *base) {
      return is_docked(dynamic_cast<Motor *>(base));
    });
    kmc_.events_.emplace_back(
        "bind_ii", motors_.p_event_.at("bind_ii").GetVal(),
        &motors_.sorted_.at("bind_ii").size_,
        &motors_.sorted_.at("bind_ii").entries_, poisson,
        [&](Object *base) {
          return weight_bind_ii(dynamic_cast<CatalyticHead *>(base));
        },
        [&](Object *base) {
          exe_bind_ii(dynamic_cast<CatalyticHead *>(base), &motors_,
                      filaments_);
        });
  }
  // Unbind_II
  auto exe_unbind_ii = [](auto *head, auto *pop, auto *fil) {
    bool executed{head->Unbind()};
    if (executed) {
      pop->FlagForUpdate();
      fil->FlagForUpdate();
    }
  };
  auto weight_unbind_ii = [](auto *head) {
    return head->GetWeight_Unbind_II();
  };
  auto is_doubly_bound = [](Object *protein) -> Vec<Object *> {
    if (protein->GetNumHeadsActive() == 2) {
      return {protein->GetHeadOne(), protein->GetHeadTwo()};
    }
    return {};
  };
  if (xlinks_.crosslinking_active_) {
    xlinks_.AddPop("unbind_ii", is_doubly_bound);
    kmc_.events_.emplace_back(
        "unbind_ii", xlinks_.p_event_.at("unbind_ii").GetVal(),
        &xlinks_.sorted_.at("unbind_ii").size_,
        &xlinks_.sorted_.at("unbind_ii").entries_, poisson,
        [&](Object *base) {
          return weight_unbind_ii(dynamic_cast<BindingHead *>(base));
        },
        [&](Object *base) {
          exe_unbind_ii(dynamic_cast<BindingHead *>(base), &xlinks_,
                        filaments_);
        });
  }
  if (motors_.active_) {
    auto is_ADPP_ii_bound = [](auto *motor) -> Vec<Object *> {
      if (motor->n_heads_active_ == 2) {
        bool found_head{false};
        CatalyticHead *chosen_head{nullptr};
        if (motor->head_one_.ligand_ == CatalyticHead::Ligand::ADPP) {
          chosen_head = &motor->head_one_;
        }
        if (motor->head_two_.ligand_ == CatalyticHead::Ligand::ADPP) {
          if (chosen_head != nullptr) {
            Sys::ErrorExit("Protein_MGR::is_ADPP_ii_bound()");
          }
          chosen_head = &motor->head_two_;
        }
        if (chosen_head != nullptr) {
          return {chosen_head};
        }
      }
      return {};
    };
    motors_.AddPop("unbind_ii", [&](Object *base) {
      return is_ADPP_ii_bound(dynamic_cast<Motor *>(base));
    });
    kmc_.events_.emplace_back(
        "unbind_ii", motors_.p_event_.at("unbind_ii").GetVal(),
        &motors_.sorted_.at("unbind_ii").size_,
        &motors_.sorted_.at("unbind_ii").entries_, poisson,
        [&](Object *base) {
          return weight_unbind_ii(dynamic_cast<CatalyticHead *>(base));
        },
        [&](Object *base) {
          exe_unbind_ii(dynamic_cast<CatalyticHead *>(base), &motors_,
                        filaments_);
        });
  }
  // Unbind_I: Unbind first (singly bound) head of a protein
  auto exe_unbind_i = [](auto *head, auto *pop, auto *fil) {
    bool executed{head->Unbind()};
    if (executed) {
      head->UntetherSatellite();
      pop->RemoveFromActive(head->parent_);
      fil->FlagForUpdate();
    }
  };
  if (xlinks_.active_) {
    xlinks_.AddPop("bound_i", is_singly_bound, dim_size, i_min, get_n_neighbs);
    for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
      kmc_.events_.emplace_back(
          "unbind_i", xlinks_.p_event_.at("unbind_i").GetVal(n_neighbs),
          &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
          &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs],
          binomial, [&](Object *base) {
            exe_unbind_i(dynamic_cast<BindingHead *>(base), &xlinks_,
                         filaments_);
          });
    }
  }
  if (motors_.active_) {
    auto weight_unbind_i = [](auto *head) {
      return head->parent_->GetWeight_Unbind_I();
    };
    auto is_ADPP_i_bound = [](auto *motor) -> Vec<Object *> {
      if (motor->n_heads_active_ == 1) {
        if (motor->GetActiveHead()->ligand_ == CatalyticHead::Ligand::ADPP) {
          return {motor->GetActiveHead()};
        }
      }
      return {};
    };
    motors_.AddPop("bound_i_ADPP", [&](Object *base) {
      return is_ADPP_i_bound(dynamic_cast<Motor *>(base));
    });
    kmc_.events_.emplace_back(
        "unbind_i", motors_.p_event_.at("unbind_i").GetVal(),
        &motors_.sorted_.at("bound_i_ADPP").size_,
        &motors_.sorted_.at("bound_i_ADPP").entries_, poisson,
        [&](Object *base) {
          return weight_unbind_i(dynamic_cast<CatalyticHead *>(base));
        },
        [&](Object *base) {
          exe_unbind_i(dynamic_cast<CatalyticHead *>(base), &motors_,
                       filaments_);
        });
  }
  // Unbind_I_Teth
  // Tether_Free
  // Untether_Free

  // vv motors only vv
  // Bind_ATP
  if (motors_.active_) {
    auto exe_bind_ATP = [](auto *head, auto *pop) {
      // printf("boop\n");
      bool executed{head->parent_->Bind_ATP(head)};
      if (executed) {
        pop->FlagForUpdate();
      }
    };
    auto is_NULL_i_bound = [](auto *motor) -> Vec<Object *> {
      if (motor->n_heads_active_ == 1) {
        if (motor->GetActiveHead()->ligand_ == CatalyticHead::Ligand::NONE) {
          return {motor->GetActiveHead()};
        }
      }
      return {};
    };
    motors_.AddPop("bound_i_NULL", [&](Object *base) {
      return is_NULL_i_bound(dynamic_cast<Motor *>(base));
    });
    kmc_.events_.emplace_back(
        "bind_ATP_i", motors_.p_event_.at("bind_ATP_i").GetVal(),
        &motors_.sorted_.at("bound_i_NULL").size_,
        &motors_.sorted_.at("bound_i_NULL").entries_, binomial,
        [&](Object *base) {
          exe_bind_ATP(dynamic_cast<CatalyticHead *>(base), &motors_);
        });
    // auto is_ATP_ii_bound = [](auto *motor) {};
  }
  // Hydrolyze_ATP
  if (motors_.active_) {
    auto exe_hydrolyze = [](auto *head, auto *pop) {
      bool executed{head->parent_->Hydrolyze(head)};
      if (executed) {
        pop->FlagForUpdate();
      }
    };
    auto is_ATP_i_bound = [](auto *motor) -> Vec<Object *> {
      if (motor->n_heads_active_ == 1) {
        if (motor->GetActiveHead()->ligand_ == CatalyticHead::Ligand::ATP) {
          return {motor->GetActiveHead()};
        }
      }
      return {};
    };
    motors_.AddPop("bound_i_ATP", [&](Object *base) {
      return is_ATP_i_bound(dynamic_cast<Motor *>(base));
    });
    kmc_.events_.emplace_back(
        "hydrolyze", motors_.p_event_.at("hydrolyze").GetVal(),
        &motors_.sorted_.at("bound_i_ATP").size_,
        &motors_.sorted_.at("bound_i_ATP").entries_, binomial,
        [&](Object *base) {
          exe_hydrolyze(dynamic_cast<CatalyticHead *>(base), &motors_);
        });
  }
  // Tether_Bound
  // Untether_Bound

  // vv xlinks only vv
  // Diffusion
  auto exe_diff = [](auto *head, auto *pop, auto *fil, int dir) {
    bool executed{head->Diffuse(dir)};
    if (executed) {
      bool still_attached{head->parent_->UpdateExtension()};
      if (!still_attached) {
        pop->FlagForUpdate();
      }
      fil->FlagForUpdate();
    }
  };
  if (xlinks_.active_) {
    for (int n_neighbs{0}; n_neighbs < _n_neighbs_max; n_neighbs++) {
      kmc_.events_.emplace_back(
          "diffuse_i_fwd",
          xlinks_.p_event_.at("diffuse_i_fwd").GetVal(n_neighbs),
          &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
          &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs],
          binomial, [&](Object *base) {
            exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_,
                     1);
          });
      kmc_.events_.emplace_back(
          "diffuse_i_bck",
          xlinks_.p_event_.at("diffuse_i_bck").GetVal(n_neighbs),
          &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
          &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs],
          binomial, [&](Object *base) {
            exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_,
                     -1);
          });
    }
  }
  if (xlinks_.crosslinking_active_) {
    auto weight_diff_ii = [](auto *head, int dir) {
      return head->GetWeight_Diffuse(dir);
    };
    xlinks_.AddPop("diffuse_ii_to_rest", is_doubly_bound);
    xlinks_.AddPop("diffuse_ii_fr_rest", is_doubly_bound);
    kmc_.events_.emplace_back(
        "diffuse_ii_to_rest",
        xlinks_.p_event_.at("diffuse_ii_to_rest").GetVal(),
        &xlinks_.sorted_.at("diffuse_ii_to_rest").size_,
        &xlinks_.sorted_.at("diffuse_ii_to_rest").entries_, poisson,
        [&](Object *base) {
          return weight_diff_ii(dynamic_cast<BindingHead *>(base), 1);
        },
        [&](Object *base) {
          exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_, 1);
        });
    kmc_.events_.emplace_back(
        "diffuse_ii_fr_rest",
        xlinks_.p_event_.at("diffuse_ii_fr_rest").GetVal(),
        &xlinks_.sorted_.at("diffuse_ii_fr_rest").size_,
        &xlinks_.sorted_.at("diffuse_ii_fr_rest").entries_, poisson,
        [&](Object *base) {
          return weight_diff_ii(dynamic_cast<BindingHead *>(base), -1);
        },
        [&](Object *base) {
          exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_, -1);
        });
  }
  // Bind_II_Teth
  // Unbind_II_Teth
  for (auto const &event : kmc_.events_) {
    printf("%s: %g\n", event.name_.c_str(), event.p_occur_);
  }
}

void ProteinManager::FlagFilamentsForUpdate() { filaments_->FlagForUpdate(); }

void ProteinManager::UpdateFilaments() { filaments_->UpdateUnoccupied(); }
