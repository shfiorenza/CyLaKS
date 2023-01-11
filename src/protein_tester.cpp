#include "cylaks/protein_tester.hpp"
#include "cylaks/filament_tester.hpp"

void ProteinTester::Initialize(FilamentTester *filaments) {
  filaments_ = filaments;
  ProteinManager::filaments_ = dynamic_cast<FilamentManager *>(filaments_);
  SetTestMode();
  CheckProbabilities();
}

void ProteinTester::UpdateFilaments() {
  filaments_->UpdateUnoccupied();
  if (Sys::test_mode_ != "filament_ablation") {
    return;
  }
  if (Sys::i_step_ == Sys::ablation_step_) {
    filaments_->protofilaments_[1].pos_[0] += 200.0;
    filaments_->protofilaments_[1].ForceUpdate();
  }
}

void ProteinTester::ReportTestStatistics() {
  for (auto const &entry : test_stats_) {
    Sys::Log("For event %s:\n", entry.first.c_str());
    for (int index{0}; index < entry.second.size(); index++) {
      auto stats = entry.second[index];
      double p{double(stats.first) / stats.second};
      double ref{test_ref_.at(entry.first)[index]};
      Sys::Log("  p[%i] = %.3g (%.3g expected) [%zu / %zu events]\n", index, p,
               ref, stats.first, stats.second);
    }
  }
}

void ProteinTester::SetTestMode() {
  Sys::Log("\n");
  Sys::Log("Initializing test '%s'. Overidden parameters listed below:\n",
           Sys::test_mode_.c_str());
  Sys::Log("(Capitalization is for automated detection by .m scipts)\n");
  if (Sys::test_mode_ == "filament_ablation") {
    InitializeTest_Filament_Ablation();
  } else if (Sys::test_mode_ == "filament_separation") {
    InitializeTest_Filament_Separation();
  } else if (Sys::test_mode_ == "filament_forced_slide") {
    InitializeTest_Filament_ForcedSlide();
  } else if (Sys::test_mode_ == "hetero_tubulin") {
    InitializeTest_Filament_HeteroTubulin();
  } else if (Sys::test_mode_ == "kinesin_mutant") {
    InitializeTest_Motor_Heterodimer();
  } else if (Sys::test_mode_ == "motor_lattice_step") {
    InitializeTest_Motor_LatticeStep();
  } else if (Sys::test_mode_ == "motor_lattice_bind") {
    InitializeTest_Motor_LatticeBind();
  } else if (Sys::test_mode_ == "xlink_diffusion") {
    InitializeTest_Xlink_Diffusion();
  } else if (Sys::test_mode_ == "xlink_bind_ii") {
    InitializeTest_Xlink_Bind_II();
  } else if (Sys::test_mode_ == "shepherding") {
    InitializeTest_Shepherding();
  } else {
    Sys::ErrorExit("ProteinTester::SetTestMode()");
  }
}

void ProteinTester::InitializeTest_Filament_Ablation() {
  using namespace Params;
  Sys::OverrideParam("filaments. COUNT", &Filaments::count, 2);
  Sys::OverrideParam("filaments. N_SITES[0]", &Filaments::n_sites[0], 875);
  Sys::OverrideParam("filaments. N_SITES[1]", &Filaments::n_sites[1], 875);
  Sys::OverrideParam("filaments. POLARITY[0]", &Filaments::polarity[0], 0);
  Sys::OverrideParam("filaments. POLARITY[1]", &Filaments::polarity[1], 0);
  Sys::OverrideParam("filaments. x_initial[0]", &Filaments::x_initial[0], 0.0);
  Sys::OverrideParam("filaments. x_initial[1]", &Filaments::x_initial[1],
                     (Filaments::n_sites[0] - 1) * Filaments::site_size);
  Sys::OverrideParam("filaments. y_initial[0]", &Filaments::y_initial[0], 0.0);
  Sys::OverrideParam("filaments. y_initial[1]", &Filaments::y_initial[1], 0.0);
  // ! FIXME update if new immobile_until syntax is permanently adopted
  // Sys::Log("All filament movement has been disabled by default.\n");
  // Filaments::translation_enabled[0] = false;
  // Filaments::translation_enabled[1] = false;
  // Filaments::rotation_enabled = false;
  Sys::Log("Motor end-pausing has been activated by default.\n");
  Motors::endpausing_active = true;
  printf("Enter ablation time: ");
  Str response;
  std::getline(std::cin, response);
  double t_ablate{(double)std::stod(response)};
  Sys::ablation_step_ = size_t(std::round(t_ablate / dt));
  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  filaments_->Initialize(this);
  InitializeEvents();
}

void ProteinTester::InitializeTest_Filament_Separation() {
  using namespace Params;
  // Initialize filament environment
  if (Filaments::n_sites[0] != Filaments::n_sites[1]) {
    printf("\nError! Filaments must be the same length.\n");
    exit(1);
  }
  Sys::OverrideParam("motors: c_bulk", &Motors::c_bulk, 0.0);
  Sys::OverrideParam("xlinks: c_bulk", &Xlinks::c_bulk, 1.0);
  // ! FIXME update if new immobile_until syntax is permanently adopted
  // Sys::OverrideParam("filaments. immobile_until[0]",
  //                    &Filaments::immobile_until[0], 0.0);
  // Sys::OverrideParam("filaments. immobile_until[1]",
  //                    &Filaments::immobile_until[1], 0.0);
  // Sys::Log("Horizontal and rotational filament movement has been
  // disabled.\n"); Filaments::translation_enabled[0] = false;
  // Filaments::translation_enabled[1] = true;
  // Filaments::rotation_enabled = false;
  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  filaments_->Initialize(this);
  printf("\nRunning interactive launcher for 'filament_separation'.\n");
  printf("You can also use the following quick-launch syntax:\n\n");
  printf("  %s params.yaml sim_name filament_separation n_xlinks\n\n",
         Sys::exe_name_.c_str());
  printf("where n_xlinks is the number of doubly bound crosslinkers to "
         "be inserted.\n\n");
  int n_xlinks{Sys::n_xlinks_};
  if (n_xlinks == -1) {
    Str response;
    printf("Microtubules are %zu sites in length.\n", Filaments::n_sites[0]);
    // for (int i_pf{0}; i_pf < filaments_->protofilaments_.size(); i_pf++) {
    //   printf("%zu\n", filaments_->protofilaments_[i_pf].sites_.size());
    // }
    printf("Enter number of crosslinkers to insert: ");
    std::getline(std::cin, response);
    n_xlinks = (int)std::stoi(response);
  }
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
    BindingSite *site_one{&filaments_->protofilaments_[0].sites_[i_site]};
    BindingSite *site_two{&filaments_->protofilaments_[1].sites_[i_site]};
    bool exe_one{xlink->Bind(site_one, &xlink->head_one_)};
    bool exe_two{xlink->Bind(site_two, &xlink->head_two_)};
    if (exe_one and exe_two) {
      bool still_attached{xlink->UpdateExtension()};
      if (still_attached) {
        xlinks_.AddToActive(xlink);
        filaments_->FlagForUpdate();
      } else {
        Sys::ErrorExit("ProteinTester::InitializeTestEnvironment() [2]");
      }
    } else {
      Sys::ErrorExit("ProteinTester::InitializeTestEnvironment() [1]");
    }
  }
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
  auto weight_diff_ii = [](auto *head, int dir) {
    return head->GetWeight_Diffuse(dir);
  };
  auto exe_diff = [](auto *head, auto *pop, auto *fil, int dir) {
    bool executed{head->Diffuse(dir)};
    if (executed) {
      bool still_attached{head->parent_->UpdateExtension()};
    }
    return executed;
  };
  kmc_.events_.emplace_back(
      "diffuse_ii_to_rest", xlinks_.p_event_.at("diffuse_ii_to_rest").GetVal(),
      &xlinks_.sorted_.at("diffuse_ii_to_rest").size_,
      &xlinks_.sorted_.at("diffuse_ii_to_rest").entries_, poisson,
      [&](Object *base) {
        return weight_diff_ii(dynamic_cast<BindingHead *>(base), 1);
      },
      [&](Object *base) {
        return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_,
                        1);
      });
  kmc_.events_.emplace_back(
      "diffuse_ii_fr_rest", xlinks_.p_event_.at("diffuse_ii_fr_rest").GetVal(),
      &xlinks_.sorted_.at("diffuse_ii_fr_rest").size_,
      &xlinks_.sorted_.at("diffuse_ii_fr_rest").entries_, poisson,
      [&](Object *base) {
        return weight_diff_ii(dynamic_cast<BindingHead *>(base), -1);
      },
      [&](Object *base) {
        return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_,
                        -1);
      });
}

void ProteinTester::InitializeTest_Filament_ForcedSlide() {

  using namespace Params;
  // Initialize filament environment
  Sys::OverrideParam("motors: c_bulk", &Motors::c_bulk, 0.0);
  // Sys::OverrideParam("xlinks: c_bulk", &Xlinks::c_bulk, 1.0);
  // double pos_to_alignment{(Filaments::n_sites[0] - Filaments::n_sites[1]) *
  //                         Filaments::site_size / 2};
  // Sys::OverrideParam("filaments: x_initial[0]", &Filaments::x_initial[0],
  //                    pos_to_alignment);
  // Sys::OverrideParam("filaments. immobile_until[0]",
  //  &Filaments::immobile_until[0], 100000.0);
  // Sys::OverrideParam("filaments. immobile_until[1]",
  //                    &Filaments::immobile_until[1], 0.0);
  // Sys::Log("Rotational filament movement has been disabled.\n");
  // Filaments::translation_enabled[0] = true;
  // Filaments::translation_enabled[1] = true;
  // Filaments::rotation_enabled = false;
  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  filaments_->Initialize(this);
  printf("\nRunning interactive launcher for 'filament_forced_slide'.\n");
  printf("You can also use the following quick-launch syntax:\n\n");
  printf("  %s params.yaml sim_name filament_forced_slide "
         "n_xlinks slide_velocity mode_flag\n",
         Sys::exe_name_.c_str());
  printf("       OR\n");
  printf("  %s params.yaml sim_name filament_forced_slide "
         "n_xlinks slide_velocity mode_flag t_pause pause_duration\n",
         Sys::exe_name_.c_str());
  printf("\nwhere n_xlinks is the number of doubly bound crosslinkers to "
         "be inserted,\nslide_velocity is the imposed average velocity "
         "on microtubules, and\nmode_flag is set to 0 for a constant "
         "force assay or 1 for a constant velocity assay.\n");
  printf("t_pause and pause_duration (in seconds) are optional and allow for a "
         "temporary pause in the applied force at a specified time.\n\n");
  // Get number of crosslinkers to initialize in the starting overlap
  int n_xlinks{Sys::n_xlinks_};
  if (n_xlinks == -1) {
    Str response;
    printf("Top microtubule is %zu sites in length.\n", Filaments::n_sites[1]);
    printf("Enter number of crosslinkers to insert: ");
    std::getline(std::cin, response);
    n_xlinks = (int)std::stoi(response);
    Sys::n_xlinks_ = n_xlinks;
  }
  int n_places{(int)Filaments::n_sites[1]};
  if (n_places > (int)Filaments::n_sites[0]) {
    n_places = Filaments::n_sites[0];
  }
  if (n_xlinks > n_places) {
    printf("\nError! Too many crosslinkers for filament length used.\n");
    exit(1);
  }
  // Get desired sliding velocity for top microtubule (bottom is fixed)
  if (Sys::slide_velocity_ == -1.0) {
    Str response;
    printf("Input 0 for constant force or 1 for constant velocity mode.\n");
    std::getline(std::cin, response);
    int vel_flag{std::stoi(response)};
    if (vel_flag == 1) {
      Sys::constant_velocity_ = true;
    } else if (vel_flag == 0) {
      Sys::constant_velocity_ = false;
    } else {
      printf("Error. Flag must be 0 or 1.\n");
      exit(1);
    }
    printf("Enter constant (or average) sliding velocity (in nm/s) desired for "
           "top MT.\n");
    std::getline(std::cin, response);
    Sys::slide_velocity_ = (double)std::stod(response);
  }
  if (Sys::slide_velocity_ >= 1000) {
    printf("\nError! Please keep velocity below 1 um/s (1000 nm/s).\n");
    exit(1);
  }
  /*
  // Check if binding/unbinding of crosslinkers should be enabled
  if (Sys::binding_active_ == -1) {
    Str response;
    printf("Enable binding/unbinding to/from solution? (y/n)\n");
    std::getline(std::cin, response);
    if (response == "n" or response == "N") {
      Sys::binding_active_ = 0;
    } else if (response == "y" or response == "Y") {
      Sys::binding_active_ = 1;
    } else {
      printf("Invalid response. Please choose y (yes) or n (no).\n");
      exit(1);
    }
  }
  */
  // Check if a pause in force-clamping is desired
  if (Sys::i_pause_ == -1) {
    Str response;
    printf("Enable pause in force clamp?\n");
    std::getline(std::cin, response);
    if (response == "n" or response == "N") {
      printf("No pause scheduled.\n");
      Sys::i_pause_ = std::numeric_limits<int>::max();
      Sys::i_resume_ = std::numeric_limits<int>::max();
    } else if (response == "y" or response == "Y") {
      Str str_pause;
      printf(
          "Enter time at which force clamp should be paused. (in seconds)\n");
      std::getline(std::cin, str_pause);
      double t_pause{(double)std::stod(str_pause)};
      Str str_duration;
      printf("Enter duration of pause.(in seconds)\n");
      std::getline(std::cin, str_duration);
      double t_resume{t_pause + (double)std::stod(str_duration)};
      if (t_resume > Params::t_run) {
        printf("Error. Pause must occur and end before simulation end.\n");
        exit(1);
      }
      Sys::i_pause_ = (int)std::round(t_pause / Params::dt);
      Sys::i_resume_ = (int)std::round(t_resume / Params::dt);
      printf("i_pause is %i\n", Sys::i_pause_);
      printf("i_resume is %i\n", Sys::i_resume_);
    } else {
      printf("Invalid response. Please choose y (yes) or n (no).\n");
      exit(1);
    }
  }
  // Rescale pause/resume variables if initialized via quick-launcher
  if (Sys::rescale_times_) {
    double t_pause{(double)Sys::i_pause_};
    double t_duration{(double)Sys::i_resume_};
    double t_resume{t_pause + t_duration};
    if (t_resume > Params::t_run) {
      printf("Error. Pause must occur and end before simulation end.\n");
      exit(1);
    }
    Sys::i_pause_ = (int)std::round(t_pause / Params::dt);
    Sys::i_resume_ = (int)std::round(t_resume / Params::dt);
    Sys::rescale_times_ = false;
  }

  // Randomly place crosslinkers on filaments w/ x = 0
  int site_indices[n_places];
  for (int index{0}; index < n_places; index++) {
    site_indices[index] = index;
  }
  SysRNG::Shuffle(site_indices, n_places, sizeof(int));
  int i_mt_short{1};
  int i_mt_long{0};
  if (Filaments::n_sites[1] > Filaments::n_sites[0]) {
    i_mt_short = 0;
    i_mt_long = 1;
  }
  for (int i_xlink{0}; i_xlink < n_xlinks; i_xlink++) {
    Protein *xlink{xlinks_.GetFreeEntry()};
    int i_site{site_indices[i_xlink]};
    BindingSite *site_one{
        &filaments_->protofilaments_[i_mt_short].sites_[i_site]};
    // printf("site one: %zu\n", site_one->index_);
    BindingSite *site_two{
        filaments_->protofilaments_[i_mt_long].GetNeighb(site_one, 0)};
    // printf("site two: %zu\n", site_two->index_);
    bool exe_one{xlink->Bind(site_one, &xlink->head_one_)};
    bool exe_two{xlink->Bind(site_two, &xlink->head_two_)};
    if (exe_one and exe_two) {
      bool still_attached{xlink->UpdateExtension()};
      if (still_attached) {
        xlinks_.AddToActive(xlink);
        filaments_->FlagForUpdate();
      } else {
        Sys::ErrorExit("ProteinTester::InitializeTestEnvironment() [2]");
      }
    } else {
      Sys::ErrorExit("ProteinTester::InitializeTestEnvironment() [1]");
    }
  }
  Sys::Log("%i crosslinkers initialized.\n", n_xlinks);
  // Probability distributions
  auto poisson = [&](double p, int n) { // For events w/ time-variable probs.
    if (p > 0.0) {
      return SysRNG::SamplePoisson(p);
    } else {
      return 0;
    }
  };
  auto binomial = [&](double p, int n) { // For events w/ constant probs.
    if (n > 0) {
      return SysRNG::SampleBinomial(p, n);
    } else {
      return 0;
    }
  };
  // Binning helper functions
  Vec<size_t> dim_size{1, 1, _n_neighbs_max + 1};
  Vec<int> i_min{0, 0, 0};
  auto get_n_neighbs = [](Object *entry) {
    Vec<int> indices_vec{entry->GetNumNeighborsOccupied()};
    return indices_vec;
  };
  // Sorting criteria functions
  auto is_unocc = [](Object *base) -> Vec<Object *> {
    BindingSite *site{dynamic_cast<BindingSite *>(base)};
    if (!site->IsOccupied()) {
      return {site};
    }
    return {};
  };
  auto is_bound_i = [](Object *base) -> Vec<Object *> {
    Protein *protein{dynamic_cast<Protein *>(base)};
    if (protein->GetNumHeadsActive() == 1) {
      if (!protein->IsTethered() or protein->HasSatellite()) {
        return {protein->GetActiveHead()};
      }
    }
    return {};
  };
  auto is_bound_ii = [&](Object *protein) -> Vec<Object *> {
    if (protein->GetNumHeadsActive() == 2) {
      return {protein->GetHeadOne(), protein->GetHeadTwo()};
    }
    return {};
  };
  // Event execution functions
  auto exe_diff = [](auto *head, auto *pop, auto *fil, int dir) {
    if ((head->site_ == head->site_->filament_->plus_end_ or
         head->site_ == head->site_->filament_->minus_end_) and
        head->parent_->GetNumHeadsActive() == 2 and
        Params::Xlinks::p_diffuse_off_end > 0.0) {
      double ran{SysRNG::GetRanProb()};
      if (ran < Params::Xlinks::p_diffuse_off_end) {
        bool executed{head->Unbind()};
        return executed;
      }
      return false;
    } else {
      bool executed{head->Diffuse(dir)};
      if (executed) {
        bool still_attached{head->parent_->UpdateExtension()};
      }
      return executed;
    }
  };
  // if (Sys::binding_active_) {
  // Bind from solution (stage 0 -> stage 1)
  // Add unoccupied site tracker for crosslinkers; segregated by n_neighbs
  filaments_->AddPop("neighbs", is_unocc, dim_size, i_min, get_n_neighbs);
  auto exe_bind_i = [&](auto *site, auto *pop, auto *fil) {
    if (Sys::i_step_ < pop->step_active_) {
      return false;
    }
    auto entry{pop->GetFreeEntry()};
    if (entry == nullptr) {
      return false;
    }
    bool executed{entry->Bind(site, entry->GetHeadOne())};
    if (executed) {
      pop->AddToActive(entry);
    }
    return executed;
  };
  // Create a binomial event for each n_neighb possibility
  for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
    kmc_.events_.emplace_back(
        "bind_i_" + std::to_string(n_neighbs) + " (xlinks)",
        xlinks_.p_event_.at("bind_i").GetVal(n_neighbs),
        &filaments_->unoccupied_.at("neighbs").bin_size_[0][0][n_neighbs],
        &filaments_->unoccupied_.at("neighbs").bin_entries_[0][0][n_neighbs],
        binomial, [&](Object *base) {
          return exe_bind_i(dynamic_cast<BindingSite *>(base), &xlinks_,
                            filaments_);
        });
  }
  // Singly bound unbinding to solution (stage 1 -> stage 0)
  xlinks_.AddPop("bound_i", is_bound_i, dim_size, i_min, get_n_neighbs);
  auto exe_unbind_i = [](auto *head, auto *pop, auto *alt_pop, auto *fil) {
    bool executed{head->Unbind()};
    if (executed) {
      if (head->parent_->HasSatellite()) {
        alt_pop->RemoveFromActive(head->parent_->teth_partner_);
        bool untethered_sat{head->UntetherSatellite()};
        if (!untethered_sat) {
          Sys::ErrorExit("exe_unbind_i");
        }
      }
      if (!head->parent_->IsTethered()) {
        pop->RemoveFromActive(head->parent_);
      }
    }
    return executed;
  };
  for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
    kmc_.events_.emplace_back(
        "unbind_i_" + std::to_string(n_neighbs) + " (xlinks)",
        xlinks_.p_event_.at("unbind_i").GetVal(n_neighbs),
        &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
        &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs], binomial,
        [&](Object *base) {
          return exe_unbind_i(dynamic_cast<BindingHead *>(base), &xlinks_,
                              &motors_, filaments_);
        });
  }
  // }
  // Doubly bound binding (stage 1 -> stage 2; singly bound -> crosslinking)
  auto weight_bind_ii = [](auto *head) {
    return head->parent_->GetWeight_Bind_II();
  };
  auto exe_bind_ii = [](auto *bound_head, auto *pop, auto *fil) {
    auto head{bound_head->GetOtherHead()};
    auto site{head->parent_->GetNeighbor_Bind_II()};
    auto executed{head->parent_->Bind(site, head)};
    if (executed) {
      // FIXME check if this check is necessary / what to do
      bool still_attached{head->parent_->UpdateExtension()};
      if (!still_attached) {
        return false;
      }
    }
    return executed;
  };
  xlinks_.AddPop("singly_bound", is_bound_i);
  kmc_.events_.emplace_back(
      "bind_ii (xlinks)", xlinks_.p_event_.at("bind_ii").GetVal(),
      &xlinks_.sorted_.at("singly_bound").size_,
      &xlinks_.sorted_.at("singly_bound").entries_, poisson,
      [&](Object *base) {
        return weight_bind_ii(dynamic_cast<BindingHead *>(base));
      },
      [&](Object *base) {
        return exe_bind_ii(dynamic_cast<BindingHead *>(base), &xlinks_,
                           filaments_);
      });
  // Doubly bound unbinding (stage 2 -> stage 1; crosslinking -> singly bound)
  auto weight_unbind_ii = [](auto *head) {
    return head->GetWeight_Unbind_II();
  };
  auto exe_unbind_ii = [](auto *head, auto *pop, auto *fil) {
    bool executed{head->Unbind()};
    return executed;
  };
  xlinks_.AddPop("doubly_bound", is_bound_ii);
  kmc_.events_.emplace_back(
      "unbind_ii (xlinks)", xlinks_.p_event_.at("unbind_ii").GetVal(),
      &xlinks_.sorted_.at("doubly_bound").size_,
      &xlinks_.sorted_.at("doubly_bound").entries_, poisson,
      [&](Object *base) {
        return weight_unbind_ii(dynamic_cast<BindingHead *>(base));
      },
      [&](Object *base) {
        return exe_unbind_ii(dynamic_cast<BindingHead *>(base), &xlinks_,
                             filaments_);
      });
  // Singly bound diffusion
  xlinks_.AddPop("bound_i", is_bound_i, dim_size, i_min, get_n_neighbs);
  for (int n_neighbs{0}; n_neighbs < _n_neighbs_max; n_neighbs++) {
    kmc_.events_.emplace_back(
        "diffuse_i_" + std::to_string(n_neighbs) + "_fwd",
        xlinks_.p_event_.at("diffuse_i_fwd").GetVal(n_neighbs),
        &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
        &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs], binomial,
        [&](Object *base) {
          return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_,
                          filaments_, 1);
        });
    kmc_.events_.emplace_back(
        "diffuse_i_" + std::to_string(n_neighbs) + "_bck",
        xlinks_.p_event_.at("diffuse_i_bck").GetVal(n_neighbs),
        &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
        &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs], binomial,
        [&](Object *base) {
          return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_,
                          filaments_, -1);
        });
  }
  // Doubly bound diffusion
  xlinks_.AddPop("diffuse_ii_to_rest", is_bound_ii);
  xlinks_.AddPop("diffuse_ii_fr_rest", is_bound_ii);
  auto weight_diff_ii = [](auto *head, int dir) {
    return head->GetWeight_Diffuse(dir);
  };
  kmc_.events_.emplace_back(
      "diffuse_ii_to_rest", xlinks_.p_event_.at("diffuse_ii_to_rest").GetVal(),
      &xlinks_.sorted_.at("diffuse_ii_to_rest").size_,
      &xlinks_.sorted_.at("diffuse_ii_to_rest").entries_, poisson,
      [&](Object *base) {
        return weight_diff_ii(dynamic_cast<BindingHead *>(base), 1);
      },
      [&](Object *base) {
        return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_,
                        1);
      });
  kmc_.events_.emplace_back(
      "diffuse_ii_fr_rest", xlinks_.p_event_.at("diffuse_ii_fr_rest").GetVal(),
      &xlinks_.sorted_.at("diffuse_ii_fr_rest").size_,
      &xlinks_.sorted_.at("diffuse_ii_fr_rest").entries_, poisson,
      [&](Object *base) {
        return weight_diff_ii(dynamic_cast<BindingHead *>(base), -1);
      },
      [&](Object *base) {
        return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_,
                        -1);
      });
}

void ProteinTester::InitializeTest_Filament_HeteroTubulin() {
  // ! FIXME update if new immobile_until syntax is permanently adopted
  // Sys::Log("All filament movement has been disabled by default.\n");
  // Params::Filaments::translation_enabled[0] = false;
  // Params::Filaments::translation_enabled[1] = false;
  // Params::Filaments::rotation_enabled = false;
  double p_hetero{Sys::p_mutant_};
  if (p_hetero == -1.0) {
    printf("Enter fraction of heterogenous tubulin: ");
    Str response_one;
    std::getline(std::cin, response_one);
    p_hetero = (double)std::stod(response_one);
  }
  if (p_hetero < 0.0 or p_hetero > 1.0) {
    printf("Error. Invalid fraction!\n");
    exit(1);
  }
  double bind_aff{Sys::binding_affinity_};
  if (bind_aff == -1.0) {
    printf("Enter decrease in binding affinity ");
    printf("(e.g., 2 will cut p_bind in half): ");
    Str response_two;
    std::getline(std::cin, response_two);
    bind_aff = (double)std::stod(response_two);
  }
  if (bind_aff <= 0.0) {
    printf("Error. Fractional change must be positive!\n");
    exit(1);
  }
  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  filaments_->Initialize(this);
  int n_sites{(int)filaments_->sites_.size()};
  int n_hetero{(int)std::round(n_sites * p_hetero)};
  // Randomly place heterogeneous sites on lattice
  int site_indices[n_sites];
  for (int index{0}; index < n_sites; index++) {
    site_indices[index] = index;
  }
  SysRNG::Shuffle(site_indices, n_sites, sizeof(int));
  for (int i_hetero{0}; i_hetero < n_hetero; i_hetero++) {
    int i_site{site_indices[i_hetero]};
    filaments_->sites_[i_site]->SetBindingAffinity(bind_aff);
  }
  InitializeEvents();
}

void ProteinTester::InitializeTest_Motor_Heterodimer() {
  // ! FIXME update if new immobile_until syntax is permanently adopted
  // Sys::Log("All filament movement has been disabled by default.\n");
  // Params::Filaments::translation_enabled[0] = false;
  // Params::Filaments::translation_enabled[1] = false;
  // Params::Filaments::rotation_enabled = false;
  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  filaments_->Initialize(this);
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
  // head_one_ is catalytic; head_two_ is passive
  auto exe_bind_i = [&](auto *site, auto *pop, auto *fil) {
    if (Sys::i_step_ < pop->step_active_) {
      return false;
    }
    auto entry{pop->GetFreeEntry()};
    Sys::Log("bound motor %zu\n", entry->GetID());
    // always bind catalytic head first
    bool executed{entry->Bind(site, &entry->head_one_)};
    if (executed) {
      pop->AddToActive(entry);
      fil->FlagForUpdate();
    }
    return executed;
  };
  auto weight_bind_i = [](auto *site) { return site->GetWeight_Bind(); };
  auto is_unocc = [](Object *site) -> Vec<Object *> {
    if (!site->IsOccupied()) {
      return {site};
    }
    return {};
  };
  filaments_->AddPop("motors", is_unocc);
  kmc_.events_.emplace_back(
      "bind_i", motors_.p_event_.at("bind_i").GetVal(),
      &filaments_->unoccupied_.at("motors").size_,
      &filaments_->unoccupied_.at("motors").entries_, poisson,
      [&](Object *base) {
        return weight_bind_i(dynamic_cast<BindingSite *>(base));
      },
      [&](Object *base) {
        return exe_bind_i(dynamic_cast<BindingSite *>(base), &motors_,
                          filaments_);
      });
  // Bind_II
  auto exe_bind_ii = [](auto *bound_head, auto *pop, auto *fil) {
    auto head{bound_head->GetOtherHead()};
    auto site{head->parent_->GetNeighbor_Bind_II()};
    auto executed{head->parent_->Bind(site, head)};
    if (executed) {
      bool still_attached{head->parent_->UpdateExtension()};
      pop->FlagForUpdate();
      fil->FlagForUpdate();
    }
    return executed;
  };
  auto weight_bind_ii = [](auto *head) {
    return head->parent_->GetWeight_Bind_II();
  };
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
        return exe_bind_ii(dynamic_cast<CatalyticHead *>(base), &motors_,
                           filaments_);
      });
  // Unbind_II
  auto exe_unbind_ii = [](auto *head, auto *pop, auto *fil) {
    bool executed{head->Unbind()};
    if (executed) {
      pop->FlagForUpdate();
      fil->FlagForUpdate();
    }
    return executed;
  };
  auto weight_unbind_ii = [](auto *head) {
    return head->GetWeight_Unbind_II();
  };
  auto is_ADPP_ii_bound = [](auto *motor) -> Vec<Object *> {
    if (motor->n_heads_active_ == 2) {
      // Always unbind active head first if both are ADPP bound
      if (motor->head_one_.ligand_ == Ligand::ADPP and
          motor->head_two_.ligand_ == Ligand::ADPP) {
        return {&motor->head_one_};
      }
      bool found_head{false};
      CatalyticHead *chosen_head{nullptr};
      if (motor->head_one_.ligand_ == Ligand::ADPP) {
        chosen_head = &motor->head_one_;
      }
      if (motor->head_two_.ligand_ == Ligand::ADPP) {
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
        return exe_unbind_ii(dynamic_cast<CatalyticHead *>(base), &motors_,
                             filaments_);
      });
  // Unbind_I: Unbind first (singly bound) head of a protein
  auto exe_unbind_i = [](auto *head, auto *pop, auto *fil) {
    bool executed{head->Unbind()};
    if (executed) {
      head->UntetherSatellite();
      pop->RemoveFromActive(head->parent_);
      fil->FlagForUpdate();
    }
    return executed;
  };
  auto weight_unbind_i = [](auto *head) {
    return head->parent_->GetWeight_Unbind_I();
  };
  auto is_ADPP_i_bound = [](auto *motor) -> Vec<Object *> {
    if (motor->n_heads_active_ == 1) {
      if (motor->GetActiveHead()->ligand_ == Ligand::ADPP) {
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
        return exe_unbind_i(dynamic_cast<CatalyticHead *>(base), &motors_,
                            filaments_);
      });
  // Bind_ATP
  auto exe_bind_ATP = [](auto *head, auto *pop) {
    bool executed{head->parent_->Bind_ATP(head)};
    if (executed) {
      pop->FlagForUpdate();
    }
    return executed;
  };
  auto is_NULL_i_bound = [](auto *motor) -> Vec<Object *> {
    if (motor->n_heads_active_ == 1) {
      if (motor->GetActiveHead()->ligand_ == Ligand::NONE) {
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
        return exe_bind_ATP(dynamic_cast<CatalyticHead *>(base), &motors_);
      });
  // Hydrolyze_ATP
  if (motors_.active_) {
    auto exe_hydrolyze = [](auto *head, auto *pop) {
      bool executed{head->parent_->Hydrolyze(head)};
      if (executed) {
        pop->FlagForUpdate();
      }
      return executed;
    };
    auto is_ATP_i_bound = [](auto *motor) -> Vec<Object *> {
      if (motor->n_heads_active_ == 1) {
        if (motor->GetActiveHead()->ligand_ == Ligand::ATP) {
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
          return exe_hydrolyze(dynamic_cast<CatalyticHead *>(base), &motors_);
        });
  }
  // Diffusion
  auto exe_diff = [](auto *head, auto *pop, auto *fil, int dir) {
    bool executed{head->Diffuse(dir)};
    if (executed) {
      pop->FlagForUpdate();
      fil->FlagForUpdate();
    }
    return executed;
  };
  auto is_singly_bound = [](auto *protein) -> Vec<Object *> {
    if (protein->n_heads_active_ == 1) {
      // only head_two can diffuse
      if (protein->GetActiveHead() == &protein->head_two_) {
        return {&protein->head_two_};
      }
    }
    return {};
  };
  Vec<int> i_min{0, 0, 0};
  Vec<size_t> dim_size{1, 1, _n_neighbs_max + 1};
  auto get_n_neighbs = [](auto *entry) {
    Vec<int> indices_vec{entry->GetNumNeighborsOccupied()};
    return indices_vec;
  };
  motors_.AddPop(
      "bound_i",
      [&](Object *base) {
        return is_singly_bound(dynamic_cast<Motor *>(base));
      },
      dim_size, i_min, get_n_neighbs);
  for (int n_neighbs{0}; n_neighbs < _n_neighbs_max; n_neighbs++) {
    kmc_.events_.emplace_back(
        "diffuse_i_fwd", xlinks_.p_event_.at("diffuse_i_fwd").GetVal(n_neighbs),
        &motors_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
        &motors_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs], binomial,
        [&](Object *base) {
          return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_,
                          filaments_, 1);
        });
    kmc_.events_.emplace_back(
        "diffuse_i_bck", xlinks_.p_event_.at("diffuse_i_bck").GetVal(n_neighbs),
        &motors_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
        &motors_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs], binomial,
        [&](Object *base) {
          return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_,
                          filaments_, -1);
        });
  }
}

void ProteinTester::InitializeTest_Motor_LatticeStep() {
  using namespace Params;
  Sys::OverrideParam("t_run", &t_run, 100.0);
  Sys::OverrideParam("t_equil", &t_equil, 0.0);
  Sys::OverrideParam("dynamic_equil_window", &dynamic_equil_window, -1);
  Sys::OverrideParam("motors: c_bulk", &Motors::c_bulk, 1.0);
  Sys::OverrideParam("motors: t_active", &Motors::t_active, 0.0);
  Sys::OverrideParam("motors: n_runs_to_exit", &Motors::n_runs_to_exit,
                     1000000);
  Sys::OverrideParam("filaments: COUNT", &Filaments::count, 1);
  Sys::OverrideParam("filaments: N_SITES[0]", &Filaments::n_sites[0], 100000);
  // ! FIXME update if new immobile_until syntax is permanently adopted
  // Sys::Log("All filament movement has been disabled by default.\n");
  // Filaments::translation_enabled[0] = false;
  // Filaments::translation_enabled[1] = false;
  // Filaments::rotation_enabled = false;
  printf("Enter test delta (-1 to check against self-coop): ");
  Str response;
  std::getline(std::cin, response);
  int test_delta{(int)std::stoi(response)};
  // Initialize sim objects
  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  // Initialize filaments
  filaments_->Initialize(this);

  double test_weight_bind{1.0};
  double test_weight_unbind{1.0};
  if (test_delta > 0) {
    test_weight_bind = Sys::weight_lattice_bind_[test_delta];
    test_weight_unbind = Sys::weight_lattice_unbind_[test_delta];
  } else {
    test_delta = 10;
  }
  // Initialize statistic trackers
  Vec<Pair<size_t, size_t>> zeros(1, {0, 0});
  Vec<double> p_theory_bind_ii(1, motors_.p_event_.at("bind_ii").GetVal() *
                                      test_weight_bind);
  test_stats_.emplace("bind_ii", zeros);
  test_ref_.emplace("bind_ii", p_theory_bind_ii);
  Vec<double> p_theory_unbind_ii(1, motors_.p_event_.at("unbind_ii").GetVal() *
                                        Square(test_weight_unbind));
  test_stats_.emplace("unbind_ii", zeros);
  test_ref_.emplace("unbind_ii", p_theory_unbind_ii);
  Vec<double> p_theory_unbind_i(1, motors_.p_event_.at("unbind_i").GetVal() *
                                       test_weight_unbind);
  test_stats_.emplace("unbind_i", zeros);
  test_ref_.emplace("unbind_i", p_theory_unbind_i);

  // Place motor head on minus end of microtubule
  BindingSite *site{filaments_->protofilaments_[0].minus_end_};
  Motor *motor{motors_.GetFreeEntry()};
  bool executed{motor->Bind(site, &motor->head_one_)};
  if (executed) {
    motors_.AddToActive(motor);
    filaments_->FlagForUpdate();
  } else {
    Sys::ErrorExit("ProteinTester::InitializeTestEnvironment()");
  }
  // Place 2nd motor test_delta distance away
  int i_site{int(site->index_ + test_delta * site->filament_->dx_)};
  BindingSite *partner_site{&filaments_->protofilaments_[0].sites_[i_site]};
  Motor *partner_motor{motors_.GetFreeEntry()};
  bool pexecuted{partner_motor->Bind(partner_site, &partner_motor->head_one_)};
  if (pexecuted) {
    motors_.AddToActive(partner_motor);
    filaments_->FlagForUpdate();
  } else {
    Sys::ErrorExit("PARTNER");
  }
  motor->head_one_.test_partner_ = &partner_motor->head_one_;
  motor->head_two_.test_partner_ = &partner_motor->head_two_;
  partner_motor->head_one_.test_partner_ = &motor->head_one_;
  partner_motor->head_two_.test_partner_ = &motor->head_two_;

  auto binomial = [&](double p, int n) {
    if (n > 0) {
      int n_exe{SysRNG::SampleBinomial(p, n)};
      if (n_exe > 1) {
        n_exe = 1;
      }
      return n_exe;
    } else {
      return 0;
    }
  };
  // Bind_ATP_I
  auto exe_bind_ATP = [](auto *head, auto *pop) {
    auto partner{head->test_partner_};
    bool executed{head->parent_->Bind_ATP(head)};
    bool pexecuted{partner->parent_->Bind_ATP(partner)};
    if (executed and pexecuted) {
      pop->FlagForUpdate();
    } else {
      Sys::ErrorExit("EXE_Bind_ATP_I()");
    }
    return executed;
  };
  auto is_NULL_i_bound = [](auto *motor) -> Vec<Object *> {
    if (motor->n_heads_active_ == 1) {
      if (motor->GetActiveHead()->ligand_ == Ligand::NONE) {
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
        return exe_bind_ATP(dynamic_cast<CatalyticHead *>(base), &motors_);
      });
  // Bind_ATP_II
  auto poisson_ATP = [&](double p, int n) {
    if (p > 0.0) {
      int n_exe{SysRNG::SamplePoisson(p)};
      if (n_exe > 1) {
        n_exe = 1;
      }
      return n_exe;
    } else {
      return 0;
    }
  };
  auto exe_bind_ATP_ii = [](auto *front_head, auto *pop, auto *fil) {
    // return;
    auto *rear_head{front_head->GetOtherHead()};
    if (front_head->trailing_) {
      return false;
    }
    auto partner{front_head->test_partner_};
    if (front_head->trailing_ != partner->trailing_) {
      partner = partner->GetOtherHead();
    }

    bool unbound{rear_head->Unbind()};
    bool executed{front_head->parent_->Bind_ATP(front_head)};

    bool punbound{partner->GetOtherHead()->Unbind()};
    bool pexecuted{partner->parent_->Bind_ATP(partner)};
    if (executed and pexecuted) {
      pop->FlagForUpdate();
      fil->FlagForUpdate();
    } else {
      Sys::ErrorExit("EXE_Bind_ATP_II()");
    }
    return executed;
  };
  auto weight_bind_ATP_ii = [](auto *head) {
    return head->parent_->GetWeight_BindATP_II(head);
  };
  auto is_NULL_ii_bound = [](auto *motor) -> Vec<Object *> {
    if (motor->n_heads_active_ == 2) {
      bool found_head{false};
      CatalyticHead *chosen_head{nullptr};
      if (motor->head_one_.ligand_ == Ligand::NONE) {
        chosen_head = &motor->head_one_;
      }
      if (motor->head_two_.ligand_ == Ligand::NONE) {
        if (chosen_head != nullptr) {
          Sys::ErrorExit("Protein_MGR::is_NULL_ii_bound()");
        }
        chosen_head = &motor->head_two_;
      }
      if (chosen_head != nullptr) {
        return {chosen_head};
      }
    }
    return {};
  };
  motors_.AddPop("bound_ii_NULL", [&](Object *base) {
    return is_NULL_ii_bound(dynamic_cast<Motor *>(base));
  });
  kmc_.events_.emplace_back(
      "bind_ATP_ii", motors_.p_event_.at("bind_ATP_ii").GetVal(),
      &motors_.sorted_.at("bound_ii_NULL").size_,
      &motors_.sorted_.at("bound_ii_NULL").entries_, poisson_ATP,
      [&](Object *base) {
        return weight_bind_ATP_ii(dynamic_cast<CatalyticHead *>(base));
      },
      [&](Object *base) {
        return exe_bind_ATP_ii(dynamic_cast<CatalyticHead *>(base), &motors_,
                               filaments_);
      });
  // Hydrolyze
  auto exe_hydrolyze = [](auto *head, auto *pop) {
    auto partner{head->test_partner_};
    bool executed{head->parent_->Hydrolyze(head)};
    bool pexecuted{partner->parent_->Hydrolyze(partner)};
    if (executed and pexecuted) {
      pop->FlagForUpdate();
    } else {
      Sys::ErrorExit("EXE_Hydrolyze");
    }
    return executed;
  };
  auto is_ATP_i_bound = [](auto *motor) -> Vec<Object *> {
    if (motor->n_heads_active_ == 1) {
      if (motor->GetActiveHead()->ligand_ == Ligand::ATP) {
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
      &motors_.sorted_.at("bound_i_ATP").entries_, binomial, [&](Object *base) {
        return exe_hydrolyze(dynamic_cast<CatalyticHead *>(base), &motors_);
      });
  // Bind_II
  auto exe_bind_ii = [&](Object *base) {
    auto bound_head{dynamic_cast<CatalyticHead *>(base)};
    auto head{bound_head->GetOtherHead()};
    auto site{head->parent_->GetNeighbor_Bind_II()};
    // ! FIXME re-introduce periodic boundary conditions
    /*
    // If dock site is plus end, unbind motor and place it on minus end
    if (site == site->filament_->plus_end_) {
      bool exe1{bound_head->Unbind()};
      auto new_site{site->filament_->minus_end_};
      bool exe2{bound_head->parent_->Bind(new_site, bound_head)};
      bool exe3{bound_head->parent_->Bind_ATP(bound_head)};
      bool exe4{bound_head->parent_->Hydrolyze(bound_head)};
      site = bound_head->parent_->GetNeighbor_Bind_II();
    }
    */
    auto partner{head->test_partner_};
    auto partner_site{partner->parent_->GetNeighbor_Bind_II()};
    auto executed{head->parent_->Bind(site, head)};
    bool pexecuted{partner->parent_->Bind(partner_site, partner)};
    if (executed and pexecuted) {
      bool still_attached{head->parent_->UpdateExtension()};
      bool wut{partner->parent_->UpdateExtension()};
      motors_.FlagForUpdate();
      filaments_->FlagForUpdate();
      test_stats_.at("bind_ii")[0].first++;
    } else {
      Sys::ErrorExit("EXE_Bind_II()");
    }
    return executed;
  };
  auto weight_bind_ii = [](auto *head) {
    return head->parent_->GetWeight_Bind_II();
  };
  auto poisson_bind_ii = [&](double p, int n) {
    test_stats_.at("bind_ii")[0].second += motors_.sorted_.at("bind_ii").size_;
    if (p > 0.0) {
      int n_exe{SysRNG::SamplePoisson(p)};
      if (n_exe > 1) {
        n_exe = 1;
      }
      return n_exe;
    } else {
      return 0;
    }
  };
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
      &motors_.sorted_.at("bind_ii").entries_, poisson_bind_ii,
      [&](Object *base) {
        return weight_bind_ii(dynamic_cast<CatalyticHead *>(base));
      },
      exe_bind_ii);
  // Unbind_II
  auto exe_unbind_ii = [&](Object *base) {
    auto head{dynamic_cast<CatalyticHead *>(base)};
    auto partner{dynamic_cast<CatalyticHead *>(head->test_partner_)};
    bool executed{head->Unbind()};
    bool pexecuted{partner->Unbind()};
    if (executed and pexecuted) {
      motors_.FlagForUpdate();
      filaments_->FlagForUpdate();
      test_stats_.at("unbind_ii")[0].first++;
    } else {
      Sys::ErrorExit("EXE_Unbind_II()");
    }
    return executed;
  };
  auto poisson_unbind_ii = [&](double p, int n) {
    test_stats_.at("unbind_ii")[0].second +=
        motors_.sorted_.at("unbind_ii").size_;
    if (p > 0.0) {
      int n_exe{SysRNG::SamplePoisson(p)};
      if (n_exe > 1) {
        n_exe = 1;
      }
      return n_exe;
    } else {
      return 0;
    }
  };
  auto weight_unbind_ii = [](auto *head) {
    return head->GetWeight_Unbind_II();
  };
  auto is_ADPP_ii_bound = [](auto *motor) -> Vec<Object *> {
    if (motor->n_heads_active_ == 2) {
      bool found_head{false};
      CatalyticHead *chosen_head{nullptr};
      if (motor->head_one_.ligand_ == Ligand::ADPP) {
        chosen_head = &motor->head_one_;
      }
      if (motor->head_two_.ligand_ == Ligand::ADPP) {
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
      &motors_.sorted_.at("unbind_ii").entries_, poisson_unbind_ii,
      [&](Object *base) {
        return weight_unbind_ii(dynamic_cast<CatalyticHead *>(base));
      },
      exe_unbind_ii);
  // Unbind_I
  auto exe_unbind_i = [&](Object *base) {
    // Count stats for unbind_i but do not actually execute it
    test_stats_.at("unbind_i")[0].first++;
    filaments_->FlagForUpdate();
    return true;
  };
  auto poisson_unbind_i = [&](double p, int n) {
    test_stats_.at("unbind_i")[0].second +=
        motors_.sorted_.at("bound_i_ADPP").size_;
    if (p > 0.0) {
      int n_exe{SysRNG::SamplePoisson(p)};
      if (n_exe > 1) {
        n_exe = 1;
      }
      return n_exe;
    } else {
      return 0;
    }
  };
  auto weight_unbind_i = [](auto *head) {
    return head->parent_->GetWeight_Unbind_I();
  };
  auto is_ADPP_i_bound = [](auto *motor) -> Vec<Object *> {
    if (motor->n_heads_active_ == 1) {
      if (motor->GetActiveHead()->ligand_ == Ligand::ADPP) {
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
      &motors_.sorted_.at("bound_i_ADPP").entries_, poisson_unbind_i,
      [&](Object *base) {
        return weight_unbind_i(dynamic_cast<CatalyticHead *>(base));
      },
      exe_unbind_i);
}

void ProteinTester::InitializeTest_Motor_LatticeBind() {
  using namespace Params;
  size_t cutoff{Motors::gaussian_range};
  // Set parameters
  Sys::OverrideParam("t_run", &t_run, 10.0);
  Sys::OverrideParam("t_equil", &t_equil, 0.0);
  Sys::OverrideParam("dynamic_equil_window", &dynamic_equil_window, -1);
  Sys::OverrideParam("motors: k_on", &Motors::k_on, 1.0);
  Sys::OverrideParam("motors: c_bulk", &Motors::c_bulk, 10.0);
  Sys::OverrideParam("motors: neighb_neighb_energy",
                     &Motors::neighb_neighb_energy, 0.0);
  Sys::OverrideParam("xlinks: c_bulk", &Xlinks::c_bulk, 0.0);
  Sys::OverrideParam("xlinks: neighb_neighb_energy",
                     &Xlinks::neighb_neighb_energy, 0.0);
  Sys::OverrideParam("filaments: COUNT", &Filaments::count, 1);
  Sys::OverrideParam("filaments: N_SITES[0]", &Filaments::n_sites[0],
                     2 * cutoff + 1);
  // ! FIXME update if new immobile_until syntax is permanently adopted
  // Sys::Log("All filament movement has been disabled by default.\n");
  // Filaments::translation_enabled[0] = false;
  // Filaments::translation_enabled[1] = false;
  // Filaments::rotation_enabled = false;
  // Initialize sim objects
  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  filaments_->Initialize(this);
  // Initialize statistic trackers
  Vec<double> p_theory(cutoff + 1, motors_.p_event_.at("bind_i").GetVal());
  for (int delta{0}; delta <= Motors::gaussian_range; delta++) {
    p_theory[delta] *= Sys::weight_lattice_bind_[delta];
  }
  test_ref_.emplace("bind", p_theory);
  Vec<Pair<size_t, size_t>> zeros(cutoff + 1, {0, 0});
  test_stats_.emplace("bind", zeros);
  // Bind motor head to middle site
  int i_site{(int)Motors::gaussian_range};
  BindingSite *site{&filaments_->protofilaments_[0].sites_[i_site]};
  Motor *motor{motors_.GetFreeEntry()};
  bool executed{motor->Bind(site, &motor->head_one_)};
  if (executed) {
    motors_.AddToActive(motor);
    filaments_->FlagForUpdate();
  } else {
    Sys::ErrorExit("ProteinTester::InitializeTestEnvironment()");
  }
  auto poisson = [&](double p, int n) {
    // Each delta distance has 2 sites available to it each timestep
    // Do not count delta = 0, where 'main' motor is permanently bound to
    for (int delta{1}; delta <= Motors::gaussian_range; delta++) {
      test_stats_.at("bind")[delta].second += 2;
    }
    if (p > 0.0) {
      return SysRNG::SamplePoisson(p);
    } else {
      return 0;
    }
  };
  auto exe_bind_i = [&](Object *base) {
    BindingSite *site{dynamic_cast<BindingSite *>(base)};
    int i_site{(int)site->index_};
    // 'main' kinesin motor will always be at index_ = lattice_cutoff_
    int delta{abs(i_site - (int)Motors::gaussian_range)};
    test_stats_.at("bind")[delta].first++;
    filaments_->FlagForUpdate();
    return true;
  };
  auto weight_bind_i = [](auto *site) { return site->GetWeight_Bind(); };
  auto is_unocc = [](Object *site) -> Vec<Object *> {
    if (!site->IsOccupied()) {
      return {site};
    }
    return {};
  };
  filaments_->AddPop("motors", is_unocc);
  kmc_.events_.emplace_back(
      "bind_i", motors_.p_event_.at("bind_i").GetVal(),
      &filaments_->unoccupied_.at("motors").size_,
      &filaments_->unoccupied_.at("motors").entries_, poisson,
      [&](Object *base) {
        return weight_bind_i(dynamic_cast<BindingSite *>(base));
      },
      exe_bind_i);
}

void ProteinTester::InitializeTest_Xlink_Diffusion() {
  using namespace Params;
  // Initialize sim objects
  Sys::OverrideParam("t_run", &t_run, 100.0);
  Sys::OverrideParam("t_equil", &t_equil, 0.0);
  Sys::OverrideParam("dynamic_equil_window", &dynamic_equil_window, -1);
  Sys::OverrideParam("xlinks: c_bulk", &Xlinks::c_bulk, 1.0);
  Sys::OverrideParam("xlinks: t_active", &Xlinks::t_active, 1.0);
  Sys::OverrideParam("filaments: COUNT", &Filaments::count, 2);
  Sys::OverrideParam("filaments: N_SITES[0]", &Filaments::n_sites[0], 1000);
  Sys::OverrideParam("filaments: N_SITES[1]", &Filaments::n_sites[1], 1000);
  // ! FIXME update if new immobile_until syntax is permanently adopted
  // Sys::Log("All filament movement has been disabled by default.\n");
  // Filaments::translation_enabled[0] = false;
  // Filaments::translation_enabled[1] = false;
  // Filaments::rotation_enabled = false;
  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  // Initialize filaments
  filaments_->Initialize(this);
  // Initialize stat trackers
  double r_y{std::fabs(Filaments::y_initial[0] - Filaments::y_initial[1])};
  // Recall, Weight = exp(0.5 * E / kbT) [assume lambda = 0.5]
  double E_max{std::log(_max_weight) * Params::kbT};
  double r_rest{Params::Xlinks::r_0};
  // E = 0.5 * k * (r - r0)^2
  double r_min{r_rest - sqrt(2 * E_max / Params::Xlinks::k_spring)};
  double r_max{r_rest + sqrt(2 * E_max / Params::Xlinks::k_spring)};
  double r_x_max{sqrt(Square(r_max) - Square(r_y))};
  size_t x_max((size_t)std::ceil(r_x_max / Filaments::site_size));
  // printf("r_max = %g\n", r_max);
  // printf("r_x_max = %g\n", r_x_max);
  // printf("x_max = %zu\n", x_max);
  Vec<double> p_theory_to(x_max + 1,
                          xlinks_.p_event_.at("diffuse_ii_to_rest").GetVal());
  Vec<double> p_theory_fr(x_max + 1,
                          xlinks_.p_event_.at("diffuse_ii_fr_rest").GetVal());
  Vec<double> wt_bind(x_max + 1, 0.0);
  Vec<double> wt_unbind(x_max + 1, 0.0);
  for (int x{0}; x <= x_max; x++) {
    double r_x{x * Filaments::site_size};
    double r{sqrt(Square(r_x) + Square(r_y))};
    if (r < r_min or r > r_max) {
      wt_bind[x] = 0.0;
      wt_unbind[x] = 0.0;
      continue;
    }
    double dr{r - Params::Xlinks::r_0};
    double dE{0.5 * Params::Xlinks::k_spring * Square(dr)};
    wt_bind[x] = exp(-(1.0 - _lambda_spring) * dE / Params::kbT);
    wt_unbind[x] = exp(_lambda_spring * dE / Params::kbT);
  }
  for (int x{0}; x <= x_max; x++) {
    // Ensure index isn't negative to avoid seg-faults
    int x_to{x > 0 ? x - 1 : 0};
    // Diffusing towards rest is considered an unbinding-type event in
    // regards to Boltzmann factors, since both events let the spring relax
    // (dividing Boltzmann factors yields E(x) - E(x_to) in the exponential)
    double wt_spring_to{wt_unbind[x] / wt_unbind[x_to]};
    // Ensure index isn't out of range to avoid seg-faults
    int x_fr{x < x_max ? x + 1 : 0};
    // Diffusing away from rest is considered a binding-type event in
    // regards to Boltzmann factors, since both events stretch the spring
    // (dividing Boltzmann factors yields E(x_fr) - E(x) in the exponential)
    double wt_spring_fr{wt_bind[x_fr] / wt_bind[x]};
    p_theory_to[x] *= wt_spring_to;
    p_theory_fr[x] *= wt_spring_fr;
  }
  test_ref_.emplace("to_rest", p_theory_to);
  test_ref_.emplace("fr_rest", p_theory_fr);
  Vec<Pair<size_t, size_t>> zeros(x_max + 1, {0, 0});
  test_stats_.emplace("to_rest", zeros);
  test_stats_.emplace("fr_rest", zeros);
  /*
  for (auto &&entry : test_stats_) {
    printf("For '%s':\n", entry.first.c_str());
    for (int x{0}; x < entry.second.size(); x++) {
      printf("  [%i] = {%zu, %zu}\n", x, entry.second[x].first,
             entry.second[x].second);
    }
  }
  */
  Protein *xlink{xlinks_.GetFreeEntry()};
  int i_site{(int)std::round(Filaments::n_sites[0] / 2)};
  BindingSite *site_one{&filaments_->protofilaments_[0].sites_[i_site]};
  BindingSite *site_two{&filaments_->protofilaments_[1].sites_[i_site]};
  bool exe_one{xlink->Bind(site_one, &xlink->head_one_)};
  bool exe_two{xlink->Bind(site_two, &xlink->head_two_)};
  if (exe_one and exe_two) {
    bool still_attached{xlink->UpdateExtension()};
    if (still_attached) {
      xlinks_.AddToActive(xlink);
      filaments_->FlagForUpdate();
    } else {
      Sys::ErrorExit("ProteinTester::InitializeTestEnvironment() [2]");
    }
  } else {
    Sys::ErrorExit("ProteinTester::InitializeTestEnvironment() [1]");
  }
  auto is_doubly_bound = [](Object *protein) -> Vec<Object *> {
    if (protein->GetNumHeadsActive() == 2) {
      return {protein->GetHeadOne(), protein->GetHeadTwo()};
    }
    return {};
  };
  xlinks_.AddPop("diffuse_ii_to_rest", is_doubly_bound);
  xlinks_.AddPop("diffuse_ii_fr_rest", is_doubly_bound);
  auto exe_diff_to = [&](Object *base) {
    BindingHead *head{dynamic_cast<BindingHead *>(base)};
    double r_x{head->pos_[0] - head->GetOtherHead()->pos_[0]};
    size_t x{(size_t)std::abs(std::round(r_x / Params::Filaments::site_size))};
    bool executed{head->Diffuse(1)};
    if (executed) {
      bool still_attached{head->parent_->UpdateExtension()};
      xlinks_.FlagForUpdate();
      filaments_->FlagForUpdate();
      test_stats_.at("to_rest")[x].first++;
    }
    return executed;
  };
  auto exe_diff_fr = [&](Object *base) {
    BindingHead *head{dynamic_cast<BindingHead *>(base)};
    double r_x{head->pos_[0] - head->GetOtherHead()->pos_[0]};
    size_t x{(size_t)std::abs(std::round(r_x / Params::Filaments::site_size))};
    bool executed{head->Diffuse(-1)};
    if (executed) {
      bool still_attached{head->parent_->UpdateExtension()};
      xlinks_.FlagForUpdate();
      filaments_->FlagForUpdate();
      test_stats_.at("fr_rest")[x].first++;
    }
    return executed;
  };
  auto weight_diff_ii = [](auto *head, int dir) {
    return head->GetWeight_Diffuse(dir);
  };
  auto poisson_to = [&](double p, int n) {
    Protein *xlink{xlinks_.active_entries_[0]};
    double r_x{xlink->head_one_.pos_[0] - xlink->head_two_.pos_[0]};
    size_t x{(size_t)std::abs(std::round(r_x / Params::Filaments::site_size))};
    if (x != 0) {
      test_stats_.at("to_rest")[x].second += 2;
    }
    if (p > 0.0) {
      return SysRNG::SamplePoisson(p);
    } else {
      return 0;
    }
  };
  auto poisson_fr = [&](double p, int n) {
    Protein *xlink{xlinks_.active_entries_[0]};
    double r_x{xlink->head_one_.pos_[0] - xlink->head_two_.pos_[0]};
    size_t x{(size_t)std::abs(std::round(r_x / Params::Filaments::site_size))};
    if (x != test_stats_.at("fr_rest").size() - 1) {
      BindingSite *site_one{xlink->head_one_.site_};
      BindingSite *site_two{xlink->head_two_.site_};
      if (site_one != site_one->filament_->plus_end_ and
          site_one != site_one->filament_->minus_end_) {
        test_stats_.at("fr_rest")[x].second += 1;
      }
      if (site_two != site_two->filament_->plus_end_ and
          site_two != site_two->filament_->minus_end_) {
        test_stats_.at("fr_rest")[x].second += 1;
      }
    }
    if (p > 0.0) {
      return SysRNG::SamplePoisson(p);
    } else {
      return 0;
    }
  };
  kmc_.events_.emplace_back(
      "diffuse_ii_to_rest", xlinks_.p_event_.at("diffuse_ii_to_rest").GetVal(),
      &xlinks_.sorted_.at("diffuse_ii_to_rest").size_,
      &xlinks_.sorted_.at("diffuse_ii_to_rest").entries_, poisson_to,
      [&](Object *base) {
        return weight_diff_ii(dynamic_cast<BindingHead *>(base), 1);
      },
      exe_diff_to);
  kmc_.events_.emplace_back(
      "diffuse_ii_fr_rest", xlinks_.p_event_.at("diffuse_ii_fr_rest").GetVal(),
      &xlinks_.sorted_.at("diffuse_ii_fr_rest").size_,
      &xlinks_.sorted_.at("diffuse_ii_fr_rest").entries_, poisson_fr,
      [&](Object *base) {
        return weight_diff_ii(dynamic_cast<BindingHead *>(base), -1);
      },
      exe_diff_fr);
}

void ProteinTester::InitializeTest_Xlink_Bind_II() {
  using namespace Params;
  Sys::OverrideParam("t_run", &t_run, 1000.0);
  Sys::OverrideParam("t_equil", &t_equil, 0.0);
  Sys::OverrideParam("dynamic_equil_window", &dynamic_equil_window, -1);
  Sys::OverrideParam("xlinks: k_off_ii", &Xlinks::k_off_ii, 143);
  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  double r_y{std::fabs(Filaments::y_initial[0] - Filaments::y_initial[1])};
  // Recall, Weight = exp(0.5 * E / kbT) [assume lambda = 0.5]
  double E_max{std::log(_max_weight) * Params::kbT};
  double r_rest{Params::Xlinks::r_0};
  // E = 0.5 * k * (r - r0)^2
  double r_min{r_rest - sqrt(2 * E_max / Params::Xlinks::k_spring)};
  double r_max{r_rest + sqrt(2 * E_max / Params::Xlinks::k_spring)};
  double r_x_max{sqrt(Square(r_max) - Square(r_y))};
  int x_max((int)std::ceil(r_x_max / Filaments::site_size));
  // printf("r_max = %g\n", r_max);
  // printf("r_x_max = %g\n", r_x_max);
  // printf("x_max = %i\n", x_max);
  Sys::OverrideParam("filaments: COUNT", &Filaments::count, 2);
  Sys::OverrideParam("filaments: N_SITES[0]", &Filaments::n_sites[0],
                     2 * x_max + 1);
  Sys::OverrideParam("filaments: N_SITES[1]", &Filaments::n_sites[1],
                     2 * x_max + 1);
  // ! FIXME update if new immobile_until syntax is permanently adopted
  // Sys::Log("All filament movement has been disabled by default.\n");
  // Filaments::translation_enabled[0] = false;
  // Filaments::translation_enabled[1] = false;
  // Filaments::rotation_enabled = false;
  // Initialize filament environment
  filaments_->Initialize(this);
  Vec<double> p_bind(2 * x_max + 1, xlinks_.p_event_.at("bind_ii").GetVal());
  Vec<double> p_unbind(2 * x_max + 1,
                       xlinks_.p_event_.at("unbind_ii").GetVal());
  double offset(Filaments::x_initial[1] - Filaments::x_initial[0]);
  for (int x{-x_max}; x <= x_max; x++) {
    int x_index{x_max + x};
    double r_x{x * Filaments::site_size + offset};
    double r{sqrt(Square(r_x) + Square(r_y))};
    if (r < r_min or r > r_max) {
      p_bind[x_index] *= 0.0;
      p_unbind[x_index] *= 0.0;
      continue;
    }
    double dr{r - Params::Xlinks::r_0};
    double dE{0.5 * Params::Xlinks::k_spring * Square(dr)};
    p_bind[x_index] *= exp(-(1.0 - _lambda_spring) * dE / Params::kbT);
    p_unbind[x_index] *= exp(_lambda_spring * dE / Params::kbT);
  }
  test_ref_.emplace("bind_ii", p_bind);
  test_ref_.emplace("unbind_ii", p_unbind);
  Vec<Pair<size_t, size_t>> zeros(2 * x_max + 1, {0, 0});
  test_stats_.emplace("bind_ii", zeros);
  test_stats_.emplace("unbind_ii", zeros);
  // Place first xlink head on lower MT; remains static for entire sim
  int i_site{(int)x_max};
  BindingSite *site{&filaments_->protofilaments_[0].sites_[i_site]};
  Protein *xlink{xlinks_.GetFreeEntry()};
  bool executed{xlink->Bind(site, &xlink->head_one_)};
  if (executed) {
    xlinks_.AddToActive(xlink);
    filaments_->FlagForUpdate();
  } else {
    Sys::ErrorExit("ProteinTester::InitializeTestEnvironment()");
  }
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
      int n_entries{(int)xlinks_.sorted_.at("bind_ii").size_};
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
        return false;
      }
      double r_x{head->pos_[0] - bound_head->pos_[0]};
      double offset(Filaments::x_initial[1] - Filaments::x_initial[0]);
      int x{(int)std::round((r_x - offset) / Params::Filaments::site_size)};
      int x_max{int(test_stats_.at("unbind_ii").size() - 1) / 2};
      test_stats_.at("bind_ii")[x + x_max].first++;
    } else {
      Sys::ErrorExit("Bind_II (TEST)");
    }
    return executed;
  };
  // Construct KMC event fr Bind_II
  kmc_.events_.emplace_back("bind_ii", xlinks_.p_event_.at("bind_ii").GetVal(),
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
      double offset(Filaments::x_initial[1] - Filaments::x_initial[0]);
      int x{(int)std::round((r_x - offset) / Params::Filaments::site_size)};
      int x_max{int(test_stats_.at("unbind_ii").size() - 1) / 2};
      test_stats_.at("unbind_ii")[x + x_max].second += 1;
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
      double offset(Filaments::x_initial[1] - Filaments::x_initial[0]);
      int x{(int)std::round((r_x - offset) / Params::Filaments::site_size)};
      int x_max{int(test_stats_.at("unbind_ii").size() - 1) / 2};
      test_stats_.at("unbind_ii")[x + x_max].first++;
      // bool head->parent_->UpdateExtension();
    } else {
      Sys::ErrorExit("Unbind_II (TEST)");
    }
    return executed;
  };
  kmc_.events_.emplace_back(
      "unbind_ii", xlinks_.p_event_.at("unbind_ii").GetVal(),
      &xlinks_.sorted_.at("unbind_ii").size_,
      &xlinks_.sorted_.at("unbind_ii").entries_, poisson_unbind_ii,
      get_weight_unbind_ii, exe_unbind_ii);
}

void ProteinTester::InitializeTest_Shepherding() {
  using namespace Params;
  Sys::OverrideParam("t_run", &t_run, 10.0);
  Sys::OverrideParam("t_equil", &t_equil, 0.0);
  Sys::OverrideParam("dynamic_equil_window", &dynamic_equil_window, -1);
  Sys::OverrideParam("motors: c_bulk", &Motors::c_bulk, 1.0);
  Sys::OverrideParam("motors: t_active", &Motors::t_active, 0.0);
  Sys::OverrideParam("motors: n_runs_to_exit", &Motors::n_runs_to_exit, -1);
  Sys::OverrideParam("xlinks: c_bulk", &Xlinks::c_bulk, 1.0);
  Sys::OverrideParam("xlinks: t_active", &Xlinks::t_active, 0.0);
  Sys::OverrideParam("filaments: COUNT", &Filaments::count, 1);
  Sys::OverrideParam("filaments: N_SITES[0]", &Filaments::n_sites[0], 1000);
  Sys::OverrideParam("filaments: f_applied[0]", &Filaments::f_applied[0], 0.0);
  Sys::OverrideParam("filaments: f_applied[1]", &Filaments::f_applied[1], 0.0);

  // ! FIXME update if new immobile_until syntax is permanently adopted
  // Sys::Log("All filament movement has been disabled by default.\n");
  // Filaments::translation_enabled[0] = false;
  // Filaments::translation_enabled[1] = false;
  // Filaments::rotation_enabled = false;
  // Initialize sim objects
  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  // Initialize filaments
  filaments_->Initialize(this);
  // Place motor head on minus end of microtubule
  BindingSite *site{filaments_->protofilaments_[0].minus_end_};
  Motor *motor{motors_.GetFreeEntry()};
  bool executed{motor->Bind(site, &motor->head_one_)};
  if (executed) {
    motors_.AddToActive(motor);
    filaments_->FlagForUpdate();
  } else {
    Sys::ErrorExit("ProteinTester::InitializeTestEnvironment() [0]");
  }
  // Place crosslinker head on neighboring site
  BindingSite *neighb{site->GetNeighbor(site->filament_->dx_)};
  Protein *xlink{xlinks_.GetFreeEntry()};
  bool pexecuted{xlink->Bind(neighb, &xlink->head_one_)};
  if (pexecuted) {
    xlinks_.AddToActive(xlink);
    filaments_->FlagForUpdate();
  } else {
    Sys::ErrorExit("ProteinTester::InitializeTestEnvironment() [1]");
  }
  auto binomial = [&](double p, int n) { // For events w/ constant probs.
    if (n > 0) {
      return SysRNG::SampleBinomial(p, n);
    } else {
      return 0;
    }
  };
  auto poisson = [&](double p, int n) { // For events w/ time-variable probs.
    if (p > 0.0) {
      return SysRNG::SamplePoisson(p);
    } else {
      return 0;
    }
  };
  // Bind_ATP
  auto is_bound_i_NULL = [](Object *base) -> Vec<Object *> {
    Motor *motor = dynamic_cast<Motor *>(base);
    if (motor->GetNumHeadsActive() == 1) {
      if (!motor->IsTethered() or motor->HasSatellite()) {
        if (motor->GetActiveHead()->ligand_ == Ligand::NONE) {
          return {motor->GetActiveHead()};
        }
      }
    }
    return {};
  };
  auto exe_bind_ATP = [](auto *head, auto *pop) {
    bool executed{head->parent_->Bind_ATP(head)};
    return executed;
  };
  motors_.AddPop("bound_i_NULL", is_bound_i_NULL);
  kmc_.events_.emplace_back(
      "bind_ATP_i", motors_.p_event_.at("bind_ATP_i").GetVal(),
      &motors_.sorted_.at("bound_i_NULL").size_,
      &motors_.sorted_.at("bound_i_NULL").entries_, binomial,
      [&](Object *base) {
        return exe_bind_ATP(dynamic_cast<CatalyticHead *>(base), &motors_);
      });
  // Hydrolyze_ATP
  auto is_bound_i_ATP = [](Object *base) -> Vec<Object *> {
    Motor *motor = dynamic_cast<Motor *>(base);
    if (motor->GetNumHeadsActive() == 1) {
      if (motor->GetActiveHead()->ligand_ == Ligand::ATP) {
        return {motor->GetActiveHead()};
      }
    }
    return {};
  };
  auto exe_hydrolyze = [](auto *head, auto *pop) {
    bool executed{head->parent_->Hydrolyze(head)};
    return executed;
  };
  motors_.AddPop("bound_i_ATP", is_bound_i_ATP);
  kmc_.events_.emplace_back(
      "hydrolyze", motors_.p_event_.at("hydrolyze").GetVal(),
      &motors_.sorted_.at("bound_i_ATP").size_,
      &motors_.sorted_.at("bound_i_ATP").entries_, binomial, [&](Object *base) {
        return exe_hydrolyze(dynamic_cast<CatalyticHead *>(base), &motors_);
      });
  // Bind_II
  auto is_docked = [](Object *base) -> Vec<Object *> {
    Motor *motor = dynamic_cast<Motor *>(base);
    auto *docked_head{motor->GetDockedHead()};
    if (docked_head != nullptr) {
      return {docked_head->GetOtherHead()};
    }
    return {};
  };
  auto weight_bind_ii = [](auto *head) {
    return head->parent_->GetWeight_Bind_II();
  };
  auto exe_bind_ii = [](auto *bound_head, auto *pop, auto *fil) {
    auto head{bound_head->GetOtherHead()};
    auto site{head->parent_->GetNeighbor_Bind_II()};
    auto executed{head->parent_->Bind(site, head)};
    if (executed) {
      // FIXME check if this check is necessary / what to do
      bool still_attached{head->parent_->UpdateExtension()};
    }
    return executed;
  };
  motors_.AddPop("docked", is_docked);
  kmc_.events_.emplace_back(
      "bind_ii (motors)", motors_.p_event_.at("bind_ii").GetVal(),
      &motors_.sorted_.at("docked").size_,
      &motors_.sorted_.at("docked").entries_, poisson,
      [&](Object *base) {
        return weight_bind_ii(dynamic_cast<CatalyticHead *>(base));
      },
      [&](Object *base) {
        return exe_bind_ii(dynamic_cast<CatalyticHead *>(base), &motors_,
                           filaments_);
      });
  // Unbind_II
  auto is_bound_ii_ADPP = [](Object *base) -> Vec<Object *> {
    Motor *motor = dynamic_cast<Motor *>(base);
    if (motor->GetNumHeadsActive() == 2) {
      Object *chosen_head{nullptr};
      if (motor->GetHeadOne()->ligand_ == Ligand::ADPP) {
        chosen_head = motor->GetHeadOne();
      }
      if (motor->GetHeadTwo()->ligand_ == Ligand::ADPP) {
        if (chosen_head != nullptr) {
          Sys::ErrorExit("Protein_MNGR::is_bound_ii_ADPP()");
        }
        chosen_head = motor->GetHeadTwo();
      }
      if (chosen_head != nullptr) {
        return {chosen_head};
      }
    }
    return {};
  };
  auto weight_unbind_ii = [](auto *head) {
    return head->GetWeight_Unbind_II();
  };
  auto exe_unbind_ii = [](auto *head, auto *pop, auto *fil) {
    bool executed{head->Unbind()};
    return executed;
  };
  motors_.AddPop("doubly_bound_ADPP", is_bound_ii_ADPP);
  kmc_.events_.emplace_back(
      "unbind_ii (motors)", motors_.p_event_.at("unbind_ii").GetVal(),
      &motors_.sorted_.at("doubly_bound_ADPP").size_,
      &motors_.sorted_.at("doubly_bound_ADPP").entries_, poisson,
      [&](Object *base) {
        return weight_unbind_ii(dynamic_cast<CatalyticHead *>(base));
      },
      [&](Object *base) {
        return exe_unbind_ii(dynamic_cast<CatalyticHead *>(base), &motors_,
                             filaments_);
      });
  // Xlink diffusion
  auto is_bound_i = [](Object *base) -> Vec<Object *> {
    Protein *protein{dynamic_cast<Protein *>(base)};
    if (protein->GetNumHeadsActive() == 1) {
      if (!protein->IsTethered() or protein->HasSatellite()) {
        return {protein->GetActiveHead()};
      }
    }
    return {};
  };
  // Sorted array dimension {i,j,k}; 1-D for neighbs so we use k and pad i & j
  Vec<size_t> dim_size{1, 1, _n_neighbs_max + 1};
  // Starting indices {i, j, k} of array for neighb coop, only use k dimension
  Vec<int> i_min{0, 0, 0};
  auto get_n_neighbs = [](Object *entry) {
    Vec<int> indices_vec{entry->GetNumNeighborsOccupied()};
    return indices_vec;
  };
  xlinks_.AddPop("bound_i", is_bound_i, dim_size, i_min, get_n_neighbs);
  auto exe_diff = [](auto *head, auto *pop, auto *fil, int dir) {
    bool executed{head->Diffuse(dir)};
    if (executed) {
      bool still_attached{head->parent_->UpdateExtension()};
      // TODO do I need this check??
      if (!still_attached) {
        // printf("what\n");
      }
    }
    return executed;
  };
  for (int n_neighbs{0}; n_neighbs < _n_neighbs_max; n_neighbs++) {
    kmc_.events_.emplace_back(
        "diffuse_i_fwd", xlinks_.p_event_.at("diffuse_i_fwd").GetVal(n_neighbs),
        &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
        &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs], binomial,
        [&](Object *base) {
          return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_,
                          filaments_, 1);
        });
    kmc_.events_.emplace_back(
        "diffuse_i_bck", xlinks_.p_event_.at("diffuse_i_bck").GetVal(n_neighbs),
        &xlinks_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
        &xlinks_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs], binomial,
        [&](Object *base) {
          return exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_,
                          filaments_, -1);
        });
  }
}