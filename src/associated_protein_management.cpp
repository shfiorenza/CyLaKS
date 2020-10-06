#include "associated_protein_management.h"
#include "master_header.h"

AssociatedProteinManagement::AssociatedProteinManagement() {}

void AssociatedProteinManagement::Initialize(system_parameters *parameters,
                                             system_properties *properties) {

  wally_ = &properties->wallace;
  gsl_ = &properties->gsl;
  parameters_ = parameters;
  properties_ = properties;
  CalculateCutoffs();
  SetParameters();
  GenerateXLinks();
  InitializeLists();
  if (wally_->test_mode_ != nullptr) {
    InitializeTestEnvironment();
    InitializeTestEvents();
    return;
  }
  InitializeEvents();
}

void AssociatedProteinManagement::CalculateCutoffs() {

  double kbT{parameters_->kbT};
  double r_0{parameters_->xlinks.r_0};
  double r_y{parameters_->microtubules.y_dist};
  double k_spring{parameters_->xlinks.k_spring};
  double site_size{parameters_->microtubules.site_size};
  /* First, calculate rest_dist_ in number of sites */
  // For now, it should always come out to zero
  int approx_rest{(int)(sqrt(r_0 * r_0 - r_y * r_y) / site_size)};
  if (approx_rest == 0) {
    rest_dist_ = approx_rest;
  } else {
    wally_->ErrorExit("AP_MGMT::CalculateCutoffs()");
  }
  double force_cutoff{20}; // in pN
  double f_spring{0.0};
  /* next, calculate extension distance cutoff */
  for (int x_dist{rest_dist_}; x_dist < 1000; x_dist++) {
    int r_x = x_dist * site_size;
    double r = sqrt(r_y * r_y + r_x * r_x);
    double dr = r - r_0;
    f_spring = fabs(dr * k_spring);
    /*
    if (f_spring >= force_cutoff) {
      dist_cutoff_ = x_dist;
      break;
    }
    */
    double U = (k_spring / 2) * dr * dr;
    double boltzmann_weight = exp(U / (2 * kbT));
    if (boltzmann_weight > 1e2) {
      dist_cutoff_ = x_dist;
      break;
    }
  }
  /*
  if (dist_cutoff_ > 9) {
    wally_->Log("As of now, x_dist must be less than 10");
    wally_->Log(" (currenly %i with k_spring=%g)\n", dist_cutoff_, k_spring);
    wally_->ErrorExit("AP_MGMT::CalculateCutoffs() [2]");
  }
  */
  if (parameters_->microtubules.count > 1) {
    crosslinking_active_ = true;
    wally_->Log("\nFor crosslinkers:\n");
    wally_->Log("  rest_dist is %i\n", rest_dist_);
    wally_->Log("  dist_cutoff is %i (f = %g pN)\n\n", dist_cutoff_, f_spring);
  } else {
    wally_->Log("\nCrosslinker double-binding is inactive.\n\n");
  }
}

void AssociatedProteinManagement::SetParameters() {

  if (parameters_->xlinks.c_bulk > 0.0) {
    population_active_ = true;
  }
  if (parameters_->motors.tethers_active and parameters_->motors.c_bulk > 0.0) {
    tethering_active_ = true;
  }
  teth_cutoff_ = properties_->kinesin4.teth_cutoff_;
  comp_cutoff_ = properties_->kinesin4.teth_cutoff_;
  /*
   For events that result in a change in energy, we use Boltzmann factors to
   scale rates appropriately. Detailed balance is satisfied with the factors:
                 exp{-(1 - lambda)*(delta_E)/(kB*T)}, and
                 exp{-(lambda)*(delta_E)/(kB*T)}
   for an event and its complement (e.g., binding and unbinding), where lambda
   is a constant that ranges from 0 to 1, delta_E is the change in energy that
   results from the event, kB is Boltzmann's constant, and T is the temperature

   Three common values of lambda demonstrate its function:
       Lambda = 0.0 means all energy dependences is in binding
       Lambda = 0.5 means energy dependence is equal for binding and unbinding
       Lambda = 1.0 means all energy dependence is in unbinding
  */
  /* Boltzmann weights for internal spring mechanism */
  double lambda_spring{0.5}; // Lambda for spring energies
  double r_0{parameters_->xlinks.r_0};
  double r_y{parameters_->microtubules.y_dist};
  double k_spring{parameters_->xlinks.k_spring};
  double site_size{parameters_->microtubules.site_size};
  double kbT{parameters_->kbT};
  weight_spring_bind_.resize(dist_cutoff_ + 1);
  weight_spring_unbind_.resize(dist_cutoff_ + 1);
  for (int x{0}; x <= dist_cutoff_; x++) {
    double r_x{x * site_size};
    double r{sqrt(r_x * r_x + r_y * r_y)};
    double dr{r - r_0};
    double energy{0.5 * k_spring * dr * dr}; // ! in pN*nm
    weight_spring_bind_[x] = exp(-(1.0 - lambda_spring) * energy / kbT);
    weight_spring_unbind_[x] = exp(lambda_spring * energy / kbT);
  }
  /* Boltzmann weights for neighbor cooperativity */
  double lambda_neighb{1.0}; // Lambda for neighbor interaction energies
  if (lambda_neighb != 1.0) {
    wally_->ErrorExit("AP_MGMT::SetParameters() [lambda_neighb]");
  }
  // Array of Boltzmann weights due to neighbor-neighbor interactions
  double energy_per_neighb{-1 * parameters_->xlinks.interaction_energy}; // kbT
  weight_neighbs_bind_.resize(max_neighbs_ + 1);
  weight_neighbs_unbind_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    double energy{n_neighbs * energy_per_neighb}; // ! in kbT
    weight_neighbs_bind_[n_neighbs] = exp(-(1.0 - lambda_neighb) * energy);
    weight_neighbs_unbind_[n_neighbs] = exp(lambda_neighb * energy);
  }

  // [DIFFUSION STATISTICS FOR CROSSLINKER W/O TETH BELOW] //
  double delta_t{parameters_->delta_t};
  double x_squared{(site_size / 1000) * (site_size / 1000)}; //! in um^2
  double tau_i{x_squared / (2 * parameters_->xlinks.diffu_coeff_i)};
  p_diffuse_i_ = delta_t / tau_i;
  double tau_ii{x_squared / (2 * parameters_->xlinks.diffu_coeff_ii)};
  p_diffuse_ii_ = delta_t / tau_ii;
  // diff_fwd and diff_bck are two separate events, which effectively
  // doubles the probability to diffuse. Thus we divide p_diff by 2.
  p_diffuse_i_ /= 2.0;
  p_diffuse_ii_ /= 2.0;
  double k_on{parameters_->xlinks.k_on};
  double c_bulk{parameters_->xlinks.c_bulk};
  p_bind_i_ = k_on * c_bulk * delta_t;
  double c_eff_teth{parameters_->motors.c_eff_tether};
  p_bind_i_teth_ = k_on * c_eff_teth * delta_t;
  double c_eff_bind{parameters_->xlinks.c_eff_bind};
  p_bind_ii_ = k_on * c_eff_bind * delta_t;
  double k_off_ii{parameters_->xlinks.k_off_ii};
  p_unbind_ii_ = k_off_ii * delta_t;
  double k_off_i{parameters_->xlinks.k_off_i};
  p_unbind_i_ = k_off_i * delta_t;
  double k_tether{parameters_->motors.k_tether};
  p_tether_ = k_tether * c_bulk * delta_t;
  double k_untether = parameters_->motors.k_untether;
  p_untether_ = k_untether * delta_t;

  /*
  // Add these johnnys to our p_theory_ map; need to pad so all are 3-D
  p_theory_["diffuse_i_fwd"] = {{p_diffuse_i_fwd_}};
  p_theory_["diffuse_i_bck"] = {{p_diffuse_i_bck_}};
  if (crosslinking_active_) {
    p_theory_["diffuse_ii_to_rest"] = {p_diffuse_ii_to_rest_};
    p_theory_["diffuse_ii_fr_rest"] = {p_diffuse_ii_fr_rest_};
  }
  if (tethering_active_) {
    p_theory_["diffuse_i_to_teth_rest"] = {p_diffuse_i_to_teth_rest_};
    p_theory_["diffuse_i_fr_teth_rest"] = {p_diffuse_i_fr_teth_rest_};
    if (crosslinking_active_) {
      p_theory_["diffuse_ii_to_both"] = p_diffuse_ii_to_both_;
      p_theory_["diffuse_ii_fr_both"] = p_diffuse_ii_fr_both_;
      p_theory_["diffuse_ii_to_self_fr_teth"] = p_diffuse_ii_to_self_fr_teth_;
      p_theory_["diffuse_ii_fr_self_to_teth"] = p_diffuse_ii_fr_self_to_teth_;
    }
  }
  p_theory_["bind_i"] = {{p_bind_i_}};
  p_theory_["unbind_i"] = {{p_unbind_i_}};
  if (crosslinking_active_) {
    p_theory_["bind_ii"] = {{{p_avg_bind_ii_}}};
    p_theory_["unbind_ii"] = {{p_unbind_ii_}};
  }
  if (tethering_active_) {
    p_theory_["bind_i_teth"] = {{{p_bind_i_teth_}}};
    p_theory_["unbind_i_teth"] = {p_unbind_i_teth_};
    if (crosslinking_active_) {
      p_theory_["bind_ii_teth"] = {{{p_avg_bind_ii_}}};
      p_theory_["unbind_ii_to_teth"] = p_unbind_ii_to_teth_;
      p_theory_["unbind_ii_fr_teth"] = p_unbind_ii_fr_teth_;
    }
    p_theory_["tether_free"] = {{{p_tether_free_}}};
    p_theory_["untether_free"] = {{{p_untether_free_}}};
  }
  // Manually construct probabilities used in Poisson processes using weights
  Vec<Vec<Vec<double>>> p_bind_ii = {weight_bind_ii_};
  Vec<Vec<Vec<double>>> p_bind_i_teth = {weight_bind_i_teth_};
  Vec<Vec<Vec<double>>> p_bind_ii_to_teth = weight_bind_ii_to_teth_;
  Vec<Vec<Vec<double>>> p_bind_ii_fr_teth = weight_bind_ii_fr_teth_;
  for (int n_neighbs{0}; n_neighbs <= 2; n_neighbs++) {
    for (int x{0}; x <= dist_cutoff_; x++) {
      p_bind_ii[0][n_neighbs][x] *= p_avg_bind_ii_;
    }
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      p_bind_i_teth[0][n_neighbs][x_dub] *= p_bind_i_teth_;
      for (int x{0}; x <= dist_cutoff_; x++) {
        p_bind_ii_to_teth[n_neighbs][x_dub][x] *= p_avg_bind_ii_;
        p_bind_ii_fr_teth[n_neighbs][x_dub][x] *= p_avg_bind_ii_;
      }
    }
  }
  // p_theory_["bind_ii"] = p_bind_ii;
  if (tethering_active_) {
    p_theory_["bind_i_teth"] = p_bind_i_teth;
    p_theory_["bind_ii_to_teth"] = p_bind_ii_to_teth;
    p_theory_["bind_ii_fr_teth"] = p_bind_ii_fr_teth;
  }
  // READ EM OUT
  for (const auto &entry : p_theory_) {
    auto label = entry.first;
    auto value = entry.second;
    for (int i{0}; i < value.size(); i++) {
      for (int j{0}; j < value[i].size(); j++) {
        for (int k{0}; k < value[i][j].size(); k++) {
          std::string name{"p_" + label + "_"};
          if (value.size() > 1) {
            name += "[" + std::to_string(i) + "]";
          }
          if (value[i].size() > 1) {
            name += "[" + std::to_string(j) + "]";
          }
          if (value[i][j].size() > 1) {
            name += "[" + std::to_string(k) + "]";
          }
          if (verbosity_ >= 3) {
            wally_->Log("%s = %g\n", name.c_str(), value[i][j][k]);
          }
          if (value[i][j][k] > 1.0) {
            wally_->Log("Error! %s = %g\n", name.c_str(), value[i][j][k]);
            // wally_->ErrorExit("AP_MGMT:SetParameters()");
          }
        }
      }
    }
  }
  */
}

void AssociatedProteinManagement::GenerateXLinks() {

  // Since only one head has to be bound, the sim will at most
  // as many xlinks as sites in the bulk (all single-bound)
  int n_mts = parameters_->microtubules.count;
  for (int i_mt{0}; i_mt < n_mts; i_mt++) {
    n_xlinks_ += parameters_->microtubules.length[i_mt];
  }
  xlinks_.resize(n_xlinks_);
  for (int id{0}; id < n_xlinks_; id++) {
    xlinks_[id].Initialize(parameters_, properties_, id);
  }
}

void AssociatedProteinManagement::InitializeLists() {

  int teth_cutoff{properties_->kinesin4.teth_cutoff_};
  // Population size trackers (not a list ok bite me)
  n_bound_i_.resize(max_neighbs_ + 1);
  n_bound_i_teth_.resize(max_neighbs_ + 1);
  n_bound_ii_.resize(max_neighbs_ + 1);
  n_bound_ii_teth_same_.resize(max_neighbs_ + 1);
  n_bound_ii_teth_oppo_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    n_bound_i_[n_neighbs] = 0;
    n_bound_ii_[n_neighbs].resize(dist_cutoff_ + 1);
    for (int x_dist{0}; x_dist <= dist_cutoff_; x_dist++) {
      n_bound_ii_[n_neighbs][x_dist] = 0;
    }
    n_bound_i_teth_[n_neighbs].resize(2 * teth_cutoff + 1);
    n_bound_ii_teth_same_[n_neighbs].resize(2 * teth_cutoff + 1);
    n_bound_ii_teth_oppo_[n_neighbs].resize(2 * teth_cutoff + 1);
    for (int x_dub{0}; x_dub <= 2 * teth_cutoff; x_dub++) {
      n_bound_i_teth_[n_neighbs][x_dub] = 0;
      n_bound_ii_teth_same_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      n_bound_ii_teth_oppo_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      for (int x{0}; x <= dist_cutoff_; x++) {
        n_bound_ii_teth_same_[n_neighbs][x_dub][x] = 0;
        n_bound_ii_teth_oppo_[n_neighbs][x_dub][x] = 0;
      }
    }
  }
  // Lists
  active_.resize(n_xlinks_);
  satellites_.resize(n_xlinks_);
  bound_unteth_.resize(n_xlinks_);
  bind_ii_candidates_.resize(n_xlinks_);
  bind_i_teth_candidates_.resize(n_xlinks_);
  bind_ii_teth_candidates_.resize(n_xlinks_);
  bound_i_.resize(max_neighbs_ + 1);
  bound_ii_.resize(max_neighbs_ + 1);
  bound_i_teth_.resize(max_neighbs_ + 1);
  bound_ii_teth_oppo_.resize(max_neighbs_ + 1);
  bound_ii_teth_same_.resize(max_neighbs_ + 1);
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    bound_i_[n_neighbs].resize(n_xlinks_);
    bound_ii_[n_neighbs].resize(dist_cutoff_ + 1);
    for (int x{0}; x <= dist_cutoff_; x++) {
      bound_ii_[n_neighbs][x].resize(n_xlinks_);
    }
    bound_i_teth_[n_neighbs].resize(2 * teth_cutoff + 1);
    bound_ii_teth_oppo_[n_neighbs].resize(2 * teth_cutoff + 1);
    bound_ii_teth_same_[n_neighbs].resize(2 * teth_cutoff + 1);
    for (int x_dub{0}; x_dub <= 2 * teth_cutoff; x_dub++) {
      bound_i_teth_[n_neighbs][x_dub].resize(n_xlinks_);
      bound_ii_teth_oppo_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      bound_ii_teth_same_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
      for (int x{0}; x <= dist_cutoff_; x++) {
        bound_ii_teth_oppo_[n_neighbs][x_dub][x].resize(n_xlinks_);
        bound_ii_teth_same_[n_neighbs][x_dub][x].resize(n_xlinks_);
      }
    }
  }
}

void AssociatedProteinManagement::InitializeTestEnvironment() {

  if (strcmp(wally_->test_mode_, "xlink_bind_ii") == 0) {
    wally_->Log("Initializing test %s\n", wally_->test_mode_);
    properties_->microtubules.InitializeTestEnvironment();
    // Stat tracker for events: 1st entry is n_executed, 2nd is n_opportunities
    Vec<std::pair<size_t, size_t>> zeros(dist_cutoff_ + 1, {0, 0});
    test_stats_.emplace("bind_ii", zeros);
    // Test reference, i.e., the theoretical probabilities we expect to observe
    Vec<double> p_theory(dist_cutoff_ + 1, p_bind_ii_);
    for (int x{0}; x <= dist_cutoff_; x++) {
      p_theory[x] *= weight_spring_bind_[x];
    }
    test_ref_.emplace("bind_ii", p_theory);
    // Bind first head of xlink to MT; static for entirety of test
    int i_site{dist_cutoff_};
    Tubulin *site{&properties_->microtubules.mt_list_[0].lattice_[i_site]};
    AssociatedProtein *xlink{GetFreeXlink()};
    site->xlink_head_ = &xlink->head_one_;
    site->occupied_ = true;
    xlink->head_one_.site_ = site;
    xlink->heads_active_++;
    AddToActive(xlink);
    properties_->microtubules.FlagForUpdate();
  }
}

void AssociatedProteinManagement::InitializeTestEvents() {

  std::string event_name; // scratch space to construct each event name
  //  Binomial probabilitiy distribution; sampled to predict most events
  auto binomial = [&](double p, int n) {
    if (n > 0) {
      return properties_->gsl.SampleBinomialDist(p, n);
    } else {
      return 0;
    }
  };
  // Poisson distribution; sampled to predict events w/ variable probabilities
  auto poisson = [&](double p, int n) {
    if (p > 0.0) {
      return gsl_->SamplePoissonDist(p);
    } else {
      return 0;
    }
  };
  // Function that returns a random probability in the range [0.0, 1.0)
  auto get_ran_prob = [&](void) { return gsl_->GetRanProb(); };
  //  Function that sets n random indices from the range [0, m)
  auto set_ran_indices = [&](int *indices, int n, int m) {
    if (m > 1) {
      properties_->gsl.SetRanIndices(indices, n, m);
    } else {
      indices[0] = 0;
    }
  };
  /*
  if (strcmp(wally_->test_mode_, "xlink_bind_ii") == 0) {
    // Bind II: binds the second crosslinker head to adjacent microtubule
    event_name = "bind_ii";
    auto exe_bind_ii = [&](ENTRY_T target) {
      POP_T *head{nullptr};
      try {
        head = std::get<POP_T *>(target);
      } catch (...) {
        wally_->ErrorExit("AP_MGMT::Bind_II()");
      }
      Bind_II(head);
      head->xlink_->UpdateExtension();
      test_stats_["bind_ii"][head->xlink_->x_dist_].first++;
    };
    auto weight_bind_ii = [&](ENTRY_T target) { return std::get < };
    // Poisson dist. is used w/ partition function for E-dependent binding
    auto poisson_ii = [&](double p, int n) {
      for (int x{0}; x <= dist_cutoff_; x++) {
        if (x == 0) {
          // Only one site (directly above) can be bound to for this config
          test_stats_["bind_ii"][x].second += n_bind_ii_candidates_;
        } else {
          // Two sites (to the left or right of above) are eligible for these
          // Include each as a separate opportunity
          test_stats_["bind_ii"][x].second += 2 * n_bind_ii_candidates_;
        }
      }
      double weight{GetWeight_Bind_II()};
      if (weight == 0.0) {
        return 0;
      }
      int n_expected{properties_->gsl.SamplePoissonDist(p * weight)};
      if (n_expected > 0) {
        if (n_expected > n) {
          wally_->Log("Rescaled Bind_II exp from %i to %i\n", n_expected, n);
          n_expected = n;
        }
        int n_removed{Set_Bind_II_Candidates(n_expected)};
        n_expected -= n_removed;
      }
      return n_expected;
    };
    events_.emplace_back(event_name, p_avg_bind_ii_, &n_bind_ii_candidates_,
                         &bind_ii_candidates_, exe_bind_ii, poisson_ii,
                         set_ran_indices);

    // Unbind II: unbinds a head of doubly-bound crosslinkers
    auto exe_unbind_ii = [&](ENTRY_T target) {
      Unbind_II(std::get<POP_T *>(target));
    };
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
      for (int x{0}; x <= dist_cutoff_; x++) {
        std::string X_DIST{"_" + std::to_string(x)};
        event_name = "unbind_ii" + N_NEIGHBS + X_DIST;
        events_.emplace_back(
            event_name, p_unbind_ii_[n_neighbs][x], &n_bound_ii_[n_neighbs][x],
            &bound_ii_[n_neighbs][x], exe_unbind_ii, binomial, set_ran_indices);
      }
    }
  }
  */
  /*
  if (strcmp(wally_->test_mode_, "xlink_bind_i_teth") == 0) {
    // Bind_I_Teth: same as above but for tethered unbound xlinks (satellites)
    event_name = "bind_i_teth";
    auto exe_bind_i_teth = [&](ENTRY_T target) {
      Bind_I_Teth(std::get<POP_T *>(target));
    };
    // Poisson dist. is used w/ partition function for E-dependent binding
    auto poisson_i_teth = [&](double p, int n) {
      double weight{GetWeight_Bind_I_Teth()};
      if (weight == 0.0) {
        return 0;
      }
      int n_expected{properties_->gsl.SamplePoissonDist(p * weight)};
      if (n_expected > 0) {
        if (n_expected > n) {
          wally_->Log("Rescaled n_bind_i_teth from %i to %i\n", n_expected, n);
          n_expected = n;
        }
        int n_removed{Set_Bind_I_Teth_Candidates(n_expected)};
        n_expected -= n_removed;
      }
      return n_expected;
    };
    events_.emplace_back(event_name, p_bind_i_teth_, &n_bind_i_teth_candidates_,
                         &bind_i_teth_candidates_, exe_bind_i_teth,
                         poisson_i_teth, set_ran_indices);
    // Unbind_I_Teth: same as above but for tethered populations
    auto exe_unbind_i_teth = [&](ENTRY_T target) {
      Unbind_I(std::get<POP_T *>(target));
    };
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
      for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
        std::string X_DUB{"_" + std::to_string(x_dub)};
        event_name = "unbind_i_teth" + N_NEIGHBS + X_DUB;
        events_.emplace_back(event_name, p_unbind_i_teth_[n_neighbs][x_dub],
                             &n_bound_i_teth_[n_neighbs][x_dub],
                             &bound_i_teth_[n_neighbs][x_dub],
                             exe_unbind_i_teth, binomial, set_ran_indices);
      }
    }
  }
  */
}

void AssociatedProteinManagement::InitializeEvents() {

  //  Binomial probabilitiy distribution; sampled to predict most events
  auto binomial = [&](double p, int n) {
    if (n > 0) {
      return properties_->gsl.SampleBinomialDist(p, n);
    } else {
      return 0;
    }
  };
  // Poisson distribution; sampled to predict events w/ variable probabilities
  auto poisson = [&](double p, int n) {
    if (p > 0.0) {
      return gsl_->SamplePoissonDist(p);
    } else {
      return 0;
    }
  };
  // Function that returns a random probability in the range [0.0, 1.0)
  auto get_ran_prob = [&](void) { return gsl_->GetRanProb(); };
  //  Function that sets n random indices from the range [0, m)
  auto set_ran_indices = [&](int *indices, int n, int m) {
    if (m > 1) {
      properties_->gsl.SetRanIndices(indices, n, m);
    } else {
      indices[0] = 0;
    }
  };
  // Event entries
  std::string event_name; // scratch space to construct each event name
  // Diffuse: steps head one site to the left or right
  auto exe_diffuse_fwd = [&](ENTRY_T target) {
    Diffuse(std::get<POP_T *>(target), 1);
  };
  auto exe_diffuse_bck = [&](ENTRY_T target) {
    Diffuse(std::get<POP_T *>(target), -1);
  };
  // Only loop up to max_neighbs - 1; p_diffuse is 0.0 with 2 neighbors
  for (int n_neighbs{0}; n_neighbs < max_neighbs_; n_neighbs++) {
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    // Currently we only consider neighb E-dep for diffusing AWAY from neighbs
    double p_diffuse_i{p_diffuse_i_ * weight_neighbs_unbind_[n_neighbs]};
    // Diffuse_i_fwd
    event_name = "diffuse_i_fwd" + N_NEIGHBS;
    events_.emplace_back(event_name, p_diffuse_i, &n_bound_i_[n_neighbs],
                         &bound_i_[n_neighbs], exe_diffuse_fwd, binomial,
                         set_ran_indices);
    // Diffuse_i_bck
    event_name = "diffuse_i_bck" + N_NEIGHBS;
    events_.emplace_back(event_name, p_diffuse_i, &n_bound_i_[n_neighbs],
                         &bound_i_[n_neighbs], exe_diffuse_bck, binomial,
                         set_ran_indices);
    if (!crosslinking_active_) {
      continue;
    }
    for (int x{0}; x <= dist_cutoff_; x++) {
      std::string X_DIST{"_" + std::to_string(x)};
      // Diffusing towards rest is considered an unbinding-type event in
      // regards to Boltzmann factors, since both events let the spring relax
      double p_diffuse_ii_to{p_diffuse_ii_ * weight_neighbs_unbind_[n_neighbs]};
      // Ensure index isn't negative to avoid seg-faults
      int x_to{x > 0 ? x - 1 : 0};
      // Dividing Boltzmann factors yields E(x_to) - E(x) in the exponential
      p_diffuse_ii_to *= weight_spring_unbind_[x_to] / weight_spring_unbind_[x];
      // Diffusing away from rest is considered a binding-type event in
      // regards to Boltzmann factors, since both events stretch the spring out
      double p_diffuse_ii_fr{p_diffuse_ii_ * weight_neighbs_unbind_[n_neighbs]};
      // Ensure index isn't out of range to avoid seg-faults
      int x_fr{x < dist_cutoff_ ? x + 1 : 0};
      p_diffuse_ii_fr *= weight_spring_bind_[x_fr] / weight_spring_bind_[x];
      // Can't diffuse to rest if we're at x = 0; both directions are fr rest
      if (x == 0) {
        p_diffuse_ii_to = 0.0;
        p_diffuse_ii_fr *= 2.0;
      }
      // Can't diffuse away from rest if we're at the cutoff distance
      if (x == dist_cutoff_) {
        p_diffuse_ii_fr = 0.0;
      }
      // Diffuse_ii_to_rest
      event_name = "diffuse_ii_to_rest" + N_NEIGHBS + X_DIST;
      events_.emplace_back(event_name, p_diffuse_ii_to,
                           &n_bound_ii_[n_neighbs][x], &bound_ii_[n_neighbs][x],
                           exe_diffuse_fwd, binomial, set_ran_indices);
      // Diffuse_ii_fr_rest
      event_name = "diffuse_ii_fr_rest" + N_NEIGHBS + X_DIST;
      events_.emplace_back(event_name, p_diffuse_ii_fr,
                           &n_bound_ii_[n_neighbs][x], &bound_ii_[n_neighbs][x],
                           exe_diffuse_bck, binomial, set_ran_indices);
    }
  }
  /*
  // Diffuse_teth: same as above but for tethered populations
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    if (!tethering_active_) {
      continue;
    }
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      std::string X_DUB{"_" + std::to_string(x_dub)};
      // Diffuse_i_to_teth rest
      event_name = "diffuse_i_to_teth_rest" + N_NEIGHBS + X_DUB;
      events_.emplace_back(
          event_name, p_diffuse_i_to_teth_rest_[n_neighbs][x_dub],
          &n_bound_i_teth_[n_neighbs][x_dub], &bound_i_teth_[n_neighbs][x_dub],
          exe_diffuse_fwd, binomial, set_ran_indices);
      // Diffuse_i_fr_teth_rest
      event_name = "diffuse_i_fr_teth_rest" + N_NEIGHBS + X_DUB;
      events_.emplace_back(
          event_name, p_diffuse_i_fr_teth_rest_[n_neighbs][x_dub],
          &n_bound_i_teth_[n_neighbs][x_dub], &bound_i_teth_[n_neighbs][x_dub],
          exe_diffuse_bck, binomial, set_ran_indices);
      if (!crosslinking_active_) {
        continue;
      }
      for (int x{0}; x <= dist_cutoff_; x++) {
        std::string X_DIST{"_" + std::to_string(x)};
        // Diffuse_ii_to_both_rests
        event_name = "diffuse_ii_to_both" + N_NEIGHBS + X_DUB + X_DIST;
        events_.emplace_back(event_name,
                             p_diffuse_ii_to_both_[n_neighbs][x_dub][x],
                             &n_bound_ii_teth_same_[n_neighbs][x_dub][x],
                             &bound_ii_teth_same_[n_neighbs][x_dub][x],
                             exe_diffuse_fwd, binomial, set_ran_indices);
        // Diffuse_ii_fr_both_rests
        event_name = "diffuse_ii_fr_both" + N_NEIGHBS + X_DUB + X_DIST;
        events_.emplace_back(event_name,
                             p_diffuse_ii_fr_both_[n_neighbs][x_dub][x],
                             &n_bound_ii_teth_same_[n_neighbs][x_dub][x],
                             &bound_ii_teth_same_[n_neighbs][x_dub][x],
                             exe_diffuse_bck, binomial, set_ran_indices);
        // Diffuse_ii_to_self_fr_teth_rest
        event_name = "diffuse_ii_to_self_fr_teth" + N_NEIGHBS + X_DUB + X_DIST;
        events_.emplace_back(event_name,
                             p_diffuse_ii_to_self_fr_teth_[n_neighbs][x_dub][x],
                             &n_bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                             &bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                             exe_diffuse_fwd, binomial, set_ran_indices);
        // Diffuse_ii_fr_self_to_teth_rest
        event_name = "diffuse_ii_fr_self_to_teth" + N_NEIGHBS + X_DUB + X_DIST;
        events_.emplace_back(event_name,
                             p_diffuse_ii_fr_self_to_teth_[n_neighbs][x_dub][x],
                             &n_bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                             &bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                             exe_diffuse_bck, binomial, set_ran_indices);
      }
    }
  }
  */
  // Bind_I: binds first crosslinker head to the microtubule
  auto exe_bind_i = [&](ENTRY_T target) { Bind_I(std::get<SITE_T *>(target)); };
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    event_name = "bind_i" + N_NEIGHBS;
    events_.emplace_back(event_name,
                         p_bind_i_ * weight_neighbs_bind_[n_neighbs],
                         &properties_->microtubules.n_unocc_xlink_[n_neighbs],
                         &properties_->microtubules.unocc_xlink_[n_neighbs],
                         exe_bind_i, binomial, set_ran_indices);
  }
  /*
  // Bind_I_Teth: same as above but for tethered unbound xlinks (satellites)
  if (tethering_active_) {
    event_name = "bind_i_teth";
    auto exe_bind_i_teth = [&](ENTRY_T target) {
      Bind_I_Teth(std::get<POP_T *>(target));
    };
    // Poisson dist. is used w/ partition function for E-dependent binding
    auto poisson_i_teth = [&](double p, int n) {
      double weight{GetWeight_Bind_I_Teth()};
      if (weight == 0.0) {
        return 0;
      }
      int n_expected{properties_->gsl.SamplePoissonDist(p * weight)};
      if (n_expected > 0) {
        if (n_expected > n) {
          n_expected = n;
        }
        Set_Bind_I_Teth_Candidates(n_expected);
      }
      return n_expected;
    };
    events_.emplace_back(event_name, p_bind_i_teth_, &n_bind_i_teth_candidates_,
                         &bind_i_teth_candidates_, exe_bind_i_teth,
                         poisson_i_teth, set_ran_indices);
  }
  */
  // Bind II: binds the second crosslinker head to adjacent microtubule
  if (crosslinking_active_) {
    event_name = "bind_ii";
    auto exe_bind_ii = [&](ENTRY_T target) {
      try {
        Bind_II(std::get<POP_T *>(target));
      } catch (...) {
        wally_->ErrorExit("AP_MGMT::Bind_II()");
      }
    };
    auto weight_bind_ii = [&](ENTRY_T target) {
      return std::get<POP_T *>(target)->xlink_->GetTotalWeight_Bind_II();
    };
    events_.emplace_back(event_name, p_bind_ii_, &n_bind_ii_candidates_,
                         &bind_ii_candidates_, exe_bind_ii, weight_bind_ii,
                         poisson, get_ran_prob, set_ran_indices);
  }
  /*
  // Bind_II_Teth: same as above but for tethered population
  if (tethering_active_ and crosslinking_active_) {
    event_name = "bind_ii_teth";
    auto exe_bind_ii_teth = [&](ENTRY_T target) {
      Bind_II(std::get<POP_T *>(target));
    };
    // Poisson dist. is used w/ partition function for E-dependent binding
    auto poisson_ii_teth = [&](double p, int n) {
      double weight{GetWeight_Bind_II_Teth()};
      if (weight == 0.0) {
        return 0;
      }
      int n_expected{properties_->gsl.SamplePoissonDist(p * weight)};
      if (n_expected > 0) {
        if (n_expected > n) {
          n_expected = n;
        }
        Set_Bind_II_Teth_Candidates(n_expected);
      }
      return n_expected;
    };
    events_.emplace_back(event_name, p_avg_bind_ii_,
                         &n_bind_ii_teth_candidates_, &bind_ii_teth_candidates_,
                         exe_bind_ii_teth, poisson_ii_teth, set_ran_indices);
  }
  */
  // Unbind II: unbinds a head of doubly-bound crosslinkers
  auto exe_unbind_ii = [&](ENTRY_T target) {
    Unbind_II(std::get<POP_T *>(target));
  };
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    if (!crosslinking_active_) {
      continue;
    }
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    for (int x{0}; x <= dist_cutoff_; x++) {
      std::string X_DIST{"_" + std::to_string(x)};
      double p_unbind_ii{p_unbind_ii_ * weight_neighbs_unbind_[n_neighbs] *
                         weight_spring_unbind_[x]};
      event_name = "unbind_ii" + N_NEIGHBS + X_DIST;
      events_.emplace_back(event_name, p_unbind_ii, &n_bound_ii_[n_neighbs][x],
                           &bound_ii_[n_neighbs][x], exe_unbind_ii, binomial,
                           set_ran_indices);
    }
  }
  /*
  // Unbind_II_Teth: same as above but for tethered populations
  auto exe_unbind_ii_teth = [&](ENTRY_T target) {
    Unbind_II(std::get<POP_T *>(target));
  };
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    if (!crosslinking_active_ or !tethering_active_) {
      continue;
    }
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      std::string X_DUB{"_" + std::to_string(x_dub)};
      for (int x{0}; x <= dist_cutoff_; x++) {
        std::string X_DIST{"_" + std::to_string(x)};
        event_name = "unbind_ii_to_teth" + N_NEIGHBS + X_DUB + X_DIST;
        events_.emplace_back(event_name,
                             p_unbind_ii_to_teth_[n_neighbs][x_dub][x],
                             &n_bound_ii_teth_same_[n_neighbs][x_dub][x],
                             &bound_ii_teth_same_[n_neighbs][x_dub][x],
                             exe_unbind_ii_teth, binomial, set_ran_indices);
        event_name = "unbind_ii_fr_teth" + N_NEIGHBS + X_DUB + X_DIST;
        events_.emplace_back(event_name,
                             p_unbind_ii_fr_teth_[n_neighbs][x_dub][x],
                             &n_bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                             &bound_ii_teth_oppo_[n_neighbs][x_dub][x],
                             exe_unbind_ii_teth, binomial, set_ran_indices);
      }
    }
  }
  */
  // Unbind_I: unbinds a head of singly-bound crosslinkres
  auto exe_unbind_i = [&](ENTRY_T target) {
    Unbind_I(std::get<POP_T *>(target));
  };
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    event_name = "unbind_i" + N_NEIGHBS;
    double p_unbind_i{p_unbind_i_ * weight_neighbs_unbind_[n_neighbs]};
    events_.emplace_back(event_name, p_unbind_i, &n_bound_i_[n_neighbs],
                         &bound_i_[n_neighbs], exe_unbind_i, binomial,
                         set_ran_indices);
  }
  /*
  // Unbind_I_Teth: same as above but for tethered populations
  auto exe_unbind_i_teth = [&](ENTRY_T target) {
    Unbind_I(std::get<POP_T *>(target));
  };
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    if (!tethering_active_) {
      continue;
    }
    std::string N_NEIGHBS{"_" + std::to_string(n_neighbs)};
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      std::string X_DUB{"_" + std::to_string(x_dub)};
      event_name = "unbind_i_teth" + N_NEIGHBS + X_DUB;
      events_.emplace_back(event_name, p_unbind_i_teth_[n_neighbs][x_dub],
                           &n_bound_i_teth_[n_neighbs][x_dub],
                           &bound_i_teth_[n_neighbs][x_dub], exe_unbind_i_teth,
                           binomial, set_ran_indices);
    }
  }
  */
  /*
if (tethering_active_) {
  // Tether_free: tethers free xlinks to bound motors
  event_name = "tether_free";
  auto exe_tether_free = [&](ENTRY_T target) {
    Tether_Free(std::get<ALT_T *>(target));
  };
  events_.emplace_back(event_name, p_tether_free_,
                       &properties_->kinesin4.n_bound_unteth_,
                       &properties_->kinesin4.bound_unteth_, exe_tether_free,
                       binomial, set_ran_indices);
  // Untether_satellite: untethers satellite xlinks from bound motors
  event_name = "untether_satellite";
  auto exe_untether_sat = [&](ENTRY_T target) {
    Untether(std::get<POP_T *>(target));
  };
  events_.emplace_back(event_name, p_untether_free_, &n_free_teth_,
                       &free_teth_, exe_untether_sat, binomial,
                       set_ran_indices);
}
                       */
}

void AssociatedProteinManagement::ReportProbabilities() {

  if (!population_active_) {
    return;
  }
  for (const auto &entry : p_theory_) {
    auto label = entry.first;
    auto value = entry.second;
    for (int i{0}; i < value.size(); i++) {
      for (int j{0}; j < value[i].size(); j++) {
        for (int k{0}; k < value[i][j].size(); k++) {
          std::string name{label};
          if (value.size() > 1) {
            name += "_" + std::to_string(i);
          }
          if (value[i].size() > 1) {
            name += "_" + std::to_string(j);
          }
          if (value[i][j].size() > 1) {
            name += "_" + std::to_string(k);
          }
          auto event = std::find_if(
              events_.begin(), events_.end(),
              [&name](EVENT_T const &event) { return event.name_ == name; });
          if (event->name_ != name) {
            wally_->Log("Couldn't find %s\n", name.c_str());
            continue;
          }
          if (event->n_opportunities_tot_ == 0) {
            wally_->Log("No statistics for %s\n", name.c_str());
            continue;
          }
          double n_exe_tot{(double)event->n_executed_tot_};
          double n_opp_tot{(double)event->n_opportunities_tot_};
          if (n_opp_tot < 0.0) {
            wally_->Log("WHAT?? n_opp = %g\n", n_opp_tot);
            exit(1);
          }
          wally_->Log("For AP event %s:\n", event->name_.c_str());
          wally_->Log("   p_theory = %g\n", value[i][j][k]);
          wally_->Log("   p_actual = %g", n_exe_tot / n_opp_tot);
          wally_->Log(" (n_exe = %i)\n", event->n_executed_tot_);
        }
      }
    }
  }
  /*
    if (wally_->test_mode_ != nullptr) {
      if (strcmp(wally_->test_mode_, "xlink_bind_ii") == 0) {
        wally_->Log("For bind_ii:\n");
        for (int x{0}; x <= dist_cutoff_; x++) {
          double p{(double)bind_ii_stats_[x].first / bind_ii_stats_[x].second};
          wally_->Log("p_theory_[%i] = %g\n", x,
                      p_avg_bind_ii_ * weight_bind_ii_[0][x]);
          wally_->Log("p_actual_[%i] = %g (n_exe = %i)\n", x, p,
                      bind_ii_stats_[x].first);
        }
      }
      if (strcmp(wally_->test_mode_, "xlink_bind_i_teth") == 0) {
        wally_->Log("For bind_i_teth:\n");
        for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
          double p{(double)bind_i_teth_stats_[x_dub].first /
                   bind_i_teth_stats_[x_dub].second};
          wally_->Log("p_theory_[%i] = %g\n", x_dub,
                      p_bind_i_teth_ * weight_bind_i_teth_[0][x_dub]);
          wally_->Log("p_actual_[%i] = %g (n_exe = %i)\n", x_dub, p,
                      bind_i_teth_stats_[x_dub].first);
        }
      }
      if (strcmp(wally_->test_mode_, "xlink_bind_ii_teth") == 0) {
        wally_->Log("For bind_ii_to_teth:\n");
        for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
          for (int x{0}; x <= dist_cutoff_; x++) {
            double p_to{(double)bind_ii_to_teth_stats_[x_dub][x].first /
                        bind_ii_to_teth_stats_[x_dub][x].second};
            wally_->Log("p_theory_[%i][%i] = %g\n", x_dub, x,
                        p_avg_bind_ii_ * weight_bind_ii_to_teth_[0][x_dub][x]);
            wally_->Log("p_actual_[%i][%i] = %g (n_exe = %i)\n", x_dub, x, p_to,
                        bind_ii_to_teth_stats_[x_dub][x].first);
          }
        }
        wally_->Log("For bind_ii_fr_teth:\n");
        for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
          for (int x{0}; x <= dist_cutoff_; x++) {
            double p_fr{(double)bind_ii_fr_teth_stats_[x_dub][x].first /
                        bind_ii_fr_teth_stats_[x_dub][x].second};
            wally_->Log("p_theory_[%i][%i] = %g\n", x_dub, x,
                        p_avg_bind_ii_ * weight_bind_ii_fr_teth_[0][x_dub][x]);
            wally_->Log("p_actual_[%i][%i] = %g (n_exe = %i)\n", x_dub, x, p_fr,
                        bind_ii_fr_teth_stats_[x_dub][x].first);
          }
        }
      }
    }
    */
}
AssociatedProtein *AssociatedProteinManagement::GetFreeXlink() {

  // Randomly choose an unbound xlink
  int i_xlink = properties_->gsl.GetRanInt(n_xlinks_);
  AssociatedProtein *xlink = &xlinks_[i_xlink];
  int attempts = 0;
  while (xlink->heads_active_ > 0 || xlink->tethered_) {
    i_xlink++;
    if (i_xlink == n_xlinks_)
      i_xlink = 0;
    xlink = &xlinks_[i_xlink];
    attempts++;
    if (attempts > n_xlinks_) {
      wally_->ErrorExit("AP_MGMT::GetFreeXlink()");
    }
  }
  return xlink;
}

void AssociatedProteinManagement::AddToActive(AssociatedProtein *xlink) {

  active_[n_active_] = xlink;
  xlink->active_index_ = n_active_;
  n_active_++;
}

void AssociatedProteinManagement::RemoveFromActive(AssociatedProtein *xlink) {

  AssociatedProtein *last_entry = active_[n_active_ - 1];
  int this_index = xlink->active_index_;
  active_[this_index] = last_entry;
  last_entry->active_index_ = this_index;
  n_active_--;
}

void AssociatedProteinManagement::FlagForUpdate() { lists_up_to_date_ = false; }

void AssociatedProteinManagement::UpdateExtensions() {

  if (!crosslinking_active_ and !tethering_active_) {
    return;
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    active_[i_entry]->UpdateExtension();
  }
}

void AssociatedProteinManagement::UpdateList_Bound_I() {

  wally_->Log(1, "Starting AP_MGMT::Update_Bound_I()\n");
  n_bind_ii_candidates_ = 0;
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    n_bound_i_[n_neighbs] = 0;
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    AssociatedProtein *xlink = active_[i_entry];
    if (xlink->heads_active_ == 1) {
      if (!xlink->tethered_ or xlink->HasSatellite()) {
        AssociatedProtein::Monomer *head = xlink->GetActiveHead();
        int n_neighbs = head->GetPRC1NeighborCount();
        bound_i_[n_neighbs][n_bound_i_[n_neighbs]++] = head;
        bind_ii_candidates_[n_bind_ii_candidates_++] = head;
      }
    }
  }
  if (wally_->test_mode_ != nullptr) {
    /*
    if (strcmp(wally_->test_mode_, "bind_ii") == 0) {
      // for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      //   for (int i_entry{0}; i_entry < n_bound_i_[n_neighbs]; i_entry++) {
      //     POP_T *head{std::get<POP_T *>(bound_i_[n_neighbs][i_entry])};
      //     if (head->site_->mt_ == 0) {
      //       int i_last{--n_bound_i_[n_neighbs]};
      //       bound_i_[n_neighbs][i_entry] = bound_i_[n_neighbs][i_last];
      //     }
      //   }
      // }
      for (int x{0}; x <= dist_cutoff_; x++) {
        if (x == 0) {
          // Only one site (directly above) can be bound to for this config
          bind_ii_stats_[x].second += n_bind_ii_candidates_;
        } else {
          // Two sites (to the left or right of above) are eligible for these
          // Include each as a separate opportunity
          bind_ii_stats_[x].second += 2 * n_bind_ii_candidates_;
        }
      }
    }
    */
  }
  /*
  if (verbosity_ >= 1) {
    wally_->Log(" n_bind_ii_candidates = %i\n", n_bind_ii_candidates_);
    for (int i_entry{0}; i_entry < n_bind_ii_candidates_; i_entry++) {
      POP_T *head{std::get<POP_T *>(bind_ii_candidates_[i_entry])};
      wally_->Log("   - xlink %i\n", head->xlink_->id_);
    }
    for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
      wally_->Log(" n_bound_i_[%i] = %i\n", n_neighbs, n_bound_i_[n_neighbs]);
      for (int i_entry{0}; i_entry < n_bound_i_[n_neighbs]; i_entry++) {
        POP_T *head{std::get<POP_T *>(bound_i_[n_neighbs][i_entry])};
        wally_->Log("   - xlink %i\n", head->xlink_->id_);
      }
    }
  }
  */
}

void AssociatedProteinManagement::UpdateList_Bound_I_Teth() {

  wally_->Log(1, "Starting AP_MGMT::Update_Bound_I_Teth()\n");
  n_bind_ii_teth_candidates_ = 0;
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      n_bound_i_teth_[n_neighbs][x_dub] = 0;
    }
  }
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    AssociatedProtein *xlink = active_[i_entry];
    if (xlink->heads_active_ != 1) {
      continue;
    }
    if (xlink->tethered_ and !xlink->HasSatellite()) {
      int x_dub = xlink->motor_->x_dist_doubled_;
      AssociatedProtein::Monomer *head = xlink->GetActiveHead();
      int n_neighbs = head->GetPRC1NeighborCount();
      int index{n_bound_i_teth_[n_neighbs][x_dub]++};
      bound_i_teth_[n_neighbs][x_dub][index] = head;
      bind_ii_teth_candidates_[n_bind_ii_teth_candidates_++] = head;
    }
  }
}

void AssociatedProteinManagement::UpdateList_Bound_II() {

  wally_->Log(1, "Starting AP_MGMT::Update_Bound_II()\n");
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    for (int x{0}; x <= dist_cutoff_; x++) {
      n_bound_ii_[n_neighbs][x] = 0;
    }
  }
  for (int i_xlink{0}; i_xlink < n_active_; i_xlink++) {
    AssociatedProtein *xlink = active_[i_xlink];
    if (xlink->heads_active_ != 2) {
      continue;
    }
    if (!xlink->tethered_ or xlink->HasSatellite()) {
      int x = xlink->x_dist_;
      // Head one
      if (wally_->test_mode_ == nullptr) {
        int neighbs_one = xlink->head_one_.GetPRC1NeighborCount();
        bound_ii_[neighbs_one][x][n_bound_ii_[neighbs_one][x]++] =
            &xlink->head_one_;
      }
      // Head two
      int neighbs_two = xlink->head_two_.GetPRC1NeighborCount();
      bound_ii_[neighbs_two][x][n_bound_ii_[neighbs_two][x]++] =
          &xlink->head_two_;
    }
  }
}

void AssociatedProteinManagement::UpdateList_Bound_II_Teth() {

  wally_->Log(1, "Starting AP_MGMT::Update_Bound_II_Teth()\n");
  for (int n_neighbs{0}; n_neighbs <= max_neighbs_; n_neighbs++) {
    for (int x_dub{2 * comp_cutoff_}; x_dub <= 2 * teth_cutoff_; x_dub++) {
      for (int x{0}; x <= dist_cutoff_; x++) {
        n_bound_ii_teth_same_[n_neighbs][x_dub][x] = 0;
        n_bound_ii_teth_oppo_[n_neighbs][x_dub][x] = 0;
      }
    }
  }
  for (int i_xlink{0}; i_xlink < n_active_; i_xlink++) {
    AssociatedProtein *xlink = active_[i_xlink];
    if (xlink->heads_active_ == 2 and xlink->tethered_) {
      // Dont count xlinks tethered to satellite motors
      if (xlink->motor_->heads_active_ > 0) {
        int x = xlink->x_dist_;
        int x_dub = xlink->motor_->x_dist_doubled_;
        // Site one
        Tubulin *site_one = xlink->head_one_.site_;
        int n_one = site_one->GetPRC1NeighborCount();
        if (x == rest_dist_) {
          int i_one = n_bound_ii_teth_same_[n_one][x_dub][x];
          bound_ii_teth_same_[n_one][x_dub][x][i_one] = &xlink->head_one_;
          n_bound_ii_teth_same_[n_one][x_dub][x]++;
          int i_two = n_bound_ii_teth_oppo_[n_one][x_dub][x];
          bound_ii_teth_oppo_[n_one][x_dub][x][i_two] = &xlink->head_one_;
          n_bound_ii_teth_oppo_[n_one][x_dub][x]++;
        } else if (site_one->EquilibriumInSameDirection()) {
          int i = n_bound_ii_teth_same_[n_one][x_dub][x];
          bound_ii_teth_same_[n_one][x_dub][x][i] = &xlink->head_one_;
          n_bound_ii_teth_same_[n_one][x_dub][x]++;
        } else {
          int i = n_bound_ii_teth_oppo_[n_one][x_dub][x];
          bound_ii_teth_oppo_[n_one][x_dub][x][i] = &xlink->head_one_;
          n_bound_ii_teth_oppo_[n_one][x_dub][x]++;
        }
        // Site two
        Tubulin *site_two = xlink->head_two_.site_;
        int n_two = site_two->GetPRC1NeighborCount();
        if (x == rest_dist_) {
          int i_one = n_bound_ii_teth_same_[n_two][x_dub][x];
          bound_ii_teth_same_[n_two][x_dub][x][i_one] = &xlink->head_two_;
          n_bound_ii_teth_same_[n_two][x_dub][x]++;
          int i_two = n_bound_ii_teth_oppo_[n_two][x_dub][x];
          bound_ii_teth_oppo_[n_two][x_dub][x][i_two] = &xlink->head_two_;
          n_bound_ii_teth_oppo_[n_two][x_dub][x]++;
        } else if (site_two->EquilibriumInSameDirection()) {
          int i = n_bound_ii_teth_same_[n_two][x_dub][x];
          bound_ii_teth_same_[n_two][x_dub][x][i] = &xlink->head_two_;
          n_bound_ii_teth_same_[n_two][x_dub][x]++;
        } else {
          int i = n_bound_ii_teth_oppo_[n_two][x_dub][x];
          bound_ii_teth_oppo_[n_two][x_dub][x][i] = &xlink->head_two_;
          n_bound_ii_teth_oppo_[n_two][x_dub][x]++;
        }
      }
    }
  }
}

void AssociatedProteinManagement::UpdateList_Bound_Unteth() {

  wally_->Log(1, "Starting AP_MGMT::Update_Bound_Unteth()\n");
  n_bound_unteth_ = 0;
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    AssociatedProtein *xlink = active_[i_entry];
    if (xlink->tethered_) {
      continue;
    }
    if (xlink->heads_active_ == 1) {
      bound_unteth_[n_bound_unteth_++] = xlink->GetActiveHead();
    } else {
      bound_unteth_[n_bound_unteth_++] = &xlink->head_one_;
    }
  }
}

void AssociatedProteinManagement::UpdateList_Free_Teth() {

  wally_->Log(1, "Starting AP_MGMT::Update_Free_Teth()\n");
  n_satellites_ = 0;
  n_bind_i_teth_candidates_ = 0;
  for (int i_entry{0}; i_entry < n_active_; i_entry++) {
    AssociatedProtein *xlink = active_[i_entry];
    if (xlink->heads_active_ == 0 and xlink->tethered_) {
      if (xlink->motor_->heads_active_ == 0) {
        wally_->ErrorExit("AP_MGMT::Update_Free_Teth()\n");
      }
      satellites_[n_satellites_++] = &xlink->head_one_;
      bind_i_teth_candidates_[n_bind_i_teth_candidates_++] = &xlink->head_one_;
    }
  }
}

void AssociatedProteinManagement::RunKMC() {

  if (!population_active_) {
    return;
  }
  UpdateLists();
  SampleEventStatistics();
  GenerateExecutionSequence();
  ExecuteEvents();
}

void AssociatedProteinManagement::UpdateLists() {

  /*
  if (lists_up_to_date_) {
    return;
  }
  */
  wally_->Log(1, "\nStarting AP_MGMT::UpdateLists()\n");
  // lists_up_to_date_ = true;
  UpdateExtensions();
  properties_->microtubules.UpdateUnoccupied();
  UpdateList_Bound_I();
  if (crosslinking_active_) {
    UpdateList_Bound_II();
  }
  if (tethering_active_) {
    UpdateList_Free_Teth();
    UpdateList_Bound_Unteth();
    UpdateList_Bound_I_Teth();
    if (crosslinking_active_) {
      UpdateList_Bound_II_Teth();
    }
    properties_->kinesin4.UpdateList_Bound_Teth();
  }
}

void AssociatedProteinManagement::SampleEventStatistics() {

  wally_->Log(1, "Starting AP_MGMT::SampleEventStatistics()\n");
  n_events_to_exe_ = 0;
  for (auto &&event : events_) {
    wally_->Log(2, "Sampling statistics for %s\n", event.name_.c_str());
    n_events_to_exe_ += event.SampleStatistics();
    if (event.n_expected_ > 0) {
      // Ensure no funny business occurs with statistics
      if (event.n_expected_ > *event.n_avail_) {
        wally_->Log("Error; %i events expected but only %i available for %s\n",
                    event.n_expected_, *event.n_avail_, event.name_.c_str());
        wally_->ErrorExit("AP_MGMT::SampleEventStatistics()");
      }
      wally_->Log(2, " %i events expected for %s (%i avail)\n",
                  event.n_expected_, event.name_.c_str(), *event.n_avail_);
      for (int i_entry{0}; i_entry < event.n_expected_; i_entry++) {
        POP_T *head{nullptr};
        try {
          head = std::get<POP_T *>(event.targets_[i_entry]);
        } catch (...) {
          unsigned long index{event.targets_[i_entry].index()};
          wally_->Log(2, "    -target %i: INDEX %lu ?!\n", i_entry, index);
          continue;
        }
        wally_->Log(2, "    -target %i: xlink %i\n", i_entry,
                    head->xlink_->id_);
      }
    }
  }
  wally_->Log(2, " %i events expected pre-correction\n", n_events_to_exe_);
  if (n_events_to_exe_ <= 1) {
    return;
  }
  // Put all events with >1 target into an active_ array
  std::pair<EVENT_T *, ENTRY_T> active_events[n_events_to_exe_];
  int i_active{0};
  for (auto &&event : events_) {
    // Add a ptr to the event for each target it has
    for (int i_tar{0}; i_tar < event.n_expected_; i_tar++) {
      active_events[i_active++] = std::make_pair(&event, event.targets_[i_tar]);
      wally_->Log(2, " Added %s to active_events\n", event.name_.c_str());
    }
  }
  // Scan through all active events to ensure that no two target the same
  // xlink
  for (int i_entry{0}; i_entry < n_events_to_exe_; i_entry++) {
    EVENT_T *event_i = active_events[i_entry].first;
    AssociatedProtein *xlink_i{nullptr};
    try {
      xlink_i = std::get<POP_T *>(active_events[i_entry].second)->xlink_;
    } catch (...) {
      continue;
    }
    for (int j_entry{i_entry + 1}; j_entry < n_events_to_exe_; j_entry++) {
      EVENT_T *event_j = active_events[j_entry].first;
      AssociatedProtein *xlink_j{nullptr};
      try {
        xlink_j = std::get<POP_T *>(active_events[j_entry].second)->xlink_;
      } catch (...) {
        continue;
      }
      // If event_i and event_j target different xlinks, continue
      if (xlink_i != xlink_j) {
        continue;
      }
      double p_one{event_i->p_occur_};
      double p_two{event_j->p_occur_};
      double ran{properties_->gsl.GetRanProb()};
      if (ran < p_one / (p_one + p_two)) {
        wally_->Log(2, "Removed xlink #%i from %s targets\n", xlink_i->id_,
                    event_i->name_.c_str());
        event_i->RemoveTarget(active_events[i_entry].second);
        active_events[i_entry] = active_events[n_events_to_exe_ - 1];
        i_entry--;
        n_events_to_exe_--;
        break;
      } else {
        wally_->Log(2, "Removed xlink #%i from %s targets\n", xlink_j->id_,
                    event_j->name_.c_str());
        event_j->RemoveTarget(active_events[j_entry].second);
        active_events[j_entry] = active_events[n_events_to_exe_ - 1];
        j_entry--;
        n_events_to_exe_--;
      }
    }
  }
}

void AssociatedProteinManagement::GenerateExecutionSequence() {

  wally_->Log(1, "Starting AP_MGMT::GenerateExecutionSequence()\n");
  if (n_events_to_exe_ == 0) {
    return;
  }
  int i_array{0};
  EVENT_T *pre_array[n_events_to_exe_];
  for (auto &&event : events_) {
    for (int i_entry{0}; i_entry < event.n_expected_; i_entry++) {
      pre_array[i_array++] = &event;
    }
  }
  if (i_array != n_events_to_exe_) {
    wally_->ErrorExit("AP_MGMT::GenerateExecutionSequence()");
  }
  if (n_events_to_exe_ > 1) {
    properties_->gsl.Shuffle(pre_array, n_events_to_exe_, sizeof(EVENT_T *));
  }
  if (n_events_to_exe_ > events_to_exe_.size()) {
    events_to_exe_.resize(n_events_to_exe_);
  }
  for (int i_event{0}; i_event < n_events_to_exe_; i_event++) {
    events_to_exe_[i_event] = pre_array[i_event];
  }
}

void AssociatedProteinManagement::ExecuteEvents() {

  wally_->Log(1, "Starting AP_MGMT::ExecuteEvents()\n");
  for (int i_event{0}; i_event < n_events_to_exe_; i_event++) {
    events_to_exe_[i_event]->Execute();
    FlagForUpdate();
  }
}

void AssociatedProteinManagement::Diffuse(POP_T *head, int dir) {

  wally_->Log(1, "Executing AP_MGMT::Diffuse()\n");
  if (head->site_ == nullptr) {
    wally_->Log("woah buddy; u cant diffuse on a NULL site!!!\n");
    return;
  }
  Tubulin *old_site = head->site_;
  int dx{dir * head->GetDirectionTowardRest()};
  // Cannot step off them MTs
  if ((old_site->index_ == 0 and dx == -1) or
      (old_site->index_ == old_site->mt_->n_sites_ - 1 and dx == 1)) {
    return;
  }
  Tubulin *new_site = &old_site->mt_->lattice_[old_site->index_ + dx];
  if (!new_site->occupied_) {
    old_site->xlink_head_ = nullptr;
    old_site->occupied_ = false;
    new_site->xlink_head_ = head;
    new_site->occupied_ = true;
    head->site_ = new_site;
    properties_->microtubules.FlagForUpdate();
  }
}

void AssociatedProteinManagement::Bind_I(SITE_T *site) {

  wally_->Log(1, "Executing AP_MGMT::Bind_I()\n");
  AssociatedProtein *xlink{GetFreeXlink()};
  site->xlink_head_ = &xlink->head_one_;
  site->occupied_ = true;
  xlink->heads_active_++;
  xlink->head_one_.site_ = site;
  AddToActive(xlink);
  properties_->microtubules.FlagForUpdate();
}

void AssociatedProteinManagement::Bind_I_Teth(POP_T *satellite_head) {

  wally_->Log(1, "Executing AP_MGMT::Bind_I_Teth()\n");
  AssociatedProtein *xlink{satellite_head->xlink_};
  Tubulin *site{xlink->GetWeightedSite_Bind_I_Teth()};
  if (site == nullptr) {
    wally_->Log("failed to Bind_I_Free_Tethered (XLINK)\n");
    return;
  }
  site->xlink_head_ = satellite_head;
  site->occupied_ = true;
  satellite_head->site_ = site;
  xlink->heads_active_++;
  properties_->microtubules.FlagForUpdate();
}

void AssociatedProteinManagement::Bind_II(POP_T *bound_head) {

  wally_->Log(1, "Executing AP_MGMT::Bind_II()\n");
  POP_T *unbound_head{bound_head->GetOtherHead()};
  AssociatedProtein *xlink{bound_head->xlink_};
  Tubulin *site{nullptr};
  bool has_satellite{false};
  if (xlink->tethered_) {
    if (xlink->motor_->heads_active_ == 0) {
      has_satellite = true;
    }
  }
  if (!xlink->tethered_ or has_satellite) {
    site = bound_head->xlink_->GetWeightedSite_Bind_II();
  } else {
    site = bound_head->xlink_->GetWeightedSite_Bind_II_Teth();
  }
  // In case another crosslinker binds and takes up the only available neighbor
  if (site == nullptr) {
    return;
    // wally_->ErrorExit("AP_MGMT::Bind_II()");
  }
  site->xlink_head_ = unbound_head;
  site->occupied_ = true;
  unbound_head->site_ = site;
  xlink->heads_active_++;
  properties_->microtubules.FlagForUpdate();
}

void AssociatedProteinManagement::Unbind_II(POP_T *head) {

  wally_->Log(1, "Executing AP_MGMT::Unbind_II()\n");
  AssociatedProtein *xlink{head->xlink_};
  Tubulin *site{head->site_};
  site->xlink_head_ = nullptr;
  site->occupied_ = false;
  head->site_ = nullptr;
  xlink->heads_active_--;
  properties_->microtubules.FlagForUpdate();
}

void AssociatedProteinManagement::Unbind_I(POP_T *head) {

  wally_->Log(1, "Executing AP_MGMT::Unbind_I()\n");
  AssociatedProtein *xlink{head->xlink_};
  Tubulin *site{head->site_};
  site->xlink_head_ = nullptr;
  site->occupied_ = false;
  head->site_ = nullptr;
  xlink->heads_active_--;
  xlink->UntetherSatellite();
  if (!xlink->tethered_) {
    RemoveFromActive(xlink);
  }
  properties_->microtubules.FlagForUpdate();
}

void AssociatedProteinManagement::Tether_Free(ALT_T *untethered_head) {

  wally_->Log(1, "Executing AP_MGMT::Tether_Free()\n");
  AssociatedProtein *xlink{GetFreeXlink()};
  Kinesin *motor{untethered_head->motor_};
  xlink->motor_ = motor;
  xlink->tethered_ = true;
  motor->xlink_ = xlink;
  motor->tethered_ = true;
  AddToActive(xlink);
  properties_->kinesin4.FlagForUpdate();
}

void AssociatedProteinManagement::Untether(POP_T *satellite_head) {

  wally_->Log(1, "Executing AP_MGMT::Untether()\n");
  AssociatedProtein *xlink{satellite_head->xlink_};
  Kinesin *motor{xlink->motor_};
  xlink->motor_ = nullptr;
  xlink->tethered_ = false;
  motor->xlink_ = nullptr;
  motor->tethered_ = false;
  RemoveFromActive(xlink);
  properties_->kinesin4.FlagForUpdate();
}
