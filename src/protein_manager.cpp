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
  // Array of MAXIMUM possible weights including neighbor interactions and
  // lattice deformation effects from MULTIPLE motor effects stacking
  double lattice_E_max_tot{-1 * Params::Motors::gaussian_ceiling_bulk};
  double wt_lattice_bind_max{exp(-(1.0 - _lambda_lattice) * lattice_E_max_tot)};
  double wt_lattice_unbind_max{exp(_lambda_lattice * lattice_E_max_tot)};
  Sys::weight_lattice_bind_max_.resize(_n_neighbs_max + 1);
  Sys::weight_lattice_unbind_max_.resize(_n_neighbs_max + 1);
  for (int n_neighbs{0}; n_neighbs <= _n_neighbs_max; n_neighbs++) {
    Sys::weight_lattice_bind_max_[n_neighbs] =
        wt_lattice_bind_max * Sys::weight_neighb_bind_[n_neighbs];
    Sys::weight_lattice_unbind_max_[n_neighbs] =
        wt_lattice_unbind_max * Sys::weight_neighb_unbind_[n_neighbs];
  }
  // Calculate lattice_alpha based on input range
  Sys::lattice_cutoff_ = Params::Motors::gaussian_range;
  double dx_cutoff{(Sys::lattice_cutoff_ + 1) * Params::Filaments::site_size};
  double lattice_E_0_solo{-1 * Params::Motors::gaussian_amp_solo};
  double lattice_alpha{-1 * lattice_E_0_solo / (dx_cutoff * dx_cutoff)};
  if (Sys::lattice_cutoff_ > 0) {
    Sys::Log("  lattice_cutoff_ is %i\n", Sys::lattice_cutoff_);
    Sys::Log("  lattice_alpha_ is %g\n", lattice_alpha);
  } else {
    Sys::Log("  Lattice cooperativity is disabled.\n");
  }
  // Array of binding/unbinding weights due to lattice deformation from a SINGLE
  // motor; multiplied together to get total weight for any arrangement
  // Array of binding/unbinding weights due to lattice deformation from a SINGLE
  // motor; multiplied together to get total weight for any arrangement
  Sys::weight_lattice_bind_.resize(Sys::lattice_cutoff_ + 1);
  Sys::weight_lattice_unbind_.resize(Sys::lattice_cutoff_ + 1);
  for (int delta{0}; delta <= Sys::lattice_cutoff_; delta++) {
    double dx{delta * Params::Filaments::site_size};
    double energy{lattice_alpha * dx * dx + lattice_E_0_solo}; // in kbT
    Sys::weight_lattice_bind_[delta] = exp(-(1.0 - _lambda_lattice) * energy);
    Sys::weight_lattice_unbind_[delta] = exp(_lambda_lattice * energy);
    // printf("weight = %#.3g (%#.3g)\n", Sys::weight_lattice_bind_[delta],
    //        Sys::weight_lattice_unbind_[delta]);
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
  double wt_unbind_ii{exp(Motors::internal_force * Motors::sigma_off_ii / kbT)};
  motors_.AddProb("unbind_ii", Motors::k_off_ii * dt * wt_unbind_ii);
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

  Sys::Log("\n");
  Sys::Log("Initializing test '%s'\n", Sys::test_mode_.c_str());
  using namespace Params;
  if (Sys::test_mode_ == "xlink_bind_ii") {
    Xlinks::k_off_ii = 14.3;
    GenerateReservoirs();
    InitializeWeights();
    SetParameters();
    double r_y{std::fabs(Filaments::y_initial[0] - Filaments::y_initial[1])};
    double r_max{xlinks_.r_max_};
    printf("r_max = %g\n", r_max);
    double r_x_max{sqrt(Square(r_max) - Square(r_y))};
    printf("r_x_max = %g\n", r_x_max);
    size_t x_max((size_t)std::ceil(r_x_max / Filaments::site_size));
    printf("x_max = %zu\n", x_max);
    // Initialize filament environment
    Filaments::count = 2;
    Filaments::n_sites[0] = Filaments::n_sites[1] = 2 * x_max + 1;
    Sys::Log("  N_SITES[0] = %i\n", Filaments::n_sites[0]);
    Sys::Log("  N_SITES[1] = %i\n", Filaments::n_sites[1]);
    Filaments::translation_enabled[0] = false;
    Filaments::translation_enabled[1] = false;
    Filaments::rotation_enabled = false;
    filaments_->Initialize(this);
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
      p_unbind[x] *= exp(lambda * dE / Params::kbT);
    }
    test_ref_.emplace("bind_ii", p_bind);
    test_ref_.emplace("unbind_ii", p_unbind);
    Vec<Pair<size_t, size_t>> zeros(x_max + 1, {0, 0});
    test_stats_.emplace("bind_ii", zeros);
    test_stats_.emplace("unbind_ii", zeros);
    // Place first xlink head on lower MT; remains static for entire sim
    int i_site{(int)x_max};
    BindingSite *site{&filaments_->proto_[0].sites_[i_site]};
    Protein *xlink{xlinks_.GetFreeEntry()};
    bool executed{xlink->Bind(site, &xlink->head_one_)};
    if (executed) {
      xlinks_.AddToActive(xlink);
      filaments_->FlagForUpdate();
    } else {
      Sys::ErrorExit("ProteinManager::InitializeTestEnvironment()");
    }
  } else if (Sys::test_mode_ == "xlink_diffusion") {
    // Initialize sim objects
    Xlinks::c_bulk = 1.0;
    Xlinks::t_active = 0.0;
    GenerateReservoirs();
    InitializeWeights();
    SetParameters();
    // Initialize filaments
    Filaments::count = 2;
    Filaments::n_sites[0] = Filaments::n_sites[1] = 100;
    Sys::Log("  N_SITES[0] = %i\n", Filaments::n_sites[0]);
    Sys::Log("  N_SITES[1] = %i\n", Filaments::n_sites[1]);
    filaments_->Initialize(this);
    // Initialize stat trackers
    double r_y{std::fabs(Filaments::y_initial[0] - Filaments::y_initial[1])};
    double r_max{xlinks_.r_max_};
    printf("r_max = %g\n", r_max);
    double r_x_max{sqrt(Square(r_max) - Square(r_y))};
    printf("r_x_max = %g\n", r_x_max);
    size_t x_max((size_t)std::ceil(r_x_max / Filaments::site_size));
    printf("x_max = %zu\n", x_max);
    Vec<double> p_theory_to(x_max + 1,
                            xlinks_.p_event_.at("diffuse_ii_to_rest").GetVal());
    Vec<double> p_theory_fr(x_max + 1,
                            xlinks_.p_event_.at("diffuse_ii_fr_rest").GetVal());
    Vec<double> wt_bind(x_max + 1, 0.0);
    Vec<double> wt_unbind(x_max + 1, 0.0);
    for (int x{0}; x <= x_max; x++) {
      double r_x{x * Filaments::site_size};
      double r{sqrt(Square(r_x) + Square(r_y))};
      if (r < xlinks_.r_min_ or r > xlinks_.r_max_) {
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
    // /*
    for (auto &&entry : test_stats_) {
      printf("For '%s':\n", entry.first.c_str());
      for (int x{0}; x < entry.second.size(); x++) {
        printf("  [%i] = {%zu, %zu}\n", x, entry.second[x].first,
               entry.second[x].second);
      }
    }
    // */
    Protein *xlink{xlinks_.GetFreeEntry()};
    int i_site{(int)std::round(Filaments::n_sites[0] / 2)};
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
  } else if (Sys::test_mode_ == "motor_lattice_bind") {
    size_t cutoff{Motors::gaussian_range};
    // Set parameters
    Xlinks::c_bulk = 0.0;
    Motors::k_on = 1.0;
    Motors::c_bulk = 10.0;
    Motors::neighb_neighb_energy = 0.0;
    Xlinks::neighb_neighb_energy = 0.0;
    Filaments::count = 1;
    Filaments::n_sites[0] = 2 * cutoff + 1;
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
    BindingSite *site{&filaments_->proto_[0].sites_[i_site]};
    Motor *motor{motors_.GetFreeEntry()};
    bool executed{motor->Bind(site, &motor->head_one_)};
    if (executed) {
      motors_.AddToActive(motor);
      filaments_->FlagForUpdate();
    } else {
      Sys::ErrorExit("ProteinManager::InitializeTestEnvironment()");
    }
  } else if (Sys::test_mode_ == "motor_lattice_step") {
    Motors::c_bulk = 1.0;
    Motors::t_active = 0.0;
    // Initialize sim objects
    GenerateReservoirs();
    InitializeWeights();
    SetParameters();
    // Initialize filaments
    Filaments::count = 1;
    Filaments::n_sites[0] = 1000;
    Filaments::translation_enabled[0] = false;
    Filaments::translation_enabled[1] = false;
    Filaments::rotation_enabled = false;
    Sys::Log("  N_SITES[0] = %i\n", Filaments::n_sites[0]);
    filaments_->Initialize(this);
    printf("Enter test delta (-1 to check against self-coop): ");
    Str response;
    std::getline(std::cin, response);
    int test_delta{(int)std::stoi(response)};
    double test_weight_bind{1.0};
    double test_weight_unbind{1.0};
    if (test_delta > 0) {
      test_weight_bind = Sys::weight_lattice_bind_[test_delta];
      test_weight_unbind = Sys::weight_lattice_unbind_[test_delta];
      for (int delta{0}; delta <= Motors::gaussian_range; delta++) {
        Sys::weight_lattice_bind_[delta] = test_weight_bind;
        Sys::weight_lattice_unbind_[delta] = test_weight_unbind;
      }
    }
    for (auto &&site : filaments_->sites_) {
      site->SetWeight_Bind(test_weight_bind);
      site->SetWeight_Unbind(test_weight_unbind);
    }
    // Initialize statistic trackers
    Vec<Pair<size_t, size_t>> zeros(1, {0, 0});
    Vec<double> p_theory_bind_ii(1, motors_.p_event_.at("bind_ii").GetVal() *
                                        test_weight_bind);
    test_stats_.emplace("bind_ii", zeros);
    test_ref_.emplace("bind_ii", p_theory_bind_ii);
    Vec<double> p_theory_unbind_ii(1,
                                   motors_.p_event_.at("unbind_ii").GetVal() *
                                       Square(test_weight_unbind));
    test_stats_.emplace("unbind_ii", zeros);
    test_ref_.emplace("unbind_ii", p_theory_unbind_ii);
    Vec<double> p_theory_unbind_i(1, motors_.p_event_.at("unbind_i").GetVal() *
                                         test_weight_unbind);
    test_stats_.emplace("unbind_i", zeros);
    test_ref_.emplace("unbind_i", p_theory_unbind_i);
    // Place motor head on minus end of microtubule
    BindingSite *site{filaments_->proto_[0].minus_end_};
    Motor *motor{motors_.GetFreeEntry()};
    bool executed{motor->Bind(site, &motor->head_one_)};
    if (executed) {
      motors_.AddToActive(motor);
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
    Filaments::translation_enabled[0] = false;
    Filaments::translation_enabled[1] = true;
    Filaments::rotation_enabled = false; // true;
    GenerateReservoirs();
    InitializeWeights();
    SetParameters();
    filaments_->Initialize(this);
    int n_xlinks{Sys::n_xlinks_};
    if (n_xlinks == -1) {
      Str response;
      printf("Enter number of crosslinkers: ");
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
  } else if (Sys::test_mode_ == "filament_ablation") {
    Filaments::count = 2;
    Sys::Log("  COUNT = 2\n");
    Filaments::n_sites[0] = Filaments::n_sites[1] = 125;
    Filaments::polarity[0] = Filaments::polarity[1] = 0;
    Sys::Log("  N_SITES[0] = %i\n", Filaments::n_sites[0]);
    Sys::Log("  N_SITES[1] = %i\n", Filaments::n_sites[1]);
    Filaments::translation_enabled[0] = false;
    Filaments::translation_enabled[1] = false;
    Filaments::rotation_enabled = false;
    Filaments::x_initial[0] = 0.0;
    Filaments::x_initial[1] =
        (Filaments::n_sites[0] - 1) * Filaments::site_size;
    Filaments::y_initial[0] = Filaments::y_initial[1] = 0.0;
    Motors::endpausing_active = false;
    printf("Enter ablation time: ");
    Str response;
    std::getline(std::cin, response);
    double t_ablate{(double)std::stod(response)};
    Sys::ablation_step_ = size_t(std::round(t_ablate / dt));
    GenerateReservoirs();
    InitializeWeights();
    SetParameters();
    filaments_->Initialize(this);
  } else if (Sys::test_mode_ == "hetero_tubulin") {
    printf("Enter fraction of heterogenous tubulin: ");
    Str response_one;
    std::getline(std::cin, response_one);
    double p_hetero{(double)std::stod(response_one)};
    if (p_hetero < 0.0 or p_hetero > 1.0) {
      printf("Invalid fraction, ya dingus!\n");
      exit(1);
    }
    printf("Enter decrease in binding affinity ");
    printf("(e.g., 2 will cut p_bind in half): ");
    Str response_two;
    std::getline(std::cin, response_two);
    double bind_aff{(double)std::stod(response_two)};
    if (bind_aff < 0.0) {
      printf("Error. Fractional change must be positive!\n");
      exit(1);
    }
    if (bind_aff == 0.0) {
      printf("You tryna start a damn singularity?!\n");
      exit(1);
    }
    GenerateReservoirs();
    InitializeWeights();
    SetParameters();
    filaments_->Initialize(this);
    int n_sites{(int)filaments_->sites_.size()};
    int n_hetero{(int)std::round(n_sites * p_hetero)};
    // printf("n_hetero = %i\n", n_hetero);
    // Randomly place heterogeneous sites on lattice
    int site_indices[n_sites];
    for (int index{0}; index < n_sites; index++) {
      site_indices[index] = index;
      // printf("i = %i\n", index);
    }
    SysRNG::Shuffle(site_indices, n_sites, sizeof(int));
    for (int i_hetero{0}; i_hetero < n_hetero; i_hetero++) {
      int i_site{site_indices[i_hetero]};
      // printf("i_site = %i\n", i_site);
      filaments_->sites_[i_site]->SetBindingAffinity(bind_aff);
      // printf("SITE %zu IS A HETERO\n", filaments_->sites_[i_site]->index_);
    }
    // exit(1);
  } else if (Sys::test_mode_ == "kinesin_mutant") {
    GenerateReservoirs();
    InitializeWeights();
    SetParameters();
    filaments_->Initialize(this);
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
        int n_entries{c * (int)xlinks_.sorted_.at("bind_ii").size_};
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
        size_t x{
            (size_t)std::abs(std::round(r_x / Params::Filaments::site_size))};
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
        size_t x{
            (size_t)std::abs(std::round(r_x / Params::Filaments::site_size))};
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
        size_t x{
            (size_t)std::abs(std::round(r_x / Params::Filaments::site_size))};
        test_stats_.at("unbind_ii")[x].first++;
        // bool head->parent_->UpdateExtension();
      } else {
        Sys::ErrorExit("Unbind_II (TEST)");
      }
    };
    kmc_.events_.emplace_back(
        "unbind_ii", xlinks_.p_event_.at("unbind_ii").GetVal(),
        &xlinks_.sorted_.at("unbind_ii").size_,
        &xlinks_.sorted_.at("unbind_ii").entries_, poisson_unbind_ii,
        get_weight_unbind_ii, exe_unbind_ii);
  } else if (Sys::test_mode_ == "xlink_diffusion") {
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
      size_t x{
          (size_t)std::abs(std::round(r_x / Params::Filaments::site_size))};
      bool executed{head->Diffuse(1)};
      if (executed) {
        bool still_attached{head->parent_->UpdateExtension()};
        xlinks_.FlagForUpdate();
        filaments_->FlagForUpdate();
        test_stats_.at("to_rest")[x].first++;
      }
    };
    auto exe_diff_fr = [&](Object *base) {
      BindingHead *head{dynamic_cast<BindingHead *>(base)};
      double r_x{head->pos_[0] - head->GetOtherHead()->pos_[0]};
      size_t x{
          (size_t)std::abs(std::round(r_x / Params::Filaments::site_size))};
      bool executed{head->Diffuse(-1)};
      if (executed) {
        bool still_attached{head->parent_->UpdateExtension()};
        xlinks_.FlagForUpdate();
        filaments_->FlagForUpdate();
        test_stats_.at("fr_rest")[x].first++;
      }
    };
    auto weight_diff_ii = [](auto *head, int dir) {
      return head->GetWeight_Diffuse(dir);
    };
    auto poisson_to = [&](double p, int n) {
      Protein *xlink{xlinks_.active_entries_[0]};
      double r_x{xlink->head_one_.pos_[0] - xlink->head_two_.pos_[0]};
      size_t x{
          (size_t)std::abs(std::round(r_x / Params::Filaments::site_size))};
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
      size_t x{
          (size_t)std::abs(std::round(r_x / Params::Filaments::site_size))};
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
        "diffuse_ii_to_rest",
        xlinks_.p_event_.at("diffuse_ii_to_rest").GetVal(),
        &xlinks_.sorted_.at("diffuse_ii_to_rest").size_,
        &xlinks_.sorted_.at("diffuse_ii_to_rest").entries_, poisson_to,
        [&](Object *base) {
          return weight_diff_ii(dynamic_cast<BindingHead *>(base), 1);
        },
        exe_diff_to);
    kmc_.events_.emplace_back(
        "diffuse_ii_fr_rest",
        xlinks_.p_event_.at("diffuse_ii_fr_rest").GetVal(),
        &xlinks_.sorted_.at("diffuse_ii_fr_rest").size_,
        &xlinks_.sorted_.at("diffuse_ii_fr_rest").entries_, poisson_fr,
        [&](Object *base) {
          return weight_diff_ii(dynamic_cast<BindingHead *>(base), -1);
        },
        exe_diff_fr);
  } else if (Sys::test_mode_ == "motor_lattice_bind") {
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
  } else if (Sys::test_mode_ == "motor_lattice_step") {
    auto binomial = [&](double p, int n) {
      if (n > 0) {
        return SysRNG::SampleBinomial(p, n);
      } else {
        return 0;
      }
    };
    // Bind_ATP_I
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
    // Bind_ATP_II
    auto poisson_ATP = [&](double p, int n) {
      if (p > 0.0) {
        return SysRNG::SamplePoisson(p);
      } else {
        return 0;
      }
    };
    auto exe_bind_ATP_ii = [](auto *front_head, auto *pop, auto *fil) {
      auto *rear_head{front_head->GetOtherHead()};
      if (front_head->trailing_) {
        return;
      }
      bool unbound{rear_head->Unbind()};
      bool executed{front_head->parent_->Bind_ATP(front_head)};
      if (executed) {
        pop->FlagForUpdate();
        fil->FlagForUpdate();
      }
    };
    auto weight_bind_ATP_ii = [](auto *head) {
      return head->parent_->GetWeight_BindATP_II(head);
    };
    auto is_NULL_ii_bound = [](auto *motor) -> Vec<Object *> {
      if (motor->n_heads_active_ == 2) {
        bool found_head{false};
        CatalyticHead *chosen_head{nullptr};
        if (motor->head_one_.ligand_ == CatalyticHead::Ligand::NONE) {
          chosen_head = &motor->head_one_;
        }
        if (motor->head_two_.ligand_ == CatalyticHead::Ligand::NONE) {
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
          exe_bind_ATP_ii(dynamic_cast<CatalyticHead *>(base), &motors_,
                          filaments_);
        });
    // Hydrolyze
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
    // Bind_II
    auto exe_bind_ii = [&](Object *base) {
      auto bound_head{dynamic_cast<CatalyticHead *>(base)};
      auto head{bound_head->GetOtherHead()};
      auto site{head->parent_->GetNeighbor_Bind_II()};
      // If dock site is plus end, unbind motor and place it on minus end
      if (site == site->filament_->plus_end_) {
        bool exe1{bound_head->Unbind()};
        auto new_site{site->filament_->minus_end_};
        bool exe2{bound_head->parent_->Bind(new_site, bound_head)};
        bool exe3{bound_head->parent_->Bind_ATP(bound_head)};
        bool exe4{bound_head->parent_->Hydrolyze(bound_head)};
        site = bound_head->parent_->GetNeighbor_Bind_II();
      }
      auto executed{head->parent_->Bind(site, head)};
      if (executed) {
        bool still_attached{head->parent_->UpdateExtension()};
        motors_.FlagForUpdate();
        filaments_->FlagForUpdate();
        test_stats_.at("bind_ii")[0].first++;
      }
    };
    auto weight_bind_ii = [](auto *head) {
      return head->parent_->GetWeight_Bind_II();
    };
    auto poisson_bind_ii = [&](double p, int n) {
      test_stats_.at("bind_ii")[0].second +=
          motors_.sorted_.at("bind_ii").size_;
      if (p > 0.0) {
        return SysRNG::SamplePoisson(p);
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
      bool executed{head->Unbind()};
      if (executed) {
        motors_.FlagForUpdate();
        filaments_->FlagForUpdate();
        test_stats_.at("unbind_ii")[0].first++;
      }
    };
    auto poisson_unbind_ii = [&](double p, int n) {
      test_stats_.at("unbind_ii")[0].second +=
          motors_.sorted_.at("unbind_ii").size_;
      // printf("sz = %zu\n", motors_.sorted_.at("unbind_ii").size_);
      if (p > 0.0) {
        return SysRNG::SamplePoisson(p);
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
    };
    auto poisson_unbind_i = [&](double p, int n) {
      test_stats_.at("unbind_i")[0].second +=
          motors_.sorted_.at("bound_i_ADPP").size_;
      if (p > 0.0) {
        return SysRNG::SamplePoisson(p);
      } else {
        return 0;
      }
    };
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
        &motors_.sorted_.at("bound_i_ADPP").entries_, poisson_unbind_i,
        [&](Object *base) {
          return weight_unbind_i(dynamic_cast<CatalyticHead *>(base));
        },
        exe_unbind_i);
  } else if (Sys::test_mode_ == "filament_separation") {
    // Poisson distribution; sampled to predict events w/ variable probabilities
    auto poisson = [&](double p, int n) {
      // printf("TO: %zu | FR: %zu\n",
      //        xlinks_.sorted_.at("diffuse_ii_to_rest").size_,
      //        xlinks_.sorted_.at("diffuse_ii_fr_rest").size_);
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
          printf("WUT\n");
        }
        // FIXME had to move this from if statement above -- why ?
      }
      filaments_->FlagForUpdate();
      xlinks_.FlagForUpdate();
    };
    auto exe_diffuse_bck = [&](Object *base) {
      auto head{dynamic_cast<BindingHead *>(base)};
      bool executed{head->Diffuse(-1)};
      if (executed) {
        bool still_attached{head->parent_->UpdateExtension()};
        if (!still_attached) {
        }
        // FIXME had to move this from if statement above -- why ?
      }
      filaments_->FlagForUpdate();
      xlinks_.FlagForUpdate();
    };
    auto get_weight_diff_ii_to = [](Object *base) {
      // printf("HI\n");
      auto head{dynamic_cast<BindingHead *>(base)};
      return head->GetWeight_Diffuse(1);
    };
    auto get_weight_diff_ii_fr = [](Object *base) {
      // printf("HII\n");
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
  } else if (Sys::test_mode_ == "filament_ablation") {
    InitializeEvents();
  } else if (Sys::test_mode_ == "hetero_tubulin") {
    InitializeEvents();
  } else if (Sys::test_mode_ == "kinesin_mutant") {
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
        return;
      }
      auto entry{pop->GetFreeEntry()};
      // always bind catalytic head first
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
    // need two events: one to bind catalytic head; one to bind passive

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
    auto is_ADPP_ii_bound = [](auto *motor) -> Vec<Object *> {
      if (motor->n_heads_active_ == 2) {
        if (motor->head_one_.ligand_ == CatalyticHead::Ligand::ADPP and
            motor->head_two_.ligand_ == CatalyticHead::Ligand::ADPP) {
          return {&motor->head_one_};
        }
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
    // Unbind_I: Unbind first (singly bound) head of a protein
    auto exe_unbind_i = [](auto *head, auto *pop, auto *fil) {
      bool executed{head->Unbind()};
      if (executed) {
        head->UntetherSatellite();
        pop->RemoveFromActive(head->parent_);
        fil->FlagForUpdate();
      }
    };
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
    // Bind_ATP
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
    // Diffusion
    auto exe_diff = [](auto *head, auto *pop, auto *fil, int dir) {
      bool executed{head->Diffuse(dir)};
      if (executed) {
        pop->FlagForUpdate();
        fil->FlagForUpdate();
      }
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
          "diffuse_i_fwd",
          xlinks_.p_event_.at("diffuse_i_fwd").GetVal(n_neighbs),
          &motors_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
          &motors_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs],
          binomial, [&](Object *base) {
            exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_,
                     1);
          });
      kmc_.events_.emplace_back(
          "diffuse_i_bck",
          xlinks_.p_event_.at("diffuse_i_bck").GetVal(n_neighbs),
          &motors_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
          &motors_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs],
          binomial, [&](Object *base) {
            exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_,
                     -1);
          });
    }
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
    Sys::Log(3, "hello\n");
    auto entry{pop->GetFreeEntry()};
    if (entry == nullptr) {
      return;
    }
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
    auto exe_bind_ATP_ii = [](auto *front_head, auto *pop, auto *fil) {
      auto *rear_head{front_head->GetOtherHead()};
      if (front_head->trailing_) {
        return;
      }
      bool unbound{rear_head->Unbind()};
      bool executed{front_head->parent_->Bind_ATP(front_head)};
      if (executed) {
        pop->FlagForUpdate();
        fil->FlagForUpdate();
      }
    };
    auto weight_bind_ATP_ii = [](auto *head) {
      return head->parent_->GetWeight_BindATP_II(head);
    };
    auto is_NULL_ii_bound = [](auto *motor) -> Vec<Object *> {
      if (motor->n_heads_active_ == 2) {
        bool found_head{false};
        CatalyticHead *chosen_head{nullptr};
        if (motor->head_one_.ligand_ == CatalyticHead::Ligand::NONE) {
          chosen_head = &motor->head_one_;
        }
        if (motor->head_two_.ligand_ == CatalyticHead::Ligand::NONE) {
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
        &motors_.sorted_.at("bound_ii_NULL").entries_, poisson,
        [&](Object *base) {
          return weight_bind_ATP_ii(dynamic_cast<CatalyticHead *>(base));
        },
        [&](Object *base) {
          exe_bind_ATP_ii(dynamic_cast<CatalyticHead *>(base), &motors_,
                          filaments_);
        });
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
        // printf("what\n");
      }
      pop->FlagForUpdate();
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

void ProteinManager::UpdateFilaments() {
  filaments_->UpdateUnoccupied();
  if (Sys::test_mode_.empty()) {
    return;
  }
  if (Sys::test_mode_ != "filament_ablation") {
    return;
  }
  if (Sys::i_step_ == Sys::ablation_step_) {
    filaments_->proto_[1].pos_[0] += 200.0;
    filaments_->proto_[1].ForceUpdate();
    // printf("HELLO\n");
  }
}
