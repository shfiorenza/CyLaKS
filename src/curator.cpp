#include "cylaks/curator.hpp"
#include "yaml-cpp/parser.h"
#include "yaml-cpp/yaml.h"
#include <iostream>

void Curator::CheckArgs(int argc, char *argv[]) {

  if (argc < 3 or argc > 6) {
    printf("\nError! Incorrect number of command-line arguments\n");
    printf("Correct format: %s parameters.yaml sim_name (required) ", argv[0]);
    printf("test_mode (optional)\n");
    printf("Currently-implemented test modes are:\n");
    for (auto const &mode : test_modes_) {
      printf("   %s\n", mode.c_str());
    }
    exit(1);
  }
  Sys::yaml_file_ = argv[1];
  Sys::sim_name_ = argv[2];
  if (argc == 4 or argc == 5 or argc == 6) {
    Sys::test_mode_ = argv[3];
    bool valid_mode{false};
    for (auto const &mode : test_modes_) {
      if (Sys::test_mode_ == mode) {
        if (argc == 5) {
          if (Sys::test_mode_ == "filament_separation") {
            Sys::n_xlinks_ = std::stoi(argv[4]);
          } else {
            break;
          }
        } else if (argc == 6) {
          if (Sys::test_mode_ == "hetero_tubulin") {
            Sys::p_mutant_ = std::stod(argv[4]);
            Sys::binding_affinity_ = std::stod(argv[5]);
          } else {
            break;
          }
        }
        valid_mode = true;
      }
    }
    if (!valid_mode) {
      printf("\nError! Invalid test mode.\n");
      printf("Currently-implemented test modes are:\n");
      for (auto const &mode : test_modes_) {
        printf("   %s\n", mode.c_str());
      }
      exit(1);
    }
  }
}

void Curator::GenerateLog() {

  char log_name[256];
  sprintf(log_name, "%s.log", Sys::sim_name_.c_str());
  // Check to see if sim files already exist
  if (std::filesystem::exists(log_name)) {
    printf("Log file with this name already exists!\n");
    printf("Do you wish to overwrite these data? y/n\n");
    Str response;
    bool response_unacceptable{true};
    size_t n_responses{0};
    while (response_unacceptable) {
      std::getline(std::cin, response);
      if (response == "n" or response == "N") {
        printf("Simulation terminated.\n");
        exit(1);
      } else if (response == "y" or response == "Y") {
        printf("Very well. Overwriting data for ");
        printf("simulation '%s'\n\n", Sys::sim_name_.c_str());
        response_unacceptable = false;
      } else {
        printf("Invalid response. Please choose yes or no.\n");
      }
      n_responses++;
      if (n_responses > 5) {
        printf("Too many incorrect inputs; terminating. ");
        exit(1);
      }
    }
  }
  Sys::log_file_ = fopen(log_name, "w");
  if (Sys::log_file_ == nullptr) {
    printf("Error; cannot open log file '%s'\n", log_name);
    exit(1);
  }
  // Daisy-chain some functs to get current date/time in a formatted string
  auto now{std::chrono::system_clock::now()};
  std::time_t now_c{std::chrono::system_clock::to_time_t(now)};
  std::tm now_tm{*std::localtime(&now_c)};
  char now_str[256];
  strftime(now_str, sizeof now_str, "%c", &now_tm);
  fprintf(Sys::log_file_, "[Log file auto-generated for simulation");
  fprintf(Sys::log_file_, " '%s' on %s]\n\n", Sys::sim_name_.c_str(), now_str);
}

void Curator::ParseParameters() {

  using namespace Sys;
  // Check to make sure param file actually exists
  if (!std::filesystem::exists(yaml_file_)) {
    Log("  Error: param file does not exist; aborting \n");
    exit(1);
  }
  // Open parameter file
  Log("Reading parameters from file '%s':\n", yaml_file_.c_str());
  YAML::Node input{YAML::LoadFile(yaml_file_)};
  // Construct function to get values from yaml file and log them
  auto ParseYAML = [&]<typename DATA_T>(DATA_T *param, Str name, Str units) {
    // Unparsed parameter value from input yaml node
    YAML::Node val;
    // The character "." is used to separate group and value labels
    // If present in the param name, split it into group & name strings
    if (name.find(".") < name.length()) {
      Str group{name.substr(0, name.find("."))};
      name = name.substr(name.find(".") + 1, name.length());
      val = input[group][name];
    }
    // Otherwise, look up value label directly
    else {
      val = input[name];
    }
    // Sometimes, size_t variables use "e" notation & must be treated as doubles
    try {
      *param = val.as<DATA_T>();
    } catch (const YAML::BadConversion err) {
      *param = (DATA_T)val.as<double>();
    }
    // Convert parameter value into a string that we can easily log
    // (Since we don't know data types a priori, we can't use Log() aka printf)
    Str val_str;
    if (val.Type() == YAML::NodeType::Scalar) {
      Log("   %s = %s %s\n", name.c_str(), val.as<Str>().c_str(),
          units.c_str());
    } else if (val.Type() == YAML::NodeType::Sequence) {
      for (int i_val{0}; i_val < val.size(); i_val++) {
        Log("   %s[%i] = %s %s\n", name.c_str(), i_val,
            val[i_val].as<Str>().c_str(), units.c_str());
      }
    } else {
      ErrorExit("Curator::ParseParameters() -- parser");
    }
  };
  // Parse thru parameters
  using namespace Params;
  Log(" General parameters:\n");
  ParseYAML(&seed, "seed", "");
  ParseYAML(&kbT, "kbT", "pN*nm");
  ParseYAML(&eta, "eta", "pN*s/um^2");
  ParseYAML(&dt, "dt", "s");
  ParseYAML(&t_run, "t_run", "s");
  ParseYAML(&t_equil, "t_equil", "s");
  ParseYAML(&t_snapshot, "t_snapshot", "s");
  ParseYAML(&dynamic_equil_window, "dynamic_equil_window", "s");
  ParseYAML(&verbosity, "verbosity", "");
  Log(" Filament parameters:\n");
  ParseYAML(&Filaments::count, "filaments.count", "filaments");
  ParseYAML(&Filaments::radius, "filaments.radius", "nm");
  ParseYAML(&Filaments::site_size, "filaments.site_size", "nm");
  ParseYAML(&Filaments::n_bd_per_kmc, "filaments.n_bd_per_kmc", "");
  // ParseYAML(&Filaments::t_ablate, "filaments.t_ablate", "s");
  ParseYAML(&Filaments::n_sites, "filaments.n_sites", "sites");
  ParseYAML(&Filaments::polarity, "filaments.polarity", "");
  ParseYAML(&Filaments::x_initial, "filaments.x_initial", "nm");
  ParseYAML(&Filaments::y_initial, "filaments.y_initial", "nm");
  ParseYAML(&Filaments::f_applied, "filaments.f_applied", "pN");
  ParseYAML(&Filaments::immobile_until, "filaments.immobile_until", "s");
  ParseYAML(&Filaments::translation_enabled, "filaments.translation_enabled",
            "");
  ParseYAML(&Filaments::rotation_enabled, "filaments.rotation_enabled", "");
  // Check to make sure there are enough vector entries for given MT count
  if (Filaments::count > Filaments::n_sites.size() or
      Filaments::count > Filaments::polarity.size() or
      Filaments::count > Filaments::x_initial.size() or
      Filaments::count > Filaments::y_initial.size() or
      Filaments::count > Filaments::immobile_until.size() or
      _n_dims_max > Filaments::translation_enabled.size()) {
    Log("Error! Incorrect number of filament parameters provided.\n");
    exit(1);
  }
  for (int i_fil{0}; i_fil < Filaments::count; i_fil++) {
    if (Filaments::polarity[i_fil] != 0 and Filaments::polarity[i_fil] != 1) {
      Log("Error! Polarity must be either 0 or 1.\n");
      exit(1);
    }
  }
  Log(" Kinesin (motor) parameters:\n");
  ParseYAML(&Motors::n_runs_to_exit, "motors.n_runs_to_exit", "runs");
  ParseYAML(&Motors::gaussian_range, "motors.gaussian_range", "sites");
  ParseYAML(&Motors::gaussian_amp_solo, "motors.gaussian_amp_solo", "kbT");
  ParseYAML(&Motors::gaussian_ceiling_bulk, "motors.gaussian_ceiling_bulk",
            "kbT");
  ParseYAML(&Motors::neighb_neighb_energy, "motors.neighb_neighb_energy",
            "kbT");
  ParseYAML(&Motors::t_active, "motors.t_active", "s");
  ParseYAML(&Motors::k_on, "motors.k_on", "1/nM*s");
  ParseYAML(&Motors::c_bulk, "motors.c_bulk", "nM");
  ParseYAML(&Motors::c_eff_bind, "motors.c_eff_bind", "nM");
  ParseYAML(&Motors::k_on_ATP, "motors.k_on_ATP", "1/s");
  ParseYAML(&Motors::c_ATP, "motors.c_ATP", "mM");
  ParseYAML(&Motors::k_hydrolyze, "motors.k_hydrolyze", "1/s");
  ParseYAML(&Motors::k_off_i, "motors.k_off_i", "1/s");
  ParseYAML(&Motors::k_off_ii, "motors.k_off_ii", "1/s");
  ParseYAML(&Motors::applied_force, "motors.applied_force", "pN");
  ParseYAML(&Motors::internal_force, "motors.internal_force", "pN");
  ParseYAML(&Motors::sigma_off_i, "motors.sigma_off_i", "nm");
  ParseYAML(&Motors::sigma_off_ii, "motors.sigma_off_ii", "nm");
  ParseYAML(&Motors::sigma_ATP, "motors.sigma_ATP", "nm");
  ParseYAML(&Motors::endpausing_active, "motors.endpausing_active", "");
  ParseYAML(&Motors::tethers_active, "motors.tethers_active", "");
  ParseYAML(&Motors::k_tether, "motors.k_tether", "1/nM*s");
  ParseYAML(&Motors::c_eff_tether, "motors.c_eff_tether", "nM");
  ParseYAML(&Motors::k_untether, "motors.k_untether", "1/s");
  ParseYAML(&Motors::r_0, "motors.r_0", "nm");
  ParseYAML(&Motors::k_spring, "motors.k_spring", "pN/nm");
  ParseYAML(&Motors::k_slack, "motors.k_slack", "pN/nm");
  Sys::Log("  Crosslinker (xlink) parameters:\n");
  ParseYAML(&Xlinks::t_active, "xlinks.t_active", "s");
  ParseYAML(&Xlinks::k_on, "xlinks.k_on", "1/nM*s");
  ParseYAML(&Xlinks::c_bulk, "xlinks.c_bulk", "nM");
  ParseYAML(&Xlinks::c_eff_bind, "xlinks.c_eff_bind", "nM");
  ParseYAML(&Xlinks::k_off_i, "xlinks.k_off_i", "1/s");
  ParseYAML(&Xlinks::k_off_ii, "xlinks.k_off_ii", "1/s");
  ParseYAML(&Xlinks::d_i, "xlinks.d_i", "um^2/s");
  ParseYAML(&Xlinks::d_ii, "xlinks.d_ii", "um^2/s");
  ParseYAML(&Xlinks::r_0, "xlinks.r_0", "nm");
  ParseYAML(&Xlinks::k_spring, "xlinks.k_spring", "pN/nm");
  ParseYAML(&Xlinks::theta_0, "xlinks.theta_0", "degrees");
  ParseYAML(&Xlinks::k_rot, "xlinks.k_rot", "pN*nm/rad");
}

void Curator::InitializeSimulation() {

  using namespace Params;
  using namespace Sys;
  // Calculate local parameters
  start_time_ = SysClock::now();
  n_steps_per_snapshot_ = (size_t)std::round(t_snapshot / dt);
  // Calculate system parameters
  verbosity_ = verbosity;
  n_steps_pre_equil_ = (size_t)std::round(t_equil / dt);
  n_steps_equil_ = n_steps_pre_equil_;
  if (n_steps_pre_equil_ == 0) {
    equilibrating_ = false;
  }
  n_steps_run_ = (size_t)std::round(t_run / dt);
  // Log parameters
  Log("\n");
  Log("  System variables calculated post-initialization:\n");
  Log("   n_steps_run = %zu\n", n_steps_run_);
  Log("   n_steps_equil = %zu\n", n_steps_equil_);
  Log("   n_steps_per_snapshot = %zu\n", n_steps_per_snapshot_);
  Log("   n_datapoints = %zu\n", n_steps_run_ / n_steps_per_snapshot_);
  Log("\n");
  // Initialize sim objects
  SysRNG::Initialize(seed);
  // If we're running a test, let proteins initialize filament environment
  if (Sys::test_mode_.empty()) {
    filaments_.Initialize(&proteins_);
  }
  proteins_.Initialize(&filaments_);
  for (auto const &pf : filaments_.proto_) {
    if (pf.sites_.size() > n_sites_max_) {
      n_sites_max_ = pf.sites_.size();
    }
  }
  if (Sys::test_mode_.empty()) {
    Log("\n");
  }
}

void Curator::GenerateDataFiles() {

  auto AddDataFile = [&](Str name) {
    data_files_.emplace(name, DataFile(name));
  };
  // Open filament pos file, which stores the N-dim coordinates of the two
  // endpoints of each filament every datapoint
  AddDataFile("filament_pos");
  // if (Params::Filaments::t_ablate > 0.0) {
  //   AddDataFile("filament_pos_postSplit");
  // }
  if (proteins_.motors_.active_ or proteins_.xlinks_.active_) {
    // Open occupancy file, which stores the species ID of each occupant
    // (or -1 for none) for all MT sites during data collection (DC)
    AddDataFile("occupancy");
    // if (Params::Filaments::t_ablate > 0.0) {
    //   AddDataFile("occupancy_postSplit");
    // }
    // Open protein ID file, which stores the unique ID of all bound proteins
    // (unbound not tracked) at their respective site indices during DC
    AddDataFile("protein_id");
    // if (Params::Filaments::t_ablate > 0.0) {
    //   AddDataFile("protein_id_postSplit");
    // }
    if (proteins_.xlinks_.crosslinking_active_) {
      AddDataFile("partner_index");
    }
  }
  if (proteins_.motors_.active_) {
    // bool; simply says if motor head is trailing or not
    AddDataFile("motor_head_trailing");
    // if (Params::Filaments::t_ablate > 0.0) {
    //   AddDataFile("motor_head_trailing_postSplit");
    // }
    if (proteins_.motors_.tethering_active_) {
      // Open tether coord file, which stores the coordinates
      // of the anchor points of tethered motors
      AddDataFile("tether_anchor_pos");
      /*
      // Open motor extension file, which stores the number of motors
      // with a certain tether extension for all possible extensions
      AddDataFile("motor_ext");
      // Open motor force file, which stores the sum
      // of forces coming from motor tether extensions
      AddDataFile("motor_force");
      */
    }
  }
  if (proteins_.xlinks_.active_ and proteins_.xlinks_.crosslinking_active_) {
    /*
    // Open xlink extension file, which stores the number of stage-2
    // xlinks at a certain extension for all possible extensions
    AddDataFile("xlink_ext");
    // Open xlink force file, which stores the sum
    // of forces coming from xlink extensions
    AddDataFile("xlink_force");
    if (proteins_.motors_.tethering_active_) {
      // Open total force file, which stores the sum of ALL
      // forces coming from xlink and motor tether extensions
      AddDataFile("total_force");
    }
    */
  }
}

void Curator::CheckPrintProgress() {

  // FIXME report t_sim & t_elapsed_irl each milestone; not step #
  using namespace Sys;
  using namespace Params;
  // Percent milestone; controls report frequency
  int p_report{10};
  // Advance simulation forward one step (or 1 dt in real time)
  i_step_++;
  // If still equilibrating, report progress and check protein equil. status
  if (equilibrating_) {
    if (i_step_ == 1) {
      Log("Pre-equilibration is 0%% complete.\n");
    }
    if (i_step_ % (n_steps_pre_equil_ / (100 / p_report)) == 0 and
        i_step_ <= n_steps_pre_equil_) {
      Log("Pre-equilibration is %g%% complete. (step #%zu | t = %g s)\n",
          double(i_step_) / n_steps_pre_equil_ * 100, i_step_, i_step_ * dt);
    }
    if (proteins_.motors_.equilibrated_ and proteins_.xlinks_.equilibrated_ and
        i_step_ >= n_steps_pre_equil_) {
      n_steps_equil_ = i_step_;
      equilibrating_ = false;
      if (dynamic_equil_window > 0.0) {
        Log("Dynamic equilibration is complete. (t = %g s)\n", i_step_ * dt);
        Log("   N_STEPS_EQUIL = %zu\n", n_steps_equil_);
      }
    }
  }
  // Otherwise, data collection must be active; simply report on that
  else {
    size_t n_steps_so_far{i_step_ - n_steps_equil_};
    if (n_steps_so_far == 1) {
      Log("Data collection is 0%% complete.\n");
    }
    if (n_steps_so_far % (n_steps_run_ / (100 / p_report)) == 0) {
      Log("Data collection is %g%% complete. (step #%zu | t = %g s)\n",
          double(n_steps_so_far) / n_steps_run_ * 100, i_step_, i_step_ * dt);
    }
  }
  // Terminate simulation once a sufficient number of steps has been taken
  if (i_step_ >= n_steps_run_ + n_steps_equil_) {
    running_ = false;
    long clock_ticks{(SysClock::now() - start_time_).count()};
    size_t ticks_per_second{SysClock::period::den};
    Log("Simulation complete. Total time to execute: %.2f s\n",
        double(clock_ticks) / ticks_per_second);
  }
}

void Curator::OutputData() {

  if (Sys::equilibrating_ or Sys::i_step_ % n_steps_per_snapshot_ != 0) {
    return;
  }
  Sys::i_datapoint_++;
  for (auto &&pf : filaments_.proto_) {
    double coord1[_n_dims_max];
    double coord2[_n_dims_max];
    for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
      coord1[i_dim] = (double)pf.plus_end_->pos_[i_dim];
      coord2[i_dim] = (double)pf.minus_end_->pos_[i_dim];
    }
    // printf("wrote plus_end = (%g, %g)\n", coord1[0], coord1[1]);
    // printf("wrote minus_end = (%g, %g)\n", coord2[0], coord2[1]);
    data_files_.at("filament_pos").Write(coord1, _n_dims_max);
    data_files_.at("filament_pos").Write(coord2, _n_dims_max);
    if (!proteins_.motors_.active_ and !proteins_.xlinks_.active_) {
      continue;
    }
    int occupancy[n_sites_max_];
    int protein_id[n_sites_max_];
    int partner_index[n_sites_max_];
    bool motor_trailing[n_sites_max_];
    double tether_anchor_pos[n_sites_max_];
    for (int i_site{0}; i_site < n_sites_max_; i_site++) {
      occupancy[i_site] = _id_site;
      protein_id[i_site] = -1;
      partner_index[i_site] = -1;
      motor_trailing[i_site] = false;
      tether_anchor_pos[i_site] = -1.0;
    }
    for (auto const &site : pf.sites_) {
      if (site.occupant_ == nullptr) {
        continue;
      }
      const size_t species_id{site.occupant_->GetSpeciesID()};
      occupancy[site.index_] = species_id;
      protein_id[site.index_] = site.occupant_->GetID();
      if (species_id == _id_xlink) {
        if (site.occupant_->parent_->n_heads_active_ == 2) {
          partner_index[site.index_] =
              site.occupant_->GetOtherHead()->site_->index_;
        }
      } else if (species_id == _id_motor) {
        motor_trailing[site.index_] = site.occupant_->Trailing();
        /*
        if (site.occupant_->parent_->tethered_) {
          auto partner{site.occupant_->parent_->partner_};
          if (partner->n_heads_active_ > 0) {
            // double anchor_coord{partner->GetAnchorCoordinate()};
            // tether_anchor_pos[site.index_] = anchor_coord;
          }
        }
        */
      }
    }
    data_files_.at("occupancy").Write(occupancy, n_sites_max_);
    data_files_.at("protein_id").Write(protein_id, n_sites_max_);
    if (proteins_.xlinks_.crosslinking_active_) {
      data_files_.at("partner_index").Write(partner_index, n_sites_max_);
    }
    if (!proteins_.motors_.active_) {
      continue;
    }
    data_files_.at("motor_head_trailing").Write(motor_trailing, n_sites_max_);
    if (!proteins_.motors_.tethering_active_) {
      continue;
    }
    data_files_.at("tether_anchor_pos").Write(tether_anchor_pos, n_sites_max_);
  }
  if (Sys::early_exit_triggered_) {
    Sys::EarlyExit();
  }
}
