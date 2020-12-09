#include "curator.hpp"
#include "yaml-cpp/yaml.h"
#include <iostream>

void Curator::CheckArgs(int argc, char *argv[]) {

  if (argc < 3 or argc > 4) {
    printf("\nError! Incorrect number of command-line arguments\n");
    printf("Correct format: %s parameters.yaml sim_name (required) ", argv[0]);
    printf("test_mode (optional)\n");
    printf("Currently-implemented test modes are:\n");
    printf("    NONE!\n");
    // printf("    xlink_bind_ii\n");
    // printf("    motor_lattice_bind\n");
    // printf("    motor_lattice_step\n");
    exit(1);
  }
  Sys::yaml_file_ = argv[1];
  Sys::name_ = argv[2];
  if (argc == 4) {
    Sys::test_mode_ = argv[3];
  }
}

void Curator::GenerateLog() {

  char log_name[256];
  sprintf(log_name, "%s.log", Sys::name_.c_str());
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
        printf("Very well. ");
        printf("Overwriting data for simulation '%s'\n\n", Sys::name_.c_str());
        response_unacceptable = false;
      } else {
        printf("ayo I said y or n. try again plz\n");
      }
      n_responses++;
      if (n_responses > 5) {
        printf("aight if u gonna be like that,");
        printf("let's just cancel this whole thing\n");
        exit(1);
      }
    }
  }
  Sys::log_file_ = fopen(log_name, "w");
  if (Sys::log_file_ == nullptr) {
    printf("Error; cannot open log file '%s'\n", log_name);
    exit(1);
  }
  // Daisy-chain c/c++ functs to get current date/time in a formatted string
  auto now{std::chrono::system_clock::now()};
  std::time_t now_c{std::chrono::system_clock::to_time_t(now)};
  std::tm now_tm{*std::localtime(&now_c)};
  char now_str[256];
  strftime(now_str, sizeof now_str, "%c", &now_tm);
  fprintf(Sys::log_file_, "[Log file auto-generated for simulation");
  fprintf(Sys::log_file_, " '%s' on %s]\n\n", Sys::name_.c_str(), now_str);
}

void Curator::ParseParameters() {

  // Check to make sure param file actually exists
  if (!std::filesystem::exists(Sys::yaml_file_)) {
    Sys::Log("  Error: param file does not exist; aborting \n");
    exit(1);
  }
  // Parse parameter file into a YAML node
  YAML::Node input = YAML::LoadFile(Sys::yaml_file_);
  // Transfer values from input param node to system_parameters structure
  try {
    params_.seed = input["seed"].as<size_t>();
  } catch (const YAML::BadConversion error) {
    params_.seed = (size_t)(input["seed"].as<double>());
  }
  params_.kbT = input["kbT"].as<double>();
  params_.eta = input["eta"].as<double>();
  params_.dt = input["dt"].as<double>();
  params_.t_run = input["t_run"].as<double>();
  params_.t_equil = input["t_equil"].as<double>();
  params_.t_snapshot = input["t_snapshot"].as<double>();
  params_.dynamic_equil = input["dynamic_equil"].as<bool>();
  params_.dynamic_equil_window = input["dynamic_equil_window"].as<double>();
  params_.verbosity = input["verbosity"].as<size_t>();
  /* Microtubule parameters below */
  YAML::Node mts = input["filaments"];
  params_.filaments.count = mts["count"].as<size_t>();
  params_.filaments.length = mts["length"].as<Vec<size_t>>();
  params_.filaments.start_coord = mts["start_coord"].as<Vec<double>>();
  params_.filaments.immobile_until = mts["immobile_until"].as<Vec<double>>();
  params_.filaments.applied_force = mts["applied_force"].as<Vec<double>>();
  params_.filaments.dim_enabled = mts["dim_enabled"].as<Vec<bool>>();
  params_.filaments.y_dist_init = mts["y_dist_init"].as<double>();
  params_.filaments.site_size = mts["site_size"].as<double>();
  params_.filaments.radius = mts["radius"].as<double>();
  /* Motor parameters below */
  YAML::Node motors = input["motors"];
  params_.motors.n_runs_desired = motors["n_runs_desired"].as<size_t>();
  try {
    params_.motors.lattice_coop_range =
        motors["lattice_coop_range"].as<size_t>();
  } catch (const YAML::BadConversion error) {
    params_.motors.lattice_coop_range =
        (size_t)std::round(motors["lattice_coop_range"].as<double>());
  }
  params_.motors.lattice_coop_Emax_solo =
      motors["lattice_coop_Emax_solo"].as<double>();
  params_.motors.lattice_coop_Emax_bulk =
      motors["lattice_coop_Emax_bulk"].as<double>();
  params_.motors.interaction_energy = motors["interaction_energy"].as<double>();
  params_.motors.t_active = motors["t_active"].as<double>();
  params_.motors.k_on = motors["k_on"].as<double>();
  params_.motors.c_bulk = motors["c_bulk"].as<double>();
  params_.motors.c_eff_bind = motors["c_eff_bind"].as<double>();
  params_.motors.k_on_ATP = motors["k_on_ATP"].as<double>();
  params_.motors.c_ATP = motors["c_ATP"].as<double>();
  params_.motors.k_hydrolyze = motors["k_hydrolyze"].as<double>();
  params_.motors.k_off_i = motors["k_off_i"].as<double>();
  params_.motors.k_off_ii = motors["k_off_ii"].as<double>();
  params_.motors.applied_force = motors["applied_force"].as<double>();
  params_.motors.internal_force = motors["internal_force"].as<double>();
  params_.motors.sigma_off_i = motors["sigma_off_i"].as<double>();
  params_.motors.sigma_off_ii = motors["sigma_off_ii"].as<double>();
  params_.motors.sigma_ATP = motors["sigma_ATP"].as<double>();
  params_.motors.k_tether = motors["k_tether"].as<double>();
  params_.motors.c_eff_tether = motors["c_eff_tether"].as<double>();
  params_.motors.k_untether = motors["k_untether"].as<double>();
  params_.motors.r_0 = motors["r_0"].as<double>();
  params_.motors.k_spring = motors["k_spring"].as<double>();
  params_.motors.k_slack = motors["k_slack"].as<double>();
  params_.motors.tethers_active = motors["tethers_active"].as<bool>();
  params_.motors.endpausing_active = motors["endpausing_active"].as<bool>();
  /* Xlink parameters below */
  YAML::Node xlinks = input["xlinks"];
  params_.xlinks.interaction_energy = xlinks["interaction_energy"].as<double>();
  params_.motors.t_active = xlinks["t_active"].as<double>();
  params_.xlinks.k_on = xlinks["k_on"].as<double>();
  params_.xlinks.c_bulk = xlinks["c_bulk"].as<double>();
  params_.xlinks.c_eff_bind = xlinks["c_eff_bind"].as<double>();
  params_.xlinks.k_off_i = xlinks["k_off_i"].as<double>();
  params_.xlinks.k_off_ii = xlinks["k_off_ii"].as<double>();
  params_.xlinks.d_i = xlinks["d_i"].as<double>();
  params_.xlinks.d_ii = xlinks["d_ii"].as<double>();
  params_.xlinks.r_0 = xlinks["r_0"].as<double>();
  params_.xlinks.k_spring = xlinks["k_spring"].as<double>();
  // Store params pointer as parameters_ in Curator
  Sys::Log("Reading parameters from '%s':\n", Sys::yaml_file_);
  Sys::Log("  General parameters:\n");
  Sys::Log("    seed = %lu\n", params_.seed);
  Sys::Log("    kbT = %g pN*nm\n", params_.kbT);
  Sys::Log("    eta = %g (pN*s)/um^2\n", params_.eta);
  Sys::Log("    dt = %g s\n", params_.dt);
  Sys::Log("    t_run = %g s\n", params_.t_run);
  Sys::Log("    t_equil = %g s\n", params_.t_equil);
  Sys::Log("    t_snapshot = %g s\n", params_.t_snapshot);
  Sys::Log("    dynamic_equil = %s\n",
           params_.dynamic_equil ? "true" : "false");
  Sys::Log("    dynamic_equil_window = %g s\n", params_.dynamic_equil_window);
  Sys::Log("    verbosity = %zu\n", params_.verbosity);
  Sys::Log("  Filament parameters:\n");
  Sys::Log("    count = %i\n", params_.filaments.count);
  // Check to make sure there are enough vector entries for given MT count
  if (params_.filaments.count > params_.filaments.length.size() or
      params_.filaments.count > params_.filaments.start_coord.size() or
      params_.filaments.count > params_.filaments.applied_force.size() or
      params_.filaments.count > params_.filaments.immobile_until.size() or
      _n_dims_max > params_.filaments.dim_enabled.size()) {
    Sys::ErrorExit("Curator::ParseParameters()");
  }
  for (int i_fil{0}; i_fil < params_.filaments.count; i_fil++) {
    Sys::Log("    length[%i] = %i sites\n", i_fil,
             params_.filaments.length[i_fil]);
    Sys::Log("    start_coord[%i] = %g sites\n", i_fil,
             params_.filaments.start_coord[i_fil]);
    Sys::Log("    applied_force[%i] = %g pN\n", i_fil,
             params_.filaments.applied_force[i_fil]);
    Sys::Log("    immobile_until[%i] = %g s\n", i_fil,
             params_.filaments.immobile_until[i_fil]);
  }
  for (int i_dim{0}; i_dim < _n_dims_max; i_dim++) {
    Sys::Log("    dim_enabled[%i] = %s\n", i_dim,
             params_.filaments.dim_enabled[i_dim] ? "true" : "false");
  }
  Sys::Log("    y_dist_init = %g nm between MTs\n",
           params_.filaments.y_dist_init);
  Sys::Log("    site_size = %g nm\n", params_.filaments.site_size);
  Sys::Log("    radius = %g nm\n", params_.filaments.radius);
  Sys::Log("  Kinesin (motor) parameters:\n");
  Sys::Log("    n_runs_desired = %zu\n", params_.motors.n_runs_desired);
  Sys::Log("    lattice_coop_range = %i\n", params_.motors.lattice_coop_range);
  Sys::Log("    lattice_coop_Emax_solo = -%g kbT\n",
           params_.motors.lattice_coop_Emax_solo);
  Sys::Log("    lattice_coop_Emax_bulk = -%g kbT\n",
           params_.motors.lattice_coop_Emax_bulk);
  Sys::Log("    interaction_energy = -%g kbT\n",
           params_.motors.interaction_energy);
  Sys::Log("    t_active = %g seconds\n", params_.motors.t_active);
  Sys::Log("    k_on = %g /(nM*s)\n", params_.motors.k_on);
  Sys::Log("    c_bulk = %g nM\n", params_.motors.c_bulk);
  Sys::Log("    c_eff_bind = %g nM\n", params_.motors.c_eff_bind);
  Sys::Log("    k_on_ATP = %g /(mM*s)\n", params_.motors.k_on_ATP);
  Sys::Log("    c_ATP = %g mM\n", params_.motors.c_ATP);
  Sys::Log("    k_hydrolyze = %g /s\n", params_.motors.k_hydrolyze);
  Sys::Log("    k_off_i = %g /s\n", params_.motors.k_off_i);
  Sys::Log("    k_off_ii = %g /s\n", params_.motors.k_off_ii);
  Sys::Log("    applied_force = %g pN\n", params_.motors.applied_force);
  Sys::Log("    internal_force = %g pN\n", params_.motors.internal_force);
  Sys::Log("    sigma_off_i = %g nm\n", params_.motors.sigma_off_i);
  Sys::Log("    sigma_off_ii = %g nm\n", params_.motors.sigma_off_ii);
  Sys::Log("    sigma_ATP = %g nm\n", params_.motors.sigma_ATP);
  Sys::Log("    k_tether = %g /(nM*s)\n", params_.motors.k_tether);
  Sys::Log("    c_eff_tether = %g nM\n", params_.motors.c_eff_tether);
  Sys::Log("    k_untether = %g /s\n", params_.motors.k_untether);
  Sys::Log("    r_0 = %g nm\n", params_.motors.r_0);
  Sys::Log("    k_spring = %g pN/nm\n", params_.motors.k_spring);
  Sys::Log("    k_slack = %g pN/nm\n", params_.motors.k_slack);
  Sys::Log("    tethers_active = %s\n",
           params_.motors.tethers_active ? "true" : "false");
  Sys::Log("    endpausing_active = %s\n",
           params_.motors.endpausing_active ? "true" : "false");
  Sys::Log("  Crosslinker (xlink) parameters:\n");
  Sys::Log("    interaction_energy = %g kbT\n",
           params_.xlinks.interaction_energy);
  Sys::Log("    t_active = %g seconds\n", params_.xlinks.t_active);
  Sys::Log("    k_on = %g /(nM*s)\n", params_.xlinks.k_on);
  Sys::Log("    c_bulk = %g nM\n", params_.xlinks.c_bulk);
  Sys::Log("    c_eff_bind = %g nM\n", params_.xlinks.c_eff_bind);
  Sys::Log("    k_off_i = %g /s\n", params_.xlinks.k_off_i);
  Sys::Log("    k_off_ii = %g /s\n", params_.xlinks.k_off_ii);
  Sys::Log("    d_i = %g um^2/s\n", params_.xlinks.d_i);
  Sys::Log("    d_ii = %g um^2/s\n", params_.xlinks.d_ii);
  Sys::Log("    r_0 = %g nm\n", params_.xlinks.r_0);
  Sys::Log("    k_spring = %g pN/nm\n", params_.xlinks.k_spring);
}

void Curator::InitializeSimulation() {

  // Calculate and/or set system parameters
  Sys::verbosity_ = params_.verbosity;
  Sys::n_steps_pre_equil_ = (size_t)std::round(params_.t_equil / params_.dt);
  Sys::n_steps_equil_ = Sys::n_steps_pre_equil_;
  Sys::n_steps_run_ = (size_t)std::round(params_.t_run / params_.dt);
  n_steps_per_snapshot_ = (size_t)std::round(params_.t_snapshot / params_.dt);
  start_time_ = SysClock::now();
  // Log system parameters
  Sys::Log("  System parameters:\n");
  Sys::Log("    n_steps_run = %zu\n", Sys::n_steps_run_);
  Sys::Log("    n_steps_equil = %zu\n", Sys::n_steps_equil_);
  Sys::Log("    n_steps_per_snapshot = %zu\n", n_steps_per_snapshot_);
  Sys::Log("    n_datapoints = %zu\n",
           Sys::n_steps_run_ / n_steps_per_snapshot_);
  Sys::Log("\n");
  // Initialize sim objects
  gsl_.Initialize(params_.seed);
  proteins_.Initialize(this, &params_);
  filaments_.Initialize(this, &params_);
}

void Curator::GenerateDataFiles() {

  // Open occupancy file, which stores the species ID of each occupant
  // (or -1 for none) for all MT sites during data collection (DC)
  AddDataFile("occupancy");
  // Motor-related files
  if (proteins_.motors_.active_) {
    // Open motor ID file, which stores the unique ID of all bound motors
    // (unbound not tracked) and their respective site indices during DC
    AddDataFile("motorID");
    // bool; simply says if motor head is trailing or not
    AddDataFile("motor_trailing");
    if (proteins_.motors_.tethering_active_) {
      // Open tether coord file, which stores the coordinates
      // of the anchor points of tethered motors
      AddDataFile("tether_coord");
      // Open motor extension file, which stores the number of motors
      // with a certain tether extension for all possible extensions
      AddDataFile("motor_dx");
      // Open motor force file, which stores the sum
      // of forces coming from motor tether extensions
      AddDataFile("motor_force");
    }
  }
  // Crosslinker-related files
  if (proteins_.xlinks_.active_) {
    // Open xlink ID file, which does the same
    // as the motor ID file but for xlinks
    AddDataFile("xlinkID");
    if (proteins_.xlinks_.crosslinking_active_) {
      // Open xlink extension file, which stores the number of stage-2
      // xlinks at a certain extension for all possible extensions
      AddDataFile("xlink_dx");
      // Open xlink force file, which stores the sum
      // of forces coming from xlink extensions
      AddDataFile("xlink_force");
    }
  }
  if (filaments_.mobile_) {
    // Open mt coord file, which stores the coordinates
    // of the left-most edge of each microtubule during DC
    AddDataFile("filament_coord");
  }
  if (proteins_.motors_.tethering_active_ and
      proteins_.xlinks_.crosslinking_active_) {
    // Open total force file, which stores the sum of ALL
    // forces coming from xlink and motor tether extensions
    AddDataFile("total_force");
  }
}

void Curator::CheckPrintProgress() {

  // Percent milestone; controls report frequency
  int p_report{10};
  double dt{params_.dt};
  // Advance simulation forward one site (or 1 dt in real time)
  Sys::i_step_++;
  // If still equilibrating, report progress and check protein equil. status
  if (Sys::equilibrating_) {
    if (Sys::i_step_ % (Sys::n_steps_pre_equil_ / (100 / p_report)) == 0 and
        Sys::i_step_ <= Sys::n_steps_pre_equil_) {
      Sys::Log("Pre-equilibration is %g%% complete. (step #%zu | t = %g s)\n",
               double(Sys::i_step_) / Sys::n_steps_pre_equil_ * 100,
               Sys::i_step_, Sys::i_step_ * dt);
    }
    if (proteins_.motors_.equilibrated_ and proteins_.xlinks_.equilibrated_ and
        Sys::i_step_ >= Sys::n_steps_pre_equil_) {
      Sys::n_steps_equil_ = Sys::i_step_;
      Sys::equilibrating_ = false;
      if (params_.dynamic_equil) {
        Sys::Log("Dynamic equilibration is complete. (t = %g s)\n",
                 Sys::i_step_ * dt);
        Sys::Log("   N_STEPS_EQUIL = %zu\n", Sys::n_steps_equil_);
      }
    }
  }
  // Otherwise data collection must be active; simply report on that
  else {
    size_t n_steps_so_far{Sys::i_step_ - Sys::n_steps_equil_};
    if (n_steps_so_far % (Sys::n_steps_run_ / (100 / p_report)) == 0) {
      Sys::Log("Data collection %g%% complete. (step #%zu | t = %g s)\n",
               double(n_steps_so_far) / Sys::n_steps_run_ * 100, Sys::i_step_,
               Sys::i_step_ * dt);
    }
  }
  // Terminate simulation once a sufficient number of steps has been taken
  if (Sys::i_step_ >= Sys::n_steps_run_ + Sys::n_steps_equil_) {
    Sys::running_ = false;
    long clock_ticks{(SysClock::now() - start_time_).count()};
    size_t ticks_per_second{SysClock::period::den};
    Sys::Log("Simulation complete. Total time to execute: %.2f s.\n",
             double(clock_ticks) / ticks_per_second);
  }
}

void Curator::OutputData() {

  if (Sys::equilibrating_ or Sys::i_step_ % n_steps_per_snapshot_ != 0) {
    return;
  }
  size_t n_pfs{params_.filaments.count};
  size_t max_length{0};
  for (int i_fil{0}; i_fil < n_pfs; i_fil++) {
    if (params_.filaments.length[i_fil] > max_length)
      max_length = params_.filaments.length[i_fil];
  }
  // Create arrays to store data; ptrs to write it to file
  int pf_coords[n_pfs];
  // Run through all MTs and get data for each
  for (int i_fil{0}; i_fil < n_pfs; i_fil++) {
    Protofilament *pf{&filaments_.list_[i_fil]};
    // Create arrays & ptrs for intraMT data
    int occupancy[max_length];
    int motor_IDs[max_length];
    int xlink_IDs[max_length];
    double tether_coords[max_length];
    bool motor_head_status[max_length];
    for (int i_site{0}; i_site < max_length; i_site++) {
      motor_head_status[i_site] = false;
    }
    // Run through all sites on this particular MT
    for (int i_site{0}; i_site < pf->n_sites_; i_site++) {
      BindingSite *site{&pf->sites_[i_site]};
      // If unoccupied, store the speciesID of tubulin to occupancy
      // and an ID of -1 (null) to motor/xlink ID files
      if (site->occupant_ == nullptr) {
        occupancy[i_site] = site->GetSpeciesID();
        motor_IDs[i_site] = -1;
        xlink_IDs[i_site] = -1;
        tether_coords[i_site] = -1;
      } else {
        occupancy[i_site] = site->occupant_->GetSpeciesID();
        switch (site->occupant_->GetSpeciesID()) {
        // If occupied by motor, store its species ID to occupancy_file,
        // its unique ID to the motor ID file, and -1 to xlink ID file
        case _id_motor:
          motor_IDs[i_site] = site->occupant_->GetID();
          xlink_IDs[i_site] = -1;
          tether_coords[i_site] = -1;
          motor_head_status[i_site] = site->HeadTrailing();
          if (site->occupant_->parent_->tethered_) {
            auto partner{site->occupant_->parent_->partner_};
            if (partner->n_heads_active_ > 0) {
              double anchor_coord{partner->GetAnchorCoordinate()};
              tether_coords[i_site] = anchor_coord;
            }
          }
          break;
        // If occupied by xlink, store its species ID to occupancy_file,
        // its unique ID to the xlink ID file, and -1 to motor ID file
        case _id_xlink:
          motor_IDs[i_site] = -1;
          xlink_IDs[i_site] = site->occupant_->GetID();
          tether_coords[i_site] = -1;
          break;
        }
      }
    }
    // Pad written data with NULL entries (-1) for shorter MTs
    for (size_t i_site{pf->n_sites_}; i_site < max_length; i_site++) {
      occupancy[i_site] = -1;
      motor_IDs[i_site] = -1;
      xlink_IDs[i_site] = -1;
      tether_coords[i_site] = -1;
    }
    pf_coords[i_fil] = pf->GetPos(0);
    /*
    // Write the data to respective files one microtubule at a time
    data_files_["occupancy"].WriteData(occupancy, max_length);
    if (proteins_.motors_.active_) {
      data_files_["motorID"].WriteData(motor_IDs, max_length);
      data_files_["motor_trailing"].WriteData(motor_head_status, max_length);
      if (proteins_.motors_.tethering_active_) {
        data_files_["tether_coord"].WriteData(tether_coords, max_length);
      }
    }
    if (proteins_.xlinks_.active_) {
      data_files_["xlinkID"].WriteData(xlink_IDs, max_length);
    }
  */
  }
  /*
  if (filaments_.mobile_) {
    data_files_["filament_coord"].WriteData(pf_coords, n_pfs);
  }
  */
}