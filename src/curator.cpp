#include "curator.hpp"
#include "yaml-cpp/yaml.h"
#include <filesystem>
#include <iostream>

void Curator::CheckArgs(char *argv[]) {

  param_file_ = argv[1];
  sim_name_ = argv[2];
  test_mode_ = argv[3];
  // This seems terrible but I cannot find a better way
  int argc{0};
  while (true) {
    if (argv[argc] == nullptr) {
      break;
    }
    argc++;
  }
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
}

void Curator::ParseParameters() {

  // Check to make sure param file actually exists
  if (!std::filesystem::exists(param_file_)) {
    Log("  Error: parameter file does not exist; aborting\n");
    exit(1);
  }
  // Parse parameter file into a YAML node
  YAML::Node input = YAML::LoadFile(param_file_);
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
  Log("Reading params from %s:\n\n", param_file_);
  Log("  General simulation parameters:\n");
  Log("    seed = %lu\n", params_.seed);
  Log("    kbT = %g pN*nm\n", params_.kbT);
  Log("    eta = %g (pN*s)/um^2\n", params_.eta);
  Log("    dt = %g s\n", params_.dt);
  Log("    t_run = %g s\n", params_.t_run);
  Log("    t_equil = %g s\n", params_.t_equil);
  Log("    t_snapshot = %g s\n", params_.t_snapshot);
  Log("    (%zu datapoints expected)\n", params_.t_run / params_.t_snapshot);
  Log("    dynamic_equil = %s\n", params_.dynamic_equil ? "true" : "false");
  Log("    verbosity = %zu\n", params_.verbosity);
  Log("\n  Filament parameters:\n");
  Log("    count = %i\n", params_.filaments.count);
  // Check to make sure there are enough vector entries for given MT count
  if (params_.filaments.count > params_.filaments.length.size() or
      params_.filaments.count > params_.filaments.start_coord.size() or
      params_.filaments.count > params_.filaments.applied_force.size() or
      params_.filaments.count > params_.filaments.immobile_until.size() or
      Sys::_n_dims_max > params_.filaments.dim_enabled.size()) {
    Log("\nToo few parameters input for filaments\n");
    ErrorExit("Curator::ParseParameters()");
  }
  for (int i_pf{0}; i_pf < params_.filaments.count; i_pf++) {
    Log("    length[%i] = %i sites\n", i_pf, params_.filaments.length[i_pf]);
    Log("    start_coord[%i] = %i sites\n", i_pf,
        params_.filaments.start_coord[i_pf]);
    Log("    applied_force[%i] = %i pN\n", i_pf,
        params_.filaments.applied_force[i_pf]);
    Log("    immobile_until[%i] = %g s\n", i_pf,
        params_.filaments.immobile_until[i_pf]);
  }
  for (int i_dim{0}; i_dim < Sys::_n_dims_max; i_dim++) {
    Log("    dim_enabled[%i] = %s\n", i_dim,
        params_.filaments.dim_enabled[i_dim] ? "true" : "false");
  }
  Log("    y_dist_init = %g nm between MTs\n", params_.filaments.y_dist_init);
  Log("    site_size = %g nm\n", params_.filaments.site_size);
  Log("    radius = %g nm\n", params_.filaments.radius);
  Log("\n  Kinesin (motor) parameters:\n");
  Log("    n_runs_desired = %zu\n", params_.motors.n_runs_desired);
  Log("    lattice_coop_range = %i\n", params_.motors.lattice_coop_range);
  Log("    lattice_coop_Emax_solo = -%g kbT\n",
      params_.motors.lattice_coop_Emax_solo);
  Log("    lattice_coop_Emax_bulk = -%g kbT\n",
      params_.motors.lattice_coop_Emax_bulk);
  Log("    interaction_energy = -%g kbT\n", params_.motors.interaction_energy);
  Log("    t_active = %g seconds\n", params_.motors.t_active);
  Log("    k_on = %g /(nM*s)\n", params_.motors.k_on);
  Log("    c_bulk = %g nM\n", params_.motors.c_bulk);
  Log("    c_eff_bind = %g nM\n", params_.motors.c_eff_bind);
  Log("    k_on_ATP = %g /(mM*s)\n", params_.motors.k_on_ATP);
  Log("    c_ATP = %g mM\n", params_.motors.c_ATP);
  Log("    k_hydrolyze = %g /s\n", params_.motors.k_hydrolyze);
  Log("    k_off_i = %g /s\n", params_.motors.k_off_i);
  Log("    k_off_ii = %g /s\n", params_.motors.k_off_ii);
  Log("    applied_force = %g pN\n", params_.motors.applied_force);
  Log("    internal_force = %g pN\n", params_.motors.internal_force);
  Log("    sigma_off_i = %g nm\n", params_.motors.sigma_off_i);
  Log("    sigma_off_ii = %g nm\n", params_.motors.sigma_off_ii);
  Log("    sigma_ATP = %g nm\n", params_.motors.sigma_ATP);
  Log("    k_tether = %g /(nM*s)\n", params_.motors.k_tether);
  Log("    c_eff_tether = %g nM\n", params_.motors.c_eff_tether);
  Log("    k_untether = %g /s\n", params_.motors.k_untether);
  Log("    r_0 = %g nm\n", params_.motors.r_0);
  Log("    k_spring = %g pN/nm\n", params_.motors.k_spring);
  Log("    k_slack = %g pN/nm\n", params_.motors.k_slack);
  Log("    tethers_active = %s\n",
      params_.motors.tethers_active ? "true" : "false");
  Log("    endpausing_active = %s\n",
      params_.motors.endpausing_active ? "true" : "false");
  Log("\n  Crosslinker (xlink) parameters:\n");
  Log("    interaction_energy = %g kbT\n", params_.xlinks.interaction_energy);
  Log("    t_active = %g seconds\n", params_.xlinks.t_active);
  Log("    k_on = %g /(nM*s)\n", params_.xlinks.k_on);
  Log("    c_bulk = %g nM\n", params_.xlinks.c_bulk);
  Log("    c_eff_bind = %g nM\n", params_.xlinks.c_eff_bind);
  Log("    k_off_i = %g /s\n", params_.xlinks.k_off_i);
  Log("    k_off_ii = %g /s\n", params_.xlinks.k_off_ii);
  Log("    r_0 = %g nm\n", params_.xlinks.r_0);
  Log("    k_spring = %g pN/nm\n", params_.xlinks.k_spring);
  Log("    d_i = %g um^2/s\n", params_.xlinks.d_i);
  Log("    d_ii = %g um^2/s\n", params_.xlinks.d_ii);
}

void Curator::InitializeSimulation() {

  double dt{params_.dt};
  n_steps_pre_equil_ = (size_t)std::round(params_.t_equil / dt);
  n_steps_equil_ = n_steps_pre_equil_;
  n_steps_run_ = (size_t)std::round(params_.t_run / dt);
  n_steps_snapshot_ = (size_t)std::round(params_.t_snapshot / dt);
  verbosity_ = params_.verbosity;

  gsl_.Initialize(params_.seed);
  proteins_.Initialize(this, &params_);
  filaments_.Initialize(this, &params_);
}

void Curator::GenerateDataFiles() {

  // Open occupancy file, which stores the species ID of each occupant
  // (or -1 for none) for all MT sites during data collection (DC)
  files_.AddDataFile(sim_name_, "occupancy");
  // Motor-related files
  if (proteins_.motors_.active_) {
    // Open motor ID file, which stores the unique ID of all bound motors
    // (unbound not tracked) and their respective site indices during DC
    files_.AddDataFile(sim_name_, "motorID");
    // bool; simply says if motor head is trailing or not
    files_.AddDataFile(sim_name_, "motor_trailing");
    if (proteins_.motors_.tethering_active_) {
      // Open tether coord file, which stores the coordinates
      // of the anchor points of tethered motors
      files_.AddDataFile(sim_name_, "tether_coord");
      // Open motor extension file, which stores the number of motors
      // with a certain tether extension for all possible extensions
      // // files_.AddDataFile(sim_name_, "motor_dx");
      // Open motor force file, which stores the sum
      // of forces coming from motor tether extensions
      // // files_.AddDataFile(sim_name_, "motor_force");
    }
  }
  // Crosslinker-related files
  if (proteins_.xlinks_.active_) {
    // Open xlink ID file, which does the same
    // as the motor ID file but for xlinks
    files_.AddDataFile(sim_name_, "xlinkID");
    if (proteins_.xlinks_.crosslinking_active_) {
      // Open xlink extension file, which stores the number of stage-2
      // xlinks at a certain extension for all possible extensions
      // // files_.AddDataFile(sim_name_, "xlink_dx");
      // Open xlink force file, which stores the sum
      // of forces coming from xlink extensions
      // // files_.AddDataFile(sim_name_, "xlink_force");
    }
  }
  if (filaments_.mobile_) {
    // Open mt coord file, which stores the coordinates
    // of the left-most edge of each microtubule during DC
    files_.AddDataFile(sim_name_, "filament_coord");
  }
  if (proteins_.motors_.tethering_active_ and
      proteins_.xlinks_.crosslinking_active_) {
    // Open total force file, which stores the sum of ALL
    // forces coming from xlink and motor tether extensions
    // // files_.AddDataFile(sim_name_, "total_force");
  }
}

void Curator::CheckPrintProgress() {

  int p_report{10};
  if (sim_equilibrating_) {
    if (i_step_ % (n_steps_pre_equil_ / (100 / p_report)) == 0 and
        i_step_ <= n_steps_pre_equil_) {
      Log("Run '%s' pre-equilibration %zu%% complete. (step %zu)\n", sim_name_,
          i_step_ / n_steps_pre_equil_ * 100, i_step_);
    }
    if (proteins_.motors_.equilibrated_ and proteins_.xlinks_.equilibrated_) {
      n_steps_equil_ = i_step_;
      sim_equilibrating_ = false;
      if (params_.dynamic_equil) {
        Log("Run '%s' dynamic equilibration complete.\n", sim_name_);
        Log("N_STEPS_EQUIL = %zu\n", n_steps_equil_);
      }
    }
  } else {
    size_t n_steps_so_far{i_step_ - n_steps_equil_};
    if (n_steps_so_far % (n_steps_run_ / (100 / p_report)) == 0) {
      Log("Run '%s' data collection %zu%% complete. (step %zu)\n", sim_name_,
          n_steps_so_far / n_steps_run_ * 100, i_step_);
    }
  }
  if (i_step_ >= n_steps_run_ + n_steps_equil_) {
    sim_running_ = false;
    Log("Run '%s' complete.\n", sim_name_);
    double sim_dur{(double)(SysClock::now() - start_time_).count()};
    Log("\nTime to execute run '%s': %.2f seconds.\n", sim_name_,
        sim_dur / SysClock::period::den);
  }
}

void Curator::OutputData() {

  if (sim_equilibrating_) {
    return;
  }
  if (i_step_ % n_steps_snapshot_ != 0) {
    return;
  }

  size_t n_pfs{params_.filaments.count};
  size_t max_length{0};
  for (int i_pf{0}; i_pf < n_pfs; i_pf++) {
    if (params_.filaments.length[i_pf] > max_length)
      max_length = params_.filaments.length[i_pf];
  }
  // Create arrays to store data; ptrs to write it to file
  int pf_coords[n_pfs];
  // Run through all MTs and get data for each
  for (int i_pf{0}; i_pf < n_pfs; i_pf++) {
    Protofilament *pf{&filaments_.list_[i_pf]};
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
        case Sys::_id_motor:
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
        case Sys::_id_xlink:
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
    pf_coords[i_pf] = pf->GetPos(0);
    // Write the data to respective files one microtubule at a time
    files_.data_["occupancy"].WriteData(occupancy, max_length);
    if (proteins_.motors_.active_) {
      files_.data_["motorID"].WriteData(motor_IDs, max_length);
      files_.data_["motor_trailing"].WriteData(motor_head_status, max_length);
      if (proteins_.motors_.tethering_active_) {
        files_.data_["tether_coord"].WriteData(tether_coords, max_length);
      }
    }
    if (proteins_.xlinks_.active_) {
      files_.data_["xlinkID"].WriteData(xlink_IDs, max_length);
    }
  }
  if (filaments_.mobile_) {
    files_.data_["filament_coord"].WriteData(pf_coords, n_pfs);
  }
}

void Curator::EvolveSimulation() {

  printf("hello walter\n");
  // proteins_.RunKMC();
  // filaments_.RunBD();
  // CheckPrintProgress();
  // OutputData();
}

void Curator::TerminateSimulation() {

  sim_running_ = false;
  Log("Sim '%s' terminated after sufficient data collection\n", sim_name_);
  Log("N_STEPS = %zu\n", i_step_ - n_steps_equil_);
  Log("N_DATAPOINTS = %zu\n", i_datapoint_);
}