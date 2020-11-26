#include "curator.hpp"

void Curator::CheckArgs(char *argv[]) {

  param_file_ = argv[1];
  sim_name_ = argv[2];
  test_mode_ = argv[3];
  int argc{0};
  while (true) {
    if (argv[argc] == nullptr) {
      break;
    }
    argc++;
  }
  if (argc <= 2 or argc >= 5) {
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

void Curator::GenerateLogFile() {

  char log_filename[256];
  sprintf(log_filename, "%s.log", sim_name_);
  // Check to see if sim files already exist
  if (FileExists(log_filename)) {
    printf("Simulation log file with this name already exists!\n");
    printf("Do you wish to overwrite these data? y/n\n");
    std::string response;
    int n_responses{0};
    bool response_unacceptable{true};
    while (response_unacceptable) {
      std::getline(std::cin, response);
      if (response == "n" or response == "N") {
        printf("Simulation terminated.\n");
        exit(1);
      } else if (response == "y" or response == "Y") {
        printf("Very well. ");
        printf("Overwriting data for sim '%s'\n\n", sim_name_);
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
  log_file_ = fopen(log_filename, "w");
  if (log_file_ == nullptr) {
    printf("Error; cannot open log file '%s'\n", log_file_name);
  }
}

void Curator::GenerateDataFiles() {

  // Open occupancy file, which stores the species ID of each occupant
  // (or -1 for none) for all MT sites during data collection (DC)
  data_files_.emplace("occupancy", SysFile(sim_name_, "occupancy"));
  // Motor-related files
  if (properties_->kinesin4.step_active_ < params_.n_steps) {
    // Open motor ID file, which stores the unique ID of all bound motors
    // (unbound not tracked) and their respective site indices during DC
    data_files_.emplace("motorID", SysFile(sim_name_, "motorID"));
    // bool; simply says if motor head is trailing or not
    data_files_.emplace("motor_trailing", SysFile(sim_name_, "motor_trailing"));
    if (properties_->kinesin4.tethering_active_) {
      // Open tether coord file, which stores the coordinates
      // of the anchor points of tethered motors
      data_files_.emplace("tether_coord", SysFile(sim_name_, "tether_coord"));
      // Open motor extension file, which stores the number of motors
      // with a certain tether extension for all possible extensions
      data_files_.emplace("motor_dx", SysFile(sim_name_, "motor_dx"));
      // Open motor force file, which stores the sum
      // of forces coming from motor tether extensions
      data_files_.emplace("motor_force", SysFile(sim_name_, "motor_force"));
    }
  }
  // Crosslinker-related files
  if (properties_->prc1.population_active_) {
    // Open xlink ID file, which does the same
    // as the motor ID file but for xlinks
    data_files_.emplace("xlinkID", SysFile(sim_name_, "xlinkID"));
    if (properties_->prc1.crosslinking_active_) {
      // Open xlink extension file, which stores the number of stage-2
      // xlinks at a certain extension for all possible extensions
      data_files_.emplace("xlink_dx", SysFile(sim_name_, "xlink_dx"));
      // Open xlink force file, which stores the sum
      // of forces coming from xlink extensions
      data_files_.emplace("xlink_force", SysFile(sim_name_, "xlink_force"));
    }
  }
  if (params_.filaments.diffusion_on) {
    // Open mt coord file, which stores the coordinates
    // of the left-most edge of each microtubule during DC
    data_files_.emplace("mt_coord", SysFile(sim_name_, "mt_coord"));
  }
  if (properties_->kinesin4.tethering_active_ and
      properties_->prc1.crosslinking_active_) {
    // Open total force file, which stores the sum of ALL
    // forces coming from xlink and motor tether extensions
    data_files_.emplace("total_force", SysFile(sim_name_, "total_force"));
  }
}

void Curator::ParseParameters() {

  // Check to make sure param file actually exists
  if (!FileExists(param_file_)) {
    Log("  Error: parameter file does not exist; aborting\n");
    exit(1);
  }
  /*
    // Parse parameter file into a YAML node
    YAML::Node input = YAML::LoadFile(param_file_);
    for (auto it : input) {
      Str key{it->first.as<Str>()};
      auto val =
    }
    */

  // Transfer values from input param node to system_parameters structure
  try {
    params_.seed = input["seed"].as<unsigned long>();
  } catch (const YAML::BadConversion error) {
    params_.seed = (unsigned long)(input["seed"].as<double>());
  }
  try {
    params_.n_steps = input["n_steps"].as<unsigned long>();
  } catch (const YAML::BadConversion error) {
    try {
      params_.n_steps = (unsigned long)input["n_steps"].as<double>();
    } catch (const YAML::BadConversion error) {
      params_.n_steps = (unsigned long)input["n_steps"].as<int>();
    }
  }
  params_.n_datapoints = input["n_datapoints"].as<int>();
  params_.data_threshold = input["data_threshold"].as<int>();
  params_.delta_t = input["delta_t"].as<double>();
  params_.kbT = input["kbT"].as<double>();
  params_.eta = input["eta"].as<double>();
  /* Microtubule parameters below */
  YAML::Node mts = input["filaments"];
  params_.filaments.count = mts["count"].as<int>();
  params_.filaments.length = mts["length"].as<std::vector<int>>();
  params_.filaments.y_dist = mts["y_dist"].as<double>();
  params_.filaments.site_size = mts["site_size"].as<double>();
  params_.filaments.radius = mts["radius"].as<double>();
  params_.filaments.elevation = mts["elevation"].as<double>();
  params_.filaments.start_coord = mts["start_coord"].as<std::vector<double>>();
  params_.filaments.immobile_until =
      mts["immobile_until"].as<std::vector<double>>();
  params_.filaments.applied_force = mts["applied_force"].as<double>();
  params_.filaments.printout_on = mts["printout_on"].as<bool>();
  params_.filaments.diffusion_on = mts["diffusion_on"].as<bool>();
  /* Motor parameters below */
  YAML::Node motors = input["motors"];
  params_.motors.n_runs_desired = motors["n_runs_desired"].as<size_t>();
  try {
    params_.motors.lattice_coop_range = motors["lattice_coop_range"].as<int>();
  } catch (const YAML::BadConversion error) {
    params_.motors.lattice_coop_range =
        (int)std::round(motors["lattice_coop_range"].as<double>());
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
  params_.motors.endpausing_active = motors["endpausing_active"].as<bool>();
  params_.motors.tethers_active = motors["tethers_active"].as<bool>();
  /* Xlink parameters below */
  YAML::Node xlinks = input["xlinks"];
  params_.xlinks.k_on = xlinks["k_on"].as<double>();
  params_.xlinks.c_bulk = xlinks["c_bulk"].as<double>();
  params_.xlinks.c_eff_bind = xlinks["c_eff_bind"].as<double>();
  params_.xlinks.k_off_i = xlinks["k_off_i"].as<double>();
  params_.xlinks.k_off_ii = xlinks["k_off_ii"].as<double>();
  params_.xlinks.r_0 = xlinks["r_0"].as<double>();
  params_.xlinks.k_spring = xlinks["k_spring"].as<double>();
  params_.xlinks.diffu_coeff_i = xlinks["diffu_coeff_i"].as<double>();
  params_.xlinks.diffu_coeff_ii = xlinks["diffu_coeff_ii"].as<double>();
  params_.xlinks.interaction_energy = xlinks["interaction_energy"].as<double>();
  // Store params pointer as parameters_ in Curator
  unsigned long n_steps{params_.n_steps};
  double delta_t{params_.delta_t};
  Log("Reading params from %s:\n\n", param_file_);
  Log("  General simulation parameters:\n");
  Log("    seed = %lu\n", params_.seed);
  Log("    n_steps = %lu\n", params_.n_steps);
  Log("    n_datapoints = %i\n", params_.n_datapoints);
  Log("    data_threshold = %i steps\n", params_.data_threshold);
  Log("    delta_t = %g s\n", params_.delta_t);
  Log("    kbT = %g pN*nm\n", params_.kbT);
  Log("    eta = %g (pN*s)/um^2\n", params_.eta);
  Log("\n  Microtubule (mt) parameters:\n");
  Log("    count = %i\n", params_.filaments.count);
  // Check to make sure there are enough vector entries for given MT count
  int n_lengths = input["filaments"]["length"].size();
  int n_start_coords = input["filaments"]["start_coord"].size();
  int n_immo = input["filaments"]["immobile_until"].size();
  if (params_.filaments.count > n_lengths or
      params_.filaments.count > n_start_coords or
      params_.filaments.count > n_immo) {
    Log("\nToo few parameters input for filaments\n");
    ErrorExit("Curator::ParseParameters()");
  }
  for (int i_mt = 0; i_mt < n_lengths; i_mt++) {
    Log("    length = %i sites for mt %i\n", params_.filaments.length[i_mt],
        i_mt);
  }
  Log("    y_dist = %g nm between MTs\n", params_.filaments.y_dist);
  Log("    site_size = %g nm\n", params_.filaments.site_size);
  Log("    radius = %g nm\n", params_.filaments.radius);
  Log("    elevation = %g nm above surface\n", params_.filaments.elevation);
  for (int i_mt = 0; i_mt < n_start_coords; i_mt++) {
    double start_coord = params_.filaments.start_coord[i_mt];
    Log("    start_coord = %g sites for mt %i\n", start_coord, i_mt);
  }
  for (int i_mt = 0; i_mt < n_immo; i_mt++) {
    double immo = params_.filaments.immobile_until[i_mt];
    Log("    immobile until = %g s for mt %i\n", immo, i_mt);
  }
  Log("    applied_force = %g pN\n", params_.filaments.applied_force);
  Log("    printout_on = %s\n",
      params_.filaments.printout_on ? "true" : "false");
  Log("    diffusion_on = %s\n",
      params_.filaments.diffusion_on ? "true" : "false");
  Log("\n  Kinesin (motor) parameters:\n");
  // Log("    lattice_coop_alpha = %g\n",
  // params_.motors.lattice_coop_alpha);
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
  Log("    k_on = %g /(nM*s)\n", params_.xlinks.k_on);
  Log("    c_bulk = %g nM\n", params_.xlinks.c_bulk);
  Log("    c_eff_bind = %g nM\n", params_.xlinks.c_eff_bind);
  Log("    k_off_i = %g /s\n", params_.xlinks.k_off_i);
  Log("    k_off_ii = %g /s\n", params_.xlinks.k_off_ii);
  Log("    r_0 = %g nm\n", params_.xlinks.r_0);
  Log("    k_spring = %g pN/nm\n", params_.xlinks.k_spring);
  Log("    diffu_coeff_i = %g um^2/s\n", params_.xlinks.diffu_coeff_i);
  Log("    diffu_coeff_ii = %g um^2/s\n", params_.xlinks.diffu_coeff_ii);
  Log("    interaction_energy = %g kbT\n", params_.xlinks.interaction_energy);
  Log("\nTotal simulation duration: %g seconds\n", delta_t * n_steps);
}

void Curator::SetLocalParameters() {

  size_t n_steps{params_.n_steps};
  int n_datapoints{params_.n_datapoints};
  data_threshold_ = params_.data_threshold;
  n_steps_recorded_ = n_steps - data_threshold_;
  n_steps_per_output_ = n_steps_recorded_ / n_datapoints;
  equil_milestone_ = data_threshold_ / 10;
  data_milestone_ = n_steps_recorded_ / 10;
  for (int i = 0; i < 5; i++) {
    t_motors_[i] = 0;
    t_xlinks_[i] = 0;
    t_MTs_[i] = 0;
  }
  t_motors_[5] = 0;
}

void Curator::InitializeSimObjects() {

  gsl.Initialize(params);
  filaments.Initialize(this, params;
  proteins.Initialize(this, params)
}

void Curator::CheckEquilibration() {

  if (!equilibrating_) {
    return;
  }
  if (i_step_ >= n_steps_equil_) {
    if (proteins.motors_.equilibrated_ and proteins.xlinks_.equilibrated_) {
      equilibrated_ = true;
      n_steps_equil_ = i_step;
    }
  }
}

void Curator::CheckPrintProgress() {

  int p_report{10};
  if (equilibrating_) {
    if (i_step_ <= n_steps_preequil_) {
      if (i_step_ % (n_steps_preequil_ / (100 / p_report)) == 0)
        Log("Run '%s' pre-equilibration %lu%% complete. (step %lu)\n",
            sim_name_, i_step_ / n_steps_equil_ * 100, i_step_);
    }
  } else {
    size_t delta{i_step_ - n_steps_equil_};
    size_t delta_tot{n_steps_tot_ - n_steps_equil_};
    if (delta % (delta_t / (100 / p_report)) == 0) {
      Log("Run '%s' data collection %lu%% complete. (step %lu)\n", sim_name_,
          delta / delta_tot * 100, i_step_);
    }
  }
  if (i_step_ >= n_steps_tot_) {
    sim_running_ = false;
  }
}

void Curator::OutputData() {

  int n_mts{params_.filaments.count};
  int max_length{0};
  for (int i_mt{0}; i_mt < n_mts; i_mt++) {
    if (params_.filaments.length[i_mt] > max_length)
      max_length = params_.filaments.length[i_mt];
  }
  // Create arrays to store data; ptrs to write it to file
  int mt_coords[n_mts];
  // For extension statistics, data is on a per-extension basis
  int motor_cutoff = properties_->kinesin4.teth_cutoff_;
  int motor_extensions[2 * motor_cutoff + 1];
  int xlink_cutoff = properties_->prc1.dist_cutoff_;
  int xlink_extensions[xlink_cutoff + 1];
  // Back to normal per-MT array format
  double motor_forces[n_mts];
  double xlink_forces[n_mts];
  double total_forces[n_mts];
  // Run through all MTs and get data for each
  for (int i_mt{0}; i_mt < n_mts; i_mt++) {
    int mt_length{params_.filaments.length[i_mt]};
    Microtubule *mt = &properties_->filaments.mt_list_[i_mt];
    // Create arrays & ptrs for intraMT data
    int motor_IDs[max_length];
    int xlink_IDs[max_length];
    int occupancy[max_length];
    double tether_coords[max_length];
    bool motor_head_status[max_length];
    for (int i_site{0}; i_site < max_length; i_site++) {
      motor_head_status[i_site] = false;
    }
    // Run through all sites on this particular MT
    for (int i_site{0}; i_site < mt_length; i_site++) {
      Tubulin *site{&mt->lattice_[i_site]};
      // If unoccupied, store the speciesID of tubulin to occupancy
      // and an ID of -1 (null) to motor/xlink ID files
      if (!site->occupied_) {
        occupancy[i_site] = site->species_id_;
        motor_IDs[i_site] = -1;
        xlink_IDs[i_site] = -1;
        tether_coords[i_site] = -1;
      }
      // If occupied by xlink, store its species ID to occupancy_file,
      // its unique ID to the xlink ID file, and -1 to motor ID file
      else if (site->xlink_head_ != nullptr) {
        AssociatedProtein *xlink{site->xlink_head_->xlink_};
        occupancy[i_site] = xlink->species_id_;
        motor_IDs[i_site] = -1;
        xlink_IDs[i_site] = xlink->id_;
        tether_coords[i_site] = -1;
      }
      // If occupied by motor, store its species ID to occupancy_file,
      // its unique ID to the motor ID file, and -1 to xlink ID file
      else if (site->motor_head_ != nullptr) {
        Kinesin *motor{site->motor_head_->motor_};
        occupancy[i_site] = motor->species_id_;
        motor_IDs[i_site] = motor->id_;
        xlink_IDs[i_site] = -1;
        tether_coords[i_site] = -1;
        motor_head_status[i_site] = site->motor_head_->trailing_;
        if (motor->tethered_) {
          AssociatedProtein *xlink{motor->xlink_};
          if (xlink->heads_active_ > 0) {
            double anchor_coord = xlink->GetAnchorCoordinate();
            // printf("wrote teth coord %g\n", anchor_coord);
            tether_coords[i_site] = anchor_coord;
            double stalk_coord = motor->GetStalkCoordinate();
            double teth_dist = fabs(anchor_coord - stalk_coord);
            if (teth_dist > motor->teth_cutoff_) {
              Log("woah, teth dist is %g\n", teth_dist);
            }
          }
        }
      }
    }
    // Pad written data with NULL entries (-1) for shorter MTs
    for (int i_site{mt_length}; i_site < max_length; i_site++) {
      occupancy[i_site] = -1;
      motor_IDs[i_site] = -1;
      xlink_IDs[i_site] = -1;
      tether_coords[i_site] = -1;
    }
    mt_coords[i_mt] = mt->coord_;
    motor_forces[i_mt] = mt->GetNetForce_Motors();
    xlink_forces[i_mt] = mt->GetNetForce_Xlinks();
    total_forces[i_mt] = mt->GetNetForce();
    // Write the data to respective files one microtubule at a time
    fwrite(occupancy, sizeof(int), max_length, properties_->occupancy_file_);
    if (properties_->kinesin4.step_active_ < params_.n_steps) {
      fwrite(motor_IDs, sizeof(int), max_length, properties_->motor_ID_file_);
      fwrite(motor_head_status, sizeof(bool), max_length,
             properties_->motor_head_status_file_);
      if (properties_->kinesin4.tethering_active_) {
        fwrite(tether_coords, sizeof(double), max_length,
               properties_->tether_coord_file_);
      }
    }
    if (properties_->prc1.population_active_) {
      fwrite(xlink_IDs, sizeof(int), max_length, properties_->xlink_ID_file_);
    }
  }
  // Scan through kinesin4/prc1 statistics to get extension occupancies
  KinesinManagement *kinesin4{&properties_->kinesin4};
  for (int i_ext{0}; i_ext <= 2 * motor_cutoff; i_ext++) {
    motor_extensions[i_ext] = kinesin4->n_bound_teth_[i_ext];
  }
  AssociatedProteinManagement *prc1{&properties_->prc1};
  for (int i_ext{0}; i_ext <= xlink_cutoff; i_ext++) {
    xlink_extensions[i_ext] = 0;
    for (int n_neighbs{0}; n_neighbs <= prc1->max_neighbs_; n_neighbs++) {
      xlink_extensions[i_ext] += prc1->n_bound_ii_[n_neighbs][i_ext];
    }
  }
  if (properties_->kinesin4.step_active_ > params_.n_steps) {
    if (properties_->kinesin4.tethering_active_) {
      fwrite(motor_forces, sizeof(double), n_mts,
             properties_->motor_force_file_);
      fwrite(motor_extensions, sizeof(int), 2 * motor_cutoff + 1,
             properties_->motor_extension_file_);
    }
  }
  if (properties_->prc1.population_active_) {
    if (properties_->prc1.crosslinking_active_) {
      fwrite(xlink_extensions, sizeof(int), xlink_cutoff + 1,
             properties_->xlink_extension_file_);
      fwrite(xlink_forces, sizeof(double), n_mts,
             properties_->xlink_force_file_);
    }
  }
  if (params_.filaments.diffusion_on) {
    fwrite(mt_coords, sizeof(int), n_mts, properties_->mt_coord_file_);
  }
  if (properties_->kinesin4.tethering_active_ and
      properties_->prc1.crosslinking_active_) {
    fwrite(total_forces, sizeof(double), n_mts, properties_->total_force_file_);
  }
}

void Curator::OutputSimDuration() {

  int n_per_sec{sys_clock::period::den};
  finish_ = sys_clock::now();
  double sim_duration{(double)(finish_ - start_).count()};
  Log("\nTime to execute sim '%s': %.2f seconds.\n", sim_name_,
      sim_duration / n_per_sec);
  Log(" (%i datapoints recorded)\n", n_datapoints_recorded_);
  Log("   -KinesinManagement::RunKMC(): %.2f\n", t_motors_[0] / n_per_sec);
  Log("      -CheckEquilibration(): %.2f\n", t_motors_[1] / n_per_sec);
  Log("      -UpdateLists(): %.2f\n", t_motors_[2] / n_per_sec);
  Log("      -SampleEventStatistics(): %.2f\n", t_motors_[3] / n_per_sec);
  Log("      -GenerateExecutionSequence(): %.2f\n", t_motors_[4] / n_per_sec);
  Log("      -ExecuteEvents() %.2f\n", t_motors_[4] / n_per_sec);
  Log("   -AssociatedProteinManagement::RunKMC(): %.2f\n",
      t_xlinks_[0] / n_per_sec);
  Log("      -UpdateLists(): %.2f\n", t_xlinks_[1] / n_per_sec);
  Log("      -SampleEventStatistics(): %.2f\n", t_xlinks_[2] / n_per_sec);
  Log("      -GenerateExecutionSequence(): %.2f\n", t_xlinks_[3] / n_per_sec);
  Log("      -ExecuteEvents(): %.2f\n", t_xlinks_[4] / n_per_sec);
  /*
  Log("   -MTs: %.2f\n", t_MTs_[0] / n_per_sec);
  Log("      -Summing forces: %.2f\n", t_MTs_[1] / n_per_sec);
  Log("      -Calculating displacement: %.2f\n", t_MTs_[2] / n_per_sec);
  Log("      -Execution: %.2f\n", t_MTs_[3] / n_per_sec);
  */
}

void Curator::CloseDataFiles() {

  fclose(log_file_);
  fclose(properties_->occupancy_file_);
  if (properties_->kinesin4.step_active_ < params_.n_steps) {
    fclose(properties_->motor_ID_file_);
    fclose(properties_->motor_head_status_file_);
    if (properties_->kinesin4.tethering_active_) {
      fclose(properties_->tether_coord_file_);
      fclose(properties_->motor_extension_file_);
      fclose(properties_->motor_force_file_);
    }
  }
  if (properties_->prc1.population_active_) {
    fclose(properties_->xlink_ID_file_);
    if (properties_->prc1.crosslinking_active_) {
      fclose(properties_->xlink_extension_file_);
      fclose(properties_->xlink_force_file_);
    }
  }
  if (params_.filaments.diffusion_on) {
    fclose(properties_->mt_coord_file_);
  }
  if (properties_->kinesin4.tethering_active_ and
      properties_->prc1.crosslinking_active_) {
    fclose(properties_->total_force_file_);
  }
}

void Curator::ErrorExit(const char *function_name) {

  Log("\nFatal error in %s\n", sim_name_);
  Log("Function name: %s\n", function_name);
  Log("Step no: #%i\n", properties_->current_step_);
  Log(" *** EXITING ***\n");
  exit(1);
}

void Curator::EvolveSimulation() {

  proteins.RunKMC();
  filaments.RunBD();
  CheckPrintProgress();
  OutputData();
  /*
  if (params_.filaments.printout_on and i_step % 1000 == 0) {
    PrintMicrotubules(0);
  }
  */
}

void Curator::PrintMicrotubules() {

  /*
  int n_mts = params_.filaments.count;
  // Figure out which MT is the farthest left
  int leftmost_coord = 0;
  for (int i_mt = 0; i_mt < n_mts; i_mt++) {
    int mt_coord = properties_->filaments.mt_list_[i_mt].coord_;
    if (i_mt == 0)
      leftmost_coord = mt_coord;
    else if (mt_coord < leftmost_coord)
      leftmost_coord = mt_coord;
  }
  // Print out MTs
  for (int i_mt = n_mts - 1; i_mt >= 0; i_mt--) {
    int mt_length = params_.filaments.length[i_mt];
    Microtubule *mt = &properties_->filaments.mt_list_[i_mt];
    int mt_coord = mt->coord_;
    int delta = mt_coord - leftmost_coord;
    if (delta < 0) {
      Log("wat\n");
      exit(1);
    }
    for (int i_entry = 0; i_entry < delta; i_entry++) {
      printf(" ");
    }
    for (int i_site = 0; i_site < mt_length; i_site++) {
      Tubulin *site = &mt->lattice_[i_site];
      if (site->occupied_ == false)
        printf("=");
      else if (site->xlink_head_ != nullptr) {
        AssociatedProtein *xlink = site->xlink_head_->xlink_;
        if (xlink->heads_active_ == 1) {
          if (xlink->tethered_ == false)
            printf("i");
          else if (xlink->motor_->heads_active_ == 0)
            printf("F");
          else
            printf("I");
        } else if (xlink->heads_active_ == 2)
          if (xlink->tethered_ == false) {
            printf("x");
            // printf("%i",
            // xlink->x_dist_);
          } else
            printf("X");
        else {
          printf("no sunny. look in wallace's print\n");
          exit(1);
        }
      } else if (site->motor_head_ != nullptr) {
        Kinesin *motor = site->motor_head_->motor_;
        if (motor->heads_active_ == 1) {
          if (motor->tethered_ == false)
            std::cout << site->motor_head_->ligand_;
          //	printf("m");
          else {
            printf("[");
            std::cout << site->motor_head_->ligand_;
            printf("]");
          }
        } else if (motor->heads_active_ == 2) {
          int i_front = motor->head_one_.site_->index_;
          int i_rear = motor->head_two_.site_->index_;
          if (i_front > i_rear) {
            if (i_site == i_rear) {
              if (motor->tethered_ == false) {
                printf("(");
                std::cout << motor->head_two_.ligand_;
                printf("|");
              } else
                printf("%i", motor->x_dist_doubled_ / 10);
              //	printf("[");
            }
            if (i_site == i_front) {
              if (motor->tethered_ == false) {
                std::cout << motor->head_one_.ligand_;
                printf(")");
              } else
                // printf("]");
                printf("%i", motor->x_dist_doubled_ % 10);
            }
          } else if (i_front < i_rear) {
            if (i_site == i_front) {
              if (motor->tethered_ == false) {
                printf("(");
                std::cout << motor->head_one_.ligand_;
                printf("|");
              } else
                printf("%i", motor->x_dist_doubled_ / 10);
              //	printf("[");
            }
            if (i_site == i_rear) {
              if (motor->tethered_ == false) {
                std::cout << motor->head_two_.ligand_;
                printf(")");
              } else
                // printf("]");
                printf("%i", motor->x_dist_doubled_ % 10);
            }
          } else {
            printf("error in print MT\n");
            exit(1);
          }
        }
      }
    }
    printf(" %i\n", mt->polarity_);
  }
  printf("\n");
  */

  /* site coordinate printout below */
  /*
     int mt1_coord = properties_->filaments.mt_list_[0].coord_;
     int mt2_coord = properties_->filaments.mt_list_[1].coord_;
     int greater_coord = 0;
     if(mt1_coord > mt2_coord)
     greater_coord = mt1_coord;
     else
     greater_coord = mt2_coord;
     int extra_digits = 0;
     for(int i_site = 0; i_site < mt_length + greater_coord; i_site++){
     if(extra_digits > 0)
     extra_digits--;
     else if(i_site%5 == 0){
     printf("%i", i_site);
     if(i_site < 10)
     extra_digits = 0;
     else if(i_site < 100)
     extra_digits = 1;
     else if(i_site < 1000)
     extra_digits = 2;
     else if(i_site < 10000)
     extra_digits = 3;
     else{
     printf("what the fuck are u. why do you need >10,000 sites??\n");
     exit(1);
     }
     }
     else if(i_site == (mt_length + greater_coord) - 1)
     printf("%i", i_site);
     else
     printf(" ");
     }
    */
}

void Curator::PrintMicrotubules(double pause_duration) {

  PrintMicrotubules();
  if (pause_duration > 0)
    PauseSim(pause_duration);
}

void Curator::PauseSim(double duration) {

  // Duration should be input in seconds
  pause_dur_.tv_sec = (int)duration;
  pause_dur_.tv_nsec = (duration - (int)duration) * 100000000;
  nanosleep(&pause_dur_, NULL);
}

void Curator::StartDataCollection() {

  size_t n_steps{params_.n_steps};
  int n_datapoints{params_.n_datapoints};
  data_threshold_ = properties_->current_step_;
  n_steps_recorded_ = n_steps - data_threshold_;
  n_steps_per_output_ = n_steps_recorded_ / n_datapoints;
  data_milestone_ = n_steps_recorded_ / 10;
  Log("Sim '%s' data collection has begun.\n", sim_name_);
  Log("DATA_THRESHOLD = %zu\n", data_threshold_);
}

void Curator::TerminateSimulation() {

  properties_->sim_running_ = false;
  size_t steps_recorded{properties_->current_step_ - data_threshold_};
  size_t n_datapoints{steps_recorded / n_steps_per_output_};
  Log("Sim '%s' data collection terminated early after sufficient unjammed "
      "kinesin unbinding events\n",
      sim_name_, properties_->kinesin4.n_runs_recorded_);
  Log("N_STEPS = %zu\n", properties_->current_step_);
  Log("N_DATAPOINTS = %zu\n", n_datapoints);
}

void Curator::CleanUp() {

  OutputSimDuration();
  CloseDataFiles();
  properties_->gsl.CleanUp();
}
