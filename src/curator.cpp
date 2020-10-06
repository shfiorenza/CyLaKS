#include "curator.h"
#include "master_header.h"

Curator::Curator() {}

void Curator::InitializeSimulation(char *argv[], system_properties *properties,
                                   system_parameters *parameters) {

  properties_ = properties;
  parameters_ = parameters;
  // Check that user-input arguments are valid and parse them if they are
  CheckArgs(argv);
  // Log file saves all sim outputs to terminal
  GenerateLogFile();
  // Parameters from YAML file are transferred to local parameter structs
  ParseParameters();
  // Once parameters are parsed, the Curator's own parameters can be set
  SetLocalParameters();
  // Microtubules, motors, and crosslinkers are all initialized at once
  InitializeSimObjects();
  // Data files save all pertinent info: occupancy, coords, extensions, etc.
  GenerateDataFiles();
}

FILE *Curator::OpenFile(const char *file_name, const char *type) {

  FILE *file_ptr;
  if ((file_ptr = fopen(file_name, type)) == NULL) {
    fprintf(stderr, "Cannot open %s\n", file_name);
    exit(1);
  }
  return file_ptr;
}

bool Curator::FileExists(std::string file_name) {

  struct stat buffer;
  return (stat(file_name.c_str(), &buffer) != -1);
}

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
    printf("    xlink_bind_ii\n");
    printf("    motor_lattice_bind\n");
    printf("    motor_lattice_step\n");
    exit(1);
  }
}

void Curator::GenerateLogFile() {

  char log_file[256];
  sprintf(log_file, "%s.log", sim_name_);
  // Check to see if sim files already exist
  if (FileExists(log_file)) {
    printf("Simulation log file with this name already exists!\n");
    printf("Do you wish to overwrite these data? y/n\n");
    std::string response;
    bool response_unacceptable{true};
    int n_responses{0};
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
  log_file_ = OpenFile(log_file, "w");
}

void Curator::GenerateDataFiles() {

  char occupancy_file[256], motor_ID_file[256], xlink_ID_file[256],
      tether_coord_file[256], mt_coord_file[256], motor_extension_file[256],
      xlink_extension_file[256], motor_force_file[256], xlink_force_file[256],
      total_force_file[256], motor_head_status_file[256];
  // Generate names of output files based on the input simulation name
  sprintf(occupancy_file, "%s_occupancy.file", sim_name_);
  sprintf(motor_ID_file, "%s_motorID.file", sim_name_);
  sprintf(xlink_ID_file, "%s_xlinkID.file", sim_name_);
  sprintf(tether_coord_file, "%s_tether_coord.file", sim_name_);
  sprintf(mt_coord_file, "%s_mt_coord.file", sim_name_);
  sprintf(motor_extension_file, "%s_motor_extension.file", sim_name_);
  sprintf(xlink_extension_file, "%s_xlink_extension.file", sim_name_);
  sprintf(motor_force_file, "%s_motor_force.file", sim_name_);
  sprintf(xlink_force_file, "%s_xlink_force.file", sim_name_);
  sprintf(total_force_file, "%s_total_force.file", sim_name_);
  sprintf(motor_head_status_file, "%s_motor_head_status.file", sim_name_);
  // Open occupancy file, which stores the species ID of each occupant
  // (or -1 for none) for all MT sites during data collection (DC)
  properties_->occupancy_file_ = OpenFile(occupancy_file, "w");
  // Motor-related files
  if (properties_->kinesin4.step_active_ < parameters_->n_steps) {
    // Open motor ID file, which stores the unique ID of all bound motors
    // (unbound not tracked) and their respective site indices during DC
    properties_->motor_ID_file_ = OpenFile(motor_ID_file, "w");
    // bool; simply says if motor head is trailing or not
    properties_->motor_head_status_file_ =
        OpenFile(motor_head_status_file, "w");
    if (properties_->kinesin4.tethering_active_) {
      // Open tether coord file, which stores the coordinates
      // of the anchor points of tethered motors
      properties_->tether_coord_file_ = OpenFile(tether_coord_file, "w");
      // Open motor extension file, which stores the number of motors
      // with a certain tether extension for all possible extensions
      properties_->motor_extension_file_ = OpenFile(motor_extension_file, "w");
      // Open motor force file, which stores the sum
      // of forces coming from motor tether extensions
      properties_->motor_force_file_ = OpenFile(motor_force_file, "w");
    }
  }
  // Crosslinker-related files
  if (properties_->prc1.population_active_) {
    // Open xlink ID file, which does the same
    // as the motor ID file but for xlinks
    properties_->xlink_ID_file_ = OpenFile(xlink_ID_file, "w");
    if (properties_->prc1.crosslinking_active_) {
      // Open xlink extension file, which stores the number of stage-2
      // xlinks at a certain extension for all possible extensions
      properties_->xlink_extension_file_ = OpenFile(xlink_extension_file, "w");
    }
    // Open xlink force file, which stores the sum
    // of forces coming from xlink extensions
    properties_->xlink_force_file_ = OpenFile(xlink_force_file, "w");
  }
  if (parameters_->microtubules.diffusion_on) {
    // Open mt coord file, which stores the coordinates
    // of the left-most edge of each microtubule during DC
    properties_->mt_coord_file_ = OpenFile(mt_coord_file, "w");
  }
  if (properties_->kinesin4.tethering_active_ and
      properties_->prc1.crosslinking_active_) {
    // Open total force file, which stores the sum of ALL
    // forces coming from xlink and motor tether extensions
    properties_->total_force_file_ = OpenFile(total_force_file, "w");
  }
}

void Curator::ParseParameters() {

  // Check to make sure param file actually exists
  if (!FileExists(param_file_)) {
    Log("  Error: parameter file does not exist; aborting\n");
    exit(1);
  }

  // Parse parameter file into a YAML node
  YAML::Node input = YAML::LoadFile(param_file_);
  // Transfer values from input param node to system_parameters structure
  try {
    parameters_->seed = input["seed"].as<unsigned long>();
  } catch (const YAML::BadConversion error) {
    parameters_->seed = (unsigned long)(input["seed"].as<double>());
  }
  try {
    parameters_->n_steps = input["n_steps"].as<unsigned long>();
  } catch (const YAML::BadConversion error) {
    try {
      parameters_->n_steps = (unsigned long)input["n_steps"].as<double>();
    } catch (const YAML::BadConversion error) {
      parameters_->n_steps = (unsigned long)input["n_steps"].as<int>();
    }
  }
  parameters_->n_datapoints = input["n_datapoints"].as<int>();
  parameters_->data_threshold = input["data_threshold"].as<int>();
  parameters_->delta_t = input["delta_t"].as<double>();
  parameters_->kbT = input["kbT"].as<double>();
  parameters_->eta = input["eta"].as<double>();
  /* Microtubule parameters below */
  YAML::Node mts = input["microtubules"];
  parameters_->microtubules.count = mts["count"].as<int>();
  parameters_->microtubules.length = mts["length"].as<std::vector<int>>();
  parameters_->microtubules.y_dist = mts["y_dist"].as<double>();
  parameters_->microtubules.site_size = mts["site_size"].as<double>();
  parameters_->microtubules.radius = mts["radius"].as<double>();
  parameters_->microtubules.elevation = mts["elevation"].as<double>();
  parameters_->microtubules.start_coord =
      mts["start_coord"].as<std::vector<double>>();
  parameters_->microtubules.immobile_until =
      mts["immobile_until"].as<std::vector<double>>();
  parameters_->microtubules.applied_force = mts["applied_force"].as<double>();
  parameters_->microtubules.printout_on = mts["printout_on"].as<bool>();
  parameters_->microtubules.diffusion_on = mts["diffusion_on"].as<bool>();
  /* Motor parameters below */
  YAML::Node motors = input["motors"];
  parameters_->motors.n_runs_desired = motors["n_runs_desired"].as<size_t>();
  try {
    parameters_->motors.lattice_coop_range =
        motors["lattice_coop_range"].as<int>();
  } catch (const YAML::BadConversion error) {
    parameters_->motors.lattice_coop_range =
        (int)std::round(motors["lattice_coop_range"].as<double>());
  }
  parameters_->motors.lattice_coop_Emax_solo =
      motors["lattice_coop_Emax_solo"].as<double>();
  parameters_->motors.lattice_coop_Emax_bulk =
      motors["lattice_coop_Emax_bulk"].as<double>();
  parameters_->motors.interaction_energy =
      motors["interaction_energy"].as<double>();
  parameters_->motors.t_active = motors["t_active"].as<double>();
  parameters_->motors.k_on = motors["k_on"].as<double>();
  parameters_->motors.c_bulk = motors["c_bulk"].as<double>();
  parameters_->motors.c_eff_bind = motors["c_eff_bind"].as<double>();
  parameters_->motors.k_on_ATP = motors["k_on_ATP"].as<double>();
  parameters_->motors.c_ATP = motors["c_ATP"].as<double>();
  parameters_->motors.k_hydrolyze = motors["k_hydrolyze"].as<double>();
  parameters_->motors.k_off_i = motors["k_off_i"].as<double>();
  parameters_->motors.k_off_ii = motors["k_off_ii"].as<double>();
  parameters_->motors.applied_force = motors["applied_force"].as<double>();
  parameters_->motors.internal_force = motors["internal_force"].as<double>();
  parameters_->motors.sigma_off_i = motors["sigma_off_i"].as<double>();
  parameters_->motors.sigma_off_ii = motors["sigma_off_ii"].as<double>();
  parameters_->motors.sigma_ATP = motors["sigma_ATP"].as<double>();
  parameters_->motors.k_tether = motors["k_tether"].as<double>();
  parameters_->motors.c_eff_tether = motors["c_eff_tether"].as<double>();
  parameters_->motors.k_untether = motors["k_untether"].as<double>();
  parameters_->motors.r_0 = motors["r_0"].as<double>();
  parameters_->motors.k_spring = motors["k_spring"].as<double>();
  parameters_->motors.k_slack = motors["k_slack"].as<double>();
  parameters_->motors.endpausing_active =
      motors["endpausing_active"].as<bool>();
  parameters_->motors.tethers_active = motors["tethers_active"].as<bool>();
  /* Xlink parameters below */
  YAML::Node xlinks = input["xlinks"];
  parameters_->xlinks.k_on = xlinks["k_on"].as<double>();
  parameters_->xlinks.c_bulk = xlinks["c_bulk"].as<double>();
  parameters_->xlinks.c_eff_bind = xlinks["c_eff_bind"].as<double>();
  parameters_->xlinks.k_off_i = xlinks["k_off_i"].as<double>();
  parameters_->xlinks.k_off_ii = xlinks["k_off_ii"].as<double>();
  parameters_->xlinks.r_0 = xlinks["r_0"].as<double>();
  parameters_->xlinks.k_spring = xlinks["k_spring"].as<double>();
  parameters_->xlinks.diffu_coeff_i = xlinks["diffu_coeff_i"].as<double>();
  parameters_->xlinks.diffu_coeff_ii = xlinks["diffu_coeff_ii"].as<double>();
  parameters_->xlinks.interaction_energy =
      xlinks["interaction_energy"].as<double>();
  // Store params pointer as parameters_ in Curator
  unsigned long n_steps{parameters_->n_steps};
  double delta_t{parameters_->delta_t};
  Log("Reading params from %s:\n\n", param_file_);
  Log("  General simulation parameters:\n");
  Log("    seed = %lu\n", parameters_->seed);
  Log("    n_steps = %lu\n", parameters_->n_steps);
  Log("    n_datapoints = %i\n", parameters_->n_datapoints);
  Log("    data_threshold = %i steps\n", parameters_->data_threshold);
  Log("    delta_t = %g s\n", parameters_->delta_t);
  Log("    kbT = %g pN*nm\n", parameters_->kbT);
  Log("    eta = %g (pN*s)/um^2\n", parameters_->eta);
  Log("\n  Microtubule (mt) parameters:\n");
  Log("    count = %i\n", parameters_->microtubules.count);
  // Check to make sure there are enough vector entries for given MT count
  int n_lengths = input["microtubules"]["length"].size();
  int n_start_coords = input["microtubules"]["start_coord"].size();
  int n_immo = input["microtubules"]["immobile_until"].size();
  if (parameters_->microtubules.count > n_lengths or
      parameters_->microtubules.count > n_start_coords or
      parameters_->microtubules.count > n_immo) {
    Log("\nToo few parameters input for microtubules\n");
    ErrorExit("Curator::ParseParameters()");
  }
  for (int i_mt = 0; i_mt < n_lengths; i_mt++) {
    Log("    length = %i sites for mt %i\n",
        parameters_->microtubules.length[i_mt], i_mt);
  }
  Log("    y_dist = %g nm between MTs\n", parameters_->microtubules.y_dist);
  Log("    site_size = %g nm\n", parameters_->microtubules.site_size);
  Log("    radius = %g nm\n", parameters_->microtubules.radius);
  Log("    elevation = %g nm above surface\n",
      parameters_->microtubules.elevation);
  for (int i_mt = 0; i_mt < n_start_coords; i_mt++) {
    double start_coord = parameters_->microtubules.start_coord[i_mt];
    Log("    start_coord = %g sites for mt %i\n", start_coord, i_mt);
  }
  for (int i_mt = 0; i_mt < n_immo; i_mt++) {
    double immo = parameters_->microtubules.immobile_until[i_mt];
    Log("    immobile until = %g s for mt %i\n", immo, i_mt);
  }
  Log("    applied_force = %g pN\n", parameters_->microtubules.applied_force);
  Log("    printout_on = %s\n",
      parameters_->microtubules.printout_on ? "true" : "false");
  Log("    diffusion_on = %s\n",
      parameters_->microtubules.diffusion_on ? "true" : "false");
  Log("\n  Kinesin (motor) parameters:\n");
  // Log("    lattice_coop_alpha = %g\n",
  // parameters_->motors.lattice_coop_alpha);
  Log("    n_runs_desired = %zu\n", parameters_->motors.n_runs_desired);
  Log("    lattice_coop_range = %i\n", parameters_->motors.lattice_coop_range);
  Log("    lattice_coop_Emax_solo = -%g kbT\n",
      parameters_->motors.lattice_coop_Emax_solo);
  Log("    lattice_coop_Emax_bulk = -%g kbT\n",
      parameters_->motors.lattice_coop_Emax_bulk);
  Log("    interaction_energy = -%g kbT\n",
      parameters_->motors.interaction_energy);
  Log("    t_active = %g seconds\n", parameters_->motors.t_active);
  Log("    k_on = %g /(nM*s)\n", parameters_->motors.k_on);
  Log("    c_bulk = %g nM\n", parameters_->motors.c_bulk);
  Log("    c_eff_bind = %g nM\n", parameters_->motors.c_eff_bind);
  Log("    k_on_ATP = %g /(mM*s)\n", parameters_->motors.k_on_ATP);
  Log("    c_ATP = %g mM\n", parameters_->motors.c_ATP);
  Log("    k_hydrolyze = %g /s\n", parameters_->motors.k_hydrolyze);
  Log("    k_off_i = %g /s\n", parameters_->motors.k_off_i);
  Log("    k_off_ii = %g /s\n", parameters_->motors.k_off_ii);
  Log("    applied_force = %g pN\n", parameters_->motors.applied_force);
  Log("    internal_force = %g pN\n", parameters_->motors.internal_force);
  Log("    sigma_off_i = %g nm\n", parameters_->motors.sigma_off_i);
  Log("    sigma_off_ii = %g nm\n", parameters_->motors.sigma_off_ii);
  Log("    sigma_ATP = %g nm\n", parameters_->motors.sigma_ATP);
  Log("    k_tether = %g /(nM*s)\n", parameters_->motors.k_tether);
  Log("    c_eff_tether = %g nM\n", parameters_->motors.c_eff_tether);
  Log("    k_untether = %g /s\n", parameters_->motors.k_untether);
  Log("    r_0 = %g nm\n", parameters_->motors.r_0);
  Log("    k_spring = %g pN/nm\n", parameters_->motors.k_spring);
  Log("    k_slack = %g pN/nm\n", parameters_->motors.k_slack);
  Log("    tethers_active = %s\n",
      parameters_->motors.tethers_active ? "true" : "false");
  Log("    endpausing_active = %s\n",
      parameters_->motors.endpausing_active ? "true" : "false");
  Log("\n  Crosslinker (xlink) parameters:\n");
  Log("    k_on = %g /(nM*s)\n", parameters_->xlinks.k_on);
  Log("    c_bulk = %g nM\n", parameters_->xlinks.c_bulk);
  Log("    c_eff_bind = %g nM\n", parameters_->xlinks.c_eff_bind);
  Log("    k_off_i = %g /s\n", parameters_->xlinks.k_off_i);
  Log("    k_off_ii = %g /s\n", parameters_->xlinks.k_off_ii);
  Log("    r_0 = %g nm\n", parameters_->xlinks.r_0);
  Log("    k_spring = %g pN/nm\n", parameters_->xlinks.k_spring);
  Log("    diffu_coeff_i = %g um^2/s\n", parameters_->xlinks.diffu_coeff_i);
  Log("    diffu_coeff_ii = %g um^2/s\n", parameters_->xlinks.diffu_coeff_ii);
  Log("    interaction_energy = %g kbT\n",
      parameters_->xlinks.interaction_energy);
  Log("\nTotal simulation duration: %g seconds\n", delta_t * n_steps);
}

void Curator::SetLocalParameters() {

  size_t n_steps{parameters_->n_steps};
  int n_datapoints{parameters_->n_datapoints};
  data_threshold_ = parameters_->data_threshold;
  n_steps_recorded_ = n_steps - data_threshold_;
  n_steps_per_output_ = n_steps_recorded_ / n_datapoints;
  equil_milestone_ = data_threshold_ / 10;
  data_milestone_ = n_steps_recorded_ / 10;
  for (int i = 0; i < 4; i++) {
    t_motors_[i] = 0;
    t_xlinks_[i] = 0;
    t_MTs_[i] = 0;
  }
}

void Curator::InitializeSimObjects() {

  // Gsl: wrapper class for GSL library; manages random number generation
  properties_->gsl.Initialize(parameters_);
  // Microtubules: discretized 1-D lattice that proteins bind to
  properties_->microtubules.Initialize(parameters_, properties_);
  // Kinesin4: active motors that step towards plus-end of MTs
  properties_->kinesin4.Initialize(parameters_, properties_);
  // PRC1: passive crosslinkers that can doubly-bind to overlapping MTs
  properties_->prc1.Initialize(parameters_, properties_);
}

void Curator::OutputData() {

  int n_mts{parameters_->microtubules.count};
  int max_length{0};
  for (int i_mt{0}; i_mt < n_mts; i_mt++) {
    if (parameters_->microtubules.length[i_mt] > max_length)
      max_length = parameters_->microtubules.length[i_mt];
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
    int mt_length{parameters_->microtubules.length[i_mt]};
    Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
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
            tether_coords[i_site] = anchor_coord;
            double stalk_coord = motor->GetStalkCoordinate();
            double teth_dist = abs(anchor_coord - stalk_coord);
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
    if (properties_->kinesin4.step_active_ < parameters_->n_steps) {
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
  if (properties_->kinesin4.step_active_ > parameters_->n_steps) {
    if (properties_->kinesin4.tethering_active_) {
      fwrite(motor_forces, sizeof(double), n_mts,
             properties_->motor_force_file_);
      fwrite(motor_extensions, sizeof(int), 2 * motor_cutoff + 1,
             properties_->motor_extension_file_);
    }
  }
  if (properties_->prc1.population_active_) {
    fwrite(xlink_forces, sizeof(double), n_mts, properties_->xlink_force_file_);
    if (properties_->prc1.crosslinking_active_) {
      fwrite(xlink_extensions, sizeof(int), xlink_cutoff + 1,
             properties_->xlink_extension_file_);
    }
  }
  if (parameters_->microtubules.diffusion_on) {
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
  /*
  Log("   -Motors: %.2f\n", t_motors_[0] / n_per_sec);
  Log("      -Calculating stats: %.2f\n", t_motors_[1] / n_per_sec);
  Log("      -Constructing list: %.2f\n", t_motors_[2] / n_per_sec);
  Log("      -Execution: %.2f\n", t_motors_[3] / n_per_sec);
  Log("   -Xlinks: %.2f\n", t_xlinks_[0] / n_per_sec);
  Log("      -Updating lists: %.2f\n", t_xlinks_[1] / n_per_sec);
  Log("           -Extensions: %.2f\n", t_xlinks_[5] / n_per_sec);
  Log("           -MT unnocupied sites: %.2f\n", t_xlinks_[6] / n_per_sec);
  Log("           -Bound_I: %.2f\n", t_xlinks_[7] / n_per_sec);
  Log("      -Refreshing populations: %.2f\n", t_xlinks_[2] / n_per_sec);
  Log("      -Generating sequence: %.2f\n", t_xlinks_[3] / n_per_sec);
  Log("      -Executing events: %.2f\n", t_xlinks_[4] / n_per_sec);
  Log("   -MTs: %.2f\n", t_MTs_[0] / n_per_sec);
  Log("      -Summing forces: %.2f\n", t_MTs_[1] / n_per_sec);
  Log("      -Calculating displacement: %.2f\n", t_MTs_[2] / n_per_sec);
  Log("      -Execution: %.2f\n", t_MTs_[3] / n_per_sec);
  */
}

void Curator::CloseDataFiles() {

  fclose(log_file_);
  fclose(properties_->occupancy_file_);
  if (properties_->kinesin4.step_active_ < parameters_->n_steps) {
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
  if (parameters_->microtubules.diffusion_on) {
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

void Curator::UpdateTimestep() {

  // Record current step number and then increment it by one
  size_t i_step{properties_->current_step_++};
  // Start global system clock on very first step
  if (i_step == 0) {
    start_ = sys_clock::now();
  }
  if (properties_->sim_equilibrating_) {
    if (i_step >= data_threshold_) {
      properties_->sim_equilibrating_ = false;
    }
    if (i_step % equil_milestone_ == 0) {
      int p{int(i_step / equil_milestone_) * 10};
      Log("Sim '%s' pre-equilibration is %lu%% complete. (step %lu)\n",
          sim_name_, p, i_step);
    }
    return;
  }
  if (!properties_->kinesin4.equilibrated_) {
    return;
  }
  size_t steps_past_threshold{i_step - data_threshold_};
  // Give updates on status of data collection (every 10 percent)
  if (steps_past_threshold % data_milestone_ == 0) {
    unsigned long p{(steps_past_threshold / data_milestone_) * 10};
    Log("Sim '%s' data collection is %lu%% complete. (step # %lu)\n", sim_name_,
        p, i_step);
  }
  if (i_step >= parameters_->n_steps) {
    properties_->sim_running_ = false;
    return;
  }
  // Collect data every n_pickup timesteps
  if (steps_past_threshold % n_steps_per_output_ == 0) {
    OutputData();
    n_datapoints_recorded_++;
  }
  /*
  if (parameters_->microtubules.printout_on and i_step % 1000 == 0) {
    PrintMicrotubules(0);
  }
  */
}

void Curator::PrintMicrotubules() {

  /*
  int n_mts = parameters_->microtubules.count;
  // Figure out which MT is the farthest left
  int leftmost_coord = 0;
  for (int i_mt = 0; i_mt < n_mts; i_mt++) {
    int mt_coord = properties_->microtubules.mt_list_[i_mt].coord_;
    if (i_mt == 0)
      leftmost_coord = mt_coord;
    else if (mt_coord < leftmost_coord)
      leftmost_coord = mt_coord;
  }
  // Print out MTs
  for (int i_mt = n_mts - 1; i_mt >= 0; i_mt--) {
    int mt_length = parameters_->microtubules.length[i_mt];
    Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
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
     int mt1_coord = properties_->microtubules.mt_list_[0].coord_;
     int mt2_coord = properties_->microtubules.mt_list_[1].coord_;
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

  size_t n_steps{parameters_->n_steps};
  int n_datapoints{parameters_->n_datapoints};
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