#include "master_header.h"
#include "curator.h"

Curator::Curator(){
}

void Curator::InitializeSimulation(char *exe_name, char *param_file, 
		char *sim_name, int argc, system_properties *properties, 
		system_parameters *parameters){

	properties_ = properties;
	// Initiate log file -- records all output of the sim
	OpenLog(sim_name);
	// Check that input has correct number of arguments
	CheckArguments(exe_name, argc);
	// Parse parameters from input YAML file to local parameter structure
	ParseParameters(parameters, param_file);
	// Set local parameters
	SetParameters();
	// Initialize MTs, kinesin, etc. 
	SetExperimentalStage();
	// Generate files which data will be written to 
	GenerateDataFiles(sim_name);
}

void Curator::OpenLog(char *sim_name){

	char log_file[160];
	sprintf(log_file, "%s.log", sim_name);
	// Check to see if sim files already exist
	if(FileExists(log_file)){
		printf("Simulation log file with this name already exists!\n");
		printf("Do you wish to overwrite these data? y/n\n");
		std::string response; 
		bool response_unacceptable = true;
		while(response_unacceptable){
			std::getline(std::cin, response);
			if(response == "n"){
				printf("Simulation terminated.\n");
				exit(1); 
			}
			else if(response == "y"){
				printf("Very well.");
			   	printf(" Overwriting data for sim '%s'\n\n", sim_name);
				response_unacceptable = false; 
			}
			else{
				printf("ayo I said y or n. try again plz\n");
			}
		}
	}
	properties_->log_file_ = OpenFile(log_file, "w");
}

void Curator::CheckArguments(char *exe_name, int argc){

	if(argc != 3){
		Log("\nWrong number of command-line arguments in main\n");
		Log("Usage: %s parameters.yaml sim_name\n\n", exe_name);
		exit(1);
	}
}

void Curator::ParseParameters(system_parameters *params, char *param_file){

	// Check to make sure param file actually exists
	if(!FileExists(param_file)){
		Log("  Error: parameter file does not exist; aborting\n"); 
		exit(1);
	}
	// Parse parameter file into a YAML node
	YAML::Node input = YAML::LoadFile(param_file);
	// Transfer values from input param node to system_parameters structure
	try{
		params->seed = input["seed"].as<long>();
	}
	catch(const YAML::BadConversion error){
		params->seed = (long)(input["seed"].as<double>());
	}   
	params->n_steps = input["n_steps"].as<int>();
	params->n_datapoints = input["n_datapoints"].as<int>();
	params->data_threshold = input["data_threshold"].as<int>();
	params->delta_t = input["delta_t"].as<double>();
	params->kbT = input["kbT"].as<double>();
	params->eta = input["eta"].as<double>();
	/* Motor parameters below */ 
	YAML::Node motors = input["motors"];
	params->motors.t_active = motors["t_active"].as<double>();
	params->motors.k_on = motors["k_on"].as<double>();
	params->motors.c_bulk = motors["c_bulk"].as<double>();
	params->motors.c_eff_bind = motors["c_eff_bind"].as<double>();
	params->motors.k_on_ATP = motors["k_on_ATP"].as<double>();
	params->motors.c_ATP = motors["c_ATP"].as<double>();
	params->motors.k_hydrolyze = motors["k_hydrolyze"].as<double>();
	params->motors.k_hydrolyze_stalled = motors["k_hydrolyze_stalled"].
		as<double>();
	params->motors.k_off_i = motors["k_off_i"].as<double>();
	params->motors.k_off_i_stalled = motors["k_off_i_stalled"].as<double>();
	params->motors.k_off_ii = motors["k_off_ii"].as<double>();
	params->motors.endpausing_active 
		= motors["endpausing_active"].as<bool>();
	params->motors.tethers_active = motors["tethers_active"].as<bool>();
	params->motors.k_tether = motors["k_tether"].as<double>();
	params->motors.c_eff_tether = motors["c_eff_tether"].as<double>();
	params->motors.k_untether = motors["k_untether"].as<double>();
	params->motors.r_0 = motors["r_0"].as<double>();
	params->motors.k_spring = motors["k_spring"].as<double>();
	params->motors.k_slack = motors["k_slack"].as<double>();
	params->motors.stall_force = motors["stall_force"].as<double>();
	params->motors.applied_force = motors["applied_force"].as<double>();
	/* Xlink parameters below */
	YAML::Node xlinks = input["xlinks"];
	params->xlinks.k_on = xlinks["k_on"].as<double>();
	params->xlinks.concentration = xlinks["concentration"].as<double>();
	params->xlinks.conc_eff_bind = xlinks["conc_eff_bind"].as<double>();
	params->xlinks.k_off_i = xlinks["k_off_i"].as<double>();
	params->xlinks.k_off_ii = xlinks["k_off_ii"].as<double>();
	params->xlinks.diffusion_const_i = 
			xlinks["diffusion_const_i"].as<double>();
	params->xlinks.diffusion_const_ii = 
			xlinks["diffusion_const_ii"].as<double>();
	params->xlinks.r_0 = xlinks["r_0"].as<double>();
	params->xlinks.k_spring = xlinks["k_spring"].as<double>();
	/* Microtubule parameters below */ 
	YAML::Node mts = input["microtubules"];
	params->microtubules.count = mts["count"].as<int>();
	params->microtubules.length = mts["length"].as<std::vector<int>>();
	params->microtubules.y_dist = mts["y_dist"].as<double>(); 
	params->microtubules.site_size = mts["site_size"].as<double>();
	params->microtubules.radius = mts["radius"].as<double>();
	params->microtubules.elevation = mts["elevation"].as<double>();
	params->microtubules.start_coord = 
			mts["start_coord"].as<std::vector<double>>();
	params->microtubules.imposed_velocity =
			mts["imposed_velocity"].as<std::vector<double>>();
	params->microtubules.immobile_until = 
			mts["immobile_until"].as<std::vector<double>>();
	// Check to make sure there are enough vector entries for given MT count
	int n_lengths = mts["length"].size();
	int n_start_coords = mts["start_coord"].size();
	int n_imp_vel = mts["imposed_velocity"].size();
	int n_immo = mts["immobile_until"].size();
	if(params->microtubules.count > n_lengths
	|| params->microtubules.count > n_start_coords
	|| params->microtubules.count > n_imp_vel
	|| params->microtubules.count > n_immo){
		Log("\nError! More MTs than given parameters; ");
		Log("check vector entries in parameter file!\n\n");
		exit(1);
	}
	params->microtubules.printout_on = mts["printout_on"].as<bool>();
	params->microtubules.diffusion_on = mts["diffusion_on"].as<bool>();
	// Store params pointer as parameters_ in Curator
	parameters_ = params;
	int n_steps = parameters_->n_steps;
	double delta_t = parameters_->delta_t;
	Log("Reading params from %s:\n\n", param_file);
	Log("  General simulation parameters:\n");
	Log("    seed = %li\n", params->seed);
	Log("    n_steps = %i\n", params->n_steps);
	Log("    n_datapoints = %i\n", params->n_datapoints);
	Log("    data_threshold = %i steps\n", params->data_threshold);
	Log("    delta_t = %g s\n", params->delta_t);
	Log("    kbT = %g pN*nm\n", params->kbT);
	Log("    eta = %g (pN*s)/um^2\n", params->eta);
	Log("\n  Kinesin (motor) parameters:\n");
	Log("    t_active = %g seconds\n", params->motors.t_active);
	Log("    k_on = %g /(nM*s)\n", params->motors.k_on);
	Log("    c_bulk = %g nM\n", params->motors.c_bulk);
	Log("    c_eff_bind = %g nM\n", params->motors.c_eff_bind);
	Log("    k_on_ATP = %g /(mM*s)\n", params->motors.k_on_ATP);
	Log("    c_ATP = %g mM\n", params->motors.c_ATP);
	Log("    k_hydrolyze = %g /s\n", params->motors.k_hydrolyze);
	Log("    k_hydrolyze_stalled = %g /s\n", 
			params->motors.k_hydrolyze_stalled);
	Log("    k_off_i = %g /s\n", params->motors.k_off_i);
	Log("    k_off_i_stalled = %g /s\n", params->motors.k_off_i_stalled); 
	Log("    k_off_ii = %g /s\n", params->motors.k_off_ii);
	Log("    k_tether = %g /(nM*s)\n", params->motors.k_tether);
	Log("    c_eff_tether = %g nM\n", params->motors.c_eff_tether);
	Log("    k_untether = %g /s\n", params->motors.k_untether);
	Log("    r_0 = %g nm\n", params->motors.r_0);
	Log("    k_spring = %g pN/nm\n", params->motors.k_spring);
	Log("    k_slack = %g pN/nm\n", params->motors.k_slack);
	Log("    stall_force = %g pN\n", params->motors.stall_force);
	Log("    applied_force = %g pN\n", params->motors.applied_force);
	Log("    tethers_active = %s\n", params->motors.tethers_active? 
			"true" : "false");	
	Log("    endpausing_active = %s\n", params->motors.endpausing_active? 
			"true" : "false");	
	Log("\n  Crosslinker (xlink) parameters:\n");
	Log("    k_on = %g /(nM*s)\n", params->xlinks.k_on);
	Log("    concentration = %g nM\n", params->xlinks.concentration);
	Log("    conc_eff_bind = %g nM\n", params->xlinks.conc_eff_bind);
	Log("    k_off_i = %g /s\n", params->xlinks.k_off_i);
	Log("    k_off_ii = %g /s\n", params->xlinks.k_off_ii);
	Log("    diffusion_constant_i = %g um^2/s\n", 
			params->xlinks.diffusion_const_i);
	Log("    diffusion_constant_ii = %g um^2/s\n", 
			params->xlinks.diffusion_const_ii);
	Log("    r_0 = %g nm\n", params->xlinks.r_0);
	Log("    k_spring = %g pN/nm\n", params->xlinks.k_spring);
	Log("\n  Microtubule (mt) parameters:\n");
	Log("    count = %i\n", params->microtubules.count);
	for(int i_mt = 0; i_mt < n_lengths; i_mt++){
		Log("    length = %i sites for mt %i\n", 
				params->microtubules.length[i_mt], i_mt);
	}
	Log("    y_dist = %g nm between MTs\n", params->microtubules.y_dist);
	Log("    site_size = %g nm\n", params->microtubules.site_size);
	Log("    radius = %g nm\n", params->microtubules.radius);
	Log("    elevation = %g nm above surface\n", 
			params->microtubules.elevation);
	for(int i_mt = 0; i_mt < n_start_coords; i_mt++){
		double start_coord = params->microtubules.start_coord[i_mt];
		Log("    start_coord = %g sites for mt %i\n", start_coord, i_mt);
	}
	for(int i_mt = 0; i_mt < n_imp_vel; i_mt++){
		double imp_vel = params->microtubules.imposed_velocity[i_mt];
		Log("    imposed_velocity = %g nm/s for mt %i\n", imp_vel, i_mt);
	}
	for(int i_mt = 0; i_mt < n_immo; i_mt++){
		double immo = params->microtubules.immobile_until[i_mt];	
		Log("    immobile until = %g s for mt %i\n", immo, i_mt);
	}
	Log("    printout_on = %s\n", params->microtubules.printout_on? 
			"true":"false");
	Log("    diffusion_on = %s\n", params->microtubules.diffusion_on? 
			"true":"false");
	Log("\nTotal simulation duration: %g seconds\n", delta_t*n_steps);
}

void Curator::SetParameters(){

	int n_steps = parameters_->n_steps;
	int n_datapoints = parameters_->n_datapoints;
	data_threshold_ = parameters_->data_threshold;
	range_of_data_ = n_steps - data_threshold_;
	n_pickup_ = range_of_data_/n_datapoints;
	equil_milestone_ = data_threshold_/10;
	data_milestone_ = range_of_data_/10;

	sim_duration_ = 0;
	for(int i = 0; i < 4; i++){
		t_motors_[i] = 0;
		t_xlinks_kmc_[i] = 0;
		t_xlinks_dif_[i] = 0;
		t_MTs_[i] = 0; 
	}
}

void Curator::SetExperimentalStage(){

	// Initialize microtubules, kinesin4, and prc1 classes 
	properties_->microtubules.Initialize(parameters_, properties_);
	properties_->kinesin4.Initialize(parameters_, properties_); 
	properties_->prc1.Initialize(parameters_, properties_);
	// Initialize the general science library (gsl) class; 
	// just an easy way of sampling distributions and referencing the RNG
	properties_->gsl.Initialize(parameters_, properties_);
}

void Curator::GenerateDataFiles(char* sim_name){

	char occupancy_file[160], 
		 motor_ID_file[160], xlink_ID_file[160], 
		 tether_coord_file[160], mt_coord_file[160], 
		 motor_extension_file[160], xlink_extension_file[160],
		 motor_force_file[160], xlink_force_file[160], total_force_file[160],
		 motor_head_status_file[160];
	// Generate names of output files based on the input simulation name
	sprintf(occupancy_file, "%s_occupancy.file", sim_name);
	sprintf(motor_ID_file, "%s_motorID.file", sim_name);
	sprintf(xlink_ID_file, "%s_xlinkID.file", sim_name);	
	sprintf(tether_coord_file, "%s_tether_coord.file", sim_name);
	sprintf(mt_coord_file, "%s_mt_coord.file", sim_name);
	sprintf(motor_extension_file, "%s_motor_extension.file", sim_name);
	sprintf(xlink_extension_file, "%s_xlink_extension.file", sim_name);
	sprintf(motor_force_file, "%s_motor_force.file", sim_name);
	sprintf(xlink_force_file, "%s_xlink_force.file", sim_name);
	sprintf(total_force_file, "%s_total_force.file", sim_name);
	sprintf(motor_head_status_file, "%s_motor_head_status.file", sim_name);
	// Open occupancy file, which stores the species ID of each occupant 
	// (or -1 for none) for all MT sites during data collection (DC)
	properties_->occupancy_file_ = OpenFile(occupancy_file, "w");
	// Open motor ID file, which stores the unique ID of all bound motors 
	// (unbound not tracked) and their respective site indices during DC
	properties_->motor_ID_file_ = OpenFile(motor_ID_file, "w");
	// Open xlink ID file, which does the same 
	// as the motor ID file but for xlinks
	properties_->xlink_ID_file_ = OpenFile(xlink_ID_file, "w");
	// Open tether coord file, which stores the coordinates 
	// of the anchor points of tethered motors
	properties_->tether_coord_file_ = OpenFile(tether_coord_file, "w");
	// Open mt coord file, which stores the coordinates 
	// of the left-most edge of each microtubule during DC
	properties_->mt_coord_file_ = OpenFile(mt_coord_file, "w");
	// Open motor extension file, which stores the number of motors 
	// with a certain tether extension for all possible extensions
	properties_->motor_extension_file_ = OpenFile(motor_extension_file, "w");
	// Open xlink extension file, which stores the number of stage-2 
	// xlinks at a certain extension for all possible extensions
	properties_->xlink_extension_file_ = OpenFile(xlink_extension_file, "w");
	// Open motor force file, which stores the sum 
	// of forces coming from motor tether extensions
	properties_->motor_force_file_ = OpenFile(motor_force_file, "w");
	// Open xlink force file, which stores the sum
	// of forces coming from xlink extensions
	properties_->xlink_force_file_ = OpenFile(xlink_force_file, "w");
	// Open total force file, which stores the sum of ALL 
	// forces coming from xlink and motor tether extensions
	properties_->total_force_file_ = OpenFile(total_force_file, "w");
	// bool; simply says if motor head is trailing or not
	properties_->motor_head_status_file_ 
		= OpenFile(motor_head_status_file, "w");
}

bool Curator::FileExists(std::string file_name){

	struct stat buffer;
	return (stat(file_name.c_str(), &buffer) != -1);
}

FILE* Curator::OpenFile(const char *file_name, const char *type){

    FILE *file_ptr;
    if ((file_ptr = fopen(file_name, type)) == NULL) {
		fprintf(stderr, "Cannot open %s\n", file_name);
		exit(1);
    }
    return file_ptr;
}

void Curator::Log(const char *msg){

	int l_return = printf("%s", msg);
	int l_freturn = fprintf(properties_->log_file_, "%s", msg);
	if(l_return < 0 || l_freturn < 0){
		Log("Error in Curator::Log_base brobro\n");
		exit(1);
	}
}

template <typename T> void Curator::Log(const char *msg, T arg){

	int l_return = printf(msg, arg);
	int l_freturn = fprintf(properties_->log_file_, msg, arg);
	if(l_return < 0 || l_freturn < 0){
		Log("Error in Curator::Log_arg brobro\n");
		exit(1);
	}
}

template <typename T1, typename T2> void Curator::Log(const char *msg, 
		T1 arg1, T2 arg2){

	int l_return = printf(msg, arg1, arg2);
	int l_freturn = fprintf(properties_->log_file_, msg, arg1, arg2);
	if(l_return < 0 || l_freturn < 0){
		Log("Error in Curator::Log_2arg brobro\n");
		exit(1);
	}
}

void Curator::UpdateTimestep(int i_step){

	properties_->current_step_ = i_step;
	if(i_step == 0) start_ = sys_clock::now(); 
	// Give updates on equilibrium process (every 10%)
	if(i_step < data_threshold_ && i_step % equil_milestone_ == 0)
		Log("Equilibration is %i percent complete (step # %i)\n", 
				(int)(i_step/equil_milestone_)*10, i_step);
	// Start data collection at appropriate step threshold
	else if(i_step >= data_threshold_){
		int delta = i_step - data_threshold_;
		// Collect data every n_pickup timesteps
		if(delta % n_pickup_ == 0) OutputData();
		// Give updates on status of simulation every 10% 
		if(delta % data_milestone_ == 0)
			Log("Data collection is %i percent complete (step # %i)\n",
					(int)(delta / data_milestone_)*10, i_step);
		// Announce when simulation is done 
		else if(delta == range_of_data_ - 1) Log("Done!\n");
	}
	if(parameters_->microtubules.printout_on && i_step % 1000 == 0)
		PrintMicrotubules(0);
}

void Curator::PrintMicrotubules(){

	int n_mts = parameters_->microtubules.count;
	// Figure out which MT is the farthest left 
	int leftmost_coord = 0;
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		int mt_coord = properties_->microtubules.mt_list_[i_mt].coord_;
		if(i_mt == 0)
			leftmost_coord = mt_coord;	
		else if(mt_coord < leftmost_coord)
			leftmost_coord = mt_coord;
	}
	// Print out MTs
	for(int i_mt = n_mts - 1; i_mt >= 0; i_mt--){
		int mt_length = parameters_->microtubules.length[i_mt]; 
		Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
		int mt_coord = mt->coord_;
		int delta = mt_coord - leftmost_coord;
		if(delta < 0){
			Log("wat\n");
			exit(1);
		}
		for(int i_entry = 0; i_entry < delta; i_entry++){
			printf(" ");
		}
		for(int i_site = 0; i_site < mt_length; i_site++){
			Tubulin *site = &mt->lattice_[i_site];
			if(site->occupied_ == false)
				printf("=");
			else if(site->xlink_ != nullptr){
				AssociatedProtein *xlink = site->xlink_;
				if(xlink->heads_active_ == 1){
					if(xlink->tethered_ == false)
						printf("i");
					else if(xlink->motor_->heads_active_ == 0)
						printf("F");
					else
						printf("I");
				}
				else if(xlink->heads_active_ == 2)
					if(xlink->tethered_ == false){
						printf("x");			
//						printf("%i", xlink->x_dist_);
					}	
					else
						printf("X");
				else{
					printf("no sunny. look in wallace's print\n");
					exit(1);
				}
			}
			else if(site->motor_head_ != nullptr){
				Kinesin *motor = site->motor_head_->motor_; 
				if(motor->heads_active_ == 1){
					if(motor->tethered_ == false)
						std::cout << site->motor_head_->ligand_;
					//	printf("m");
					else{
						printf("[");
						std::cout << site->motor_head_->ligand_;
						printf("]");
					}
				}
				else if(motor->heads_active_ == 2){
					int i_front = motor->head_one_.site_->index_;
					int i_rear = motor->head_two_.site_->index_;
					if(i_front > i_rear){
						if(i_site == i_rear){
							if(motor->tethered_ == false){
								printf("(");
								std::cout << motor->head_two_.ligand_;
								printf("|");
							}
							else
								printf("%i", motor->x_dist_doubled_ / 10);
							//	printf("[");
						}
						if(i_site == i_front){
							if(motor->tethered_ == false){
								std::cout << motor->head_one_.ligand_;
								printf(")");
							}
							else
//								printf("]");
								printf("%i", motor->x_dist_doubled_ % 10);
						}
					} 
					else if(i_front < i_rear){
						if(i_site == i_front){	
							if(motor->tethered_ == false){
								printf("(");
								std::cout << motor->head_one_.ligand_;
								printf("|");
							}
							else
								printf("%i", motor->x_dist_doubled_ / 10);
							//	printf("[");
						}
						if(i_site == i_rear){
							if(motor->tethered_ == false){
								std::cout << motor->head_two_.ligand_;
								printf(")");
							}
							else
						//		printf("]");
								printf("%i", motor->x_dist_doubled_ % 10);
						}
					}
					else{
						printf("error in print MT\n");
						exit(1);
					}
				}
			}
		}
		printf(" %i\n", mt->polarity_);
	}   
	printf("\n");

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

void Curator::PrintMicrotubules(double pause_duration){

	PrintMicrotubules();
	if(pause_duration > 0) PauseSim(pause_duration);
}

void Curator::PauseSim(double duration){

	// Duration should be input in seconds
	pause_dur_.tv_sec = (int) duration;	 
	pause_dur_.tv_nsec = (duration - (int)duration)*100000000;
	nanosleep(&pause_dur_, NULL);
}

void Curator::OutputData(){

	int n_mts = parameters_->microtubules.count; 
	int max_length = 0;
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		if(parameters_->microtubules.length[i_mt] > max_length)
			max_length = parameters_->microtubules.length[i_mt];
	}
	// Get file pointers from system properties
	FILE *occupancy_file = properties_->occupancy_file_;
	FILE *motor_ID_file = properties_->motor_ID_file_;
	FILE *xlink_ID_file = properties_->xlink_ID_file_;
	FILE *tether_coord_file = properties_->tether_coord_file_;
	FILE *mt_coord_file = properties_->mt_coord_file_;
	FILE *motor_extension_file = properties_->motor_extension_file_;
	FILE *xlink_extension_file = properties_->xlink_extension_file_;
	FILE *motor_force_file = properties_->motor_force_file_;
	FILE *xlink_force_file = properties_->xlink_force_file_;
	FILE *total_force_file = properties_->total_force_file_;
	FILE *motor_head_status_file = properties_->motor_head_status_file_;
	// Create arrays to store data; ptrs to write it to file 
	double mt_coord_array[n_mts];
	double *mt_coord_ptr = mt_coord_array;
	// For extension statistics, data is on a per-extension basis
	int motor_ext_cutoff = properties_->kinesin4.dist_cutoff_;
	int motor_extension_array[2*motor_ext_cutoff + 1];
	int *motor_extension_ptr = motor_extension_array; 
	int xlink_ext_cutoff = properties_->prc1.dist_cutoff_; 
	int xlink_extension_array[xlink_ext_cutoff + 1];
	int *xlink_extension_ptr = xlink_extension_array;
	// Back to normal per-MT array format
	double motor_force_array[n_mts];
	double *motor_force_ptr = motor_force_array;
	double xlink_force_array[n_mts];
	double *xlink_force_ptr = xlink_force_array;
	double total_force_array[n_mts];
	double *total_force_ptr = total_force_array;	
	// Run through all MTs and get data for each
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		int mt_length = parameters_->microtubules.length[i_mt]; 
		Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
		// Create arrays & ptrs for intraMT data 
		int motor_ID_array[max_length],
		xlink_ID_array[max_length], 
		occupancy_array[max_length];
		int	*motor_ID_ptr = motor_ID_array, 
			*xlink_ID_ptr = xlink_ID_array, 
			*occupancy_ptr = occupancy_array;
		double tether_coord_array[max_length];
		double *teth_coord_ptr = tether_coord_array;
		bool motor_head_status_array[max_length];
		for(int i_site = 0; i_site < max_length; i_site++){   
			motor_head_status_array[i_site] = false;
		};
		bool *motor_head_status_ptr = motor_head_status_array;
		// Run through all sites on this particular MT
		for(int i_site = 0; i_site < mt_length; i_site++){
			Tubulin *site = &mt->lattice_[i_site];
			// If unoccupied, store the speciesID of tubulin to occupancy
			// and an ID of -1 (null) to motor/xlink ID files 
			if(site->occupied_ == false){
				occupancy_array[i_site] = site->speciesID_;
				motor_ID_array[i_site] = -1;
				xlink_ID_array[i_site] = -1;
				tether_coord_array[i_site] = -1;
			}
			// If occupied by xlink, store its species ID to occupancy_file,
			// its unique ID to the xlink ID file, and -1 to motor ID file
			else if(site->xlink_ != nullptr){
				occupancy_array[i_site] = site->xlink_->speciesID_;
				motor_ID_array[i_site] = -1;
				xlink_ID_array[i_site] = site->xlink_->ID_;
				tether_coord_array[i_site] = -1;
			}
			// If occupied by motor, store its species ID to occupancy_file, 
			// its unique ID to the motor ID file, and -1 to xlink ID file
			else if(site->motor_head_ != nullptr){
				Kinesin *motor = site->motor_head_->motor_;
				occupancy_array[i_site] = motor->speciesID_;
				motor_ID_array[i_site] = motor->ID_;
				xlink_ID_array[i_site] = -1;
				motor_head_status_array[i_site] 
					= site->motor_head_->trailing_;
				if(motor->tethered_ == true){
					if(motor->xlink_->heads_active_ > 0){
						AssociatedProtein* xlink = motor->xlink_;
						double anchor_coord = xlink->GetAnchorCoordinate();
						tether_coord_array[i_site] = anchor_coord; 
						double stalk_coord = motor->GetStalkCoordinate();
						double teth_dist = abs(anchor_coord - stalk_coord);
						if(teth_dist > motor->dist_cutoff_){
							Log("woah, teth dist is %g\n", teth_dist);
						}
					}
					else tether_coord_array[i_site] = -1;
				}
				else tether_coord_array[i_site] = -1;
			}
		}
		// Pad written data with zeros for shorter MTs
		for(int i_site = mt_length; i_site < max_length; i_site++){
				occupancy_array[i_site] = -1;
				motor_ID_array[i_site] = -1;
				xlink_ID_array[i_site] = -1;
				tether_coord_array[i_site] = -1;
		}
		mt_coord_array[i_mt] = mt->coord_; 
		motor_force_array[i_mt] = mt->GetNetForce_Motors();
		xlink_force_array[i_mt] = mt->GetNetForce_Xlinks();
		total_force_array[i_mt] = mt->GetNetForce();
		// Write the data to respective files one microtubule at a time
		fwrite(occupancy_ptr, sizeof(int), max_length, occupancy_file);
		fwrite(motor_ID_ptr, sizeof(int), max_length, motor_ID_file);
		fwrite(xlink_ID_ptr, sizeof(int), max_length, xlink_ID_file);
		fwrite(teth_coord_ptr, sizeof(double), max_length, 
				tether_coord_file);
		fwrite(motor_head_status_ptr, sizeof(bool), max_length, 
				motor_head_status_file);
	}	
	// Scan through kinesin4/prc1 statistics to get extension occupancies 
	for(int i_ext = 0; i_ext <= 2*motor_ext_cutoff; i_ext++){
		KinesinManagement *kinesin4 = &properties_->kinesin4; 
		motor_extension_array[i_ext] = kinesin4->n_bound_tethered_[i_ext];
	}
	for(int i_ext = 0; i_ext <= xlink_ext_cutoff; i_ext++){
		AssociatedProteinManagement *prc1 = &properties_->prc1; 
		int max = prc1->max_neighbs_;
		xlink_extension_array[i_ext] = prc1->n_sites_ii_[max+1][i_ext]; 
	}
	// Write the data to respective files one timestep at a time 
	fwrite(mt_coord_ptr, sizeof(double), n_mts, mt_coord_file);
	fwrite(motor_force_ptr, sizeof(double), n_mts, motor_force_file);
	fwrite(xlink_force_ptr, sizeof(double), n_mts, xlink_force_file);
	fwrite(total_force_ptr, sizeof(double), n_mts, total_force_file);
	fwrite(motor_extension_ptr, sizeof(int), 2*motor_ext_cutoff + 1, 
			motor_extension_file);
	fwrite(xlink_extension_ptr, sizeof(int), xlink_ext_cutoff + 1, 
			xlink_extension_file);
	fwrite(motor_force_ptr, sizeof(double), n_mts, motor_force_file);
	fwrite(xlink_force_ptr, sizeof(double), n_mts, xlink_force_file);
	fwrite(total_force_ptr, sizeof(double), n_mts, total_force_file);
}

void Curator::OutputSimDuration(){

	finish_ = sys_clock::now();
	auto elapsed = std::chrono::duration_cast<t_unit>(finish_ - start_);
	sim_duration_ = elapsed.count();
	Log("\nTime to execute sim: %.2f seconds.\n", sim_duration_/n_per_sec_);
	Log("   -Motors: %.2f\n", t_motors_[0]/n_per_sec_);
	Log("      -Calculating stats: %.2f\n", t_motors_[1]/n_per_sec_);
	Log("      -Constructing list: %.2f\n", t_motors_[2]/n_per_sec_);
	Log("      -Execution: %.2f\n", t_motors_[3]/n_per_sec_);
	Log("   -Xlinks (KMC): %.2f\n", t_xlinks_kmc_[0]/n_per_sec_);
	Log("      -Calculating stats: %.2f\n", t_xlinks_kmc_[1]/n_per_sec_);
	Log("      -Constructing list: %.2f\n", t_xlinks_kmc_[2]/n_per_sec_);
	Log("      -Execution: %.2f\n", t_xlinks_kmc_[3]/n_per_sec_);
	Log("   -Xlinks (Diffusion): %.2f\n", t_xlinks_dif_[0]/n_per_sec_);
	Log("      -Calculating stats: %.2f\n", t_xlinks_dif_[1]/n_per_sec_);
	Log("      -Constructing list: %.2f\n", t_xlinks_dif_[2]/n_per_sec_);
	Log("      -Execution: %.2f\n", t_xlinks_dif_[3]/n_per_sec_);
	Log("   -MTs: %.2f\n", t_MTs_[0]/n_per_sec_);
	Log("      -Summing forces: %.2f\n", t_MTs_[1]/n_per_sec_);
	Log("      -Calculating displacement: %.2f\n", t_MTs_[2]/n_per_sec_);
	Log("      -Execution: %.2f\n", t_MTs_[3]/n_per_sec_);
}

void Curator::CloseDataFiles(){

	fclose(properties_->log_file_);
	fclose(properties_->occupancy_file_);
	fclose(properties_->motor_ID_file_);
	fclose(properties_->xlink_ID_file_);
	fclose(properties_->tether_coord_file_);
	fclose(properties_->mt_coord_file_);
	fclose(properties_->motor_extension_file_);
	fclose(properties_->xlink_extension_file_);
	fclose(properties_->motor_force_file_);
	fclose(properties_->xlink_force_file_);
	fclose(properties_->total_force_file_);
	fclose(properties_->motor_head_status_file_);
}

void Curator::CleanUp(){

	OutputSimDuration();
	CloseDataFiles();
	properties_->gsl.CleanUp();
}
