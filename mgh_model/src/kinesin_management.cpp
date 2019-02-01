#include "master_header.h"
#include "kinesin_management.h"

KinesinManagement::KinesinManagement(){
}

void KinesinManagement::Initialize(system_parameters *parameters, 
        system_properties *properties){

    parameters_ = parameters;
    properties_ = properties;

    GenerateMotors();	
    SetParameters();	
	InitializeLists();
	InitializeEvents();
}

void KinesinManagement::GenerateMotors(){

    int n_mts = parameters_->microtubules.count;
    int n_sites = parameters_->microtubules.length;
    // Since only one head has to be bound, the most that will ever
    // be needed (all single-bound) is the total number of sites 
    n_motors_ = n_mts*n_sites;
    motors_.resize(n_motors_);
    for(int ID = 0; ID < n_motors_; ID++){
        motors_[ID].Initialize(parameters_, properties_, ID);
    }
}

void KinesinManagement::SetParameters(){

	// Simulation constants
    double delta_t = parameters_->delta_t;
	double site_size = parameters_->microtubules.site_size;
	double kbT = parameters_->kbT;
	double r_0 = motors_[0].r_0_;
	double k_spring = motors_[0].k_spring_;
	double k_slack = motors_[0].k_slack_;
	double r_y = parameters_->microtubules.y_dist / 2;

	// Non-tethered statistics
    double k_on = parameters_->motors.k_on;
    double c_motor = parameters_->motors.c_bulk;
	p_bind_i_ = k_on * c_motor * delta_t;
	double k_on_ATP = parameters_->motors.k_on_ATP;
	double c_ATP = parameters_->motors.c_ATP;
	p_bind_ATP_ = k_on_ATP * c_ATP * delta_t;
	double k_hydrolyze = parameters_->motors.k_hydrolyze;
	p_hydrolyze_ = k_hydrolyze * delta_t;
	double c_eff = parameters_->motors.c_eff_bind;
	p_bind_ii_ = k_on * c_eff * delta_t; 
	double k_off_ii = parameters_->motors.k_off_ii;
	p_unbind_ii_ = k_off_ii * delta_t; 
	double k_off_i = parameters_->motors.k_off_i;
	p_unbind_i_ = k_off_i * delta_t; 
	// Sound the alarm if our timestep is too large
	if(p_bind_i_ > 1) 
		printf("WARNING: p_bind_i=%g for motors\n", p_bind_i_);
	if(p_bind_ATP_ > 1)
		printf("WARNING: p_bind_ATP=%g for motors\n", p_bind_ATP_);	
	if(p_hydrolyze_ > 1)
		printf("WARNING: p_phos=%g for motors\n", p_hydrolyze_);	
	if(p_bind_ii_ > 1)
		printf("WARNING: p_bind_ii=%g for motors\n", p_bind_ii_);
	if(p_unbind_ii_ > 1)
		printf("WARNING: p_unbind_ii=%g for motors\n", p_unbind_ii_);
	if(p_unbind_i_ > 1)
		printf("WARNING: p_unbind_i=%g for motors\n", p_unbind_i_);

	// Get tether information from motors
	rest_dist_ = motors_[0].rest_dist_;
	comp_cutoff_ = motors_[0].comp_cutoff_;
	dist_cutoff_ = motors_[0].dist_cutoff_;
	if(parameters_->motors.tethers_active){
		printf("\nFor motors:\n");
		printf("  rest_dist is %g\n", rest_dist_);
		printf("  comp_cutoff is %i\n", comp_cutoff_);
		printf("  dist_cutoff is %i\n", dist_cutoff_);
	}
	// Tethered statistics
	double k_tether = parameters_->motors.k_tether;
	double k_untether= parameters_->motors.k_untether;
	double c_eff_teth = parameters_->motors.c_eff_tether;
	if(!parameters_->motors.tethers_active){
		k_tether = 0;
		k_untether = 0;
		c_eff_teth = 0;
	}
	p_bind_i_tethered_ = k_on * c_eff_teth * delta_t;
	p_tether_free_ = k_tether * c_motor * delta_t;
	p_tether_bound_ = k_tether * c_eff_teth * delta_t;
	p_untether_free_ = k_untether * delta_t;
	// Sound the alarm if our timestep is too large
	if(p_bind_i_tethered_ > 1)
		printf("WARNING: p_bind_i_teth=%g for motors\n", 
				p_bind_i_tethered_);
	if(p_tether_free_ > 1) 
		printf("WARNING: p_teth_free=%g for motors\n", p_tether_free_);
	if(p_tether_bound_ > 1)
		printf("WARNING: p_teth_bound=%g for motors\n", p_tether_bound_);
	if(p_untether_free_ > 1)
		printf("WARNING: p_unteth_free=%g for motors\n", p_untether_free_);	
	// For events that depend on tether stretch, each different extension 
	// has its own rate; "base" refers to when the tether is unstretched 
	double p_bind_base = k_on * c_eff_teth * delta_t; 
	double p_unbind_base = p_unbind_i_;
	double p_teth_base = k_tether * c_eff_teth * delta_t;
	double p_unteth_base = p_untether_free_;
	p_unbind_i_tethered_.resize(2*dist_cutoff_ + 1);
	p_untether_bound_.resize(2*dist_cutoff_ + 1);
	for(int x_dub = 0; x_dub <= 2*dist_cutoff_; x_dub++){
		// Calculate tether length for this x_dist 
		double r_x = (double)x_dub * site_size / 2;
		double r = sqrt(r_y*r_y + r_x*r_x);
		// Calculate extension of tether for given x_dub
		double dr = r - r_0; 
		// Calculate potential energy of this tether extension
		double U_teth;
		if(dr > 0)
			U_teth = k_spring * dr * dr / 2;
		else
			U_teth = k_slack * dr * dr / 2;
		// Weight of creating this configuration
		double weight_create = exp(-U_teth/(2*kbT));
		// Weight of annihilating this configuration
		double weight_annihilate = exp(U_teth/(2*kbT));
		// If tethering is disabled, all weights are automatically zero
		if(!parameters_->motors.tethers_active){
			weight_create = 0;
			weight_annihilate = 0; 
		}
		// Otherwise, only weights below comp_cutoff_ are zero
		else if(x_dub < 2*comp_cutoff_){
			weight_create = 0;
			weight_annihilate = 0; 
		}
		// Calculate appropriately-weighted probabilities
		p_unbind_i_tethered_[x_dub] = weight_annihilate * p_unbind_base;
		p_untether_bound_[x_dub] = weight_annihilate * p_unteth_base;
		// Sound the alarm if our timestep is too large
		if(p_unbind_i_tethered_[x_dub] > 1)
			printf("WARNING: p_unbind_i_teth=%g for 2x=%i\n", 
					p_unbind_i_tethered_[x_dub], x_dub);
		if(p_untether_bound_[x_dub] > 1)
			printf("WARNING: p_unteth_bound=%g for 2x=%i\n", 
					p_untether_bound_[x_dub], x_dub);
	}
}

void KinesinManagement::InitializeLists(){

	// One dimensional stuff
	active_.resize(n_motors_);
	free_tethered_.resize(n_motors_);
	bound_untethered_.resize(n_motors_);
	docked_.resize(n_motors_); 
	bound_NULL_.resize(n_motors_);
	bound_ATP_.resize(n_motors_);
	bound_ADPP_i_.resize(n_motors_);
	bound_ADPP_ii_.resize(n_motors_);
	// Two dimensional stuff
	n_bound_tethered_.resize(2*dist_cutoff_ + 1);
	n_bound_ADPP_i_tethered_.resize(2*dist_cutoff_ + 1);
	bound_tethered_.resize(2*dist_cutoff_ + 1);
	bound_ADPP_i_tethered_.resize(2*dist_cutoff_ + 1);
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		n_bound_tethered_[x_dist_dub] = 0;
		n_bound_ADPP_i_tethered_[x_dist_dub] = 0;
		bound_tethered_[x_dist_dub].resize(n_motors_); 
		bound_ADPP_i_tethered_[x_dist_dub].resize(n_motors_);
	}
}

void KinesinManagement::InitializeEvents(){

	int tot_size = 6;
	if(parameters_->motors.tethers_active)
		tot_size = 12 + 4*dist_cutoff_;
	events_.reserve(tot_size);

	// Modular binomial lambda expression for use in event structs
	auto binomial = [&](double p, int n, int ID){
		if(n > 0){
			return properties_->gsl.SampleBinomialDist_Kinesin(p, n, ID);
		}
		else return 0;
	};

	event bind_i = event(0, 10, "bind_i", "unocc", binomial, 
			&properties_->microtubules.n_unoccupied_, p_bind_i_); 
	events_.emplace_back(std::move(bind_i)); 

	event bind_ATP = event(1, 20, "bind_ATP", "bound_NULL", binomial, 
			&n_bound_NULL_, p_bind_ATP_); 
	events_.emplace_back(std::move(bind_ATP)); 

	event hydrolyze = event(2, 30 ,"hydrolyze", "bound_ATP", binomial, 
			&n_bound_ATP_, p_hydrolyze_);
	events_.emplace_back(std::move(hydrolyze)); 

	event bind_ii = event(3, 40, "bind_ii", "docked", binomial, 
			&n_docked_, p_bind_ii_);
	events_.emplace_back(std::move(bind_ii)); 

	event unbind_ii = event(4, 50, "unbind_ii", "bound_ADPP_ii", binomial, 
			&n_bound_ADPP_ii_, p_unbind_ii_);
	events_.emplace_back(std::move(unbind_ii)); 

	event unbind_i = event(5, 60, "unbind_i", "bound_ADPP_i", binomial, 
			&n_bound_ADPP_i_, p_unbind_i_);
	events_.emplace_back(std::move(unbind_i)); 
	// If tethers ARE active, serialize all extension-based events
	if(parameters_->motors.tethers_active){

		auto poisson_bind = [&](double p, int n, int ID){
			if(n > 0){
				double wt = GetWeight_BindTethered();
				double avg = p * wt; 
				return properties_->gsl.SamplePoissonDist_Kinesin(avg, ID);
			}
			else return 0;
		};
		event bind_i_teth = event(6, 11, "bind_i_teth", "free_teth", 
				poisson_bind, &n_free_tethered_, p_tether_free_);
		events_.emplace_back(std::move(bind_i_teth));

		for(int x_dub(0); x_dub <= 2*dist_cutoff_; x_dub++){
			event unbind_i_teth = event(7 + x_dub, 600 + x_dub, 
					"unbind_i_teth", "bound_ADPP_i_teth",
					binomial, &n_bound_ADPP_i_tethered_[x_dub], 
					p_unbind_i_tethered_[x_dub]);
			events_.emplace_back(std::move(unbind_i_teth)); 
		}

		event teth_free = event(8 + 2*dist_cutoff_, 70, 
				"tether_free", "unteth_xlinks", 
				binomial, &properties_->prc1.n_bound_untethered_, 
				p_tether_free_);
		events_.emplace_back(std::move(teth_free));

		auto poisson_teth = [&](double p, int n, int ID){
			if(n > 0){
				double wt = GetWeight_TetherBound();
				double avg = p * wt;
				return properties_->gsl.SamplePoissonDist_Kinesin(avg, ID);
			}
			else return 0;
		};
		event teth_bound = event(9 + 2*dist_cutoff_, 71, 
				"tether_bound", "unteth_xlinks", 
				poisson_teth, &n_bound_untethered_, p_tether_bound_);
		events_.emplace_back(std::move(teth_bound));

		event untether_free = event(10 + 2*dist_cutoff_, 80, 
				"untether_free", "free_tethered", 
				binomial, &n_free_tethered_, p_untether_free_);
		events_.emplace_back(std::move(untether_free));

		int offset =  11 + 2*dist_cutoff_;
		for(int x_dub(0); x_dub <= 2*dist_cutoff_; x_dub++){
			event untether_bound = event(offset + x_dub, 800 + x_dub, 
					"untether_bound", "bound_tethered", 
					binomial, &n_bound_tethered_[x_dub], 
					p_untether_bound_[x_dub]);
			events_.emplace_back(std::move(untether_bound));
		}
	}
}

int KinesinManagement::GetNumBoundUntethered(){

	UpdateBoundUntethered();
	return n_bound_untethered_;
}

double KinesinManagement::GetWeight_BindTethered(){

	double weight = 0;
	for(int i_entry = 0; i_entry < n_free_tethered_; i_entry++){
		weight += free_tethered_[i_entry]->GetTotalBindingWeight();
	}
	return weight; 
}

double KinesinManagement::GetWeight_TetherBound(){

	double weight = 0;
	for(int i_entry = 0; i_entry < n_bound_untethered_; i_entry++){
		weight += bound_untethered_[i_entry]->GetTotalTetheringWeight();
	}
	return weight;
}

Kinesin* KinesinManagement::GetFreeMotor(){

	// Randomly pick a motor from the reservoir
	int i_motor = properties_->gsl.GetRanInt(n_motors_);
	Kinesin *motor = &motors_[i_motor];
	int attempts = 0;
	while(motor->heads_active_ > 0 
	|| motor->tethered_ == true){
		i_motor++;
		if(i_motor == n_motors_) i_motor = 0;
		motor = &motors_[i_motor];
		attempts++;
		if(attempts > n_motors_){
			printf("error in get_free_motor\n");
			exit(1);
		}
	}
	return motor;
}

Kinesin* KinesinManagement::GetBoundUntetheredMotor(){

	UpdateBoundUntethered();
	int i_motor = properties_->gsl.GetRanInt(n_bound_untethered_); 
	return bound_untethered_[i_motor];
}

void KinesinManagement::UpdateAllLists(){

	properties_->microtubules.UpdateUnoccupied();
	UpdateDocked();
	UpdateBoundNULL();
	UpdateBoundATP();
	UpdateBoundADPP_I();
	UpdateBoundADPP_II();
	if(parameters_->motors.tethers_active){
		properties_->prc1.UpdateBoundUntethered();
		UpdateFreeTethered();
		UpdateBoundADPP_I_Tethered();
		UpdateBoundUntethered();
		UpdateBoundTethered();
	}
}

void KinesinManagement::UpdateFreeTethered(){

	n_free_tethered_ = 0;
	for(int i_motor = 0; i_motor < n_active_; i_motor++){
		Kinesin *motor = active_[i_motor];
		if(motor->heads_active_ == 0
		&& motor->tethered_ == true){
			free_tethered_[n_free_tethered_] = motor;
			n_free_tethered_++;
		}
	}
}

void KinesinManagement::UpdateDocked(){

	n_docked_ = 0;
	for(int i_motor = 0; i_motor < n_active_; i_motor++){
		Kinesin *motor = active_[i_motor];
		if(motor->heads_active_ == 1){
			if(motor->GetActiveHead()->ligand_ == "ADPP"){
				double site_coord = motor->GetDockedCoordinate();
				int i_site = site_coord - motor->mt_->coord_;
				if(i_site >= 0 
				&& i_site <= parameters_->microtubules.length - 1){
					// Ensure site isn't occupied
					if(!motor->mt_->lattice_[i_site].occupied_){
						docked_[n_docked_] = motor->GetDockedHead(); 
						n_docked_++;
					}
				}
			}
		}
	}
}

void KinesinManagement::UpdateBoundNULL(){

	n_bound_NULL_ = 0;
	for(int i_motor = 0; i_motor < n_active_; i_motor++){
		Kinesin *motor = active_[i_motor];
		if(motor->heads_active_ == 1){
			if(motor->head_one_.site_ != nullptr
			&& motor->head_one_.ligand_ == "NULL"){
				bound_NULL_[n_bound_NULL_] = &motor->head_one_;
				n_bound_NULL_++;
			}
			if(motor->head_two_.site_ != nullptr
			&& motor->head_two_.ligand_ == "NULL"){
				bound_NULL_[n_bound_NULL_] = &motor->head_two_;
				n_bound_NULL_++;
			}
		}
		else if(motor->heads_active_ == 2){
			if(motor->head_one_.site_ != nullptr
			&& motor->head_one_.ligand_ == "NULL"
			&& motor->head_one_.trailing_){
				bound_NULL_[n_bound_NULL_] = &motor->head_one_;
				n_bound_NULL_++;
			}
			if(motor->head_two_.site_ != nullptr
			&& motor->head_two_.ligand_ == "NULL"
			&& motor->head_two_.trailing_){
				bound_NULL_[n_bound_NULL_] = &motor->head_two_;
				n_bound_NULL_++;
			}
		}
	}
}

void KinesinManagement::UpdateBoundATP(){

	n_bound_ATP_ = 0;
	for(int i_motor = 0; i_motor < n_active_; i_motor++){
		Kinesin *motor = active_[i_motor];
		if(motor->head_one_.site_ != nullptr
		&& motor->head_one_.ligand_ == "ATP"){
			bound_ATP_[n_bound_ATP_] = &motor->head_one_;
			n_bound_ATP_++;
		}
		if(motor->head_two_.site_ != nullptr
		&& motor->head_two_.ligand_ == "ATP"){
			bound_ATP_[n_bound_ATP_] = &motor->head_two_;
			n_bound_ATP_++;
		}
	}
}

void KinesinManagement::UpdateBoundADPP_I(){

	n_bound_ADPP_i_ = 0;
	for(int i_motor = 0; i_motor < n_active_; i_motor++){
		Kinesin *motor = active_[i_motor];
		if(motor->heads_active_ == 1){
			if(!motor->tethered_){
				if(motor->head_one_.site_ != nullptr
				&& motor->head_one_.ligand_ == "ADPP"){
					bound_ADPP_i_[n_bound_ADPP_i_] = &motor->head_one_;
					n_bound_ADPP_i_++;
				}
				if(motor->head_two_.site_ != nullptr
				&& motor->head_two_.ligand_ == "ADPP"){
					bound_ADPP_i_[n_bound_ADPP_i_] = &motor->head_two_;
					n_bound_ADPP_i_++;
				}
			}
			else if(motor->xlink_->heads_active_ == 0){
				if(motor->head_one_.site_ != nullptr
				&& motor->head_one_.ligand_ == "ADPP"){
					bound_ADPP_i_[n_bound_ADPP_i_] = &motor->head_one_;
					n_bound_ADPP_i_++;
				}
				if(motor->head_two_.site_ != nullptr
				&& motor->head_two_.ligand_ == "ADPP"){
					bound_ADPP_i_[n_bound_ADPP_i_] = &motor->head_two_;
					n_bound_ADPP_i_++;
				}
			}
		}
	}
}

void KinesinManagement::UpdateBoundADPP_I_Tethered(){

	for(int x_dub = 0; x_dub <= 2*dist_cutoff_; x_dub++){
		n_bound_ADPP_i_tethered_[x_dub] = 0;
	}
	for(int i_motor = 0; i_motor < n_active_; i_motor++){
		Kinesin *motor = active_[i_motor];
		if(motor->heads_active_ == 1
		&& motor->tethered_){
			if(motor->xlink_->heads_active_ > 0){
				motor->UpdateExtension(); 
				if(motor->tethered_){
					int x_dub = motor->x_dist_doubled_;
					int id = n_bound_ADPP_i_tethered_[x_dub];
					if(motor->head_one_.site_ != nullptr
					&& motor->head_one_.ligand_ == "ADPP"){
						bound_ADPP_i_tethered_[x_dub][id] 
							= &motor->head_one_;
						n_bound_ADPP_i_tethered_[x_dub]++;
						id++;
					}
					if(motor->head_two_.site_ != nullptr
					&& motor->head_two_.ligand_ == "ADPP"){
						bound_ADPP_i_tethered_[x_dub][id] 
							= &motor->head_two_;
						n_bound_ADPP_i_tethered_[x_dub]++;
						id++;
					}
				}
			}
		}
	}
}

void KinesinManagement::UpdateBoundADPP_II(){

	n_bound_ADPP_ii_ = 0;
	for(int i_motor = 0; i_motor < n_active_; i_motor++){
		Kinesin *motor = active_[i_motor];
		if(motor->heads_active_ == 2){
			if(motor->head_one_.ligand_ == "ADPP"){
				bound_ADPP_ii_[n_bound_ADPP_ii_] = &motor->head_one_;
				n_bound_ADPP_ii_++;
			}
			if(motor->head_two_.ligand_ == "ADPP"){
				bound_ADPP_ii_[n_bound_ADPP_ii_] = &motor->head_two_;
				n_bound_ADPP_ii_++;
			}
		}
	}
}

void KinesinManagement::UpdateBoundUntethered(){

	n_bound_untethered_ = 0;
	for(int i_motor = 0; i_motor < n_active_; i_motor++){
		Kinesin *motor = active_[i_motor];
		if(motor->tethered_ == false
		&& motor->heads_active_ > 0){
			bound_untethered_[n_bound_untethered_] = motor;
			n_bound_untethered_++; 
		}
	}
}

void KinesinManagement::UpdateBoundTethered(){

	for(int x_dub = 0; x_dub <= 2*dist_cutoff_; x_dub++){
		n_bound_tethered_[x_dub] = 0;
	}
	for(int i_motor = 0; i_motor < n_active_; i_motor++){
		Kinesin *motor = active_[i_motor];
		if(motor->heads_active_ > 0
		&& motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				motor->UpdateExtension(); 
				if(motor->tethered_ == true){
					int x_dub = motor->x_dist_doubled_; 
					int index = n_bound_tethered_[x_dub];
					bound_tethered_[x_dub][index] = motor;
					n_bound_tethered_[x_dub]++; 
				}
			}
		}
	}
}

void KinesinManagement::GenerateKMCList(){

	sys_time start = sys_clock::now();
	// Update population lists and predicted events
	UpdateAllLists();
	UpdateEvents();
	// Track the time it takes to update lists and sample statistics
	sys_time finish = sys_clock::now();
	auto elapsed = std::chrono::duration_cast<t_unit>(finish - start);
	properties_->wallace.t_motors_[1] += elapsed.count();
	start = sys_clock::now();
	// currently hardcoded ... FIXME eventually 
	// Stat correcton; ensure there aren't more events than population size
	while(events_[3].n_expected_ + events_[5].n_expected_ 
	> n_bound_ADPP_i_){
		double ran = properties_->gsl.GetRanProb();
		double p_tot = events_[3].p_occur_ 
					 + events_[5].p_occur_;
		if(ran < events_[3].p_occur_/p_tot
		&& events_[3].n_expected_ > 0)
			events_[3].n_expected_--;
		else if(events_[5].n_expected_ > 0) 
			events_[5].n_expected_--;
	}

	// XXX ADD TETHERED STAT CORRECTIONS XXX
	int n_events = 0;
	int pre_array[5*events_.size()];
	// Scan over all events; record those with >0 expected in this timestep
	for(int i_entry = 0; i_entry < events_.size(); i_entry++){
		for(int i = 0; i < events_[i_entry].n_expected_; i++){
			pre_array[n_events] = events_[i_entry].kmc_code_;
//			printf("bruh: %i\n", events_[i_entry].kmc_code_);
//			std::cout << events_[i_entry].label_ << std::endl;
			n_events++; 
			// Make sure we don't bamboozle ourselves here
			if(n_events > 5 * events_.size()){
				printf("Error in GenerateKMCList for motors!!\n");
				exit(1);
			}
		}
	}
	// If total expected events is greater than 0, construct kmc_list_
	if(n_events > 0){
		// Trim array to appropriate size/data
		int reduced_array[n_events];
		for(int i_entry = 0; i_entry < n_events; i_entry++){
			reduced_array[i_entry] = pre_array[i_entry];
		}
		// If there are more than 1 expected events, shuffle their order
		if(n_events > 1) gsl_ran_shuffle(properties_->gsl.rng_, 
										 reduced_array, n_events, 
										 sizeof(int));
		// Transfer shuffled array into kmc_list_ structure
        kmc_list_.resize(n_events);
        for(int i_entry = 0; i_entry < n_events; i_entry++){
            kmc_list_[i_entry] = reduced_array[i_entry];
        }
	}
	// Otherwise, simply clear kmc_list_
    else{
        kmc_list_.clear();
    }
	// Conserve memory
	kmc_list_.shrink_to_fit();
	// Track the time it takes to construct KMC list
	finish = sys_clock::now();
	elapsed = std::chrono::duration_cast<t_unit>(finish - start);
	properties_->wallace.t_motors_[2] += elapsed.count();
}

void KinesinManagement::UpdateEvents(){

	for(int i_entry = 0; i_entry < events_.size(); i_entry++){
		/*
		event *entry = &events_[i_entry];
		printf("event \"");
		std::cout << entry->label_;
		printf("\" targets population \"");
		std::cout << entry->target_pop_;
		printf("\", which has a current size of %i\n", *entry->pop_ptr_);
		*/
		events_[i_entry].SampleStatistics(); 
	}
}

void KinesinManagement::RunKMC(){

	sys_time start1 = sys_clock::now();
    GenerateKMCList();
	sys_time start2 = sys_clock::now();
	if(!kmc_list_.empty()){
		int x_dub = 0;
//		printf("\nStart of Kinesin KMC cycle\n");
		for(int i_entry = 0; i_entry < kmc_list_.size(); i_entry++){
//			printf("executing event: ");
//			std::cout << kmc_list_[i_entry]->label_ << std::endl;
			if(kmc_list_[i_entry] >= 600 && kmc_list_[i_entry] < 700){
				x_dub = kmc_list_[i_entry] % 100;
				kmc_list_[i_entry] = 61;
			}
			else if(kmc_list_[i_entry] >= 700 && kmc_list_[i_entry] < 800){
				x_dub = kmc_list_[i_entry] % 100;
				kmc_list_[i_entry] = 81;
			}
			switch(kmc_list_[i_entry]){
				case 10:
		//			printf("Bind_I\n");
					KMC_Bind_I();
					break;
				case 11:
//					printf("Bind_I_Tethered\n");
					KMC_Bind_I_Tethered();
					break;
				case 20:
//					printf("Bind_ATP\n");
					KMC_Bind_ATP();
					break;
				case 30:
//					printf("Hydrolyze\n");
					KMC_Hydrolyze();
					break;
				case 40:
//					printf("Bind_II\n");
					KMC_Bind_II();
					break;
				case 50:
//					printf("Unbind_II\n");
					KMC_Unbind_II();
					break;
				case 60:
//					printf("Unbind_I\n");
					KMC_Unbind_I();
					break;
				case 61:
//					printf("Unbind_I_Teth\n");
					KMC_Unbind_I_Tethered(x_dub);
					break;
				case 70:
//					printf("Tether_Free\n");
					KMC_Tether_Free();
					break;
				case 71:
//					printf("Tether_Bound\n");
					KMC_Tether_Bound();
					break;
				case 80:	
//					printf("Untether_Free\n");
					KMC_Untether_Free();
					break;
				case 81:
//					printf("Untether_Bound\n");
					KMC_Untether_Bound(x_dub);
					break;
			}
		}
	}
	// Track the total time elapsed during kinesin4.RunKMC
	sys_time finish = sys_clock::now();
	auto elapsed = std::chrono::duration_cast<t_unit>(finish - start1);
	properties_->wallace.t_motors_[0] += elapsed.count();
	// Track the time spent executing events
	elapsed = std::chrono::duration_cast<t_unit>(finish - start2);
	properties_->wallace.t_motors_[3] += elapsed.count();
}

void KinesinManagement::KMC_Bind_I(){

	properties_->microtubules.UpdateUnoccupied();
	if(properties_->microtubules.n_unoccupied_ > 0){
		// Get random free motor
		Kinesin *motor = GetFreeMotor();
		// Get random unoccupied site
		Tubulin *site = properties_->microtubules.GetUnoccupiedSite();
		// Update site details
		site->motor_head_ = &motor->head_one_;
		site->occupied_ = true;
		// Update motor details
		motor->mt_ = site->mt_;
		motor->head_one_.site_ = site;
		motor->head_one_.ligand_ = "NULL";
		motor->head_one_.trailing_ = false;
		motor->head_two_.trailing_ = true;
		motor->heads_active_++;
		// Update active_ list
		active_[n_active_] = motor;
		motor->active_index_ = n_active_;
		n_active_++;

	}
	else{
		printf("Failed to Bind_I (motors): no unoccupied sites.\n");
	}
}

void KinesinManagement::KMC_Bind_I_Tethered(){

	UpdateFreeTethered();
	properties_->microtubules.UpdateUnoccupied();
	if(n_free_tethered_ > 0
	&& properties_->microtubules.n_unoccupied_ > 0){
		// Pick a random tethered_free motor to bind
		int i_motor = properties_->gsl.GetRanInt(n_free_tethered_);
		Kinesin *motor = free_tethered_[i_motor];
		motor->UpdateNeighborSites();
		int attempts = 0;
		bool neighbors_exist = true;
		while(motor->n_neighbor_sites_ == 0){
			if(attempts > 10*n_free_tethered_){
				neighbors_exist = false;
				break;
			}
			i_motor = properties_->gsl.GetRanInt(n_free_tethered_);
			motor = free_tethered_[i_motor];
			motor->UpdateNeighborSites();
			attempts++; 
		}
		if(neighbors_exist){
			Tubulin *site = motor->GetWeightedNeighborSite();
			// Update site details
			site->motor_head_ = &motor->head_one_;
			site->occupied_ = true;
			// Update motor details
			motor->mt_ = site->mt_;
			motor->head_one_.site_ = site;
			motor->head_one_.ligand_ = "NULL";
			motor->head_one_.trailing_ = false;
			motor->head_two_.trailing_ = true;
			motor->heads_active_++;
			motor->UpdateExtension();
		}
		else{
			printf("Failed to bind_i_teth\n");
			exit(1);
		}
	}
	else if(properties_->microtubules.n_unoccupied_ == 0){
		printf("Error in Bind_I_Tethered: no unoccupied sites\n");
		//      exit(1);
	}
	else{
		printf("Error in Bind_I_Tethered: no tethered free motors\n");
		//      exit(1);
	}
}

void KinesinManagement::KMC_Bind_ATP(){

	UpdateBoundNULL();
	if(n_bound_NULL_ > 0){
		// Get a random bound_NULL motor
		int i_entry = properties_->gsl.GetRanInt(n_bound_NULL_);
		Kinesin::head *head = bound_NULL_[i_entry];
		// Verify correctness
		if(head->ligand_ != "NULL"
		|| head->site_ == nullptr){
			printf("Error in KMC_Bind_ATP()\n");
			exit(1);
		}
		// Update motor head
		head->ligand_ = "ATP";
		Microtubule *mt = head->site_->mt_;
		int i_site = head->site_->index_;
		int dx = mt->delta_x_;
		int mt_length = parameters_->microtubules.length - 1;
		if(parameters_->motors.endpausing_active){
			if(!(i_site == 0 && dx == -1) 
			&& !(i_site == mt_length && dx == 1)){
				if(!mt->lattice_[i_site + dx].occupied_){
					head->motor_->ChangeConformation();
				}
			}
		}
		else head->motor_->ChangeConformation(); 
	}
	else{
		printf("Failed to Bind_ATP: no bound_NULL motors.\n");
	}
}

void KinesinManagement::KMC_Hydrolyze(){

	UpdateBoundATP();
	if(n_bound_ATP_ > 0){
		// Get a random bound_ATP motor
		int i_entry = properties_->gsl.GetRanInt(n_bound_ATP_);
		Kinesin::head *head = bound_ATP_[i_entry];
		// Verify correctness
		if(head->ligand_ != "ATP"
		|| head->site_ == nullptr){
			printf("Error in KMC_Hydrolyze()\n");
			exit(1);
		}
		head->ligand_ = "ADPP";
	}
	else{
		printf("Failed to Hydrolyze: no bound_ATP motors.\n");
	}
}

void KinesinManagement::KMC_Bind_II(){

	UpdateDocked();
	if(n_docked_ > 0){
		// Get a random docked motor
		int i_entry = properties_->gsl.GetRanInt(n_docked_);
		Kinesin::head *docked_head = docked_[i_entry];
		Kinesin *motor = docked_head->motor_; 
		// Verify correctness
		if(docked_head->ligand_ != "ADP"
		|| motor->GetActiveHead()->ligand_ != "ADPP"
		|| motor->frustrated_){
			printf("Error in KMC_Bind_II()\n");
			exit(1);
		}
		// Verify that proposed site is unoccupied
		int i_dock = motor->GetDockedCoordinate() - motor->mt_->coord_;
		Tubulin *dock_site = &motor->mt_->lattice_[i_dock];
		if(dock_site->occupied_){
			printf("Error in KMC_Bind_II(): Dock is occupied!!\n");
			exit(1);
		}	
		// Update site
		dock_site->motor_head_ = docked_head; 
		dock_site->occupied_ = true;
		// Update motor
		docked_head->site_ = dock_site;
		docked_head->ligand_ = "NULL";
		motor->heads_active_++;
	}
	else{
		printf("Failed to Bind_II: no docked motors.\n");
//		properties_->wallace.PrintMicrotubules(2);
	}
}

void KinesinManagement::KMC_Unbind_II(){

	UpdateBoundADPP_II();
	if(n_bound_ADPP_ii_ > 0){
		// Get a random bound_ADPP motor head
		int i_entry = properties_->gsl.GetRanInt(n_bound_ADPP_ii_);
		Kinesin::head *head = bound_ADPP_ii_[i_entry];
		if(head->site_ == nullptr
		|| head->ligand_ != "ADPP"){
			printf("Error in KMC_Unbind (motors): ");
			std::cout << head->ligand_;
			printf(" bound to head\n");
			exit(1);
		}
		if(head->motor_->heads_active_ == 1){
			printf("Error TWO in KMC_Unbind (motors).\n");
			exit(1);
		}
		// Update site
		head->site_->occupied_ = false;
		head->site_->motor_head_ = nullptr;
		// Update motor
		head->site_ = nullptr;
		head->ligand_ = "ADP";
		head->motor_->heads_active_--;
		if(head->motor_->frustrated_){
			if(head->trailing_){
//				printf("STEPPED TRAILING HEAD:\n");
				head->motor_->ChangeConformation();
				if(head->motor_->frustrated_){
					printf("Error THREE in KMC_Unbind().\n");
					exit(1);
				}
			}
			else{
//				printf("STEPPED LEADING HEAD:\n");
				head->motor_->frustrated_ = false; 
			}
		}
	}
	else{
		printf("Failed to KMC_Unbind_II: no bound_ADPP motors.\n");
	}
}

void KinesinManagement::KMC_Unbind_I(){

	UpdateBoundADPP_I();
	if(n_bound_ADPP_i_ > 0){
		// Get a random bound_ADPP motor head
		int i_entry = properties_->gsl.GetRanInt(n_bound_ADPP_i_);
		Kinesin::head *head = bound_ADPP_i_[i_entry];
		if(head->site_ == nullptr
		|| head->ligand_ != "ADPP"){
			printf("Error in KMC_Unbind (motors): ");
			std::cout << head->ligand_;
			printf(" bound to head\n");
			exit(1);
		}
		if(head->motor_->heads_active_ == 2){
			printf("Error TWO in KMC_Unbind (motors).\n");
			exit(1);
		}
		// Update site
		head->site_->occupied_ = false;
		head->site_->motor_head_ = nullptr;
		// Update motor
		head->site_ = nullptr;
		head->ligand_ = "ADP";
		head->motor_->heads_active_--;
		head->motor_->mt_ = nullptr; 
		// If this motor has a satellite xlink, untether it
		if(head->motor_->xlink_ != nullptr){
			if(head->motor_->xlink_->heads_active_ == 0){
				head->motor_->UntetherSatellite();
			}
			else{
				printf("Error in KMC_Unbind_I: tethed to bound xlink?\n");
				exit(1);
			}
		}
		// Remove this motor from active_, replace with last entry
		Kinesin *last_entry = active_[n_active_ - 1];
		int this_index = head->motor_->active_index_; 
		if(this_index != n_active_ - 1){
			active_[this_index] = last_entry; 
			last_entry->active_index_ = this_index; 
		}
		n_active_--;
	}
	else{
		printf("Failed to KMC_Unbind_I: no bound_ADPP motors.\n");
	}
}

void KinesinManagement::KMC_Unbind_I_Tethered(int x_dub){

	UpdateBoundADPP_I_Tethered();
	int n_bound = n_bound_ADPP_i_tethered_[x_dub]; 
	if(n_bound > 0){
		// Randomly pick a motor
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin::head *head = bound_ADPP_i_tethered_[x_dub][i_entry];
		if(head->motor_->x_dist_doubled_ != x_dub){
			printf("\"not a fan\" -- unbind_i_teth (motors)\n");
			exit(1);
		}
		// Update site
		head->site_->motor_head_ = nullptr;
		head->site_->occupied_ = false;
		// Update motor details
		head->site_ = nullptr;
		head->ligand_ = "ADP";
		head->motor_->heads_active_--;
		head->motor_->mt_ = nullptr;
	}
	else{
		printf("Error in Unbind_I_Tethered: no pseudo bound motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Tether_Free(){

	properties_->prc1.UpdateBoundUntethered();
	if(properties_->prc1.n_bound_untethered_ > 0){
		Kinesin *motor = GetFreeMotor();
		// Randomly pick an xlink
		AssociatedProtein *xlink 
			= properties_->prc1.GetBoundUntetheredXlink();
		// Update motor and xlink details
		motor->xlink_ = xlink;
		motor->tethered_ = true;
		xlink->tethered_ = true;
		xlink->motor_ = motor; 	
		// Update active_ list
		active_[n_active_] = motor;
		motor->active_index_ = n_active_;
		n_active_++;
	}
	else{
		printf("Error in Tether_Free: no bound untethered xlinks\n");
//		exit(1);		
	}
}

void KinesinManagement::KMC_Tether_Bound(){

	UpdateBoundUntethered();
	if(n_bound_untethered_ > 0){
		// Randomly pick a motor
		int i_motor = properties_->gsl.GetRanInt(n_bound_untethered_);
		Kinesin *motor = bound_untethered_[i_motor];
		// Update motor's neighbor xlinks (i.e., those it can teth to)
		motor->UpdateNeighborXlinks(); 
		int attempts = 0; 
		bool neighbors_exist = true; 
		// Ensure we get a motor that has eligible neighbors
		while(motor->n_neighbor_xlinks_ == 0){
			if(attempts > 10*n_bound_untethered_){
				neighbors_exist = false;
				break;
			}
			i_motor = properties_->gsl.GetRanInt(n_bound_untethered_);
			motor = bound_untethered_[i_motor];
			motor->UpdateNeighborXlinks();
			attempts++; 
		}
		if(neighbors_exist){
			AssociatedProtein* xlink = motor->GetWeightedNeighborXlink();
			// Update motor and xlink details
			motor->xlink_ = xlink;
			motor->tethered_ = true;
			xlink->tethered_ = true;
			xlink->motor_ = motor;
			motor->UpdateExtension();
			if(motor->tethered_ == false){
				printf("what in tether_bound nation\n");
				exit(1);
			}
		}
		else{
			printf("Failed to tether bound motor\n");
		}
	}
	else{
		printf("Error in Tether_Bound: no bound untethered motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Untether_Free(){

	UpdateFreeTethered();
	if(n_free_tethered_ > 0){
		// Randomly pick a motor
		int i_entry = properties_->gsl.GetRanInt(n_free_tethered_);
		Kinesin *motor = free_tethered_[i_entry];
		AssociatedProtein *xlink = motor->xlink_;
		// Update motor and xlink detail
		xlink->motor_ = nullptr; 
		xlink->tethered_ = false; 
		motor->xlink_ = nullptr;
		motor->tethered_ = false;
		// Remove this motor from active_, replace with last entry
		// (only if there are more than 1 active motor)
		if(n_active_ > 1){
			int this_index = motor->active_index_; 
			Kinesin *last_entry = active_[n_active_ - 1];
			active_[n_active_ - 1] = nullptr; 
			active_[this_index] = last_entry; 
		}
		n_active_--;
	}
	else{
		printf("Error in Untether_Free: no free tethered motors!\n");
	}
}

void KinesinManagement::KMC_Untether_Bound(int x_dub){

	UpdateBoundTethered();
	int n_bound_tethered = n_bound_tethered_[x_dub];
	if(n_bound_tethered  > 0){
		// Randomly pick a motor
		int i_entry = properties_->gsl.GetRanInt(n_bound_tethered);
		Kinesin *motor = bound_tethered_[x_dub][i_entry];
		if(motor->x_dist_doubled_ != x_dub){
			printf("error in Untether_Bound (motor) x_dub: ");
			printf("%i received; %i in motor\n", x_dub, 
					motor->x_dist_doubled_);
			exit(1);
		}
		AssociatedProtein *xlink = motor->xlink_;
		// Update motor and xlink details
		xlink->motor_ = nullptr;
		xlink->tethered_ = false;
		motor->xlink_ = nullptr;
		motor->tethered_ = false; 
		// Update statistics
		motor->UpdateExtension();
	}
	else{
		printf("Error in Untether_Bound: no bound tethered motors!\n");
//		exit(1);
	}
}
