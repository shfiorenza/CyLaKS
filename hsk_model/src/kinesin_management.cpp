#include "master_header.h"
#include "kinesin_management.h"

KinesinManagement::KinesinManagement(){
}

void KinesinManagement::Initialize(system_parameters *parameters, 
		system_properties *properties){

	parameters_ = parameters;
	properties_ = properties;

	SetParameters();	
	GenerateActiveMotors();	
	GenerateUnboundList();
}

void KinesinManagement::SetParameters(){

	n_tot_ = 0;
	n_bound_ = 0;
	n_unbound_ = 0;

	double k_on = parameters_->k_on;
	double c_motor = parameters_->c_motor;
	double k_off = parameters_->k_off;
	double switch_rate = parameters_->switch_rate;
	double motor_speed = parameters_->motor_speed;
	double delta_t = parameters_->delta_t;
	p_bind_ = k_on*c_motor*delta_t;
	p_unbind_ = k_off*delta_t;
	p_switch_ = switch_rate*delta_t;
	p_step_ = motor_speed*delta_t;

	alpha_ = parameters_->alpha;
	beta_ = parameters_->beta;
}

void KinesinManagement::GenerateActiveMotors(){

	int n_mts = parameters_->n_microtubules;
	int n_sites = parameters_->length_of_microtubule;
	n_tot_ = n_mts*n_sites;
	active_motors.resize(n_tot_);
	for(int ID = 0; ID < n_tot_; ID++){
		active_motors[ID].Initialize(parameters_, properties_, ID);
	}
}

void KinesinManagement::GenerateUnboundList(){

	for(int i_entry; i_entry < n_tot_; i_entry++){
		if(active_motors[i_entry].site_ == nullptr){
			Kinesin *motor = &active_motors[i_entry];
			unbound_list.push_back(motor);
			n_unbound_++;
		}
		else{
			printf("An error occured while generating unbound_list\n");
			exit(1);
		}
	}
}

void KinesinManagement::UnboundCheck(int ID){

	if(active_motors[ID].site_ != nullptr){
		printf("Error with motor #%i; bound but in unbound_list\n", ID);
		exit(1);
	}
}

void KinesinManagement::BoundCheck(int ID){

	if(active_motors[ID].site_ == nullptr){
		printf("Error with motor #%i; unbound but in bound_list\n", ID);
		exit(1);
	}
}

bool KinesinManagement::BoundaryCheck(Kinesin *motor){

	bool boundary_status; 
	if(motor->site_->index_ == motor->mt_->plus_end_ ||
	   motor->site_->index_ == motor->mt_->minus_end_){
		boundary_status = true;
	}
	else{
		boundary_status = false;
	}
	return boundary_status;
}

void KinesinManagement::GenerateKMCList(){

	int n_to_bind = GetNumToBind();
	int n_to_unbind = GetNumToUnbind();
	int n_to_switch = GetNumToSwitch();	
	int n_to_step = GetNumToStep();
	int n_events = n_to_bind + n_to_unbind + n_to_switch + n_to_step;
	if(n_events > 0){
		int pre_list[n_events];
		int i_event = 0;
		for(i_event; i_event < n_to_bind; i_event++){
			pre_list[i_event] = 0;
		}
		for(i_event; i_event < n_to_bind + n_to_unbind; i_event++){
			pre_list[i_event] = 1;
		}
		for(i_event; i_event < n_events - n_to_step; i_event++){
			pre_list[i_event] = 2;
		}
		for(i_event; i_event < n_events; i_event++){
			pre_list[i_event] = 3;
		}
		gsl_ran_shuffle(properties_->gsl.rng, pre_list, n_events, sizeof(int));
		kmc_list.resize(n_events);
		for(int i_event = 0; i_event < n_events; i_event++){
			kmc_list[i_event] = pre_list[i_event];
		}
	}
	else{
		kmc_list.clear();
	}
}

int KinesinManagement::GetNumToBind(){

	if(properties_->microtubules.n_unoccupied_ !=
	   properties_->microtubules.unoccupied_list.size()){
		printf("bro what: THE PREQUEL");
		exit(1);
	}
	int n_to_bind;
	double n_avg = p_bind_*properties_->microtubules.n_unoccupied_;
	n_to_bind = properties_->gsl.SamplePoissonDist(n_avg);
	return n_to_bind;
}

int KinesinManagement::GetNumToUnbind(){

	if(n_bound_ != bound_list.size()){
		printf("bro what\n");
		exit(1);
	}
	int n_to_unbind;
	n_to_unbind = properties_->gsl.SampleBinomialDist(p_unbind_, n_bound_);
	return n_to_unbind;
}

int KinesinManagement::GetNumToSwitch(){

	int n_to_switch;
	int n_switchable = 0;
	for(int i_motor = 0; i_motor < n_bound_; i_motor++){
		int i_mt = bound_list[i_motor]->mt_->index_;
		Microtubule *mt = &properties_->microtubules.active_mts[i_mt];
		int plus_end = mt->plus_end_;
		int minus_end = mt->minus_end_;	
		int i_mt_adj = mt->mt_index_adj_;
		Microtubule *mt_adj = &properties_->microtubules.active_mts[i_mt_adj];
		int i_site = bound_list[i_motor]->site_->index_;
		// Exclude boundary sites from statistics (no switching there)
		if(i_site != plus_end && i_site != minus_end){
			if(mt_adj->lattice[i_site].occupant == nullptr){
				n_switchable++;
			}
		}
	}
	n_to_switch = properties_->gsl.SampleBinomialDist(p_switch_, n_switchable);
	return n_to_switch;
}

int KinesinManagement::GetNumToStep(){

	int n_to_step; 
	int n_steppable = 0;
	for(int i_motor = 0; i_motor < n_bound_; i_motor++){
		int i_mt = bound_list[i_motor]->mt_->index_;
		Microtubule *mt = &properties_->microtubules.active_mts[i_mt];
		int plus_end = mt->plus_end_;
		int delta_x = mt->delta_x_;
		int i_site = bound_list[i_motor]->site_->index_;
		// Exclude plus_end from statistics because of end-pausing
		if(i_site != plus_end){
			if(mt->lattice[i_site + delta_x].occupant == nullptr){
				n_steppable++; 
			}
		}
	}
	n_to_step = properties_->gsl.SampleBinomialDist(p_step_, n_steppable);
	return n_to_step;
}

void KinesinManagement::RunKMC(){

	GenerateKMCList();
	if(kmc_list.empty() != true){
		int n_events = kmc_list.size();
		for(int i_event; i_event < n_events; i_event++){
			int kmc_event = kmc_list[i_event];
			switch(kmc_event){
				case 0:	RunKMC_Bind();
						break;
				case 1:	RunKMC_Unbind();
						break;
				case 2: RunKMC_Switch();
						break;
				case 3: RunKMC_Step();
						break;
			}
			RunKMC_Boundaries(n_events);
		}
	}
	else{
		RunKMC_Boundaries(1);
	}
}

void KinesinManagement::RunKMC_Bind(){

	if(unbound_list.empty() != true){
		// Get unbound motor
		int i_entry = properties_->gsl.GetRanInt(unbound_list.size());
		Kinesin *motor = unbound_list[i_entry];
		UnboundCheck(motor->ID_);
		// Find unoccupied site to bind to
		Tubulin *site = properties_->microtubules.GetUnoccupiedSite();
		Microtubule *mt = site->parent;
		// Place motor on site; remove site from unoccupied_list
		site->occupant = motor;
		properties_->microtubules.RemoveFromUnoccupiedList(site);
		// Update motor details
		motor->mt_ = mt;
		motor->site_ = site;
		// Remove motor from unbound_list; place it in bound_list
		unbound_list.erase(unbound_list.begin() + i_entry);
		n_unbound_--;
		bound_list.emplace_back(motor);
		n_bound_++;
	}
	else{
		printf("Error in RunKMC_Bind(): unbound_list is empty\n");
		exit(1);
	}
}

void KinesinManagement::RunKMC_Unbind(){

	if(bound_list.empty() != true){
		int i_entry, n_attempts = 0;
		Kinesin* motor;
		bool failure = false;
		// Attempt to find a bound motor; exclude boundary sites
		do{	if(n_attempts > 2*bound_list.size()){
				failure = true;
				break;
			}
			i_entry = properties_->gsl.GetRanInt(bound_list.size());
			motor = bound_list[i_entry];
			BoundCheck(motor->ID_);
			n_attempts++;
		}while(BoundaryCheck(motor) == true);
		// Only attempt to unbind if an eligible motor was found
		if(failure != true){
			// Remove motor from site; add site to unoccupied_list
			motor->site_->occupant = nullptr;
			properties_->microtubules.AddToUnoccupiedList(motor->site_);
			// Update motor details
			motor->mt_ = nullptr;
			motor->site_ = nullptr;
			// Remove motor from bound_list; place it in unbound_list
			bound_list.erase(bound_list.begin() + i_entry);
			n_bound_--;
			unbound_list.emplace_back(motor);
			n_unbound_++;
		}
		else{
			printf("Failed to unbind.\n");
		}
	}
	else{
		printf("Error in RunKMC_Unbind(): bound_list is empty\n");
		exit(1);
	}
}

void KinesinManagement::RunKMC_Switch(){

	if(bound_list.empty() != true){
		int i_motor, i_site, i_mt_adj, n_attempts = 0;
		Kinesin *motor;
		Tubulin *new_site;
		Microtubule *adjacent_mt;
		bool failure = false;
		// Attempt to find a motor capable of stepping
		do{	if(n_attempts > 2*bound_list.size()){
				failure = true;
				break;
			}
			i_motor = properties_->gsl.GetRanInt(bound_list.size());
			motor = bound_list[i_motor];
			BoundCheck(motor->ID_);
			i_site = motor->site_->index_;
			i_mt_adj = motor->mt_->mt_index_adj_;
			adjacent_mt = &properties_->microtubules.active_mts[i_mt_adj];
			new_site = &adjacent_mt->lattice[i_site];
			n_attempts++;
		}while(BoundaryCheck(motor) == true || new_site->occupant != nullptr);
		// Only attempt to switch if an eligible motor was found
		if(failure != true){
			Tubulin *old_site = &motor->mt_->lattice[i_site];
			properties_->microtubules.OccupiedCheck(old_site);
			properties_->microtubules.UnoccupiedCheck(new_site);
			// Remove motor from old site; Place it on new site
			old_site->occupant = nullptr;
			new_site->occupant = motor;
			// Update motor details
			motor->site_ = new_site;
			motor->mt_ = adjacent_mt;
			// Update MT lists 
			properties_->microtubules.AddToUnoccupiedList(old_site);	
			properties_->microtubules.RemoveFromUnoccupiedList(new_site);
		}	
		else{
			int ID = motor->ID_;
			int i_mt = motor->mt_->index_;
			printf("Motor %i failed to switch @ %i_%i\n", ID, i_mt, i_site);
		}
	}
	else{
		printf("Error in RunKMC_Switch(): bound_list is empty\n");
		exit(1);
	}
}

void KinesinManagement::RunKMC_Step(){

	if(bound_list.empty() != true){
		int i_motor, i_site, plus_end, minus_end, delta_x, n_attempts = 0;
		Kinesin *motor;
		bool failure = false;
		// Attempt to find a motor capable of stepping
		do{	if(n_attempts > 2*bound_list.size()){
				failure = true;
				break;
			}
			i_motor = properties_->gsl.GetRanInt(bound_list.size());
			motor = bound_list[i_motor];
			BoundCheck(motor->ID_);
			i_site = motor->site_->index_;
			plus_end = motor->mt_->plus_end_;
			minus_end = motor->mt_->minus_end_;
			delta_x = motor->mt_->delta_x_;
			n_attempts++;
			// Ensure chosen motor is not on plus_end (enforce end-pausing)
			while(motor->site_->index_ == plus_end){
				if(n_attempts > 2*bound_list.size()){
					failure = true;
					break;
				}
				i_motor = properties_->gsl.GetRanInt(bound_list.size());
				motor = bound_list[i_motor];
				BoundCheck(motor->ID_);
				plus_end = motor->mt_->plus_end_;
				minus_end = motor->mt_->minus_end_;
				delta_x = motor->mt_->delta_x_;
				i_site = motor->site_->index_;
				n_attempts++;
			}
		}while(motor->mt_->lattice[i_site + delta_x].occupant != nullptr);
		// Only attempt to step if an eligible motor was found
		if(failure != true){
			Tubulin *old_site = &motor->mt_->lattice[i_site];
			properties_->microtubules.OccupiedCheck(old_site);
			Tubulin *new_site = &motor->mt_->lattice[i_site + delta_x];
			properties_->microtubules.UnoccupiedCheck(new_site);
			// Remove motor from old site; Place it on new site
			old_site->occupant = nullptr;
			new_site->occupant = motor;
			// Update motor details
			motor->site_ = new_site;
			// If stepping from minus_end, supress old site from unoccupied_list
			if(i_site != minus_end){
				properties_->microtubules.AddToUnoccupiedList(old_site);	
			}
			// If stepping to plus_end, site is absent from unoccupied_list
			if(i_site + delta_x != plus_end){
				properties_->microtubules.RemoveFromUnoccupiedList(new_site);
			}
		}
		else{
			int ID = motor->ID_;
			int i_mt = motor->mt_->index_;
			printf("Motor %i failed to step @ %i_%i\n", ID, i_mt, i_site);
		}
	}
	else{
		printf("Error in RunKMC_Step(): bound_list is empty\n");
		exit(1);
	}
}

void KinesinManagement::RunKMC_Boundaries(int n_events){

	int n_mts = parameters_->n_microtubules;
	double alpha_eff = (alpha_*p_step_)/n_events;
	double beta_eff = (beta_*p_step_)/n_events;
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		Microtubule *mt = &properties_->microtubules.active_mts[i_mt];
		Tubulin *minus_end = &mt->lattice[mt->minus_end_];
		Tubulin *plus_end = &mt->lattice[mt->plus_end_];
		double ran1 = properties_->gsl.GetRanProb();
		double ran2 = properties_->gsl.GetRanProb();
		// Insert motors into lattice with probability alpha_eff	
		if(minus_end->occupant == nullptr && ran1 < alpha_eff){
			// Get unbound motor
			int i_entry = properties_->gsl.GetRanInt(unbound_list.size());
			Kinesin *motor = unbound_list[i_entry];
			UnboundCheck(motor->ID_);
			// Place motor on site
			minus_end->occupant = motor;
			// Update motor details
			motor->mt_ = minus_end->parent;
			motor->site_ = minus_end;
			// Remove motor from unbound_list; place it into bound_list
			unbound_list.erase(unbound_list.begin() + i_entry);
			n_unbound_--;
			bound_list.emplace_back(motor);
			n_bound_++;
		}
		// Remove motors from minus_end with probability beta_eff
		if(plus_end->occupant != nullptr && ran2 < beta_eff){
			// Get motor bound to minus_end
			Kinesin *motor = plus_end->occupant;
			BoundCheck(motor->ID_);
			// Remove motor from site
			plus_end->occupant = nullptr;
			// Update motor details
			motor->mt_ = nullptr;
			motor->site_ = nullptr;
			// Generate iterator pointing to this site's entry in bound_list
			auto motor_entry = std::find(bound_list.begin(), bound_list.end(), motor);
			// Get numerical index (i.e. distance from vector's start) of the iterator
			int i_entry = std::distance(bound_list.begin(), motor_entry);
			// Erase entry in bound_list that corresponds to this motor
			bound_list.erase(bound_list.begin() + i_entry);
			n_bound_--;
			// Place motor in unbound_list 
			unbound_list.emplace_back(motor);
			n_unbound_++;
		} 
	}
}
