#include "master_header.h"
#include "kinesin_management.h"

KinesinManagement::KinesinManagement(){
}

void KinesinManagement::Initialize(system_parameters *parameters, 
		system_properties *properties){

	parameters_ = parameters;
	properties_ = properties;

	SetParameters();	
	GenerateMotors();	
}

void KinesinManagement::SetParameters(){

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

void KinesinManagement::GenerateMotors(){

	int n_mts = parameters_->n_microtubules;
	int n_sites = parameters_->length_of_microtubule;
	// Since each kinesin occupies two sites, the most that will ever
	// be needed (all sites occupied) is half the total number of sites 
	n_motors_tot_ = n_mts*n_sites; // /2; FIXME
	motor_list_.resize(n_motors_tot_);
	for(int ID = 0; ID < n_motors_tot_; ID++){
		motor_list_[ID].Initialize(parameters_, properties_, ID);
	}
	// Also resize bound_list_ here since there's no other good place to do it
	bound_list_.resize(n_motors_tot_);
}

void KinesinManagement::UnboundCheck(Kinesin *motor){

	// Check motor
	if(motor->front_site_ != nullptr){// || motor->rear_site_ != nullptr){
		printf("Error with motor #%i: classified as unbound, but", motor->ID_);
		printf(" at least one of its heads is attached to tubulin\n");
		exit(1);
	}
	// Check motor_list
	int ID = motor->ID_;
	if(motor_list_[ID].front_site_ != nullptr){
//	|| motor_list_[ID].rear_site_ != nullptr){
		printf("Error: motor #%i is out of sync with motor_list\n", ID);
		exit(1);
	}
}

void KinesinManagement::BoundCheck(Kinesin *motor){

	// Check motor
	if(motor->front_site_ == nullptr){// || motor->rear_site_ == nullptr){
		printf("Error with motor #%i: classified as bound, but" , motor->ID_);
		printf(" at least one of its heads is not attached to tubulin\n");
		exit(1);
	}
	// Check motor_list
	int ID = motor->ID_;
	if(motor_list_[ID].front_site_ == nullptr){
 //	|| motor_list_[ID].rear_site_ == nullptr){
		printf("Error: motor #%i is out of sync with motor_list\n", ID);
		exit(1);
	}
}

bool KinesinManagement::OnBoundarySite(Kinesin *motor){

	bool boundary_status; 
	if(motor->front_site_->index_ == motor->mt_->plus_end_
	|| motor->front_site_->index_ == motor->mt_->minus_end_){
//	|| motor->rear_site_->index_ == motor->mt_->minus_end_){
		boundary_status = true;
	}
	else{
		boundary_status = false;
	}
	return boundary_status;
}

void KinesinManagement::UpdateBoundList(){
	
	int i_bound = 0;
	for(int i_motor = 0; i_motor < n_motors_tot_; i_motor++){
		Kinesin *motor = &motor_list_[i_motor];
		if(motor->bound_ == true){
			bound_list_[i_bound] = motor;
			i_bound++;		
		}
	}
	if(i_bound != n_bound_){
		printf("something awful in update_bound_list bruh\n");
		exit(1);
	}
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
		kmc_list_.resize(n_events);
		for(int i_event = 0; i_event < n_events; i_event++){
			kmc_list_[i_event] = pre_list[i_event];
		}
	}
	else{
		kmc_list_.clear();
	}
}

int KinesinManagement::GetNumToBind(){

//	properties_->microtubules.UpdateNumUnoccupiedPairs();
//	int n_unoccupied_pairs = properties_->microtubules.n_unoccupied_pairs_;
//	double n_avg = p_bind_*n_unoccupied_pairs;
	properties_->microtubules.UpdateNumUnoccupied();
	int n_unoccupied = properties_->microtubules.n_unoccupied_;
	double n_avg = p_bind_*n_unoccupied;
	int n_to_bind = properties_->gsl.SamplePoissonDist(n_avg);
	if(n_unoccupied > 0){
//	if(n_unoccupied_pairs > 0){
//		double p_bind = n_to_bind/(double)n_unoccupied_pairs;
		double p_bind = n_to_bind/(double)n_unoccupied;
		properties_->p_bind_cum_ += p_bind;
		properties_->n_binds_++;
	}
	return n_to_bind;
}

int KinesinManagement::GetNumToUnbind(){

	int n_to_unbind = properties_->gsl.SampleBinomialDist(p_unbind_, n_bound_);
	if(n_bound_ > 0){
		double p_unbind = n_to_unbind/(double)n_bound_;
		properties_->p_unbind_cum_ += p_unbind;
		properties_->n_unbinds_++;
	}
	return n_to_unbind;
}

int KinesinManagement::GetNumToSwitch(){

	int n_to_switch;
	int n_switchable = 0;
	int n_unbound = 0;
	for(int i_motor = 0; i_motor < n_motors_tot_; i_motor++){
		if(motor_list_[i_motor].bound_ == true){
			Kinesin *motor = &motor_list_[i_motor];
			Microtubule *mt = motor->mt_;
			int plus_end = mt->plus_end_;
			int minus_end = mt->minus_end_;	
			int i_adj = mt->mt_index_adj_;
			Microtubule *mt_adj = &properties_->microtubules.mt_list_[i_adj];
			int i_front_site = motor->front_site_->index_;
//			int i_rear_site = motor->rear_site_->index_;
			// Exclude boundary sites from statistics (no switching there)
			if(i_front_site != plus_end){// && i_rear_site != minus_end){
				if(mt_adj->lattice_[i_front_site].occupied_ == false){
//				&& mt_adj->lattice_[i_rear_site].occupied_ == false){
					n_switchable++;
				}
			}
		}
		else{
			n_unbound++;
		}
	}
	if(n_unbound != n_motors_tot_ - n_bound_){
		printf("Something went wrong in Kinesin MGMT!\n");
		exit(1);
	}
	n_to_switch = properties_->gsl.SampleBinomialDist(p_switch_, n_switchable);
	return n_to_switch;
}

int KinesinManagement::GetNumToStep(){

	int n_to_step; 
	int n_steppable = 0;
	int n_unbound = 0;
	for(int i_motor = 0; i_motor < n_motors_tot_; i_motor++){
		if(motor_list_[i_motor].bound_ == true){
			Kinesin *motor = &motor_list_[i_motor];
			Microtubule *mt = motor->mt_;
			int plus_end = mt->plus_end_;
			int delta_x = mt->delta_x_;
			int i_front_site = motor->front_site_->index_;
			// Exclude plus_end from statistics because of end-pausing
			if(i_front_site != plus_end){
				if(mt->lattice_[i_front_site + delta_x].occupied_ == false){
					n_steppable++; 
				}
			}
		}
		else{
			n_unbound++;
		}
	}
	if(n_unbound != n_motors_tot_ - n_bound_){
		printf("Something went wrong in Kinesin MGMT (2):\n");
		printf("  %i != %i - %i\n", n_unbound, n_motors_tot_, n_bound_); 
		exit(1);
	}
	n_to_step = properties_->gsl.SampleBinomialDist(p_step_, n_steppable);
	return n_to_step;
}

void KinesinManagement::RunKMC(){

	GenerateKMCList();
	if(kmc_list_.empty() == false){
		int n_events = kmc_list_.size();
		for(int i_event; i_event < n_events; i_event++){
			int kmc_event = kmc_list_[i_event];
			switch(kmc_event){
				case 0:	RunKMC_Bind();
						break;
				case 1:	RunKMC_Unbind();
						break;
//				case 2: RunKMC_Switch();
//						break;
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

	// Make sure that at least one unbound motor exists
	if(n_bound_ != n_motors_tot_){
        // Randomly pick a motor from the reservoir
        int i_motor = properties_->gsl.GetRanInt(n_motors_tot_);
        Kinesin *motor = &motor_list_[i_motor];
        while(motor->bound_ == true){
            i_motor++;
            if(i_motor == n_motors_tot_)
				i_motor = 0;
			motor = &motor_list_[i_motor];
		}
		UnboundCheck(motor);
		properties_->microtubules.UpdateUnoccupiedList();
		if(properties_->microtubules.n_unoccupied_ > 0){
//		properties_->microtubules.UpdateUnoccupiedPairsList();
//		if(properties_->microtubules.n_unoccupied_pairs_ > 0){
			// Find unoccupied pair of sites to bind to
			Tubulin *front_site = properties_->microtubules.GetUnoccupiedSite();
//			Tubulin *rear_site = properties_->microtubules.GetUnoccupiedPair_1();
//			Tubulin *front_site = properties_->microtubules.GetUnoccupiedPair_2(rear_site);
			Microtubule *mt = front_site->mt_;
			// Place motor on sites
//			rear_site->motor_ = motor;
//			rear_site->occupied_ = true;
			front_site->motor_ = motor;
			front_site->occupied_ = true;
			// Update motor details
			motor->bound_ = true;
			motor->mt_ = mt;
//			motor->rear_site_ = rear_site;
			motor->front_site_ =  front_site;
			// Update statistics
			n_bound_++;
		}
		else{
			printf("Failed to bind\n");
		}
	}
	else{
		printf("Error in RunKMC_Bind: no unbound motors\n");
		exit(1);
	}
}

void KinesinManagement::RunKMC_Unbind(){

	// Make sure that at least one bound motor exists
	if(n_bound_ > 0){
		UpdateBoundList();
		// Randomly pick a bound motor
		int i_entry = properties_->gsl.GetRanInt(n_bound_);
		Kinesin* motor = bound_list_[i_entry];
		BoundCheck(motor);
		int n_attempts = 0;
		bool failure = false;
		// Ensure that motors on boundary sites are excluded
		while(OnBoundarySite(motor) == true){
			if(n_attempts > 2*n_bound_){
				failure = true;
				break;
			}
			i_entry = properties_->gsl.GetRanInt(n_bound_);
			motor = bound_list_[i_entry];
			BoundCheck(motor);
			n_attempts++;
		}	
		// Make sure an eligible motor was found
		if(failure == false){
			// Remove motor from sites
			motor->front_site_->motor_ = nullptr;
			motor->front_site_->occupied_ = false;
//			motor->rear_site_->motor_ = nullptr;
//			motor->rear_site_->occupied_ = false;
			// Update motor details
			motor->bound_ = false;
			motor->mt_ = nullptr;
			motor->front_site_ = nullptr;
//			motor->rear_site_ = nullptr;
			// Update statistics
			n_bound_--;
		}
		else{
			printf("Failed to unbind\n");
		}
	}
	else{
		printf("Error in RunKMC_Unbind: no bound motors\n");
		exit(1);
	}
}
/*
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
*/
void KinesinManagement::RunKMC_Step(){

	// Make sure there is at least one bound motor
	if(n_bound_ > 0){
		UpdateBoundList();
		int i_entry, i_plus_end, delta_x;
		Kinesin *motor;
		Tubulin *front_site, *rear_site;
		int n_attempts = 0;
		bool failure = false;
		// Ensure that blocked motors and motors on the plus end are excluded
		do{	
			if(n_attempts > 2*n_bound_){
				failure = true;
				break;
			}
			i_entry = properties_->gsl.GetRanInt(n_bound_);
			motor = bound_list_[i_entry];
			BoundCheck(motor);
			front_site = motor->front_site_;
//			rear_site = motor->rear_site_;
			i_plus_end = motor->mt_->plus_end_;
			delta_x = motor->mt_->delta_x_;
			n_attempts++;
			// Ensure chosen motor is not on plus_end (enforce end-pausing)
			while(front_site->index_ == i_plus_end){
				if(n_attempts > 2*n_bound_){
					failure = true;
					break;
				}
				i_entry = properties_->gsl.GetRanInt(n_bound_);
				motor = bound_list_[i_entry];
				BoundCheck(motor);
				front_site = motor->front_site_;
//				rear_site = motor->rear_site_;
				i_plus_end = motor->mt_->plus_end_;
				delta_x = motor->mt_->delta_x_;
				n_attempts++;
			}
		}while(motor->mt_->lattice_[front_site->index_ + delta_x].occupied_ == true);
		// Make sure an eligible motor was found
		if(failure == false){
			Tubulin *old_front_site = front_site;
//			Tubulin *old_rear_site = rear_site;
			properties_->microtubules.OccupiedCheck(old_front_site);
//			properties_->microtubules.OccupiedCheck(old_rear_site);
			Tubulin *new_front_site = &motor->mt_->lattice_[old_front_site->index_ + delta_x];
//			Tubulin *new_rear_site = old_front_site;
			properties_->microtubules.UnoccupiedCheck(new_front_site);
			// Take rear motor head off of old rear site; place it on old front site
			old_front_site->motor_ = nullptr;
			old_front_site->occupied_ = false;
//			old_rear_site->motor_ = nullptr;
//			old_rear_site->occupied_ = false;
//			motor->rear_site_ = old_front_site;
			// Take front motor head off of old front site; place it on new front site
			motor->front_site_ = new_front_site;
			new_front_site->motor_ = motor;
			new_front_site->occupied_ = true;
		}
		else{
			int ID = motor->ID_;
			int i_mt = motor->mt_->index_;
			int i_front = motor->front_site_->index_;
//			int i_rear = motor->rear_site_->index_;
			printf("Motor %i failed to step @ %i_%i\n", ID, i_mt, i_front);//, i_rear); //FIXME
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
//	printf("%g_%g (%i)\n", alpha_eff, beta_eff, n_events);
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
		int i_plus_end = mt->plus_end_;
		int i_minus_end = mt->minus_end_;
		int delta_x = mt->delta_x_;
		Tubulin *minus_end = &mt->lattice_[i_minus_end];
//		Tubulin *minus_neighbor = &mt->lattice_[i_minus_end + delta_x];
		Tubulin *plus_end = &mt->lattice_[i_plus_end];
//		Tubulin *plus_neighbor = &mt->lattice_[i_plus_end - delta_x];
		double ran1 = properties_->gsl.GetRanProb();
		double ran2 = properties_->gsl.GetRanProb();
		// Insert motors onto minus_end with probability alpha_eff	
		if(ran1 < alpha_eff
		&& minus_end->occupied_ == false){ 
//		&& minus_neighbor->occupied_ == false){
			// Start at random point in motor_list, then scan for unbound motor
        	int i_motor = properties_->gsl.GetRanInt(n_motors_tot_);
        	Kinesin *motor = &motor_list_[i_motor];
        	while(motor->bound_ == true){
           		i_motor++;
				if(i_motor == n_motors_tot_)
					i_motor = 0;
			motor = &motor_list_[i_motor];
			}	
			UnboundCheck(motor);
			// Place motor on sites
			minus_end->motor_ = motor;
			minus_end->occupied_ = true;
//			minus_neighbor->motor_ = motor;
//			minus_neighbor->occupied_ = true;
			// Update motor details
			motor->bound_ = true;
			motor->mt_ = minus_end->mt_;
//			motor->rear_site_ = minus_end;
			motor->front_site_ = minus_end;
//			motor->front_site_ =  minus_neighbor;
			// Update statistics
			n_bound_++;
		}
		// Remove motors from plus_end with probability beta_eff
		if(ran2 < beta_eff
		&& plus_end->occupied_ == true){
//		&& plus_neighbor->occupied_ == true){
			// Get motor bound to minus_end
			Kinesin *motor = plus_end->motor_;
			BoundCheck(motor);
			// Remove motor from sites
			plus_end->motor_ = nullptr;
			plus_end->occupied_ = false;
//			plus_neighbor->motor_ = nullptr;
//			plus_neighbor->occupied_ = false;
			// Update motor details
			motor->bound_ = false;
			motor->mt_ = nullptr;
			motor->front_site_ = nullptr;
//			motor->rear_site_ = nullptr;
			// Update statistics
			n_bound_--;
		} 
	}
}

