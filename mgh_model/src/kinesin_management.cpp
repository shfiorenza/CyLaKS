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
	p_bind_tethered_ = p_bind_*2;	//FIXME
	p_unbind_ = k_off*delta_t;
	p_tether_unbound_ = p_bind_*2;  //FIXME 
	p_switch_ = switch_rate*delta_t;
	p_step_ = motor_speed*delta_t;

	alpha_ = parameters_->alpha;
	beta_ = parameters_->beta;
}

void KinesinManagement::GenerateMotors(){

	int n_mts = parameters_->n_microtubules;
	int n_sites = parameters_->length_of_microtubule;
	// Since only one head has to be bound, the most that will ever
	// be needed (all single-bound) is the total number of sites 
	n_motors_ = n_mts*n_sites;
	motor_list_.resize(n_motors_);
	for(int ID = 0; ID < n_motors_; ID++){
		motor_list_[ID].Initialize(parameters_, properties_, ID);
	}
	// Also resize bound_list_ here since there's no other good place to do it
	bound_list_.resize(n_motors_);
	tethered_unbound_list_.resize(n_motors_);
}

void KinesinManagement::UnboundCheck(Kinesin *motor){

	// Check motor
	if(motor->front_site_ != nullptr || motor->rear_site_ != nullptr){
		printf("Error with motor #%i: classified as unbound, but", motor->ID_);
		printf(" at least one of its heads is attached to tubulin\n");
		exit(1);
	}
	// Check motor_list
	int ID = motor->ID_;
	if(motor_list_[ID].front_site_ != nullptr
	|| motor_list_[ID].rear_site_ != nullptr){
		printf("Error: motor #%i is out of sync with motor_list\n", ID);
		exit(1);
	}
}

void KinesinManagement::BoundCheck(Kinesin *motor){

	// Check motor
	if(motor->front_site_ == nullptr || motor->rear_site_ == nullptr){
		printf("Error with motor #%i: classified as bound, but" , motor->ID_);
		printf(" at least one of its heads is not attached to tubulin\n");
		exit(1);
	}
	// Check motor_list
	int ID = motor->ID_;
	if(motor_list_[ID].front_site_ == nullptr
 	|| motor_list_[ID].rear_site_ == nullptr){
		printf("Error: motor #%i is out of sync with motor_list\n", ID);
		exit(1);
	}
}

bool KinesinManagement::BoundaryStatus(Kinesin *motor){

	bool boundary_status; 
	if(motor->front_site_->index_ == motor->mt_->plus_end_
	|| motor->rear_site_->index_ == motor->mt_->minus_end_){
		boundary_status = true;
	}
	else{
		boundary_status = false;
	}
	return boundary_status;
}

void KinesinManagement::UpdateLists(){
	
	int i_bound = 0;
	int i_tethered = 0;
	int n_unbound = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motor_list_[i_motor];
		if(motor->bound_ == true){
			bound_list_[i_bound] = motor;
			i_bound++;		
		}
		else if(motor->tethered_ == true){
			tethered_unbound_list_[i_tethered] = motor;
			i_tethered++;
			n_unbound++;
		}
		else
			n_unbound++;
	}
	if(i_bound != n_single_bound_ + n_double_bound_
	|| n_unbound != n_motors_ - n_single_bound_ - n_double_bound_){
		printf("something awful in update_bound_list bruh\n");
		exit(1);
	}
}

Kinesin* KinesinManagement::GetUnboundMotor(){

	// Make sure an unbound motor exists
	if(n_single_bound_ + n_double_bound_ != n_motors_){
        // Randomly pick a motor from the reservoir
        int i_motor = properties_->gsl.GetRanInt(n_motors_);
        Kinesin *motor = &motor_list_[i_motor];
        while(motor->bound_ == true){
            i_motor++;
            if(i_motor == n_motors_)
				i_motor = 0;
			motor = &motor_list_[i_motor];
		}
		UnboundCheck(motor);
		return motor;
	}
	else{
		printf("ERROR IN GET UNBOUND MOTOR\n");
		exit(1);
	}	
}


void KinesinManagement::RunDiffusion(){

	UpdateLists();
	int mt_length = parameters_->length_of_microtubule;
	// List of step directions for all motors (will be shuffled)
	int step_dir [n_double_bound_];
	if(n_double_bound_ > 0){
		int n_half = (int) (n_double_bound_/2);
		// First half steps backwards
		for(int i = 0; i < n_half; i++){
			step_dir[i] = -1;
		}
		// Second half steps forwards
		for(int i = n_half; i < n_double_bound_; i++){
			step_dir[i] = 1;
		}
		gsl_ran_shuffle(properties_->gsl.rng, step_dir, n_double_bound_, sizeof(int));
	}
	// Run through bound_list and step all motors appropriately 
	for(int i_motor = 0; i_motor < n_double_bound_; i_motor++){
		Kinesin *motor = bound_list_[i_motor];
		int dx = step_dir[i_motor] * motor->mt_->delta_x_;
		int i_front = motor->front_site_->index_;
		int i_rear = motor->rear_site_->index_;
		// Don't let motors step off boundaries and into the void
		if(!((i_front == 0 || i_rear == 0) && dx == -1)
		&& !((i_rear == mt_length - 1 || i_front == mt_length - 1) && dx == 1)){
			Tubulin *old_front_site = motor->front_site_;
			Tubulin *old_rear_site = motor->rear_site_;
			if(dx == motor->mt_->delta_x_){
				Tubulin *new_front_site = &motor->mt_->lattice_[i_front + dx];
				Tubulin *new_rear_site = old_front_site;
				if(new_front_site->occupied_ == false){
					old_rear_site->motor_ = nullptr;
					old_rear_site->occupied_ = false;
					new_front_site->motor_ = motor;
					new_front_site->occupied_ = true;
					motor->front_site_ = new_front_site;
					motor->rear_site_ = new_rear_site;
				}
			}
			else{
				Tubulin *new_front_site = old_rear_site;
				Tubulin *new_rear_site = &motor->mt_->lattice_[i_rear - dx];
				if(new_rear_site->occupied_ == false){
					old_front_site->motor_ = nullptr;
					old_front_site->occupied_ = false;
					new_rear_site->motor_ = motor;
					new_rear_site->occupied_ = true;
					motor->front_site_ = new_front_site;
					motor->rear_site_ = new_rear_site;
				}
			}
		}
	}
}

void KinesinManagement::GenerateKMCList(){

	int n_bind_i = GetNumToBind_I_Free();
	int n_bind_i_tethered = GetNumToBind_I_Tethered();
	int n_unbind_i = GetNumToUnbind_I();
	int n_tether_unbound = GetNumToTether_Unbound();
	int n_switch = GetNumToSwitch();	
	int n_step = GetNumToStep();
	int n_events = n_bind_i + n_bind_i_tethered + n_unbind_i + n_tether_unbound + n_switch + 
				   n_step;
	if(n_events > 0){
		int pre_list[n_events];
		int i_event = 0;
		for(i_event; i_event < n_bind_i; i_event++){
			pre_list[i_event] = 0;
		}
		for(i_event; i_event < n_bind_i + n_bind_i_tethered; i_event++){
			pre_list[i_event] = 1;
		}
		for(i_event; i_event < n_bind_i + n_bind_i_tethered + n_unbind_i; i_event++){
			pre_list[i_event] = 2;
		}
		for(i_event; i_event < n_events - n_step - n_switch; i_event++){
			pre_list[i_event] = 3;
		}
		for(i_event; i_event < n_events - n_step; i_event++){
			pre_list[i_event] = 4;
		}
		for(i_event; i_event < n_events; i_event++){
			pre_list[i_event] = 5;
		}
		gsl_ran_shuffle(properties_->gsl.rng, pre_list, n_events, sizeof(int));
		kmc_list_.resize(n_events);
		for(int j_event = 0; j_event < n_events; j_event++){
			kmc_list_[j_event] = pre_list[j_event];
		}
	}
	else{
		kmc_list_.clear();
	}
}

int KinesinManagement::GetNumToBind_I_Free(){

	properties_->microtubules.UpdateNumUnoccupiedPairs();
	int n_unoccupied_pairs = properties_->microtubules.n_unoccupied_pairs_;
	double n_avg = p_bind_*n_unoccupied_pairs;
	int n_to_bind = properties_->gsl.SamplePoissonDist(n_avg);
	return n_to_bind;
}

int KinesinManagement::GetNumToBind_I_Tethered(){
	
	int n_avg = 0;
	int n_to_bind_i_tethered = 0;
	for(int i_motor = 0; i_motor > n_tethered_unbound_; i_motor++){
		Kinesin *motor = tethered_unbound_list_[i_motor];
		AssociatedProtein *xlink = motor->xlink_;
//		int i_anchor = xlink->GetAnchorIndex();
		int i_anchor;
		// standard error check
		if(xlink->heads_active_ == 0){
			printf("not good in gnb_i_teth bro\n");
			exit(1);
		}
		// use index of only site if single-bound (distribution 
		// of second head's diffusion will average out)
		else if(xlink->heads_active_ == 1){
			if(xlink->site_one_ != nullptr){
				i_anchor = xlink->site_one_->index_;		
			}
			else if(xlink->site_two_ != nullptr){
				i_anchor = xlink->site_two_->index_;
			}
			else
				printf("ummm... getnumtobinditethered.\n");
		}
		// Avg index of head sites of double-bound (systematic
		// rounding-down should average out as well)
		else if(xlink->heads_active_ == 2){
			int i_first = xlink->site_one_->index_;
			int i_second = xlink->site_two_->index_;
			i_anchor = (i_first + i_second)/2; 
		}
		// Calculate this xlink's contribution to expected binds via tethers
		int n_entries = xlink->n_neighbor_motors_;
		for(int i_entry = 0; i_entry < n_entries; i_entry++){

		}
	}
	return n_to_bind_i_tethered;
}

int KinesinManagement::GetNumToUnbind_I(){
	
	int n_to_unbind = properties_->gsl.SampleBinomialDist(p_unbind_, n_double_bound_);
	return n_to_unbind;
}

int KinesinManagement::GetNumToTether_Unbound(){

	int n_untethered_xlinks = properties_->prc1.n_untethered_;
	double n_tether_avg = p_tether_unbound_*n_untethered_xlinks; 	
	int n_to_tether_unbound = properties_->gsl.SamplePoissonDist(n_tether_avg);
	return n_to_tether_unbound;
}

int KinesinManagement::GetNumToSwitch(){

	int mt_length = parameters_->length_of_microtubule;
	int n_to_switch;
	int n_switchable = 0;
	int n_unbound = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		if(motor_list_[i_motor].bound_ == true){
			Kinesin *motor = &motor_list_[i_motor];
			Microtubule *mt = motor->mt_;
			Microtubule *mt_adj = mt->neighbor_; //define neighb u idiot
			int mt_coord = mt->coord_;
			int mt_adj_coord = mt_adj->coord_;
			int offset = mt_adj_coord - mt_coord;
			int i_front_site = motor->front_site_->index_;
			int i_rear_site_adj = i_front_site - offset;
			int i_rear_site = motor->rear_site_->index_;
			int i_front_site_adj = i_rear_site - offset;
			// Exclude boundary sites from statistics (no switching there)
			// and non-overlapping microtubule segments
			if((i_front_site != mt->plus_end_ && i_rear_site != mt->minus_end_)
			&& (i_front_site_adj > 0 && i_front_site_adj < mt_length - 1)
			&& (i_rear_site_adj > 0 && i_rear_site_adj < mt_length - 1)){
				if(mt_adj->lattice_[i_front_site_adj].occupied_ == false
				&& mt_adj->lattice_[i_rear_site_adj].occupied_ == false){
					n_switchable++;
/*					printf("%i_%i/%i; %i_%i/%i\n", mt->index_, 
					i_front_site, i_rear_site, mt_adj->index_, 
					i_front_site_adj, i_rear_site_adj);
*/				}
			}

		}
		else{
			n_unbound++;
		}
	}
	if(n_unbound != n_motors_ - n_double_bound_){
		printf("Something went wrong in Kinesin MGMT!\n");
		exit(1);
	}
	n_to_switch = properties_->gsl.SampleBinomialDist(p_switch_, n_switchable);
//	printf("%i (out of %i)\n", n_to_switch, n_switchable);
	return n_to_switch;
}

int KinesinManagement::GetNumToStep(){

	int n_to_step; 
	int n_steppable = 0;
	int n_unbound = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
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
	if(n_unbound != n_motors_ - n_double_bound_){
		printf("Something went wrong in Kinesin MGMT (2):\n");
		printf("  %i != %i - %i\n", n_unbound, n_motors_, n_double_bound_); 
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
				case 0:	KMC_Bind_I_Free();
						break;
				case 1: KMC_Bind_I_Tethered();
						break;
				case 2:	KMC_Unbind_I();
						break;
				case 3: KMC_Tether_Unbound();
						break;
				case 4: KMC_Switch();
						break;
				case 5: KMC_Step();
						break;
			}
			KMC_Boundaries(n_events);
		}
	}
	else{
		KMC_Boundaries(1);
	}
}

void KinesinManagement::KMC_Bind_I_Free(){

	// Make sure that at least one unbound motor exists
	if(n_double_bound_ != n_motors_){
        Kinesin *motor = GetUnboundMotor();
		properties_->microtubules.UpdateUnoccupiedPairsList();
		if(properties_->microtubules.n_unoccupied_pairs_ > 0){
			// Find unoccupied pair of sites to bind to
			Tubulin *rear_site = properties_->microtubules.GetUnoccupiedPair_1();
			Tubulin *front_site = properties_->microtubules.GetUnoccupiedPair_2(rear_site);
			Microtubule *mt = front_site->mt_;
			// Place motor on sites
			rear_site->motor_ = motor;
			rear_site->occupied_ = true;
			front_site->motor_ = motor;
			front_site->occupied_ = true;
			// Update motor details
			motor->bound_ = true;
			motor->mt_ = mt;
			motor->rear_site_ = rear_site;
			motor->front_site_ =  front_site;
			// Update statistics
			n_double_bound_++;
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

void KinesinManagement::KMC_Bind_I_Tethered(){

}

void KinesinManagement::KMC_Unbind_I(){

	// Make sure that at least one bound motor exists
	if(n_double_bound_ > 0){
		UpdateLists();
		// Randomly pick a bound motor
		int i_entry = properties_->gsl.GetRanInt(n_double_bound_);
		Kinesin* motor = bound_list_[i_entry];
		BoundCheck(motor);
		int n_attempts = 0;
		bool failure = false;
		// Ensure that motors on boundary sites are excluded
		while(BoundaryStatus(motor) == true){
			if(n_attempts > 2*n_double_bound_){
				failure = true;
				break;
			}
			i_entry = properties_->gsl.GetRanInt(n_double_bound_);
			motor = bound_list_[i_entry];
			BoundCheck(motor);
			n_attempts++;
		}	
		// Make sure an eligible motor was found
		if(failure == false){
			// Remove motor from sites
			motor->front_site_->motor_ = nullptr;
			motor->front_site_->occupied_ = false;
			motor->rear_site_->motor_ = nullptr;
			motor->rear_site_->occupied_ = false;
			// Update motor details
			motor->bound_ = false;
			motor->mt_ = nullptr;
			motor->front_site_ = nullptr;
			motor->rear_site_ = nullptr;
			// Update statistics
			n_double_bound_--;
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

void KinesinManagement::KMC_Tether_Unbound(){

	// Make sure there is at least one unbound motor
	if(n_single_bound_ + n_double_bound_ != n_motors_){
        Kinesin *motor = GetUnboundMotor();
		if(properties_->prc1.n_untethered_ > 0){
			AssociatedProtein *xlink = properties_->prc1.GetUntetheredXlink();
			// Update motor and xlink details
			motor->xlink_ = xlink;
			motor->tethered_ = true;
			xlink->tethered_ = true;
			xlink->motor_ = motor; 	
			// Update statistics 
			n_tethered_unbound_++;
			properties_->prc1.n_tethered_++;
			properties_->prc1.n_untethered_--;
		}
	}
	else{
		printf("Error in KMC_Tether_Unbound: no unbound motors\n");
		exit(1);
	}
}

void KinesinManagement::KMC_Switch(){

	int mt_length = parameters_->length_of_microtubule;
	// Make sure there is at least one bound motor
	if(n_double_bound_ > 0){
		UpdateLists();
		int i_entry, offset, i_front_new, i_rear_new;
		Kinesin *motor;
		Tubulin *old_front_site, *old_rear_site, 
				*new_front_site, *new_rear_site;
		Microtubule *mt, *mt_adj;
		int n_attempts = 0;
		bool failure = false;
		bool flag = false;
		// Ensure that blocked motors and motors on the plus end are excluded
		do{	
			if(n_attempts > 2*n_double_bound_){
				failure = true;
				break;
			}
			i_entry = properties_->gsl.GetRanInt(n_double_bound_);
			motor = bound_list_[i_entry];
			BoundCheck(motor);
			old_front_site = motor->front_site_;
			old_rear_site = motor->rear_site_;
			mt = motor->mt_;
			mt_adj = mt->neighbor_;
			offset = mt_adj->coord_ - mt->coord_;
			i_front_new = old_rear_site->index_ - offset;
			i_rear_new = old_front_site->index_ - offset;
			while((i_front_new <= 0 || i_front_new >= mt_length - 1)
			|| (i_rear_new <= 0 || i_rear_new >= mt_length - 1)){
				if(n_attempts > 2*n_double_bound_){
					failure = true;
					break;
				}
				i_entry = properties_->gsl.GetRanInt(n_double_bound_);
				motor = bound_list_[i_entry];
				BoundCheck(motor);
				old_front_site = motor->front_site_;
				old_rear_site = motor->rear_site_;
				mt = motor->mt_;
				mt_adj = mt->neighbor_;
				offset = mt_adj->coord_ - mt->coord_;
				i_front_new = old_rear_site->index_ - offset;
				i_rear_new = old_front_site->index_ - offset;
				n_attempts++;
			}
			new_front_site = &mt_adj->lattice_[i_front_new];
			new_rear_site = &mt_adj->lattice_[i_rear_new];
			n_attempts++;
		}while(BoundaryStatus(motor) == true
		    || new_front_site->occupied_ == true
			|| new_rear_site->occupied_ == true);
		// Only attempt to switch if an eligible motor was found
		if(failure != true){
			properties_->microtubules.OccupiedCheck(old_front_site);
			properties_->microtubules.OccupiedCheck(old_rear_site);
			properties_->microtubules.UnoccupiedCheck(new_rear_site);
			properties_->microtubules.UnoccupiedCheck(new_front_site);
			// Remove motor from old sites; Place it on new sites
			old_front_site->occupied_ = false;
			old_front_site->motor_ = nullptr;
			old_rear_site->occupied_ = false;
			old_rear_site->motor_ = nullptr;
			new_front_site->occupied_ = true;
			new_front_site->motor_ = motor;
			new_rear_site->occupied_ = true;
			new_rear_site->motor_ = motor;
			// Update motor details
			motor->front_site_ = new_front_site;
			motor->rear_site_ = new_rear_site;
			motor->mt_ = mt_adj;
		}	
		else{
			int ID = motor->ID_;
			int i_mt = motor->mt_->index_;
			int i_front = motor->front_site_->index_;
			int i_rear = motor->rear_site_->index_;
			printf("Motor %i failed to switch @ %i_%i/%i\n", ID, i_mt, i_front, i_rear);
		}
	}
	else{
		printf("Error in RunKMC_Switch(): bound_list is empty\n");
		exit(1);
	}
}

void KinesinManagement::KMC_Step(){

	// Make sure there is at least one bound motor
	if(n_double_bound_ > 0){
		UpdateLists();
		int i_entry, i_plus_end, delta_x;
		Kinesin *motor;
		Tubulin *front_site, *rear_site;
		int n_attempts = 0;
		bool failure = false;
		// Ensure that blocked motors and motors on the plus end are excluded
		do{	
			if(n_attempts > 2*n_double_bound_){
				failure = true;
				break;
			}
			i_entry = properties_->gsl.GetRanInt(n_double_bound_);
			motor = bound_list_[i_entry];
			BoundCheck(motor);
			front_site = motor->front_site_;
			rear_site = motor->rear_site_;
			i_plus_end = motor->mt_->plus_end_;
			delta_x = motor->mt_->delta_x_;
			n_attempts++;
			// Ensure chosen motor is not on plus_end (enforce end-pausing)
			while(front_site->index_ == i_plus_end){
				if(n_attempts > 2*n_double_bound_){
					failure = true;
					break;
				}
				i_entry = properties_->gsl.GetRanInt(n_double_bound_);
				motor = bound_list_[i_entry];
				BoundCheck(motor);
				front_site = motor->front_site_;
				rear_site = motor->rear_site_;
				i_plus_end = motor->mt_->plus_end_;
				delta_x = motor->mt_->delta_x_;
				n_attempts++;
			}
		}while(motor->mt_->lattice_[front_site->index_ + delta_x].occupied_ == true);
		// Make sure an eligible motor was found
		if(failure == false){
			Tubulin *old_front_site = front_site;
			Tubulin *old_rear_site = rear_site;
			properties_->microtubules.OccupiedCheck(old_front_site);
			properties_->microtubules.OccupiedCheck(old_rear_site);
			Tubulin *new_front_site = &motor->mt_->lattice_[old_front_site->index_ + delta_x];
			Tubulin *new_rear_site = old_front_site;
			properties_->microtubules.UnoccupiedCheck(new_front_site);
			// Take rear motor head off of old site; place it on new site
			old_rear_site->motor_ = nullptr;
			old_rear_site->occupied_ = false;
			motor->rear_site_ = old_front_site;
			// Take front motor head off of old site; place it on new  site
			motor->front_site_ = new_front_site;
			new_front_site->motor_ = motor;
			new_front_site->occupied_ = true;
		}
		else{
			int ID = motor->ID_;
			int i_mt = motor->mt_->index_;
			int i_front = motor->front_site_->index_;
			int i_rear = motor->rear_site_->index_;
			printf("Motor %i failed to step @ %i_%i/%i\n", ID, i_mt, i_front, i_rear);
		}
	}
	else{
		printf("Error in RunKMC_Step(): bound_list is empty\n");
		exit(1);
	}
}

void KinesinManagement::KMC_Boundaries(int n_events){

	int n_mts = parameters_->n_microtubules;
	double alpha_eff = (alpha_*p_step_)/n_events;
	double beta_eff = (beta_*p_step_)/n_events;
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
		int i_plus_end = mt->plus_end_;
		int i_minus_end = mt->minus_end_;
		int delta_x = mt->delta_x_;
		Tubulin *minus_end = &mt->lattice_[i_minus_end];
		Tubulin *minus_neighbor = &mt->lattice_[i_minus_end + delta_x];
		Tubulin *plus_end = &mt->lattice_[i_plus_end];
		Tubulin *plus_neighbor = &mt->lattice_[i_plus_end - delta_x];
		double ran1 = properties_->gsl.GetRanProb();
		double ran2 = properties_->gsl.GetRanProb();
		// Insert motors onto minus_end with probability alpha_eff	
		if(ran1 < alpha_eff
		&& minus_end->occupied_ == false
		&& minus_neighbor->occupied_ == false){
        	Kinesin *motor = GetUnboundMotor();
			// Place motor on sites
			minus_end->motor_ = motor;
			minus_end->occupied_ = true;
			minus_neighbor->motor_ = motor;
			minus_neighbor->occupied_ = true;
			// Update motor details
			motor->bound_ = true;
			motor->mt_ = minus_end->mt_;
			motor->rear_site_ = minus_end;
			motor->front_site_ =  minus_neighbor;
			// Update statistics
			n_double_bound_++;
		}
		// Remove motors from plus_end with probability beta_eff
		if(ran2 < beta_eff
		&& plus_end->motor_ != nullptr
		&& plus_neighbor->motor_ != nullptr){
			// Get motor bound to minus_end
			Kinesin *motor = plus_end->motor_;
			BoundCheck(motor);
			// Remove motor from sites
			plus_end->motor_ = nullptr;
			plus_end->occupied_ = false;
			plus_neighbor->motor_ = nullptr;
			plus_neighbor->occupied_ = false;
			// Update motor details
			motor->bound_ = false;
			motor->mt_ = nullptr;
			motor->front_site_ = nullptr;
			motor->rear_site_ = nullptr;
			// Update statistics
			n_double_bound_--;
		} 
	}
}

