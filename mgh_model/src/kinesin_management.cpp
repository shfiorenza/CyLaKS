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
	InitiateLists();
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
}

void KinesinManagement::SetParameters(){

    double k_on = parameters_->k_on;
    double c_motor = parameters_->c_motor;
    double k_off = parameters_->k_off;
    double switch_rate = parameters_->switch_rate;
    double motor_speed = parameters_->motor_speed;
    double delta_t = parameters_->delta_t;
    p_bind_i_free_ = k_on*c_motor*delta_t;
    p_bind_i_tethered_ = p_bind_i_free_ * 4;
	p_bind_ii_ = 1; 					// FIXME should be like ~ tau/delta_t
    p_unbind_ = k_off * delta_t;
	p_unbind_tethered_ = p_unbind_ / 2; 
    p_tether_free_ = p_bind_i_free_;  				//FIXME
	p_untether_free_ = p_unbind_ / 2;
	p_step_untethered_ = motor_speed*delta_t;
	p_switch_ = switch_rate*delta_t;
   	// Generate untethering and stepping rates for all tether extensions	
	// Everything is 2*dist_cutoff to allow for half-integer distances, 
	// so the 3rd entry will correspond to a distance of 1.5, etc. 
	dist_cutoff_ = (int) motor_list_[0].dist_cutoff_;
	p_untether_bound_.resize(2*dist_cutoff_ + 1);
	p_step_tethered_.resize(2*dist_cutoff_ + 1);
	for(int x_dist_doubled = 0; x_dist_doubled <= 2*dist_cutoff_; x_dist_doubled++){
		double site_size = motor_list_[0].site_size_;
		double r_0 = motor_list_[0].r_0_;
		double kbT = motor_list_[0].kbT_;
		double k_spring = motor_list_[0].k_spring_;
		double r_x = x_dist_doubled*site_size/2;
		double r_y = 17.5;		// in nm; from MT to midpoint of prc1
		double r = sqrt(r_y*r_y + r_x*r_x);
		double dr = r - r_0;
		double stall_force = 5; 	// in pN
		// If extension is positive, treat tether like a spring
		if(dr >= 0){
			// Untether probability is same as xlink unbinding
			double k_untether = exp(dr*dr*k_spring/(2*kbT));
			p_untether_bound_[x_dist_doubled] = k_untether*delta_t;
			// Stepping probability has force-dependent relation
			double force = dr*k_spring;
			if(force < stall_force){
				// see sci advancements supplement
				double cosine = r_x/r;
				double coeff = 1 - cosine*(force/stall_force);
				p_step_tethered_[x_dist_doubled] = coeff*p_step_untethered_;
			}
			else{
				p_step_tethered_[x_dist_doubled] = 0;
			}
		}
		else{
			p_untether_bound_[x_dist_doubled] = p_untether_free_;
			p_step_tethered_[x_dist_doubled] = p_step_untethered_;
		}
	}
    alpha_ = parameters_->alpha;
    beta_ = parameters_->beta;
}

void KinesinManagement::InitiateLists(){

	free_tethered_list_.resize(n_motors_);
	pseudo_bound_list_.resize(n_motors_);
	bound_list_.resize(n_motors_);
	bound_tethered_list_.resize(2*dist_cutoff_ + 1);
	n_bound_tethered_.resize(2*dist_cutoff_ + 1);
	for(int x_dist_doubled; x_dist_doubled <= 2*dist_cutoff_; x_dist_doubled++){
		bound_tethered_list_[x_dist_doubled].resize(n_motors_);
		n_bound_tethered_[x_dist_doubled] = 0;
	}
	bound_untethered_list_.resize(n_motors_);
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

void KinesinManagement::UpdateFreeTetheredList(){

	int i_entry = 0;
	int n_entries = 0; 
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motor_list_[i_motor];
		if(motor->heads_active_ == 0
		&& motor->tethered_ == true){
			free_tethered_list_[i_entry] = motor;
			i_entry++;
			n_entries++;
		}
	}
	if(n_entries != n_free_tethered_){
		printf("something wrong in update_free_tethered.\n");
		exit(1);
	}
}

void KinesinManagement::UpdatePseudoBoundList(){
	
	int i_entry = 0;
	int n_entries = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motor_list_[i_motor];
		if(motor->heads_active_ == 1){
			pseudo_bound_list_[i_entry] = motor;
			i_entry++;
			n_entries++;
		}
	}
	if(n_entries != n_single_bound_){
		printf("something wrong in update_pseudo_bound (motors)\n");
		exit(1);
	}
}

void KinesinManagement::UpdateBoundList(){
	
	int i_entry = 0;
	int n_entries = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motor_list_[i_motor];
		if(motor->heads_active_ == 2){
			bound_list_[i_entry] = motor;
			i_entry++;
			n_entries++;
		}
	}
	if(n_entries != n_double_bound_){
		printf("something wrong in update_bound_list (motors)\n");
		exit(1);
	}
}

void KinesinManagement::UpdateBoundTetheredList(){

	int i_entry[2*dist_cutoff_ + 1];
	int n_entries[2*dist_cutoff_ + 1]; 
	for(int x_dist_doubled = 0; x_dist_doubled <= 2*dist_cutoff_; x_dist_doubled++){
		i_entry[x_dist_doubled] = 0;
		n_entries[x_dist_doubled] = 0;
	}
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motor_list_[i_motor];
		if(motor->heads_active_ == 2
		&& motor->tethered_ == true){
			int x_dist_doubled = motor->x_dist_doubled_; 
			int index = i_entry[x_dist_doubled]; 
			bound_tethered_list_[x_dist_doubled][index] = motor;
			i_entry[x_dist_doubled]++;
			n_entries[x_dist_doubled]++;
		}
	}
	for(int x_dist_doubled = 0; x_dist_doubled <= dist_cutoff_; x_dist_doubled++){
		if(n_entries[x_dist_doubled] != n_bound_tethered_[x_dist_doubled]){
			printf("something wrong in update_bound_tethered (motors)\n");
			exit(1);
		}
	}
}

void KinesinManagement::UpdateBoundUntetheredList(){

	int i_entry = 0;
	int n_entries = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motor_list_[i_motor];
		if(motor->heads_active_ == 2
		&& motor->tethered_ == false){
			bound_untethered_list_[i_entry] = motor;
			i_entry++;
			n_entries++;
		}
	}
	if(n_entries != n_bound_untethered_){
		printf("something wrong in update_bound_untethered *motor)\n");
		exit(1);
	}
}

void KinesinManagement::UpdateExtensions(){

	// Clear extension-based statistics
	for(int x_dist_doubled = 0; x_dist_doubled <= 2*dist_cutoff_; x_dist_doubled++){
		n_bound_tethered_[x_dist_doubled] = 0;	
	}
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motor_list_[i_motor];
		if(motor->heads_active_ == 2
		&& motor->tethered_ == true){
			AssociatedProtein *xlink = motor->xlink_;
			double stalk_coord = motor->GetStalkCoordinate();
			double anchor_coord = xlink->GetAnchorCoordinate();
			double x_dist_doubled = 2*fabs(anchor_coord - stalk_coord);
			// If new x_dist_doubled exceeds cutoff, force an untether event
			if(x_dist_doubled > 2*dist_cutoff_){
				motor->ForceUntether();
				n_bound_untethered_++;
			}
			else{
				motor->x_dist_doubled_ = x_dist_doubled;
				motor->UpdateExtension();
				n_bound_tethered_[x_dist_doubled]++; 
			}
		}
	}
}

Kinesin* KinesinManagement::GetFreeMotor(){

	// Randomly pick a motor from the reservoir
	int i_motor = properties_->gsl.GetRanInt(n_motors_);
	Kinesin *motor = &motor_list_[i_motor];
	int attempts = 0;
	while(motor->heads_active_ > 0 
			|| motor->tethered_ == true){
		i_motor++;
		if(i_motor == n_motors_)
			i_motor = 0;
		motor = &motor_list_[i_motor];
		attempts++;
		if(attempts > n_motors_){
			printf("error in get free motor\n");
			exit(1);
		}
	}
	UnboundCheck(motor);
	return motor;
}

void KinesinManagement::RunDiffusion(){

	int i_step = properties_->current_step_;
	// FIXME currently only valid for delta_t = 0.0005 FIXME
	int i_tau_bound = 12;
	if(i_step % i_tau_bound == 0)
		RunDiffusion_Bound();
	UpdateExtensions();	
}

void KinesinManagement::RunDiffusion_Bound(){

    UpdateBoundList();
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
        int dx = step_dir[i_motor];
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

	/* see xlinks analog for a much simpler case and explanation */
	int n_events = 0;
    int n_bind_i_free = GetNumToBind_I_Free();
	n_events += n_bind_i_free;
    int n_bind_i_tethered = GetNumToBind_I_Tethered();
	n_events += n_bind_i_tethered;
	int n_bind_ii = GetNumToBind_II();
	n_events += n_bind_ii;
	int n_unbind_ut = GetNumToUnbind_Untethered();
	n_events += n_unbind_ut;
	int n_unbind_t = GetNumToUnbind_Tethered();
	n_events += n_unbind_t;
	int n_tether_free = GetNumToTether_Free();
	n_events += n_tether_free;
	int n_tether_bound = GetNumToTether_Bound();
	n_events += n_tether_bound;
    int n_switch = GetNumToSwitch();	
	n_events += n_switch;
	int n_untether_free = GetNumToUntether_Free();
	n_events += n_untether_free;
	int n_step_untethered = GetNumToStep_Untethered(); 
	n_events += n_step_untethered;
	int n_untether_bound[2*dist_cutoff_ + 1];
	int n_step_tethered[2*dist_cutoff_ + 1];
	for(int x_dist_doubled = 0; x_dist_doubled <= 2*dist_cutoff_; x_dist_doubled++){
			n_untether_bound[x_dist_doubled] = GetNumToUntether_Bound(x_dist_doubled);
			n_step_tethered[x_dist_doubled] = GetNumToStep_Tethered(x_dist_doubled);
			n_events += n_untether_bound[x_dist_doubled];
			n_events += n_step_tethered[x_dist_doubled];
	}
    if(n_events > 0){
        int pre_list[n_events];
		int kmc_index = 0;
		for(int i_event = 0; i_event < n_bind_i_free; i_event++){
			pre_list[kmc_index] = 10;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_bind_i_tethered; i_event++){
			pre_list[kmc_index] = 11;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_bind_ii; i_event++){
			pre_list[kmc_index] = 12;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_unbind_ut; i_event++){
			pre_list[kmc_index] = 20;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_unbind_t; i_event++){
			pre_list[kmc_index] = 21;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_tether_free; i_event++){
			pre_list[kmc_index] = 30;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_tether_bound; i_event++){
			pre_list[kmc_index] = 31;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_switch; i_event++){
			pre_list[kmc_index] = 40;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_untether_free; i_event++){
			pre_list[kmc_index] = 50;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_step_untethered; i_event++){
			pre_list[kmc_index] = 60;
			kmc_index++;
		}
		for(int x_dist_doubled = 0; x_dist_doubled <= 2*dist_cutoff_; x_dist_doubled++){
			for(int i_event = 0; i_event < n_untether_bound[x_dist_doubled]; i_event++){
				pre_list[kmc_index] = 500 + x_dist_doubled; 	
				kmc_index++;
			}
			for(int i_event = 0; i_event < n_step_tethered[x_dist_doubled]; i_event++){
				pre_list[kmc_index] = 600 + x_dist_doubled;
				kmc_index++;
			}
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

int KinesinManagement::GetNumToBind_I_Free(){

    properties_->microtubules.UpdateNumUnoccupiedPairs();
    int n_unoccupied_pairs = properties_->microtubules.n_unoccupied_pairs_;
    double n_avg = p_bind_*n_unoccupied_pairs;
    int n_to_bind = properties_->gsl.SamplePoissonDist(n_avg);
	printf("%i\n", n_to_bind);
    return n_to_bind;
}

int KinesinManagement::GetNumToBind_I_Tethered(){

    UpdateTetheredFreeList();
    double weights_summed = 0;
    // Sum over all tethered but unbound motors
    for(int i_motor = 0; i_motor < n_tethered_free_; i_motor++){
        Kinesin *motor = tethered_free_list_[i_motor];
        if(motor->tethered_ == true)
            motor->UpdateNeighborSites();
        // Get weight of all neighbor sites
        int n_neighbors = motor->n_neighbor_sites_;
        for(int i_neighb = 0; i_neighb < n_neighbors; i_neighb++){
            Tubulin *site = motor->neighbor_sites_[i_neighb];
            double weight = motor->GetBindingWeight(site);
            weights_summed += weight;
        }
    }
    return 0;
}

int KinesinManagement::GetNumToUnbind_I(){

    int n_to_unbind = properties_->gsl.SampleBinomialDist(p_unbind_, n_double_bound_);
    return n_to_unbind;
}

int KinesinManagement::GetNumToTether_Free(){

    int n_untethered_xlinks = properties_->prc1.n_untethered_;
    double n_tether_avg = p_tether_unbound_*n_untethered_xlinks; 	
    int n_to_tether_unbound = properties_->gsl.SamplePoissonDist(n_tether_avg);
    //	printf("%i\n", n_untethered_xlinks);
    return n_to_tether_unbound;
}

int KinesinManagement::GetNumToSwitch(){

    int mt_length = parameters_->length_of_microtubule;
    int n_to_switch;
    int n_switchable = 0;
    int n_unbound = 0;
    for(int i_motor = 0; i_motor < n_motors_; i_motor++){
        if(motor_list_[i_motor].heads_active_ == 2){
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
        if(motor_list_[i_motor].heads_active_ == 2){
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
        for(int i_event = 0; i_event < n_events; i_event++){
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
	printf("hi\n");
    if(n_double_bound_ != n_motors_){
        Kinesin *motor = GetFreeMotor();
		printf("heyy\n");
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
            motor->heads_active_ = 2;
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
            motor->heads_active_ = 0;
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
        Kinesin *motor = GetFreeMotor();
        if(properties_->prc1.n_untethered_ > 0){
            AssociatedProtein *xlink = properties_->prc1.GetUntetheredXlink();
            // Update motor and xlink details
            motor->xlink_ = xlink;
            motor->tethered_ = true;
            xlink->tethered_ = true;
            xlink->motor_ = motor; 	
            // Update statistics 
            n_tethered_free_++;
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
        }while(motor->mt_->
				lattice_[front_site->index_ + delta_x].occupied_ == true);
        // Make sure an eligible motor was found
        if(failure == false){
            Tubulin *old_front_site = front_site;
            Tubulin *old_rear_site = rear_site;
            properties_->microtubules.OccupiedCheck(old_front_site);
            properties_->microtubules.OccupiedCheck(old_rear_site);
            Tubulin *new_front_site = &motor->mt_->lattice_[old_front_site->
				index_ + delta_x];
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
            printf("Motor %i failed to step @ %i_%i/%i\n", 
					ID, i_mt, i_front, i_rear);
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
            Kinesin *motor = GetFreeMotor();
            // Place motor on sites
            minus_end->motor_ = motor;
            minus_end->occupied_ = true;
            minus_neighbor->motor_ = motor;
            minus_neighbor->occupied_ = true;
            // Update motor details
            motor->heads_active_ = 2;
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
            motor->heads_active_ = 0;
            motor->mt_ = nullptr;
            motor->front_site_ = nullptr;
            motor->rear_site_ = nullptr;
            if(motor->tethered_ == true){
                motor->xlink_->motor_ = nullptr;
                motor->xlink_->tethered_ = false;
                motor->tethered_ = false;
                motor->xlink_ = nullptr;
            }
            // Update statistics
            n_double_bound_--;
        } 
    }
}

