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
	// Statistics for diffusion
	p_diffuse_fwd_untethered_ = 0.5 * delta_t / tau_;
	p_diffuse_bck_untethered_ = 0.5 * delta_t / tau_; 
	// Generate stepping rates based on extension of tether: rates are 
	// increased if stepping towards rest length, and reduced if stepping
	// away from it and increasing tether extension as a result
	dist_cutoff_ = (int) motor_list_[0].dist_cutoff_;
	p_diffuse_to_tether_.resize(2*dist_cutoff_ + 1);
	p_diffuse_from_tether_.resize(2*dist_cutoff_ + 1);
	double site_size = motor_list_[0].site_size_;
	double r_0 = motor_list_[0].r_0_;
	double kbT = motor_list_[0].kbT_;
	double k_spring = motor_list_[0].k_spring_;
	double k_eff_slack = motor_list_[0].k_eff_slack_;
	double r_y = 17.5;		// in nm; from MT to midpoint of prc1
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		// Calculate tether length for this x_dist as well as +/- 1 it
		double r_x = x_dist_dub * site_size / 2;
		double r_x_fwd = (x_dist_dub + 1) *  site_size / 2;
		double r_x_bck = (x_dist_dub - 1) * site_size / 2;
		double r = sqrt(r_y*r_y + r_x*r_x);
		double r_fwd = sqrt(r_y*r_y + r_x_fwd*r_x_fwd);
		double r_bck = sqrt(r_y*r_y + r_x_bck*r_x_bck);
		// Calculate extension of tether for given x_dist_dub
		double dr = r - r_0; 
		// Calculate extension if motor diffuses toward xlink
		double dr_toward = r_bck - r_0;
		// Calculate extension if motor diffuses away from xlink
		double dr_from = r_fwd - r_0;
		// use k_spring if extension is positive
		if(dr >= 0){
			// Get the corresponding changes in potential energy
			double dU_to = (k_spring/2)*(dr_toward*dr_toward - dr*dr);
			double dU_from = (k_spring/2)*(dr_from*dr_from - dr*dr);
			// Weights according to Lanksy et al. 
			double weight_towards = exp(-dU_to/(2*kbT));
			double weight_away = exp(-dU_from/(2*kbT));
			double p_to = weight_towards * 0.5 * delta_t / tau_;
			double p_from = weight_away * 0.5 * delta_t / tau_;
			p_diffuse_to_tether_[x_dist_dub] = p_to;
			if(x_dist_dub < 2*dist_cutoff_)
				p_diffuse_from_tether_[x_dist_dub] = p_from;
			else
				p_diffuse_from_tether_[x_dist_dub] = 0;
		}
		// otherwise, use k_eff to model 'slack' in the tether
		else{
			// Get the corresponding changes in potential energy
			double dU_to = (k_eff_slack/2)*(dr_toward*dr_toward - dr*dr);
			double dU_from = (k_eff_slack/2)*(dr_from*dr_from - dr*dr);
			// Weights according to Lanksy et al.
			double weight_towards = exp(-dU_to/(2*kbT));
			double weight_away = exp(-dU_from/(2*kbT));
			double p_to = weight_towards * 0.5 * delta_t / tau_;
			double p_from = weight_away * 0.5 * delta_t / tau_;
			p_diffuse_from_tether_[x_dist_dub] = p_from;
			p_diffuse_to_tether_[x_dist_dub] = p_to;
		}
	}
	// Statistics for KMC
	p_bind_i_free_ = k_on*c_motor*delta_t;
	p_bind_i_tethered_ = p_bind_i_free_ * 4;
	p_bind_ii_ = 0.05; 				// FIXME should be like ~ tau/delta_t
	p_unbind_untethered_ = k_off * delta_t;
	p_unbind_tethered_ = p_unbind_untethered_ / 2; 
	p_unbind_pseudo_ = 0.005;		// FIXME
	p_tether_free_ = p_bind_i_free_;  	
	p_untether_free_ = p_unbind_untethered_ / 2;
	p_step_untethered_ = motor_speed*delta_t;
	p_switch_ = switch_rate*delta_t;
	// Generate untethering and stepping rates for all tether extensions	
	// Everything is 2*dist_cutoff to allow for half-integer distances, 
	// so the 3rd entry will correspond to a distance of 1.5, etc. 
	p_untether_bound_.resize(2*dist_cutoff_ + 1);
	p_step_tethered_.resize(2*dist_cutoff_ + 1);
	double stall_force = motor_list_[0].stall_force_;		// in pN
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		double r_x = x_dist_dub*site_size/2;
		double r = sqrt(r_y*r_y + r_x*r_x);
		double dr = r - r_0;
		// If extension is positive, treat tether like a spring
		if(dr >= 0){
			// Calculate spring potential energy (PE) for this x_dist
			double U_teth = (k_spring/2)*dr*dr;
			// Get untethering weight for this PE
			double unteth_weight = exp(U_teth/(2*kbT));
			p_untether_bound_[x_dist_dub] = p_untether_free_*unteth_weight;
			// Stepping probability has force-dependent relation
			double force = dr*k_spring;
			if(force < stall_force){
				// see sci advancements supplement
				double cosine = r_x/r;
				double coeff = 1 - cosine*(force/stall_force);
				p_step_tethered_[x_dist_dub] = coeff*p_step_untethered_;
			}
			else{
				p_step_tethered_[x_dist_dub] = 0;
			}
		}
		// Otherwise, use k_eff for slack
		else{
			double U_teth = (k_eff_slack/2)*dr*dr;
			double unteth_weight = exp(U_teth/(2*kbT));
			p_untether_bound_[x_dist_dub] = p_untether_free_*unteth_weight;
			double force = dr*k_eff_slack;
			if(force < stall_force){
				double cosine = r_x/r;
				double coeff = 1 - cosine*(force/stall_force);
				p_step_tethered_[x_dist_dub] = coeff*p_step_untethered_;
			}
			else{
				p_step_tethered_[x_dist_dub] = 0;
			}
		}
	}
	alpha_ = parameters_->alpha;
	beta_ = parameters_->beta;
}
void KinesinManagement::InitiateLists(){

	// One dimensional stuff
	free_tethered_list_.resize(n_motors_);
	pseudo_bound_list_.resize(n_motors_);
	eligible_pseudo_list_.resize(n_motors_); 
	bound_list_.resize(n_motors_);
	bound_untethered_list_.resize(n_motors_);
	bound_tethered_list_.resize(n_motors_);
	switchable_list_.resize(n_motors_);
	stepable_untethered_list_.resize(n_motors_);
	// Two dimensional stuff
	n_bound_tethered_.resize(2*dist_cutoff_ + 1);
	n_stepable_tethered_.resize(2*dist_cutoff_ + 1);
	bound_tethered_table_.resize(2*dist_cutoff_ + 1);
	stepable_tethered_table_.resize(2*dist_cutoff_ + 1);
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		bound_tethered_table_[x_dist_dub].resize(n_motors_);
		stepable_tethered_table_[x_dist_dub].resize(n_motors_);
		n_bound_tethered_[x_dist_dub] = 0;
		n_stepable_tethered_[x_dist_dub] = 0;
	}
}

void KinesinManagement::UnboundCheck(Kinesin *motor){

    // Check motor
    if(motor->front_site_ != nullptr || motor->rear_site_ != nullptr){
        printf("Error with motor #%i: classified as unbound, ", motor->ID_);
        printf("but at least one of its heads is attached to tubulin\n");
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
        printf("Error with motor #%i: classified as bound, " , motor->ID_);
        printf("but at least one of its heads is not attached to tubulin\n");
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

    if(motor->front_site_->index_ == motor->mt_->plus_end_
	|| motor->rear_site_->index_ == motor->mt_->minus_end_){
        return true;
    }
    else{
        return false;
    }
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
	if(n_entries != n_pseudo_bound_){
		printf("something wrong in update_pseudo_bound (motors)\n");
		exit(1);
	}
}

void KinesinManagement::UpdateEligiblePseudoList(){

	n_eligible_pseudo_ = 0; 
	int i_entry = 0; 
	int mt_length = parameters_->length_of_microtubule; 
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motor_list_[i_motor]; 
		if(motor->heads_active_ == 1){
			Microtubule *mt = motor->mt_;
			Tubulin *bound_site = motor->GetActiveHeadSite();
			int i_site = bound_site->index_;
			int dx = mt->delta_x_;
			// Don't access lattice sites that don't exist (or boundaries)
			if(!(i_site <= 0 && dx == -1)
			&& !(i_site >= mt_length - 1 && dx == 1)){
				// Add motor if it has an unoccupied site to either side
				if(mt->lattice_[i_site + dx].occupied_ == false
				|| mt->lattice_[i_site - dx].occupied_ == false){
					eligible_pseudo_list_[i_entry] = motor;
					i_entry++;
					n_eligible_pseudo_++; 
				}
			}
		}
	}
}

/*
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
	if(n_entries != n_bound_){
		printf("something wrong in update_bound_list (motors):");
		printf(" %i vs %i\n", n_entries, n_bound_);
		exit(1);
	}
}
*/

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
		printf("something wrong in update_bound_untethered:");
		printf(" %i entries in list, %i from statistics\n", 
				n_entries, n_bound_untethered_);
		exit(1);
	}
}

void KinesinManagement::UpdateBoundTetheredList(){

	int i_entry = 0;
	int n_entries = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motor_list_[i_motor];
		if(motor->heads_active_ == 2
		&& motor->tethered_ == true){
			bound_tethered_list_[i_entry] = motor;
			i_entry++;
			n_entries++;
		}
	}
	if(n_entries != n_bound_tethered_tot_){
		printf("something wrong in update_teth_list (motors 1D)");
		printf(" %in statistics, %i entries tho\n", n_bound_tethered_tot_, 
				n_entries);
		exit(1);
	}
}

void KinesinManagement::UpdateBoundTetheredTable(){

	int i_entry[2*dist_cutoff_ + 1];
	int n_entries[2*dist_cutoff_ + 1]; 
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		i_entry[x_dist_dub] = 0;
		n_entries[x_dist_dub] = 0;
	}
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motor_list_[i_motor];
		if(motor->heads_active_ == 2
		&& motor->tethered_ == true){
			motor->UpdateExtension();
			if(motor->tethered_ == true){
				int x_dist_dub = motor->x_dist_doubled_; 
				int index = i_entry[x_dist_dub]; 
				bound_tethered_table_[x_dist_dub][index] = motor;
				i_entry[x_dist_dub]++;
				n_entries[x_dist_dub]++;
			}
		}
	}
	for(int x_dist_dub = 0; x_dist_dub <= dist_cutoff_; x_dist_dub++){
		if(n_entries[x_dist_dub] != n_bound_tethered_[x_dist_dub]){
			printf("something wrong in update_bound_tethered (motors)");
			printf("for ext %i, %i in stats but %i entries counted\n", 
					x_dist_dub, n_bound_tethered_[x_dist_dub], 
					n_entries[x_dist_dub]);
			exit(1);
		}
	}
}

void KinesinManagement::UpdateSwitchableList(){

	n_switchable_ = 0; 
	int i_entry = 0;
	int mt_length = parameters_->length_of_microtubule;
    for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motor_list_[i_motor];
		if(motor->heads_active_ == 2
		&& motor->tethered_ == true
		&& BoundaryStatus(motor) == false){
            Microtubule *mt = motor->mt_;
            Microtubule *mt_adj = mt->neighbor_; 	// FIXME define neighb
            int mt_coord = mt->coord_;
            int mt_adj_coord = mt_adj->coord_;
            int offset = mt_adj_coord - mt_coord;
            int i_front = motor->front_site_->index_;
            int i_rear = motor->rear_site_->index_;
            int i_front_adj = i_rear - offset;
            int i_rear_adj = i_front - offset;
            // Exclude non-overlapping microtubule segments
            if((i_front_adj > 0 && i_front_adj < mt_length - 1)
			&& (i_rear_adj > 0 && i_rear_adj < mt_length - 1)){
                if(mt_adj->lattice_[i_front_adj].occupied_ == false
				&& mt_adj->lattice_[i_rear_adj].occupied_ == false){
					switchable_list_[i_entry] = motor; 
					i_entry++;
                    n_switchable_++;
				}
            }
        }
    }
}

void KinesinManagement::UpdateStepableTetheredTable(){

	int i_entry[2*dist_cutoff_ + 1];
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		i_entry[x_dist_dub] = 0;
		n_stepable_tethered_[x_dist_dub] = 0; 	
	}
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motor_list_[i_motor];
		if(motor->heads_active_ == 2
		&& motor->tethered_ == true){
			motor->UpdateExtension();
			Microtubule *mt = motor->mt_;
			int plus_end = mt->plus_end_;
			int delta_x = mt->delta_x_;
			int i_front = motor->front_site_->index_;
			// end pausing
			if(i_front != plus_end){
				if(mt->lattice_[i_front + delta_x].occupied_ == false){
					int x_dist_dub = motor->x_dist_doubled_;
					int index = i_entry[x_dist_dub];
//					printf("added ext %i to list index %i\n", 
//							motor->x_dist_doubled_, x_dist_dub);
					stepable_tethered_table_[x_dist_dub][index] = motor;
					i_entry[x_dist_dub]++;
					n_stepable_tethered_[x_dist_dub]++;
				}
			}
		}
	}
}

void KinesinManagement::UpdateStepableUntetheredList(){

	n_stepable_untethered_ = 0; 
	int i_entry = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motor_list_[i_motor];
        if(motor->heads_active_ == 2
		&& motor->tethered_ == false){
            Microtubule *mt = motor->mt_;
            int plus_end = mt->plus_end_;
            int delta_x = mt->delta_x_;
            int i_front_site = motor->front_site_->index_;
            // Exclude plus_end from statistics because of end-pausing
            if(i_front_site != plus_end){
                if(mt->lattice_[i_front_site + delta_x].occupied_ == false){
					stepable_untethered_list_[i_entry] = motor;
					i_entry++; 
					n_stepable_untethered_++;
                }
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

void KinesinManagement::GenerateDiffusionList(){

	int n_events = 0; 
	// Untethered statistics
	UpdateBoundUntetheredList();
	int n_fwd_unteth = GetNumToStepForward_Unteth();
	int n_bck_unteth = GetNumToStepBackward_Unteth();
	while(n_fwd_unteth + n_bck_unteth > n_bound_untethered_){
		double ran = properties_->gsl.GetRanProb();
		if(ran < 0.5)
			n_fwd_unteth--;
		else
			n_bck_unteth--;
	}
	n_events += n_fwd_unteth;
	n_events += n_bck_unteth;
	// Tethered statistics
	UpdateBoundTetheredTable();
	int n_toward_rest[2*dist_cutoff_ + 1];
	int n_from_rest[2*dist_cutoff_ + 1];
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		n_toward_rest[x_dist_dub] = GetNumToStepTowardRest(x_dist_dub);
		n_from_rest[x_dist_dub] = GetNumToStepFromRest(x_dist_dub);
		int n_to = n_toward_rest[x_dist_dub];
		int n_from = n_from_rest[x_dist_dub];
		int n_tethered = n_bound_tethered_[x_dist_dub];
		while(n_to + n_from > n_tethered){
			double ran = properties_->gsl.GetRanProb();
			double p_to = p_diffuse_to_tether_[x_dist_dub];
			double p_from = p_diffuse_from_tether_[x_dist_dub];
			double tot_prob = p_to + p_from;  // XXX is this correct? XXX
			if(ran < p_to/tot_prob){
				n_toward_rest[x_dist_dub]--;
				n_to = n_toward_rest[x_dist_dub];
			}
			else{
				n_from_rest[x_dist_dub]--;
				n_to = n_from_rest[x_dist_dub];
			}
		}
		n_events += n_toward_rest[x_dist_dub];
		n_events += n_from_rest[x_dist_dub];
	}
	// Only generate list if we have more than zero diffusion events
	if(n_events > 0){
		int pre_list[n_events];
		int diff_index = 0;
		for(int i_event = 0; i_event < n_fwd_unteth; i_event++){
			pre_list[diff_index] = 10;
			diff_index++;
		}
		for(int i_event = 0; i_event < n_bck_unteth; i_event++){
			pre_list[diff_index] = 20;
			diff_index++;
		}
		for(int x_dist_dub = 0; x_dist_dub < 2*dist_cutoff_; x_dist_dub++){
			int n_step_to = n_toward_rest[x_dist_dub];
			for(int i_event = 0; i_event < n_step_to; i_event++){
				pre_list[diff_index] = 300 + x_dist_dub;
				diff_index++;
			}
			int n_step_from = n_from_rest[x_dist_dub];
			for(int i_event = 0; i_event < n_step_from; i_event++){
				pre_list[diff_index] = 400 + x_dist_dub;
				diff_index++;
			}
		}
		RandomNumberManagement *gsl = &properties_->gsl;
		gsl_ran_shuffle(gsl->rng, pre_list, n_events, sizeof(int));
		diffusion_list_.resize(n_events);
		for(int i_event = 0; i_event < n_events; i_event++){
			diffusion_list_[i_event] = pre_list[i_event];
		}
	}
	else
		diffusion_list_.clear();
}

int KinesinManagement::GetNumToStepForward_Unteth(){

	int n_bound = n_bound_untethered_;
	double p_step = p_diffuse_fwd_untethered_;
	double n_step_avg = n_bound * p_step;
	int n_to_step = properties_->gsl.SamplePoissonDist(n_step_avg);
	return n_to_step;
}

int KinesinManagement::GetNumToStepBackward_Unteth(){

	int n_bound = n_bound_untethered_;
	double p_step = p_diffuse_bck_untethered_;
	double n_step_avg = n_bound * p_step;
	int n_to_step = properties_->gsl.SamplePoissonDist(n_step_avg);
	return n_to_step;
}

int KinesinManagement::GetNumToStepTowardRest(int x_dist_doubled){

	int n_bound = n_bound_tethered_[x_dist_doubled];
	double p_step = p_diffuse_to_tether_[x_dist_doubled];
	double n_step_avg = n_bound * p_step;
	int n_to_step = properties_->gsl.SamplePoissonDist(n_step_avg);
	return n_to_step;
}

int KinesinManagement::GetNumToStepFromRest(int x_dist_doubled){

	int n_bound = n_bound_tethered_[x_dist_doubled];
	double p_step = p_diffuse_from_tether_[x_dist_doubled];
	double n_step_avg = n_bound * p_step;
	int n_to_step = properties_->gsl.SamplePoissonDist(n_step_avg);
	return n_to_step;
}

void KinesinManagement::RunDiffusion(){

	GenerateDiffusionList();
	if(diffusion_list_.empty() == false){
		int n_events = diffusion_list_.size();
		int x_dist_dub = 0;
		for(int i_event = 0; i_event < n_events; i_event++){
			int diff_event = diffusion_list_[i_event];
			if(diff_event >= 300 && diff_event < 400){
				x_dist_dub = diff_event	% 100;
				diff_event = 30;
			}
			if(diff_event >= 400 && diff_event < 500){
				x_dist_dub = diff_event % 100;
				diff_event = 40;
			}
			switch(diff_event){
				case 10:
					printf("unteth step fwd\n");
					RunDiffusion_Forward_Untethered();
					break;	
				case 20:
					printf("unteth step bck\n");
					RunDiffusion_Backward_Untethered();
					break;
				case 30:
					printf("teth step to\n");
					RunDiffusion_Toward_Tether(x_dist_dub);
					break;
				case 40:
					printf("teth step from\n");
					RunDiffusion_From_Tether(x_dist_dub);
					break;
			}
		}
	}
}

void KinesinManagement::RunDiffusion_Forward_Untethered(){

	UpdateBoundUntetheredList();
	int n_bound = n_bound_untethered_;
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_untethered_list_[i_entry];
		Microtubule *mt = motor->mt_;
		int i_front = motor->front_site_->index_; 
		int dx = mt->delta_x_;
		int plus_end = mt->plus_end_;
		if(i_front != plus_end){
			if(mt->lattice_[i_front + dx].occupied_ == false){
				// Get new sites
				Tubulin *new_front = &mt->lattice_[i_front + dx];
				Tubulin *new_rear = motor->front_site_;
				Tubulin *old_rear = motor->rear_site_;
				// Update new site
				new_front->motor_ = motor;
				new_front->occupied_ = true;
				// Update old site
				old_rear->motor_ = nullptr;
				old_rear->occupied_ = false;
				// Update motor
				motor->front_site_ = new_front;
				motor->rear_site_ = new_rear;

			}
			else{
				printf("oh well fwd\n");
			}
		}
		else{
			printf("cant diffuse outta this one brotha\n");
		}
	}
	else{
		printf("ya blew it. we failed to step untethered motor fwd\n");
		exit(1);
	}

}

void KinesinManagement::RunDiffusion_Backward_Untethered(){

	UpdateBoundUntetheredList();
	int n_bound = n_bound_untethered_;
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_untethered_list_[i_entry];
		Microtubule *mt = motor->mt_;
		int i_rear = motor->rear_site_->index_; 
		int dx = mt->delta_x_;
		int minus_end = mt->minus_end_;
		if(i_rear != minus_end){
			if(mt->lattice_[i_rear - dx].occupied_ == false){
				// Get new sites
				Tubulin *new_rear = &mt->lattice_[i_rear - dx];
				Tubulin *new_front = motor->rear_site_;
				Tubulin *old_front = motor->front_site_;
				// Update new site
				new_rear->motor_ = motor;
				new_rear->occupied_ = true;
				// Update old site
				old_front->motor_ = nullptr;
				old_front->occupied_ = false;
				// Update motor
				motor->front_site_ = new_front;
				motor->rear_site_ = new_rear;

			}
			else{
				printf("oh well bck\n");
			}
		}
		else{
			printf("cant diffuse outta this one brotha\n");
		}

	}
	else{
		printf("ya blew it. we failed to step untethered motor bck\n");
		exit(1);
	}
}

void KinesinManagement::RunDiffusion_Toward_Tether(int x_dist_doubled){
	
	UpdateBoundTetheredTable();
	int n_bound = n_bound_tethered_[x_dist_doubled];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_tethered_table_[x_dist_doubled][i_entry];
		Microtubule *mt = motor->mt_;
		Tubulin *near_site = motor->GetSiteCloserToXlink();
		Tubulin *far_site = motor->GetSiteFartherFromXlink();


	}
	else{
		printf("ya blew it. we failed to step tethered motor towards\n");
		exit(1);
	}
}

void KinesinManagement::RunDiffusion_From_Tether(int x_dist_doubled){

}

void KinesinManagement::GenerateKMCList(){

	/* see xlinks analog for a much simpler case and explanation */
	int n_events = 0;
    int n_bind_i_free = GetNumToBind_I_Free();
	n_events += n_bind_i_free;
	UpdateFreeTetheredList();
    int n_bind_i_tethered = GetNumToBind_I_Tethered();
	n_events += n_bind_i_tethered;
	UpdateEligiblePseudoList();
	int n_bind_ii = GetNumToBind_II();
	n_events += n_bind_ii;
	int n_unbind_ut = GetNumToUnbind_Untethered();
	n_events += n_unbind_ut;
	int n_unbind_t = GetNumToUnbind_Tethered();
	n_events += n_unbind_t;
	UpdatePseudoBoundList();
	int n_unbind_p = GetNumToUnbind_Pseudo();
	while(n_unbind_p + n_bind_ii > n_pseudo_bound_){
		if(n_bind_ii > 0){
			n_bind_ii--;
			n_events--;
		}
		else{
			n_unbind_p--;
		}
	}
	n_events += n_unbind_p;
	int n_tether_free = GetNumToTether_Free();
	n_events += n_tether_free;
	UpdateBoundUntetheredList();
	int n_tether_bound = GetNumToTether_Bound();
	n_events += n_tether_bound;
	UpdateSwitchableList();
    int n_switch = GetNumToSwitch();	
	n_events += n_switch;
	int n_untether_free = GetNumToUntether_Free();
	// Ensure we don't get too many binding events from free motors (if we
	// have 1 and roll for 1 untether, then obviously we can't bind 1 too)
	while(n_untether_free + n_bind_i_tethered > n_free_tethered_){
			if(n_bind_i_tethered > 0){
				n_bind_i_tethered--;
				n_events--;	
			}
			else{
				n_untether_free--;
			}
	}
	n_events += n_untether_free;
	UpdateStepableUntetheredList();
	int n_step_ut = GetNumToStep_Untethered(); 
	// Ensure we don't get too many events for untethered bound motors
	while(n_unbind_ut + n_tether_bound + n_step_ut > n_bound_untethered_){
		// If the # of step events is 0 or less, remove tether events
		// (unbinding tends to have priority over all over events)
		if(n_step_ut > 0){
			n_step_ut--;
		}
		else if(n_tether_bound > 0){
			n_tether_bound--;
			n_events--;
		}
		else{
			n_unbind_ut--;
			n_events--;
		}
	}
	n_events += n_step_ut;
	// Handle the statistics of differnt tether extensions separately 
	UpdateStepableTetheredTable();
	UpdateBoundTetheredTable();
	int n_untether_bound[2*dist_cutoff_ + 1];
	int n_step_tethered[2*dist_cutoff_ + 1];
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		n_untether_bound[x_dist_dub] = GetNumToUntether_Bound(x_dist_dub);
		n_step_tethered[x_dist_dub] = GetNumToStep_Tethered(x_dist_dub);
		int n_unteth = n_untether_bound[x_dist_dub];
		int n_step_teth = n_step_tethered[x_dist_dub];
		int n_tethered = n_bound_tethered_[x_dist_dub];
		// Prevent too many steps from occuring (imagine if 2 motors
		// are bound and we roll for 1 unbind but 2 steps)
		while(n_unteth + n_step_teth > n_tethered){
			printf("*** STAT CORRECTION ***\n");
			printf("to unteth: %i, to step: %i, avail: %i\n", 
					n_unteth, n_step_teth, n_tethered);
			if(n_unteth > 10000)
				exit(1);
			if(n_step_teth > 0){
				n_step_tethered[x_dist_dub]--;
				n_step_teth = n_step_tethered[x_dist_dub];
			}
			else if (n_unteth > 0){
				n_untether_bound[x_dist_dub]--;
				n_unteth = n_untether_bound[x_dist_dub];
			}
		} 
		n_events += n_untether_bound[x_dist_dub];
		n_events += n_step_tethered[x_dist_dub];
	}
//	printf("n_events: %i \n", n_events);
    if(n_events > 0){
        int pre_list[n_events];
		int kmc_index = 0;
		for(int i_event = 0; i_event < n_bind_i_free; i_event++){
			printf("10\n");
			pre_list[kmc_index] = 10;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_bind_i_tethered; i_event++){
			printf("11\n");
			pre_list[kmc_index] = 11;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_bind_ii; i_event++){
			printf("12\n");
			pre_list[kmc_index] = 12;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_unbind_ut; i_event++){
			printf("20\n");
			pre_list[kmc_index] = 20;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_unbind_t; i_event++){
			printf("21\n");
			pre_list[kmc_index] = 21;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_unbind_p; i_event++){
			printf("22\n");
			pre_list[kmc_index] = 22;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_tether_free; i_event++){
			printf("30\n");
			pre_list[kmc_index] = 30;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_tether_bound; i_event++){
			printf("31\n");
			pre_list[kmc_index] = 31;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_switch; i_event++){
			printf("40\n");
			pre_list[kmc_index] = 40;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_untether_free; i_event++){
			printf("50\n");
			pre_list[kmc_index] = 50;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_step_ut; i_event++){
			printf("60\n");
			pre_list[kmc_index] = 60;
			kmc_index++;
		}
		for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
			int n_unteth_b = n_untether_bound[x_dist_dub];
//			printf("n_unteth_b: %i (ext %i)\n", n_unteth_b, x_dist_dub);
			for(int i_event = 0; i_event < n_unteth_b; i_event++){
//				printf("n_events: %i, index: %i, 2x: %i\n", 
//						n_events, kmc_index, x_dist_dub);
//				fflush(stdout);
				pre_list[kmc_index] = 500 + x_dist_dub; 	
				printf("%i\n", 500 + x_dist_dub);
				kmc_index++;
			}
			int n_step_t = n_step_tethered[x_dist_dub];
			for(int i_event = 0; i_event < n_step_t; i_event++){
				pre_list[kmc_index] = 600 + x_dist_dub;
				printf("%i\n", 600 + x_dist_dub);
				kmc_index++;
			}
		}
		RandomNumberManagement *gsl = &properties_->gsl;
        gsl_ran_shuffle(gsl->rng, pre_list, n_events, sizeof(int));
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
    double n_avg = p_bind_i_free_*n_unoccupied_pairs;
    int n_to_bind = properties_->gsl.SamplePoissonDist(n_avg);
    return n_to_bind;
}

int KinesinManagement::GetNumToBind_I_Tethered(){

    double weights_summed = 0;
    // Sum over all tethered but unbound motors
    for(int i_motor = 0; i_motor < n_free_tethered_; i_motor++){
        Kinesin *motor = free_tethered_list_[i_motor];
		motor->UpdateNeighborSites();
        // Get weight of all neighbor sites
        int n_neighbs = motor->n_neighbor_sites_;
        for(int i_neighb = 0; i_neighb < n_neighbs; i_neighb++){
            Tubulin *site = motor->neighbor_sites_[i_neighb];
            double weight = motor->GetBindingWeight(site);
            weights_summed += weight;
        }
    }
    return weights_summed;
}

int KinesinManagement::GetNumToBind_II(){

	double p_bind = p_bind_ii_; 
	int n_able = n_eligible_pseudo_;
	int n_to_bind = properties_->gsl.SampleBinomialDist(p_bind, n_able);
	return n_to_bind;
}

int KinesinManagement::GetNumToUnbind_Untethered(){

	int n_bound = n_bound_untethered_;
	double p_unbind = p_unbind_untethered_;
    int n_to_unbind = properties_->gsl.SampleBinomialDist(p_unbind, n_bound);
    return n_to_unbind;
}

int KinesinManagement::GetNumToUnbind_Tethered(){
		
	int n_bound = n_bound_tethered_tot_;
	double p_unbind = p_unbind_tethered_;
	int n_to_unbind = properties_->gsl.SampleBinomialDist(p_unbind, n_bound);
	return n_to_unbind; 
}

int KinesinManagement::GetNumToUnbind_Pseudo(){

	int n_bound = n_pseudo_bound_; 
	double p_unbind = p_unbind_pseudo_;
	int n_to_unbind = properties_->gsl.SampleBinomialDist(p_unbind, n_bound);
	return n_to_unbind;
}

int KinesinManagement::GetNumToTether_Free(){

    int n_untethered_xlinks = properties_->prc1.n_untethered_;
	// Calculate how many free motors tether within delta_t on avg
    double n_tether_avg = p_tether_free_*n_untethered_xlinks; 	
    int n_to_tether = properties_->gsl.SamplePoissonDist(n_tether_avg);
    return n_to_tether;
}

int KinesinManagement::GetNumToTether_Bound(){

	double weights_summed = 0;
	for(int i_motor = 0; i_motor < n_bound_untethered_; i_motor++){
		Kinesin *motor = bound_untethered_list_[i_motor];
		motor->UpdateNeighborXlinks();
		int n_neighbs = motor->n_neighbor_xlinks_; 
		for(int i_neighb = 0; i_neighb < n_neighbs; i_neighb++){
			AssociatedProtein *xlink = motor->neighbor_xlinks_[i_neighb];
			double weight = motor->GetTetheringWeight(xlink);
			weights_summed += weight;
		}
	}
	return weights_summed;
}

int KinesinManagement::GetNumToSwitch(){

	double p_switch = p_switch_;
	int n_able = n_switchable_;
    int n_to_switch = properties_->gsl.SampleBinomialDist(p_switch, n_able);
    return n_to_switch;
}

int KinesinManagement::GetNumToUntether_Bound(int x_dist_doubled){

	double p_unteth = p_untether_bound_[x_dist_doubled];
	int n_teth = n_bound_tethered_[x_dist_doubled]; 
	int n_to_unteth = properties_->gsl.SampleBinomialDist(p_unteth, n_teth);
	if(n_teth < 0){
		printf("p_unteth for dist %i: %g \n", x_dist_doubled, p_unteth);
		printf("n_teth for dist %i: %i \n", x_dist_doubled, n_teth);
		printf("n to unteth for dist %i: %i\n",x_dist_doubled, n_to_unteth);
	}
	return n_to_unteth;
}


int KinesinManagement::GetNumToUntether_Free(){

	double p_unteth = p_untether_free_;
	int n_teth = n_free_tethered_;
	int n_to_unteth = properties_->gsl.SampleBinomialDist(p_unteth, n_teth);
	return n_to_unteth; 
}

int KinesinManagement::GetNumToStep_Tethered(int x_dist_doubled){

	double p_step = p_step_tethered_[x_dist_doubled];
	int n_stepable = n_stepable_tethered_[x_dist_doubled];
	int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_stepable);
	return n_to_step;
}

int KinesinManagement::GetNumToStep_Untethered(){

	double p_step = p_step_untethered_;
	int n_stepable = n_stepable_untethered_;
    int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_stepable);
    return n_to_step;
}

void KinesinManagement::RunKMC(){

//	printf(" *** START OF KMC CYCLE ***\n");
    GenerateKMCList();
    if(kmc_list_.empty() == false){
        int n_events = kmc_list_.size();
		int x_dist_doubled = 0;
        for(int i_event = 0; i_event < n_events; i_event++){
            int kmc_event = kmc_list_[i_event];
			if(kmc_event >= 500 && kmc_event < 600){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 51;
			}
			if(kmc_event >= 600 && kmc_event < 700){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 61; 
			}	
			properties_->wallace.PrintMicrotubules(0.05);
            switch(kmc_event){
				case 10:
						printf("free motor pseudo-bound\n");
						fflush(stdout);
						KMC_Bind_I_Free();
                        break;
                case 11:
						printf("tethered motor pseudo-bound\n");
						fflush(stdout);
						KMC_Bind_I_Tethered();
                        break;
				case 12:
						printf("pseudo motor bound\n");
						fflush(stdout);
						KMC_Bind_II();
						break;
				case 20:
						printf("untethered motor unbound\n");
						fflush(stdout);
						KMC_Unbind_Untethered();
                        break;
                case 21:
						printf("tethered motor unbound\n");
						fflush(stdout);
						KMC_Unbind_Tethered();
                        break;
				case 22:
						printf("pseudo-bound motor unbound\n");
						fflush(stdout);
						KMC_Unbind_Pseudo();
						break;
                case 30:
						printf("free motor tethered\n");
						fflush(stdout);
						KMC_Tether_Free();
                        break;
				case 31:
						printf("bound motor tethered\n");
						fflush(stdout);
						KMC_Tether_Bound();
						break;
				case 40:
						printf("motor switched\n");
						fflush(stdout);
						KMC_Switch();
						properties_->wallace.PrintMicrotubules(1.5);
                        break;
				case 50:
						printf("free motor untethered\n");
						fflush(stdout);
						KMC_Untether_Free();
						break;
				case 51:
						printf("bound motor (ext %i) untethered\n", 
								x_dist_doubled);
						fflush(stdout);
						KMC_Untether_Bound(x_dist_doubled);
						break;
				case 60:
						printf("untethered motor stepped\n");
						fflush(stdout);
						KMC_Step_Untethered();
						break;
				case 61:
						printf("tethered motor (ext %i) stepped\n", 
								x_dist_doubled);
						fflush(stdout);
						KMC_Step_Tethered(x_dist_doubled);
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
	properties_->microtubules.UpdateUnoccupiedList();
	if(properties_->microtubules.n_unoccupied_ > 0){	
        Kinesin *motor = GetFreeMotor();
		Tubulin *site = properties_->microtubules.GetUnoccupiedSite();
		Microtubule *mt = site->mt_;
		// Update site details
		site->motor_ = motor;
		site->occupied_ = true;
		// Update motor details
		motor->mt_ = mt;
		motor->front_site_ = site;
		motor->heads_active_ = 1;
		// Update statistics
		n_pseudo_bound_++;
	}
	else{
		printf("Error in Bind_I_Free: no unoccupied sites.\n");
		exit(1);
	}
}

void KinesinManagement::KMC_Bind_I_Tethered(){
	
	// Make sure at least one unoccupied site exists
	properties_->microtubules.UpdateUnoccupiedList();
	UpdateFreeTetheredList();
	if(properties_->microtubules.n_unoccupied_ > 0
	&& n_free_tethered_ > 0){
		// Pick a random tethered free motor to bind
		int i_motor = properties_->gsl.GetRanInt(n_free_tethered_);
		Kinesin *motor = free_tethered_list_[i_motor];
		motor->UpdateNeighborSites();
		int x_dist_dub = motor->SampleTailExtensionDoubled();
		// Ensure motor has a neighbor site for desired tether extension
		// Roll for a new extension after a certain # of tries
		int attempts = 0;
		int switches = 0;
		bool failure = false;
		while(motor->NeighborSiteExists(x_dist_dub) == false){
			if(attempts > 2*n_free_tethered_){
				if(switches > 10){
					printf("failed to bind_i_tethered\n");
					failure = true;
					break;
				}
				x_dist_dub = motor->SampleTailExtensionDoubled();
				attempts = 0;
				switches++;
			}
			i_motor = properties_->gsl.GetRanInt(n_free_tethered_);
			motor = free_tethered_list_[i_motor];
			motor->UpdateNeighborSites();
			attempts++;
		}
		if(failure == false){
			Tubulin *site = motor->GetNeighborSite(x_dist_dub); 
			Microtubule *mt = site->mt_; 
			// Update site details
			site->motor_ = motor;
			site->occupied_ = true;
			// Update motor detail
			motor->mt_ = mt;
			motor->front_site_ = site;
			motor->heads_active_ = 1;
			// Update statistics;
			n_free_tethered_--;
			n_pseudo_bound_++; 
		}
	}
	else{
		printf("Error in Bind_I_Tethered: no unoccupied sites.\n");
		exit(1);
	}
}

void KinesinManagement::KMC_Bind_II(){

	UpdateEligiblePseudoList();
	if(n_eligible_pseudo_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_eligible_pseudo_);
		Kinesin *motor = eligible_pseudo_list_[i_entry];
		Microtubule *mt = motor->mt_;
		Tubulin *bound_site = motor->GetActiveHeadSite();
		int i_site = bound_site->index_;
		int dx = mt->delta_x_;
		Tubulin *front_site,
				*rear_site; 
		if(mt->lattice_[i_site + dx].occupied_ == false){
			front_site = &mt->lattice_[i_site + dx];
			rear_site = bound_site;
		}
		else if(mt->lattice_[i_site - dx].occupied_ == false){
			front_site = bound_site;
			rear_site = &mt->lattice_[i_site - dx];
		}
		else{
			printf("Error with eligible pseudo bound motors list \n");
			exit(1); 
		}
		// Update site details
		front_site->motor_ = motor;
		front_site->occupied_ = true;
		rear_site->motor_ = motor;
		rear_site->occupied_ = true;
		// Update motor details
		motor->front_site_ = front_site;
		motor->rear_site_ = rear_site;
		motor->heads_active_ = 2;	
		// Update statistics
		n_pseudo_bound_--;
		n_eligible_pseudo_--;
		n_bound_++;
		if(motor->tethered_ == true){
			int x_dub_pre = motor->x_dist_doubled_;
			motor->UpdateExtension();
			// weird exception if we triggered a force untether
			if(motor->tethered_ == false){
				// These cancel out the subtraction in the force untether
				// since this particular xlink never contributed to stats, 
				// the subtraction just fucks with stuff
				n_bound_tethered_tot_++;
				n_bound_tethered_[x_dub_pre]++;
			}
			else{
				int x_dist_doubled = motor->x_dist_doubled_; 
				n_bound_tethered_[x_dist_doubled]++;
				n_bound_tethered_tot_++;
			}
		}
		else{
			n_bound_untethered_++;
		}
	}
	else{
		printf("Error in Bind_II: no eligible pseudo. \n");
		exit(1);
	}
}	

void KinesinManagement::KMC_Unbind_Untethered(){

    // Make sure that at least one bound motor exists
	UpdateBoundUntetheredList();
    if(n_bound_untethered_ > 0){
        // Randomly pick a bound motor (that's not on the boundary)
		int i_entry = 0;
		Kinesin *motor;
		do{
			i_entry = properties_->gsl.GetRanInt(n_bound_untethered_);
			motor = bound_untethered_list_[i_entry];
		}while(BoundaryStatus(motor) == true);
		// Update site details
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
		n_bound_untethered_--; 
		n_bound_--;
	}
	else{
		printf("Error in Unbind_Untethered: no bound untethered motors!\n");
        exit(1);
    }
}

void KinesinManagement::KMC_Unbind_Tethered(){

	UpdateBoundTetheredList();
	if(n_bound_tethered_tot_ > 0){
		int i_entry = 0;
		Kinesin *motor;
		do{	
			i_entry = properties_->gsl.GetRanInt(n_bound_tethered_tot_);
			motor = bound_tethered_list_[i_entry];
		}while(BoundaryStatus(motor) == true);
		// Update site details
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
		n_free_tethered_++;
		int x_dub_pre = motor->x_dist_doubled_;
		printf("motor that was unbound had ext of %i\n", x_dub_pre);
		motor->UpdateExtension();
		n_bound_tethered_[x_dub_pre]--;
		n_bound_tethered_tot_--;
		n_bound_--;
	}
	else{
		printf("Error in Unbind_Tethered: no bound tethered motors!\n");
		exit(1);
	}
}

void KinesinManagement::KMC_Unbind_Pseudo(){

	UpdatePseudoBoundList();
	if(n_pseudo_bound_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_pseudo_bound_);
		Kinesin *motor = pseudo_bound_list_[i_entry];
		// Update site details
		Tubulin *site = motor->GetActiveHeadSite();
		site->motor_ = nullptr;
		site->occupied_ = false;
		// Update motor details
		motor->heads_active_ = 0;
		motor->mt_ = nullptr;
		motor->front_site_ = nullptr;
		motor->rear_site_ = nullptr;
		// Update statistics
		if(motor->tethered_ == true){
			n_free_tethered_++; 
		}
		n_pseudo_bound_--;
	}
	else{
		printf("Error in Unbind:Pseudo: no pseudo bound motors!\n");
		exit(1);
	}
}

void KinesinManagement::KMC_Tether_Free(){

	// Make sure there is at least one unbound motor
	Kinesin *motor = GetFreeMotor();
	if(properties_->prc1.n_untethered_ > 0){
		AssociatedProtein *xlink = properties_->prc1.GetUntetheredXlink();
		// Update motor and xlink details
		motor->xlink_ = xlink;
		motor->tethered_ = true;
		xlink->tethered_ = true;
		xlink->motor_ = motor; 	
		// Update statistics 
		n_free_tethered_++;
		properties_->prc1.n_tethered_++;
		properties_->prc1.n_untethered_--;
	}
}

void KinesinManagement::KMC_Tether_Bound(){

	UpdateBoundUntetheredList();
	if(n_bound_untethered_ > 0){
		// Pick a random free tethered motor
		int i_motor = properties_->gsl.GetRanInt(n_bound_untethered_);
		Kinesin *motor = bound_untethered_list_[i_motor];
		motor->UpdateNeighborXlinks();
		int x_dist_dub = motor->SampleTailExtensionDoubled();
		// Ensure motor has a neighbor xlink for desired tether extension
		// Roll for a new extension after a certain # of tries
		int attempts = 0;
		int switches = 0; 
		bool failure = false;
		while(motor->NeighborXlinkExists(x_dist_dub) == false){
			if(attempts > 2*n_bound_untethered_){
				if(switches > 50){
					printf("failed to tether_bound\n");
					failure = true;
					break;
				}
				x_dist_dub = motor->SampleTailExtensionDoubled();
				attempts = 0;
				switches++;
			}
			i_motor = properties_->gsl.GetRanInt(n_bound_untethered_);
			motor = bound_untethered_list_[i_motor];
			attempts++;
		}
		if(failure == false){
			AssociatedProtein* xlink = motor->GetNeighborXlink(x_dist_dub);
			if(xlink == nullptr){
				printf("why tho (kmc tether_bound)\n");
				exit(1);
			}
			// Update motor and xlink details
			motor->xlink_ = xlink;
			motor->tethered_ = true;
			xlink->tethered_ = true;
			xlink->motor_ = motor;
			// Update statistics
			motor->UpdateExtension();
			if(motor->x_dist_doubled_ != x_dist_dub){
				printf("somethin weird in Tether_Bound: ");
				printf("%i on motor; %i desired\n", motor->x_dist_doubled_, 
						x_dist_dub);
			}
			n_bound_tethered_[x_dist_dub]++;
			n_bound_tethered_tot_++;
			n_bound_untethered_--;
			properties_->prc1.n_tethered_++;
			properties_->prc1.n_untethered_--;
		}
		else{
			printf("hu ");
		}

	}
	else{
		printf("Error in Tether_Bound: no bound untethered motors!\n");
		exit(1);
	}
}

void KinesinManagement::KMC_Switch(){

	// Make sure there is at least one switchable motor 
	UpdateSwitchableList();
	if(n_switchable_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_switchable_);
		Kinesin *motor = switchable_list_[i_entry]; 
		Tubulin *old_front_site = motor->front_site_;
		Tubulin	*old_rear_site = motor->rear_site_;
		Microtubule *mt = motor->mt_; 
		Microtubule *mt_adj = mt->neighbor_; 
		int offset = mt_adj->coord_ - mt->coord_;
		int i_front_new = old_rear_site->index_ - offset;
		int i_rear_new = old_front_site->index_ - offset;
		Tubulin *new_front_site = &mt_adj->lattice_[i_front_new];
		Tubulin *new_rear_site = &mt_adj->lattice_[i_rear_new];
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
		printf("Error in RunKMC_Switch(): no switchable motors\n");
		exit(1);
    }
}

void KinesinManagement::KMC_Untether_Bound(int x_dist_doubled){

	UpdateBoundTetheredTable();
	int n_bound_tethered = n_bound_tethered_[x_dist_doubled];
	if(n_bound_tethered  > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound_tethered);
		Kinesin *motor = bound_tethered_table_[x_dist_doubled][i_entry];
		motor->UpdateExtension();
		int x_dub_pre = motor->x_dist_doubled_;
		if(x_dub_pre != x_dist_doubled){
			printf("error in Untether_Bound (motor): ");
			printf("%i received; %i in motor\n", x_dist_doubled, x_dub_pre);
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
		n_bound_tethered_[x_dub_pre]--;
		n_bound_tethered_tot_--;
		n_bound_untethered_++;
		properties_->prc1.n_tethered_--;
		properties_->prc1.n_untethered_++;
	}
	else{
		printf("Error in Untether_Bound: no bound tethered motors!\n");
		exit(1);
	}
}

void KinesinManagement::KMC_Untether_Free(){

	UpdateFreeTetheredList();
	if(n_free_tethered_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_free_tethered_);
		Kinesin *motor = free_tethered_list_[i_entry];
		AssociatedProtein *xlink = motor->xlink_;
		// Update motor and xlink detail
		xlink->motor_ = nullptr; 
		xlink->tethered_ = false; 
		motor->xlink_ = nullptr;
		motor->tethered_ = false;
		// Update statistics
		n_free_tethered_--;
		properties_->prc1.n_tethered_--;
		properties_->prc1.n_untethered_++;

	}
	else{
		printf("Error in Untether_Free: no free tethered motors!\n");
		exit(1);
	}
}

void KinesinManagement::KMC_Step_Tethered(int x_dist_doubled){

	UpdateStepableTetheredTable();
	int n_stepable = n_stepable_tethered_[x_dist_doubled]; 
	if(n_stepable > 0){
		int i_entry = properties_->gsl.GetRanInt(n_stepable);
		Kinesin *motor = stepable_tethered_table_[x_dist_doubled][i_entry];
		int x_dub_pre = motor->x_dist_doubled_;
		if(x_dub_pre != x_dist_doubled){
			printf("error in Step_Tethered (motor): ");
			printf("%i received; %i in motor\n", x_dist_doubled, x_dub_pre);
			exit(1);
		}
		Microtubule *mt = motor->mt_;
		Tubulin *old_front_site = motor->front_site_;
		Tubulin *old_rear_site = motor->rear_site_;
		int i_old_front = old_front_site->index_;
		int dx = mt->delta_x_;
		Tubulin *new_front_site = &mt->lattice_[i_old_front + dx];
		Tubulin *new_rear_site = old_front_site; 
		// Update front head of motor
		motor->front_site_ = new_front_site;
		new_front_site->motor_ = motor;
		new_front_site->occupied_ = true;
		// Update rear head of motor 
		motor->rear_site_ = new_rear_site;
		old_rear_site->motor_ = nullptr;
		old_rear_site->occupied_ = false;
		// Update statistics
		motor->UpdateExtension();
		int x_dub_post = motor->x_dist_doubled_;
		int plus_end = motor->mt_->plus_end_; 
		int i_front = motor->front_site_->index_;
		if(x_dub_pre == x_dub_post
		&& i_front != plus_end){
			printf("error in step_teth (motor)\n");
			//exit(1);
		}
		n_bound_tethered_[x_dub_pre]--;
		n_bound_tethered_[x_dub_post]++;
	}
	else{
		printf("Error in Step_Tethered (not rly tho): no stepable motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Step_Untethered(){

    // Make sure there is at least one stepable untethered motor
	UpdateStepableUntetheredList();
	if(n_stepable_untethered_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_stepable_untethered_);
		Kinesin *motor = stepable_untethered_list_[i_entry];
		Microtubule *mt = motor->mt_; 
		Tubulin *old_front_site = motor->front_site_; 
		Tubulin *old_rear_site = motor->rear_site_;
		int i_old_front = old_front_site->index_;
		int dx = mt->delta_x_;
		Tubulin *new_front_site = &mt->lattice_[i_old_front + dx];
		Tubulin *new_rear_site = old_front_site; 
		// Update front head of motor
		motor->front_site_ = new_front_site;
		new_front_site->motor_ = motor;
		new_front_site->occupied_ = true;
		// Update rear head of motor 
		motor->rear_site_ = new_rear_site;
		old_rear_site->motor_ = nullptr;
		old_rear_site->occupied_ = false;
	}
	else{
		printf("Error in Step_Untethered: no stepable untethered motors\n");
    }
}

void KinesinManagement::KMC_Boundaries(int n_events){

    int n_mts = parameters_->n_microtubules;
    double alpha_eff = (alpha_*p_step_untethered_)/n_events;
    double beta_eff = (beta_*p_step_untethered_)/n_events;
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
            n_bound_untethered_++;
        }
        // Remove motors from plus_end with probability beta_eff
        if(ran2 < beta_eff
                && plus_end->motor_ != nullptr
                && plus_neighbor->motor_ != nullptr){
            // Get motor bound to minus_end
            Kinesin *motor = plus_end->motor_;
	//		motor->UpdateExtension();
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
			if(motor->tethered_ == true){
				int x_dist_dub = motor->x_dist_doubled_;
				n_bound_tethered_[x_dist_dub]--; 
			}
			else{
				n_bound_untethered_++;
			}
        } 
    }
}

