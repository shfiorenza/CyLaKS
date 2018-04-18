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
	printf("%i motors\n", n_motors_);
    motor_list_.resize(n_motors_);
    for(int ID = 0; ID < n_motors_; ID++){
        motor_list_[ID].Initialize(parameters_, properties_, ID);
    }
}

void KinesinManagement::SetParameters(){

//	tau_ = parameters_->tau_m;
    double delta_t = parameters_->delta_t;
	double site_size = parameters_->site_size;
	// Statistics for diffusion
	double D_const = parameters_->D_motor;
	double x_squared = (site_size/1000)*(site_size/1000); // convert to um^2
	tau_ = x_squared / (2 * D_const);
	p_diffuse_fwd_untethered_ = delta_t / tau_;
	p_diffuse_bck_untethered_ = delta_t / tau_; 
	// Generate stepping rates based on extension of tether: rates are 
	// increased if stepping towards rest length, and reduced if stepping
	// away from rest length (which increases tether extension)
	dist_cutoff_ = motor_list_[0].dist_cutoff_;
	comp_cutoff_ = motor_list_[0].comp_cutoff_;
	rest_dist_ = motor_list_[0].rest_dist_;
	p_diffuse_to_tether_rest_.resize(2*dist_cutoff_ + 1);
	p_diffuse_from_tether_rest_.resize(2*dist_cutoff_ + 1);
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
		// Calculate extension if motor diffuses toward/away from rest
		double dr_toward; 
		double dr_from; 
		// If extended (pos extension), r_bck is towards rest
		if(dr >= 0){
			dr_toward = r_bck - r_0;
			dr_from = r_fwd - r_0;
		}
		// If compressed (neg extension), r_fwd is towards rest
		else{
			dr_toward = r_fwd - r_0;
			dr_from = r_bck - r_0; 
		}
		double dU_from, 
			   dU_to, 
			   weight_to, 
			   weight_from;
		if(x_dist_dub == 2*rest_dist_){
			// Weights according to Lanksy et al. 
			dU_from = (k_eff_slack/2)*(dr_from*dr_from - dr*dr);
			dU_to = (1/2)*(k_spring*dr_toward*dr_toward - k_eff_slack*dr*dr);
			weight_to = exp(-dU_to/(2*kbT));
		   	weight_from = exp(-dU_from/(2*kbT));
        }
		// use k_spring if extension is positive
		else if(dr > 0){
			dU_to = (k_spring/2)*(dr_toward*dr_toward - dr*dr);
			dU_from = (k_spring/2)*(dr_from*dr_from - dr*dr);
			weight_to = exp(-dU_to/(2*kbT));
			if(x_dist_dub >= 2*dist_cutoff_)
				weight_from = 0;
			else
		   		weight_from = exp(-dU_from/(2*kbT));
		}
		// otherwise, use k_eff to model 'slack' in the tether
		else{
			dU_to = (k_eff_slack/2)*(dr_toward*dr_toward - dr*dr);
			dU_from	= (k_eff_slack/2)*(dr_from*dr_from - dr*dr);
			if(x_dist_dub < 2*comp_cutoff_)
				weight_to = 0;
			else
		   		weight_to = exp(-dU_to/(2*kbT));
			if(x_dist_dub <= 2*comp_cutoff_)
				weight_from = 0;
			else
		   		weight_from = exp(-dU_from/(2*kbT));
		}
		double p_to = weight_to * delta_t / tau_;
		double p_from = weight_from * delta_t / tau_;
		p_diffuse_from_tether_rest_[x_dist_dub] = p_from;
		p_diffuse_to_tether_rest_[x_dist_dub] = p_to;
	}
	// Statistics for KMC
    double k_on = parameters_->k_on_motor;
    double c_motor = parameters_->c_motor;
	p_bind_i_free_ = k_on * c_motor * delta_t;
	p_bind_i_tethered_ = p_bind_i_free_ * 4;
	double c_eff_motor_bind = parameters_->c_eff_motor_bind;
	p_bind_ii_ = k_on * c_eff_motor_bind * delta_t;
	double k_off_pseudo = parameters_->k_off_pseudo; 
	p_unbind_pseudo_ = k_off_pseudo * delta_t;		// FIXME
    double k_off = parameters_->k_off_motor;
	p_unbind_stepable_untethered_ = k_off * delta_t;
	double unbind_ratio = parameters_->k_off_ratio;
	p_unbind_stalled_untethered_ = k_off * delta_t / unbind_ratio;
	p_unbind_tethered_ = p_unbind_stepable_untethered_ / 2; 
	c_eff_ = parameters_->c_eff_motor_teth;
	double k_tether_free = parameters_->k_tether_free;
	p_tether_free_ = k_tether_free * c_motor * delta_t;
	double k_untether_free = parameters_->k_untether_free;
	p_untether_free_ = k_untether_free * delta_t;
    double motor_speed = parameters_->motor_speed;
	p_step_untethered_ = motor_speed * delta_t / site_size;
	double k_failstep = parameters_->failstep_rate; 
	p_failstep_untethered_ = k_failstep * delta_t;
    double switch_rate = parameters_->switch_rate;
	p_switch_ = switch_rate*delta_t;
	// Generate untethering and stepping rates for all tether extensions	
	// Everything is 2*dist_cutoff to allow for half-integer distances, 
	// so the 3rd entry will correspond to a distance of 1.5, etc. 
	double k_unteth = parameters_->k_untether;
	double stall_force = motor_list_[0].stall_force_;		// in pN
	p_untether_bound_.resize(2*dist_cutoff_ + 1);
	p_step_to_teth_rest_.resize(2*dist_cutoff_ + 1);
	p_step_from_teth_rest_.resize(2*dist_cutoff_ + 1);
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		double r_x = x_dist_dub * site_size / 2;
		double r = sqrt(r_y*r_y + r_x*r_x);
		double dr = r - r_0;
		double cosine = r_x / r;
		// If extension is positive, treat tether like a spring
		if(dr >= 0){
			// Calculate spring potential energy (PE) for this x_dist
			double U_teth = (k_spring/2)*dr*dr;
			// Get untethering weight for this PE
			double unteth_weight = exp(U_teth/(2*kbT));
			p_untether_bound_[x_dist_dub] = k_unteth*unteth_weight*delta_t;
			// Stepping probability has force-dependent relation
			// (see sci advi suppl)
			double force = dr*k_spring;
			double coeff_to = 1 + cosine * (force / stall_force);
			double coeff_from = 1 - cosine * (force / stall_force); 
			double p_to = coeff_to * p_step_untethered_;
			double p_from = coeff_from * p_step_untethered_; 
			if(x_dist_dub >= 2*dist_cutoff_){
				p_step_to_teth_rest_[x_dist_dub] = p_to;
				p_step_from_teth_rest_[x_dist_dub] = 0;
			}
			else if(force < stall_force){
				p_step_to_teth_rest_[x_dist_dub] = p_to;
				p_step_from_teth_rest_[x_dist_dub] = p_from;
			}
			else{
				p_step_to_teth_rest_[x_dist_dub] = p_to; 
				p_step_from_teth_rest_[x_dist_dub] = 0;
			}
		}
		// Otherwise, use k_eff for slack
		else{
			double U_teth = (k_eff_slack/2)*dr*dr;
			double unteth_weight = exp(U_teth/(2*kbT));
			p_untether_bound_[x_dist_dub] = k_unteth*unteth_weight*delta_t;
			double force = -1 * dr * k_eff_slack;
			double coeff_to = 1 + cosine * (force / stall_force);
			double coeff_from = 1 - cosine * (force / stall_force);
			double p_to = coeff_to * p_step_untethered_;
			double p_from = coeff_from * p_step_untethered_;
			if(x_dist_dub < 2*comp_cutoff_){
				p_step_to_teth_rest_[x_dist_dub] = 0;
				p_step_from_teth_rest_[x_dist_dub] = 0;
			}
			else if(x_dist_dub == 2*comp_cutoff_){
				p_step_to_teth_rest_[x_dist_dub] = p_to;
				p_step_from_teth_rest_[x_dist_dub] = 0;
			}
			else if(force < stall_force){
				p_step_to_teth_rest_[x_dist_dub] = p_to;
				p_step_from_teth_rest_[x_dist_dub] = p_from;
			}
			else{
				p_step_to_teth_rest_[x_dist_dub] = 0;
				p_step_from_teth_rest_[x_dist_dub] = 0;
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
	bound_untethered_list_.resize(n_motors_);
	stepable_untethered_list_.resize(n_motors_);
	stalled_untethered_list_.resize(n_motors_);
	bound_tethered_list_.resize(n_motors_);
	stepable_untethered_list_.resize(n_motors_);
	switchable_list_.resize(n_motors_);
	// Two dimensional stuff
	n_bound_tethered_.resize(2*dist_cutoff_ + 1);
	n_stepable_to_teth_rest_.resize(2*dist_cutoff_ + 1);
	n_stepable_from_teth_rest_.resize(2*dist_cutoff_ + 1);
	bound_tethered_table_.resize(2*dist_cutoff_ + 1);
	stepable_to_rest_table_.resize(2*dist_cutoff_ + 1);
	stepable_from_rest_table_.resize(2*dist_cutoff_ + 1);
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		bound_tethered_table_[x_dist_dub].resize(n_motors_);
		stepable_to_rest_table_[x_dist_dub].resize(n_motors_);
		stepable_from_rest_table_[x_dist_dub].resize(n_motors_);
		n_bound_tethered_[x_dist_dub] = 0;
		n_stepable_to_teth_rest_[x_dist_dub] = 0;
		n_stepable_from_teth_rest_[x_dist_dub] = 0;
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

/*
bool KinesinManagement::BoundaryStatus(Kinesin *motor){

    if(motor->front_site_->index_ == motor->mt_->plus_end_
	|| motor->rear_site_->index_ == motor->mt_->minus_end_){
        return true;
    }
    else{
        return false;
    }
}
*/

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
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motor_list_[i_motor]; 
		if(motor->heads_active_ == 1){
			Microtubule *mt = motor->mt_;
			Tubulin *bound_site = motor->GetActiveHeadSite();
			int i_site = bound_site->index_;
			int i_plus = mt->plus_end_;
			int i_minus = mt->minus_end_;
			int dx = mt->delta_x_;
			// Don't access lattice sites that don't exist (or boundaries)
			if(i_site == i_plus){
				if(mt->lattice_[i_site - dx].occupied_ == false){
					eligible_pseudo_list_[i_entry] = motor;
					i_entry++;
					n_eligible_pseudo_++;
				}
			}
			else if(i_site == i_minus){
				if(mt->lattice_[i_site + dx].occupied_ == false){
					eligible_pseudo_list_[i_entry] = motor;
					i_entry++;
					n_eligible_pseudo_++;
				}
			}
			else{
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
		printf("something wrong in update_unteth_list (motors 1D)");
		printf(" %i in statistics, %i entries tho\n", n_bound_untethered_, 
				n_entries);
		exit(1);
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
            // Exclude plus_end (can't step off of MTs)
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

void KinesinManagement::UpdateStalledUntetheredList(){

	n_stalled_untethered_ = 0; 
	int i_entry = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motor_list_[i_motor];
        if(motor->heads_active_ == 2
		&& motor->tethered_ == false){
            Microtubule *mt = motor->mt_;
            int plus_end = mt->plus_end_;
            int delta_x = mt->delta_x_;
            int i_front_site = motor->front_site_->index_;
            // Motors on plus-end are automatically considered stalled
            if(i_front_site == plus_end){
				stalled_untethered_list_[i_entry] = motor;
				i_entry++;
				n_stalled_untethered_++;
			}
			else if(mt->lattice_[i_front_site + delta_x].occupied_ == true){
				stalled_untethered_list_[i_entry] = motor;
				i_entry++; 
				n_stalled_untethered_++;
			}
        }
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

void KinesinManagement::UpdateStepableTetheredTables(){

	int i_entry_to[2*dist_cutoff_ + 1];
	int i_entry_from[2*dist_cutoff_ + 1];
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		i_entry_to[x_dist_dub] = 0;
		i_entry_from[x_dist_dub] = 0;
		n_stepable_to_teth_rest_[x_dist_dub] = 0; 	
		n_stepable_from_teth_rest_[x_dist_dub] = 0;
	}
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motor_list_[i_motor];
		if(motor->heads_active_ == 2
		&& motor->tethered_ == true){
			motor->UpdateExtension();
			Microtubule *mt = motor->mt_;
			int plus_end = mt->plus_end_;
			int delta_x = mt->delta_x_;
			int dx_rest = motor->GetDirectionTowardRest();
			int i_front = motor->front_site_->index_;
			// end pausing
			if(i_front != plus_end
			&& motor->AtCutoff() == false){
				if(mt->lattice_[i_front + delta_x].occupied_ == false){
					int x_dist_dub = motor->x_dist_doubled_;
					// if MT's dx is towards rest, add to to_rest list
					if(delta_x == dx_rest){
						int index = i_entry_to[x_dist_dub];
						stepable_to_rest_table_[x_dist_dub][index] = motor;	
						i_entry_to[x_dist_dub]++;
						n_stepable_to_teth_rest_[x_dist_dub]++;
					}
					// otherwise, add to from_rest list
					else if(delta_x == -1 * dx_rest){
						int index = i_entry_from[x_dist_dub];
						stepable_from_rest_table_[x_dist_dub][index] = motor;
						i_entry_from[x_dist_dub]++;
						n_stepable_from_teth_rest_[x_dist_dub]++;
					}
					else{
						printf("hmmm??? teth tables\n");
						exit(1);
					}
				}
			}
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
		&& motor->tethered_ == true){
//		&& BoundaryStatus(motor) == false){
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

void KinesinManagement::GenerateDiffusionList(){

	int n_events = 0; 
	// Untethered statistics
	UpdateBoundUntetheredList();
	int n_fwd_unteth = GetNumToStepForward_Unteth();
	int n_bck_unteth = GetNumToStepBackward_Unteth();
	while(n_fwd_unteth + n_bck_unteth > n_bound_untethered_){
		double ran = properties_->gsl.GetRanProb();
		if(ran < 0.5
		&& n_fwd_unteth > 0)
			n_fwd_unteth--;
		else if(n_bck_unteth > 0)
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
			double p_to = p_diffuse_to_tether_rest_[x_dist_dub];
			double p_from = p_diffuse_from_tether_rest_[x_dist_dub];
			double tot_prob = p_to + p_from;  // XXX is this correct? XXX
			if(ran < p_to/tot_prob
			&& n_to > 0){
				n_toward_rest[x_dist_dub]--;
				n_to = n_toward_rest[x_dist_dub];
			}
			else if(n_from > 0){
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
		for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
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
	int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_bound);
	return n_to_step;
}

int KinesinManagement::GetNumToStepBackward_Unteth(){

	int n_bound = n_bound_untethered_;
	double p_step = p_diffuse_bck_untethered_;
	int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_bound);
	return n_to_step;
}

int KinesinManagement::GetNumToStepTowardRest(int x_dist_doubled){

	int n_bound = n_bound_tethered_[x_dist_doubled];
	double p_step = p_diffuse_to_tether_rest_[x_dist_doubled];
	int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_bound);
	return n_to_step;
}

int KinesinManagement::GetNumToStepFromRest(int x_dist_doubled){

	int n_bound = n_bound_tethered_[x_dist_doubled];
	double p_step = p_diffuse_from_tether_rest_[x_dist_doubled];
	int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_bound);
	return n_to_step;
}

void KinesinManagement::RunDiffusion(){

//	printf("Start of kinesin diffusion cycle\n");
	GenerateDiffusionList();
	if(diffusion_list_.empty() == false){
		int n_events = diffusion_list_.size();
//		printf("%i KINESIN DIFFUSION EVENTS\n", n_events);
		//printf("%i DIFFUSION EVENTS\n", n_events);
		int x_dist_dub;
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
//			properties_->wallace.PrintMicrotubules(0.000);
			switch(diff_event){
				case 10:
//					printf("unteth step fwd\n");
					RunDiffusion_Forward_Untethered();
					break;	
				case 20:
//					printf("unteth step bck\n");
					RunDiffusion_Backward_Untethered();
					break;
				case 30:
//					printf("teth step to (%i)[%i avail]\n", x_dist_dub, 
//							n_bound_tethered_[x_dist_dub]);
					RunDiffusion_Toward_Rest(x_dist_dub);
					break;
				case 40:
//					printf("teth step from (%i)[%i avail]\n", x_dist_dub, 
//							n_bound_tethered_[x_dist_dub]);
					RunDiffusion_From_Rest(x_dist_dub);
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
		int i_plus = mt->plus_end_;
		if(i_front != i_plus){
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
	//			printf("oh well fwd\n");
			}
		}
		else{
	//		printf("cant diffuse outta this one brotha\n");
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
		int i_minus = mt->minus_end_;
		if(i_rear != i_minus){
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
	//			printf("oh well bck\n");
			}
		}
		else{
	//		printf("cant diffuse outta this one brotha\n");
		}

	}
	else{
		printf("ya blew it. we failed to step untethered motor bck\n");
		exit(1);
	}
}

void KinesinManagement::RunDiffusion_Toward_Rest(int x_dist_doubled){

	int mt_length = parameters_->length_of_microtubule;
	int mt_array_length = mt_length - 1;	
	UpdateBoundTetheredTable();
	int n_bound = n_bound_tethered_[x_dist_doubled];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_tethered_table_[x_dist_doubled][i_entry];
		Microtubule *mt = motor->mt_;
		Tubulin *near_site = motor->GetSiteCloserToRest();
		int dx = motor->GetDirectionTowardRest();
		int i_near = near_site->index_;
		if(!(i_near == mt_array_length && dx == 1)
		&& !(i_near == 0 && dx == -1)){
			if(mt->lattice_[i_near + dx].occupied_ == false){
				Tubulin *new_front, *new_rear, 
						*old_front, *old_rear;
				if(near_site == motor->front_site_){
					old_front = near_site;
					new_front = &mt->lattice_[i_near + dx];
					old_rear = motor->rear_site_;
					new_rear = motor->front_site_;
					// Update new site
					new_front->motor_ = motor;
					new_front->occupied_ = true;
					// Update old site
					old_rear->motor_ = nullptr;
					old_rear->occupied_ = false;
				}
				else if(near_site == motor->rear_site_){
					old_rear = near_site;
					new_rear = &mt->lattice_[i_near + dx];
					old_front = motor->front_site_;
					new_front = motor->rear_site_;
					// Update new site
					new_rear->motor_ = motor;
					new_rear->occupied_ = true;
					// Update old site
					old_front->motor_ = nullptr;
					old_front->occupied_ = false;
				}
				else{
					printf("woah woah woah. why is she the WRONG site?");
					printf(" (motors diffuse toward tether)\n");
					exit(1);
				}
				int x_dub_pre = motor->x_dist_doubled_;
				// Update motor
				motor->front_site_ = new_front;
				motor->rear_site_ = new_rear;
				motor->UpdateExtension();
				// Update statistics
				// Make sure an untether event wasn't forced
				if(motor->tethered_ == true){
					int x_dub_post = motor->x_dist_doubled_;
					n_bound_tethered_[x_dub_pre]--;
					n_bound_tethered_[x_dub_post]++;
					// Update prc1 site statistics
					AssociatedProtein *xlink = motor->xlink_;
					AssociatedProteinManagement *prc1 = &properties_->prc1;
					if(xlink->heads_active_ == 1){
						prc1->n_sites_i_tethered_[x_dub_pre]--;
						prc1->n_sites_i_tethered_[x_dub_post]++;
					}
					else if(xlink->heads_active_ == 2){
						int x_dist = xlink->x_dist_;
						prc1->n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
						prc1->n_sites_ii_tethered_[x_dub_post][x_dist] += 2;
					}
					else{
						printf("wat in diff_teth_to\n");
						exit(1);
					}
				}
			}
			else{
	//			printf("aww darn teth toward (%i ext)\n", x_dist_doubled);
			}
		}	
	}
	else{
		printf("ya blew it. we failed to step tethered motor towards\n");
		exit(1);
	}
}

void KinesinManagement::RunDiffusion_From_Rest(int x_dist_doubled){

	int mt_length = parameters_->length_of_microtubule;
	int mt_array_length = mt_length - 1;	
	UpdateBoundTetheredTable();
	int n_bound = n_bound_tethered_[x_dist_doubled];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_tethered_table_[x_dist_doubled][i_entry];
		Microtubule *mt = motor->mt_;
		Tubulin *far_site = motor->GetSiteFartherFromRest();
		int dx = motor->GetDirectionTowardRest();
		int i_far = far_site->index_;
//		printf("i: %i, dx: %i\n", i_far, dx);
		if(!(i_far == mt_array_length && dx == -1)
		&& !(i_far == 0 && dx == 1)){
			if(mt->lattice_[i_far - dx].occupied_ == false){
				Tubulin *new_front, *new_rear, 
						*old_front, *old_rear;
				if(far_site == motor->front_site_){
					old_front = far_site;
					new_front = &mt->lattice_[i_far - dx];
					old_rear = motor->rear_site_;
					new_rear = motor->front_site_;
					// Update new site
					new_front->motor_ = motor;
					new_front->occupied_ = true;
					// Update old site
					old_rear->motor_ = nullptr;
					old_rear->occupied_ = false;
				}
				else if(far_site == motor->rear_site_){
					old_rear = far_site;
					new_rear = &mt->lattice_[i_far - dx];
					old_front = motor->front_site_;
					new_front = motor->rear_site_;
					// Update new site  
					new_rear->motor_ = motor;
					new_rear->occupied_ = true;
					// Update old site
					old_front->motor_ = nullptr;
					old_front->occupied_ = false; 
				}
				else{
					printf("woah woah woah. why is she the WRONG site?");
					printf(" (motors diffuse away from tether)\n");
					exit(1);
				}
				int x_dub_pre = motor->x_dist_doubled_;
				// Update motor
				motor->front_site_ = new_front;
				motor->rear_site_ = new_rear;
				motor->UpdateExtension();
				// Update statistics
				// Make sure an untether event wasn't forced
				if(motor->tethered_ == true){
					int x_dub_post = motor->x_dist_doubled_;
					n_bound_tethered_[x_dub_pre]--;
					n_bound_tethered_[x_dub_post]++;
					// Update prc1 site statistics
					AssociatedProtein *xlink = motor->xlink_;
					AssociatedProteinManagement *prc1 = &properties_->prc1;
					if(xlink->heads_active_ == 1){
						prc1->n_sites_i_tethered_[x_dub_pre]--;
						prc1->n_sites_i_tethered_[x_dub_post]++;
					}
					else if(xlink->heads_active_ == 2){
						int x_dist = xlink->x_dist_;
						prc1->n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
						prc1->n_sites_ii_tethered_[x_dub_post][x_dist] += 2;
					}
					else{
						printf("wat in diff_teth_from\n");
						exit(1);
					}
				}
			}
			else{
	//			printf("aw darn teth away from (%i ext)\n", x_dist_doubled);
			}
		}
	}
	else{
		printf("ya blew it. we failed to step tethered motor away from\n");
		exit(1);
	}
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
	UpdateStepableUntetheredList();
	int n_unbind_stepable_ut = GetNumToUnbind_Stepable_Untethered();
	int n_step_ut = GetNumToStep_Untethered(); 
	while(n_unbind_stepable_ut + n_step_ut > n_stepable_untethered_){
		if(n_step_ut > 0){
			n_step_ut--;
		}
		else if(n_unbind_stepable_ut > 0){
			n_unbind_stepable_ut--;
		}
	}
	n_events += n_unbind_stepable_ut;
	n_events += n_step_ut;
	UpdateStalledUntetheredList();
	int n_unbind_stalled_ut = GetNumToUnbind_Stalled_Untethered();
	int n_failstep_ut = GetNumToFailstep_Untethered();
	while(n_unbind_stalled_ut + n_failstep_ut > n_stalled_untethered_){
		if(n_failstep_ut > 0){
			n_failstep_ut--;
		}
		else if(n_unbind_stalled_ut > 0){
			n_unbind_stalled_ut--;
		}
	}
	n_events += n_unbind_stalled_ut;
	n_events += n_failstep_ut;
	UpdateBoundTetheredList();
	int n_unbind_t = GetNumToUnbind_Tethered();
	n_events += n_unbind_t;
	UpdatePseudoBoundList();
	int n_unbind_p = GetNumToUnbind_Pseudo();
	while(n_unbind_p + n_bind_ii > n_pseudo_bound_){
		if(n_bind_ii > 0){
			n_bind_ii--;
			n_events--;
		}
		else if (n_unbind_p > 0){
			n_unbind_p--;
		}
	}
	n_events += n_unbind_p;
	// XXX is the following necessary? think about it
/*
	while(n_bind_ii + n_unbind_p > n_eligible_pseudo_
	&& n_eligible_pseudo_ > 0){
		if(n_bind_ii > 0){
			n_bind_ii--;
			n_events--;
		}
	}
*/
	int n_tether_free = GetNumToTether_Free();
	n_events += n_tether_free;
	UpdateBoundUntetheredList();
	int n_tether_bound = GetNumToTether_Bound();
	n_events += n_tether_bound;
	
// XXX HARD CODE DISABLE OF SWITCHING;XXX DO NOT FORGET ABOUT!! XXX
/*	UpdateSwitchableList();
    int n_switch = 0; //GetNumToSwitch();	
	n_events += n_switch;
*/	
	int n_untether_free = GetNumToUntether_Free();
	while(n_untether_free + n_bind_i_tethered > n_free_tethered_){
			if(n_bind_i_tethered > 0){
				n_bind_i_tethered--;
				n_events--;	
			}
			else if(n_untether_free > 0){
				n_untether_free--;
			}
	}
	n_events += n_untether_free;
	
	// XXX BELOW IS BOTCHED/INCOMPLETE; CAUTION WHEN UNCOMMENTING
/*
		double p_mobile = p_unbind_unteth_mobile_;
		double p_stalled = p_unbind_unteth_stalled_;
		double p_tot = p_mobile + p_stalled;
		double ran = properties_->gsl.GetRanProb();
		else if(n_tether_bound > 0){
			n_tether_bound--;
			n_events--;
		}
		else if(ran < p_mobile/p_tot
		&& n_bound_unteth_mobile_ > 0){
			n_unbind_ut_mobile--;
			n_events--;
		}
		else if(n_bound_unteth_stalled_ > 0){
			n_unbind_ut_stalled--;
			n_events--;
		}
		else if(n_bound_unteth_mobile_ > 0){
			n_unbind_ut_mobile--;
			n_events--;
		}
	}
*/
	// XXX END BOTCHED WARNING HERE XXX
	// Handle the statistics of differnt tether extensions separately 
	UpdateStepableTetheredTables();
	UpdateBoundTetheredTable();
	int n_untether_bound[2*dist_cutoff_ + 1];
	int n_step_to_rest[2*dist_cutoff_ + 1];
	int n_step_from_rest[2*dist_cutoff_ + 1];
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		n_untether_bound[x_dist_dub] = GetNumToUntether_Bound(x_dist_dub);
		n_step_to_rest[x_dist_dub] = GetNumToStep_ToTethRest(x_dist_dub);
		n_step_from_rest[x_dist_dub] = GetNumToStep_FromTethRest(x_dist_dub);
		int n_unteth = n_untether_bound[x_dist_dub];
		int n_step_to = n_step_to_rest[x_dist_dub];
		int n_step_from = n_step_from_rest[x_dist_dub];
		int n_tethered = n_bound_tethered_[x_dist_dub];
		// Prevent too many steps from occuring (imagine if 2 motors
		// are bound and we roll for 1 unbind but 2 steps)
		while(n_unteth + n_step_to + n_step_from > n_tethered){
//			printf("*** STAT CORRECTION ***\n");
//			printf("to unteth: %i, to step: %i, avail: %i\n", 
//					n_unteth, n_step_teth, n_tethered);
			double p_to = p_step_to_teth_rest_[x_dist_dub];
			double p_from = p_step_from_teth_rest_[x_dist_dub];
			double p_tot = p_to + p_from;
			double ran = properties_->gsl.GetRanProb();
			if(n_unteth > 10000){
				printf("lots o untetheres man\n");
				exit(1);
			}
			if(n_tethered < 0){
				printf("negative tethered mots man\n");
				exit(1);
			}
			if(ran < p_to / p_tot
			&& n_step_to > 0){
				n_step_to_rest[x_dist_dub]--;
				n_step_to = n_step_to_rest[x_dist_dub];
			}
			else if(n_step_from > 0){
				n_step_from_rest[x_dist_dub]--;
				n_step_from = n_step_from_rest[x_dist_dub];
			}
			else if(n_step_to > 0){
				n_step_to_rest[x_dist_dub]--;
				n_step_to = n_step_to_rest[x_dist_dub];
			}
			else if (n_unteth > 0){
				n_untether_bound[x_dist_dub]--;
				n_unteth = n_untether_bound[x_dist_dub];
			}
		} 
		n_events += n_untether_bound[x_dist_dub];
		n_events += n_step_to_rest[x_dist_dub];
		n_events += n_step_from_rest[x_dist_dub];
	}
	
//	printf("n_events: %i \n", n_events);
    if(n_events > 0){
        int pre_list[n_events];
		int kmc_index = 0;
		for(int i_event = 0; i_event < n_bind_i_free; i_event++){
//			printf("10\n");
			pre_list[kmc_index] = 10;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_bind_i_tethered; i_event++){
//			printf("11\n");
			pre_list[kmc_index] = 11;
			kmc_index++;
		}
		
		for(int i_event = 0; i_event < n_bind_ii; i_event++){
//			printf("12\n");
			pre_list[kmc_index] = 12;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_unbind_stepable_ut; i_event++){
//			printf("20\n");
			pre_list[kmc_index] = 20;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_unbind_stalled_ut; i_event++){
//			printf("20\n");
			pre_list[kmc_index] = 21;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_unbind_t; i_event++){
//			printf("21\n");
			pre_list[kmc_index] = 22;
			kmc_index++;
		}
		
		for(int i_event = 0; i_event < n_unbind_p; i_event++){
//			printf("22\n");
			pre_list[kmc_index] = 23;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_tether_free; i_event++){
//			printf("30\n");
			pre_list[kmc_index] = 30;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_tether_bound; i_event++){
//			printf("31\n");
			pre_list[kmc_index] = 31;
			kmc_index++;
		}
/*		for(int i_event = 0; i_event < n_switch; i_event++){
//			printf("40\n");
			pre_list[kmc_index] = 40;
			kmc_index++;
		}
*/		
		for(int i_event = 0; i_event < n_untether_free; i_event++){
//			printf("50\n");
			pre_list[kmc_index] = 50;
			kmc_index++;
		}

		for(int i_event = 0; i_event < n_step_ut; i_event++){
//			printf("60\n");
			pre_list[kmc_index] = 60;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_failstep_ut; i_event++){
			pre_list[kmc_index] = 61;
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
//				printf("%i\n", 500 + x_dist_dub);
				kmc_index++;
			}
			int n_step_to = n_step_to_rest[x_dist_dub];
			for(int i_event = 0; i_event < n_step_to; i_event++){
				pre_list[kmc_index] = 600 + x_dist_dub;
				printf("%i\n", 600 + x_dist_dub);
				kmc_index++;
			}
			int n_step_from = n_step_from_rest[x_dist_dub];
			for(int i_event = 0; i_event < n_step_from; i_event++){
				pre_list[kmc_index] = 700 + x_dist_dub;
				printf("%i\n", 700 + x_dist_dub);
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

    properties_->microtubules.UpdateNumUnoccupied();
    int n_unocc = properties_->microtubules.n_unoccupied_;
	double p_bind = p_bind_i_free_; 
//	double p_avg = p_bind * n_unocc;
//	int n_to_bind = properties_->gsl.SamplePoissonDist(p_avg);
	int n_to_bind = properties_->gsl.SampleBinomialDist(p_bind, n_unocc);
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
	double n_equil = c_eff_ * weights_summed;
	double n_avg = p_bind_i_free_ * n_equil;
	int n_to_bind = properties_->gsl.SamplePoissonDist(n_avg);
	return n_to_bind;
}

int KinesinManagement::GetNumToBind_II(){

	double p_bind = p_bind_ii_; 
	int n_able = n_eligible_pseudo_;
	int n_to_bind = properties_->gsl.SampleBinomialDist(p_bind, n_able);
	return n_to_bind;
}

int KinesinManagement::GetNumToUnbind_Pseudo(){

	int n_bound = n_pseudo_bound_; 
	double p_unbind = p_unbind_pseudo_;
	int n_to_unbind = properties_->gsl.SampleBinomialDist(p_unbind, n_bound);
	return n_to_unbind;
}

int KinesinManagement::GetNumToUnbind_Stepable_Untethered(){

	int n_bound = n_stepable_untethered_;
	double p_unbind = p_unbind_stepable_untethered_;
    int n_to_unbind = properties_->gsl.SampleBinomialDist(p_unbind, n_bound);
    return n_to_unbind;
}

int KinesinManagement::GetNumToUnbind_Stalled_Untethered(){

	int n_bound = n_stalled_untethered_;
	double p_unbind = p_unbind_stalled_untethered_;
    int n_to_unbind = properties_->gsl.SampleBinomialDist(p_unbind, n_bound);
    return n_to_unbind;
}

int KinesinManagement::GetNumToUnbind_Tethered(){
		
	int n_bound = n_bound_tethered_tot_;
	double p_unbind = p_unbind_tethered_;
	int n_to_unbind = properties_->gsl.SampleBinomialDist(p_unbind, n_bound);
	return n_to_unbind; 
}

int KinesinManagement::GetNumToTether_Free(){

    int n_untethered_xlinks = properties_->prc1.n_untethered_;
	// Calculate how many free motors tether within delta_t on avg
	int n_unteth = n_untethered_xlinks;
	double p_teth = p_tether_free_;
  	int n_to_tether = properties_->gsl.SampleBinomialDist(p_teth, n_unteth);
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
	double n_equil = c_eff_ * weights_summed;  
	double n_avg = p_tether_free_ * n_equil;
	int n_to_teth = properties_->gsl.SamplePoissonDist(n_avg);
	return n_to_teth;
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

int KinesinManagement::GetNumToStep_Untethered(){

	double p_step = p_step_untethered_;
	int n_stepable = n_stepable_untethered_;
    int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_stepable);
    return n_to_step;
}

int KinesinManagement::GetNumToFailstep_Untethered(){

	int n_bound = n_stalled_untethered_;
	double p_unbind = p_failstep_untethered_;
	int n_to_unbind = properties_->gsl.SampleBinomialDist(p_unbind, n_bound);
	return n_to_unbind;
}

int KinesinManagement::GetNumToStep_ToTethRest(int x_dist_doubled){

	double p_step = p_step_to_teth_rest_[x_dist_doubled];
	int n_stepable = n_stepable_to_teth_rest_[x_dist_doubled];
	int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_stepable);
	return n_to_step;
}

int KinesinManagement::GetNumToStep_FromTethRest(int x_dist_doubled){

	double p_step = p_step_from_teth_rest_[x_dist_doubled];
	int n_stepable = n_stepable_from_teth_rest_[x_dist_doubled];
	int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_stepable);
	return n_to_step;
}

int KinesinManagement::GetNumToSwitch(){

	double p_switch = p_switch_;
	int n_able = n_switchable_;
    int n_to_switch = properties_->gsl.SampleBinomialDist(p_switch, n_able);
    return n_to_switch;
}

void KinesinManagement::RunKMC(){

//	printf("Start of Kinesin KMC cycle\n");
    GenerateKMCList();
    if(kmc_list_.empty() == false){
        int n_events = kmc_list_.size();
//		printf("%i MOTOR KMC EVENTS\n", n_events);
		int x_dist_doubled = 0;
        for(int i_event = 0; i_event < n_events; i_event++){
            int kmc_event = kmc_list_[i_event];
			if(kmc_event >= 500 && kmc_event < 600){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 51;
			}
			if(kmc_event >= 600 && kmc_event < 700){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 62; 
			}	
			if(kmc_event >= 700 && kmc_event < 800){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 63;
			}
//			properties_->wallace.PrintMicrotubules(0.000);
            switch(kmc_event){
				case 10:
//						printf("free motor pseudo-bound\n");
						KMC_Bind_I_Free();
                        break;
                case 11:
//						printf("tethered motor pseudo-bound\n");
						KMC_Bind_I_Tethered();
                        break;
				case 12:
//						printf("pseudo motor bound\n");
						KMC_Bind_II();
						break;
			
				case 20:
//						printf("untethered stepable motor unbound\n");
						KMC_Unbind_Stepable_Untethered();
                        break;
				case 21:
//						printf("untethered stalled motor unbound\n");
						KMC_Unbind_Stalled_Untethered();
						break;
                case 22:
//						printf("tethered motor unbound\n");
						KMC_Unbind_Tethered();
                        break;
				case 23:
//						printf("pseudo-bound motor unbound\n");
						KMC_Unbind_Pseudo();
						break;
                case 30:
//						printf("free motor tethered\n");
						KMC_Tether_Free();
                        break;
				case 31:
//						printf("bound motor tethered\n");
						KMC_Tether_Bound();
						break;
				case 40:
//						printf("motor switched\n");
						KMC_Switch();
                        break;
				case 50:
//						printf("free motor untethered\n");
						KMC_Untether_Free();
						break;
				case 51:
//						printf("bound motor (ext %i) untethered\n", 
//								x_dist_doubled);
						KMC_Untether_Bound(x_dist_doubled);
						break;
				case 60:
//						printf("untethered motor stepped\n");
						KMC_Step_Untethered();
						break;
						
				case 61:
//						printf("motor failstepped\n");
						KMC_Failstep_Untethered();
						break;
						
				case 62:
//						printf("tethered motor (ext %i) stepped\n", 
//								x_dist_doubled);
						KMC_Step_ToTethRest(x_dist_doubled);
						break;
				case 63:
						KMC_Step_FromTethRest(x_dist_doubled);
						break;
            }
//          KMC_Boundaries(n_events);
        }
    }
    else{
//      KMC_Boundaries(1);
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
//		exit(1);
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
				if(switches > 50){
//					printf("failed to bind_i_tethered\n");
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
			// update sites for prc1 management
			motor->UpdateExtension();
			if(motor->tethered_ == true){
				AssociatedProtein *xlink = motor->xlink_;
				AssociatedProteinManagement *prc1 = &properties_->prc1;
				int x_dub_post = motor->x_dist_doubled_;
				if(x_dub_post != x_dist_dub){
					printf("WOAH THERE MISTER in bind_i_teth (motor)\n");
					exit(1);
				}
				if(xlink->heads_active_ == 1){
					prc1->n_sites_i_tethered_[x_dist_dub]++;	
					prc1->n_sites_i_untethered_--;
				}
				else if(xlink->heads_active_ == 2){
					int x_dist = xlink->x_dist_;
					prc1->n_sites_ii_tethered_[x_dist_dub][x_dist] += 2;
					prc1->n_sites_ii_untethered_[x_dist] -= 2;
				}
				else{
					printf("GOD DAMNIT bind_i_teth motor\n");
					exit(1);
				}
			}
			else{
				printf("NOPE in bind_i_teth\n");
				exit(1);
			}
		}
	}
	else{
		printf("Error in Bind_I_Tethered: no unoccupied sites.\n");
//		exit(1);
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
		int i_plus = mt->plus_end_;
		int i_minus = mt->minus_end_;
		Tubulin *front_site,
				*rear_site; 
		if(motor->heads_active_ != 1){
			printf("nope in pseudos\n");
			exit(1);
		}
		if(i_site == i_plus){
			if(mt->lattice_[i_site - dx].occupied_ == false){
				front_site = bound_site;
				rear_site = &mt->lattice_[i_site - dx];	
			}
			else{
				printf("something wrong in eligible pseudo list\n");
				exit(1);
			}
		}
		else if(i_site == i_minus){
			if(mt->lattice_[i_site + dx].occupied_ == false){
				front_site = &mt->lattice_[i_site + dx];
				rear_site = bound_site;
			}
			else{
				printf("something wrong in eligible pseudo list\n");
				exit(1);
			}
		}
		else if(mt->lattice_[i_site + dx].occupied_ == false
			 && mt->lattice_[i_site - dx].occupied_ == false){
			double ran = properties_->gsl.GetRanProb();
			if(ran < 0.5){
				front_site = &mt->lattice_[i_site + dx];
				rear_site = bound_site;
			}
			else{
				front_site = bound_site;
				rear_site = &mt->lattice_[i_site - dx];
			}
		}
		else if(mt->lattice_[i_site + dx].occupied_ == false){
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
		// If motor is tethered, we gotta do a whole lot of bullshit
		if(motor->tethered_ == true){
			AssociatedProtein* xlink = motor->xlink_;
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
			// Otherwise, routine statistic stuff
			else{
				// KMC stuff
				int x_dist_dub = motor->x_dist_doubled_; 
				n_bound_tethered_[x_dist_dub]++;
				n_bound_tethered_tot_++;
				// PRC1 sites for diffusion
				AssociatedProteinManagement *prc1 = &properties_->prc1;
				if(xlink->heads_active_ == 1){
					prc1->n_sites_i_tethered_[x_dub_pre]--;
					prc1->n_sites_i_tethered_[x_dist_dub]++;
				}
				else if(xlink->heads_active_ == 2){
					int x_dist = xlink->x_dist_;
					prc1->n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
					prc1->n_sites_ii_tethered_[x_dist_dub][x_dist] += 2;
				}
			}
		}
		// Otherwise, simply add to n_motors that are bound and untethered
		else{
			n_bound_untethered_++;
		}
	}
	else{
		printf("Error in Bind_II: no eligible pseudo. \n");
//		exit(1);
	}
//	properties_->wallace.PrintMicrotubules(2);
}	

void KinesinManagement::KMC_Unbind_Stepable_Untethered(){

    // Make sure that at least one bound motor exists
	UpdateStepableUntetheredList();
    if(n_stepable_untethered_ > 0){
        // Randomly pick a bound motor (that's not on the boundary)
		int i_entry = 0;
		Kinesin *motor;
//		do{
			i_entry = properties_->gsl.GetRanInt(n_stepable_untethered_);
			motor = stepable_untethered_list_[i_entry];
//		}while(BoundaryStatus(motor) == true);
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
	}
	else{
		printf("Error in Unbind_Untethered MOBILE: no bound untethered motors!\n");
  //      exit(1);
    }
}

void KinesinManagement::KMC_Unbind_Stalled_Untethered(){

    // Make sure that at least one bound motor exists
	UpdateStalledUntetheredList();
    if(n_stalled_untethered_ > 0){
        // Randomly pick a bound motor (that's not on the boundary)
		int i_entry = 0;
		Kinesin *motor;
//		do{
			i_entry = properties_->gsl.GetRanInt(n_stalled_untethered_);
			motor = stalled_untethered_list_[i_entry];
//		}while(BoundaryStatus(motor) == true);
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
	}
	else{
		printf("Error in Unbind_Untethered STALLED: no bound untethered motors!\n");
  //      exit(1);
    }
}

void KinesinManagement::KMC_Failstep_Untethered(){

    // Make sure that at least one bound motor exists
	UpdateStalledUntetheredList();
    if(n_stalled_untethered_ > 0){
        // Randomly pick a bound motor (that's not on the boundary)
		int i_entry = properties_->gsl.GetRanInt(n_stalled_untethered_);
		Kinesin *motor = stalled_untethered_list_[i_entry];
		// Update site details
		motor->rear_site_->motor_ = nullptr;
		motor->rear_site_->occupied_ = false;
		// Update motor details
		motor->heads_active_--;
		motor->rear_site_ = nullptr;
		// Update statistics
		n_bound_untethered_--; 
		n_pseudo_bound_++;
		
	}
	else{
		printf("Error in Failstep_UT: no stalled untethered motors!\n");
		printf("%i motors\n", n_stalled_untethered_);
  //      exit(1);
    }
}

void KinesinManagement::KMC_Unbind_Tethered(){

	UpdateBoundTetheredList();
	if(n_bound_tethered_tot_ > 0){
		int i_entry = 0;
		Kinesin *motor;
//		do{	
			i_entry = properties_->gsl.GetRanInt(n_bound_tethered_tot_);
			motor = bound_tethered_list_[i_entry];
//		}while(BoundaryStatus(motor) == true);
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
//		printf("motor that was unbound had ext of %i\n", x_dub_pre);
		motor->UpdateExtension();
		n_bound_tethered_[x_dub_pre]--;
		n_bound_tethered_tot_--;
		// Update sites for prc1 management
		AssociatedProtein *xlink = motor->xlink_;
		AssociatedProteinManagement *prc1 = &properties_->prc1;
		if(xlink->heads_active_ == 1){
			prc1->n_sites_i_tethered_[x_dub_pre]--;
			// Treat xlinks tethered to free motors as untethered xlinks
			prc1->n_sites_i_untethered_++;
		}
		else{
			int x_dist = xlink->x_dist_;
			prc1->n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
			// Treat xlinks tethered to free motors as untethered xlinks
			prc1->n_sites_ii_untethered_[x_dist] += 2;
		}
	}
	else{
		printf("Error in Unbind_Tethered: no bound tethered motors!\n");
//		exit(1);
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
			int x_dist_dub = motor->x_dist_doubled_;
			// Update sites for prc1 management
			AssociatedProtein *xlink = motor->xlink_;
			AssociatedProteinManagement *prc1 = &properties_->prc1;	
			if(xlink->heads_active_ == 1){
				prc1->n_sites_i_tethered_[x_dist_dub]--;
				prc1->n_sites_i_untethered_++;
			}
			else{
				int x_dist = xlink->x_dist_;
				prc1->n_sites_ii_tethered_[x_dist_dub][x_dist] -= 2;
				prc1->n_sites_ii_untethered_[x_dist] += 2;
			}
		}
		n_pseudo_bound_--;
	}
	else{
		printf("Error in Unbind:Pseudo: no pseudo bound motors!\n");
//		exit(1);
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
		AssociatedProteinManagement *prc1 = &properties_->prc1;
		prc1->n_tethered_++;
		prc1->n_untethered_--;
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
//					printf("failed to tether_bound\n");
					failure = true;
					break;
				}
				int x_dub_pre = x_dist_dub;
				while(x_dub_pre == x_dist_dub){
					x_dist_dub = motor->SampleTailExtensionDoubled();
				}
//				printf("2x: %i\n", x_dist_dub);
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
			AssociatedProteinManagement *prc1 = &properties_->prc1;
			prc1->n_tethered_++;
			prc1->n_untethered_--;
			if(xlink->heads_active_ == 1){
				prc1->n_sites_i_untethered_--;
				prc1->n_sites_i_tethered_[x_dist_dub]++;
			}
			else if(xlink->heads_active_ == 2){
				int x_dist = xlink->x_dist_;
				xlink->UpdateExtension();
				int x_dist_post = xlink->x_dist_;
				if(x_dist != x_dist_post){
					printf("error in tetherbound 2....\n");
					exit(1);
				}
				prc1->n_sites_ii_untethered_[x_dist] -= 2;
				prc1->n_sites_ii_tethered_[x_dist_dub][x_dist] += 2;
			}
			else{
				printf("error in tether_bound...\n");
				exit(1);
			}
		}
	}
	else{
		printf("Error in Tether_Bound: no bound untethered motors!\n");
//		exit(1);
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
//		exit(1);
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
		// Update sites for prc1_management
		AssociatedProteinManagement *prc1 = &properties_->prc1;
		if(xlink->heads_active_ == 1){
			prc1->n_sites_i_tethered_[x_dub_pre]--;
			prc1->n_sites_i_untethered_++;
		}
		else if(xlink->heads_active_ ==2){
			int x_dist = xlink->x_dist_;
			prc1->n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
			prc1->n_sites_ii_untethered_[x_dist] += 2;
		}
		else{
			printf("god damnit in kmc_untether_bound\n");
			exit(1);
		}
	}
	else{
		printf("Error in Untether_Bound: no bound tethered motors!\n");
//		exit(1);
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
//		exit(1);
	}
}

void KinesinManagement::KMC_Step_ToTethRest(int x_dist_doubled){

	UpdateStepableTetheredTables();
	int n_stepable = n_stepable_to_teth_rest_[x_dist_doubled]; 
	if(n_stepable > 0){
		int i_entry = properties_->gsl.GetRanInt(n_stepable);
		Kinesin *motor = stepable_to_rest_table_[x_dist_doubled][i_entry];
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
		int dx_rest = motor->GetDirectionTowardRest();
		if(dx != dx_rest){
			printf("wat in step to rest motor KMC\n");
			exit(1);
		}
//		printf("coord is %g, ext is %i, dx is %i\n", 
//				motor->GetStalkCoordinate(), x_dist_doubled, dx);
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
		// Update site statistics for prc1
		AssociatedProtein *xlink = motor->xlink_;
		AssociatedProteinManagement *prc1 = &properties_->prc1;
		if(xlink->heads_active_ == 1){
			prc1->n_sites_i_tethered_[x_dub_pre]--;
			prc1->n_sites_i_tethered_[x_dub_post]++;
		}
		else if(xlink->heads_active_ == 2){
			int x_dist = xlink->x_dist_;
			prc1->n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
			prc1->n_sites_ii_tethered_[x_dub_post][x_dist] += 2;
		}
		else{
			printf("nope. tethered stepe motor.\n");
			exit(1);
		}
	}
	else{
		printf("Error in Step_Tethered (not rly tho): no stepable motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Step_FromTethRest(int x_dist_doubled){

	UpdateStepableTetheredTables();
	int n_stepable = n_stepable_from_teth_rest_[x_dist_doubled]; 
	if(n_stepable > 0){
		int i_entry = properties_->gsl.GetRanInt(n_stepable);
		Kinesin *motor = stepable_from_rest_table_[x_dist_doubled][i_entry];
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
		int dx_rest = motor->GetDirectionTowardRest();
		if(dx != -1 * dx_rest){
			printf("wat in step to rest motor KMC\n");
			exit(1);
		}
//		printf("coord is %g, ext is %i, dx is %i\n", 
//				motor->GetStalkCoordinate(), x_dist_doubled, dx);
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
		// Update site statistics for prc1
		AssociatedProtein *xlink = motor->xlink_;
		AssociatedProteinManagement *prc1 = &properties_->prc1;
		if(xlink->heads_active_ == 1){
			prc1->n_sites_i_tethered_[x_dub_pre]--;
			prc1->n_sites_i_tethered_[x_dub_post]++;
		}
		else if(xlink->heads_active_ == 2){
			int x_dist = xlink->x_dist_;
			prc1->n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
			prc1->n_sites_ii_tethered_[x_dub_post][x_dist] += 2;
		}
		else{
			printf("nope. tethered stepe motor.\n");
			exit(1);
		}
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
