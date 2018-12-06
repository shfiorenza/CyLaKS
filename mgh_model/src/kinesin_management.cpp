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
	InitializeSerialPop();
	InitializeFunctionMap();
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

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    double delta_t = parameters_->delta_t;
	double site_size = parameters_->microtubules.site_size;
	/* Statistics for diffusion */
	double D_const = parameters_->motors.diffusion_const;
	double x_squared = (site_size/1000)*(site_size/1000); // convert to um^2
	tau_ = x_squared / (2 * D_const);
	p_diffuse_forward_ = delta_t / tau_;
	p_diffuse_backward_ = delta_t / tau_; 
	// Generate stepping rates based on extension of tether: rates are 
	// increased if stepping towards rest length, and reduced if stepping
	// away from rest length (which increases tether extension)
	rest_dist_ = motors_[0].rest_dist_;
	comp_cutoff_ = motors_[0].comp_cutoff_;
	dist_cutoff_ = motors_[0].dist_cutoff_;
	if(world_rank == 0){
		printf("For motors:\n");
		printf("  rest_dist is %g\n", rest_dist_);
		printf("  comp_cutoff is %i\n", comp_cutoff_);
		printf("  dist_cutoff is %i\n", dist_cutoff_);
	}
	p_diffuse_to_teth_rest_.resize(2*dist_cutoff_ + 1);
	p_diffuse_from_teth_rest_.resize(2*dist_cutoff_ + 1);
	double kbT = parameters_->kbT;
	double r_0 = motors_[0].r_0_;
	double k_spring = motors_[0].k_spring_;
	double k_eff_slack = motors_[0].k_slack_;
	double r_y = parameters_->microtubules.y_dist / 2;
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
			dU_to = (0.5)*(k_spring*dr_toward*dr_toward - k_eff_slack*dr*dr);
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
		p_diffuse_from_teth_rest_[x_dist_dub] = p_from;
		p_diffuse_to_teth_rest_[x_dist_dub] = p_to;
	}
	/* Statistics for KMC */
    double k_on_i = parameters_->motors.k_on_i;
    double c_motor = parameters_->motors.concentration;
	p_bind_i_ = k_on_i * c_motor * delta_t;
	double c_eff_teth = parameters_->motors.conc_eff_tether;
	if(!parameters_->motors.tethers_active){
		c_eff_teth = 0;
	}
	p_bind_i_tethered_ = k_on_i * c_eff_teth * delta_t;
	double k_on_ii = parameters_->motors.k_on_ii;
	double c_eff_motor_bind = parameters_->motors.conc_eff_bind;
	p_bind_ii_ = k_on_ii * c_eff_motor_bind * delta_t;
	if(world_rank == 0){
		if(p_bind_i_ > 1){
			printf("WARNING: p_bind_i=%g for motors\n", p_bind_i_);
		}
		if(p_bind_ii_ > 1){
			printf("WARNING: p_bind_ii=%g for motors\n", p_bind_ii_);
		}
	}
	double k_off_i = parameters_->motors.k_off_i; 
	p_unbind_i_ = k_off_i * delta_t;
    double k_off = parameters_->motors.k_off_ii;
	p_unbind_ii_ = k_off * delta_t;
	p_bind_ii_to_teth_.resize(2*dist_cutoff_ + 1);
	p_bind_ii_from_teth_.resize(2*dist_cutoff_ + 1);
	p_unbind_i_tethered_.resize(2*dist_cutoff_ + 1);
	p_unbind_ii_to_teth_.resize(2*dist_cutoff_ + 1);
	p_unbind_ii_from_teth_.resize(2*dist_cutoff_ + 1);
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
		double U_tot,
			   dU_to,
			   dU_from;
		if(x_dist_dub == 2*rest_dist_){
			U_tot = (k_eff_slack/2)*dr*dr;
			dU_from = (k_eff_slack/2)*(dr_from*dr_from - dr*dr);
			dU_to = (0.5)*(k_spring*dr_toward*dr_toward - k_eff_slack*dr*dr);
        }
		// use k_spring if extension is positive
		else if(dr > 0){
			U_tot = (k_spring/2)*dr*dr;
			dU_to = (k_spring/2)*(dr_toward*dr_toward - dr*dr);
			dU_from = (k_spring/2)*(dr_from*dr_from - dr*dr);
		}
		// otherwise, use k_eff to model 'slack' in the tether
		else{
			U_tot = (k_eff_slack/2)*dr*dr;
			dU_to = (k_eff_slack/2)*(dr_toward*dr_toward - dr*dr);
			dU_from	= (k_eff_slack/2)*(dr_from*dr_from - dr*dr);
		}
		double weight_at = exp(-U_tot/(2*kbT));
		double weight_to = exp(-dU_to/(2*kbT));
		double weight_from = exp(-dU_from/(2*kbT));
		if(x_dist_dub >= 2*dist_cutoff_
		|| x_dist_dub <= 2*comp_cutoff_){
			weight_from = 0;
		}
		if(x_dist_dub < 2*comp_cutoff_
		|| x_dist_dub > 2*dist_cutoff_){
			weight_to = 0;
		}
		if(!parameters_->motors.tethers_active){
			weight_at = 0;
			weight_to = 0;
			weight_from = 0;
		}
		p_bind_ii_to_teth_[x_dist_dub] = weight_to * p_bind_ii_; 
		p_bind_ii_from_teth_[x_dist_dub] = weight_from * p_bind_ii_; 
		p_unbind_i_tethered_[x_dist_dub] = weight_at * p_unbind_i_;
		p_unbind_ii_to_teth_[x_dist_dub] = weight_to * p_unbind_ii_;
		p_unbind_ii_from_teth_[x_dist_dub] = weight_from * p_unbind_ii_; 
		if(world_rank == 0){
			if(p_bind_ii_to_teth_[x_dist_dub] > 1){
				printf("WARNING: p_bind_ii_to_teth=%g for 2x=%i\n", 
						p_bind_ii_to_teth_[x_dist_dub], x_dist_dub);
			}
			if(p_bind_ii_from_teth_[x_dist_dub] > 1){
				printf("WARNING: p_bind_ii_from_teth=%g for 2x=%i\n", 
						p_bind_ii_from_teth_[x_dist_dub], x_dist_dub);
			}
			if(p_unbind_i_tethered_[x_dist_dub] > 1){
				printf("WARNING: p_unbind_i_tethered=%g for 2x=%i\n", 
						p_unbind_i_tethered_[x_dist_dub], x_dist_dub);
			}
			if(p_unbind_ii_to_teth_[x_dist_dub] > 1){
				printf("WARNING: p_unbind_ii_to_teth=%g for 2x=%i\n", 
						p_unbind_ii_to_teth_[x_dist_dub], x_dist_dub);
			}
			if(p_unbind_ii_from_teth_[x_dist_dub] > 1){
				printf("WARNING: p_unbind_ii_from_teth=%g for 2x=%i\n", 
						p_unbind_ii_from_teth_[x_dist_dub], x_dist_dub);
			}
		}
	}
	double k_tether_free = parameters_->motors.k_tether_free;
	if(!parameters_->motors.tethers_active){
		k_tether_free = 0;
	}
	p_tether_free_ = k_tether_free * c_motor * delta_t;
	p_tether_bound_ = k_tether_free * c_eff_teth * delta_t;
	double k_untether_free = parameters_->motors.k_untether_free;
	if(!parameters_->motors.tethers_active){
		k_untether_free = 0;
	}
	p_untether_free_ = k_untether_free * delta_t;
    double motor_speed = parameters_->motors.velocity;
	p_step_ = motor_speed * delta_t / site_size;
	if(p_step_ > 1
	&& world_rank == 0){
		printf("WARNING: p_step=%g for motors\n", p_step_);
	}
	// Generate untethering and stepping rates for all tether extensions	
	// Everything is 2*dist_cutoff to allow for half-integer distances, 
	// so the 3rd entry will correspond to a distance of 1.5, etc. 
	double k_unteth = parameters_->motors.k_untether;
	if(!parameters_->motors.tethers_active){
		k_unteth = 0;
	}
	double stall_force = parameters_->motors.stall_force;
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
			double p_to = coeff_to * p_step_;
			double p_from = coeff_from * p_step_; 
			if(!parameters_->motors.tethers_active){
				p_to = 0;
				p_from = 0;
			}
			if(world_rank == 0){
				if(p_to > 1){
					printf("WARNING: p_step_to=%g for 2x=%i\n", 
							p_to, x_dist_dub);
				}
				if(p_from > 1){
					printf("WARNING: p_step_from=%g for 2x=%i\n", 
							p_from, x_dist_dub);
				}
			}
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
			double p_to = coeff_to * p_step_;
			double p_from = coeff_from * p_step_;
			if(!parameters_->motors.tethers_active){
				p_to = 0;
				p_from = 0;
			}
			if(world_rank == 0){
				if(p_to > 1){
					printf("WARNING: p_step_to=%g for 2x=%i\n", 
							p_to, x_dist_dub);
				}
				if(p_from > 1){
					printf("WARNING: p_step_from=%g for 2x=%i\n", 
							p_from, x_dist_dub);
				}
			}
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
}

void KinesinManagement::InitiateLists(){

	// One dimensional stuff
	free_tethered_list_.resize(n_motors_);
	bound_i_list_.resize(n_motors_);
	bound_i_bindable_list_.resize(n_motors_); 
	bound_ii_list_.resize(n_motors_);
	bound_untethered_.resize(n_motors_);
	stepable_list_.resize(n_motors_);
	bound_ii_tethered_list_.resize(n_motors_);
	stepable_list_.resize(n_motors_);
	// Two dimensional stuff
	n_bindable_to_teth_.resize(2*dist_cutoff_ + 1); 
	n_bindable_from_teth_.resize(2*dist_cutoff_ + 1);
	n_bound_i_tethered_.resize(2*dist_cutoff_ + 1);
	n_bound_ii_tethered_.resize(2*dist_cutoff_ + 1);
	n_bound_tethered_.resize(2*dist_cutoff_ + 1);
	n_stepable_to_teth_rest_.resize(2*dist_cutoff_ + 1);
	n_stepable_from_teth_rest_.resize(2*dist_cutoff_ + 1);
	bindable_to_teth_.resize(2*dist_cutoff_ + 1);
	bindable_from_teth_.resize(2*dist_cutoff_ + 1);
	bound_i_tethered_.resize(2*dist_cutoff_ + 1);
	bound_ii_tethered_table_.resize(2*dist_cutoff_ + 1);
	bound_tethered_.resize(2*dist_cutoff_ + 1);
	stepable_to_rest_table_.resize(2*dist_cutoff_ + 1);
	stepable_from_rest_table_.resize(2*dist_cutoff_ + 1);
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		bindable_to_teth_[x_dist_dub].resize(n_motors_);
		bindable_from_teth_[x_dist_dub].resize(n_motors_);
		bound_i_tethered_[x_dist_dub].resize(n_motors_);
		bound_ii_tethered_table_[x_dist_dub].resize(n_motors_);
		bound_tethered_[x_dist_dub].resize(n_motors_); 
		stepable_to_rest_table_[x_dist_dub].resize(n_motors_);
		stepable_from_rest_table_[x_dist_dub].resize(n_motors_);
		n_bindable_to_teth_[x_dist_dub] = 0;
		n_bindable_from_teth_[x_dist_dub] = 0;
		n_bound_i_tethered_[x_dist_dub] = 0;
		n_bound_ii_tethered_[x_dist_dub] = 0;
		n_bound_tethered_[x_dist_dub] = 0;
		n_stepable_to_teth_rest_[x_dist_dub] = 0;
		n_stepable_from_teth_rest_[x_dist_dub] = 0;
	}
}

void KinesinManagement::InitializeSerialPop(){

	int tot_size = 16 + 16*dist_cutoff_;
	if(!parameters_->motors.tethers_active)
		tot_size = 5;
	serial_pop_.resize(tot_size);

	pop_t bind_i = {0, "bind_i", -1, -1};
	serial_pop_[0] = bind_i;

	// If tethers are active, serialize all for loops
	if(parameters_->motors.tethers_active){
		pop_t bind_i_tethered = {0, "bind_i_tethered", -1, -1};
		serial_pop_[1] = bind_i_tethered;
		pop_t bind_ii = {0, "bind_ii", -1, -1};
		serial_pop_[2] = bind_ii;
		for(int x_dub(0); x_dub <= 2*dist_cutoff_; x_dub++){
			pop_t bind_ii_to_teth = {0, "bind_ii_to_teth", -1, x_dub};
			pop_t bind_ii_fr_teth = {0, "bind_ii_fr_teth", -1, x_dub};
			serial_pop_[3 + x_dub] = bind_ii_to_teth;
			serial_pop_[4 + 2*dist_cutoff_ + x_dub] = bind_ii_fr_teth; 
		}
		int offset = 4 + 4*dist_cutoff_;

		pop_t unbind_i = {0, "unbind_i", -1, -1};
		serial_pop_[offset + 1] = unbind_i;
		for(int x_dub(0); x_dub <= 2*dist_cutoff_; x_dub++){
			pop_t unbind_i_tethered = {0, "unbind_i_tethered", -1, x_dub};
			serial_pop_[offset + 2 + x_dub] = unbind_i_tethered;
		}
		int offset2 = 6 * 6*dist_cutoff_; 

		pop_t unbind_ii = {0, "unbind_ii", -1, -1};
		serial_pop_[offset2 + 1] = unbind_ii;
		for(int x_dub(0); x_dub <= 2*dist_cutoff_; x_dub++){
			pop_t unbind_to_teth = {0, "unbind_ii_to_teth", -1, x_dub};
			pop_t unbind_fr_teth = {0, "unbind_ii_fr_teth", -1, x_dub};
			serial_pop_[offset2 + 2 + x_dub] = unbind_to_teth;
			serial_pop_[offset2 + 3 + 2*dist_cutoff_ + x_dub] = unbind_fr_teth;
		}
		int offset3 = 9 + 10*dist_cutoff_;

		pop_t tether_free = {0, "tether_free", -1, -1};
		serial_pop_[offset3 + 1] = tether_free;
		pop_t tether_bound = {0, "tether_bound", -1, -1};
		serial_pop_[offset3 + 2] = tether_bound;
		pop_t untether_free = {0, "untether_free", -1, -1};
		serial_pop_[offset3 + 3] = untether_free;
		for(int x_dub(0); x_dub <= 2*dist_cutoff_; x_dub++){
			pop_t untether_bound = {0, "untether_bound", -1, x_dub};
			serial_pop_[offset3 + 4 + x_dub] = untether_bound; 
		}
		int offset4 = 13 + 12*dist_cutoff_; 

		pop_t step = {0, "step", -1, -1};
		serial_pop_[offset4 + 1] = step;
		for(int x_dub(0); x_dub <= 2*dist_cutoff_; x_dub++){
			pop_t step_to_teth = {0, "step_to_teth", -1, x_dub};
			pop_t step_fr_teth = {0, "step_fr_teth", -1, x_dub};
			serial_pop_[offset4 + 2 + x_dub] = step_to_teth;
			serial_pop_[offset4 + 3 + 2*dist_cutoff_ + x_dub] = step_fr_teth; 
		}
	}
	// If tethering is disabled, populations are much simpler
	else{
		pop_t bind_ii = {0, "bind_ii", -1, -1};
		serial_pop_[1] = bind_ii;

		pop_t unbind_i = {0, "unbind_i", -1, -1};
		serial_pop_[2] = unbind_i;

		pop_t unbind_ii = {0, "unbind_ii", -1, -1};
		serial_pop_[3] = unbind_ii;

		pop_t step = {0, "step", -1, -1};
		serial_pop_[4] = step;
	}
	// Copy into serial_kmc_ vector
	serial_kmc_ = serial_pop_;
}

void KinesinManagement::InitializeFunctionMap(){

	// Use lambda expressions to create functions that sample statistical
	// distributions with the appropriate p & n values for each population,
	// then bind them to a string key via a std::make_pair and store in map

	// Function that gets num to bind_i
	auto bind_i = [&](int x_dub){
		if(properties_->microtubules.n_unoccupied_ > 0){
			return properties_->gsl.SampleBinomialDist(p_bind_i_, 
					properties_->microtubules.n_unoccupied_, 0);
		}
		else return 0;
	};

	// If tethering isn't enabled, population types are simple
	if(!parameters_->motors.tethers_active){
		sampling_functs.insert(std::make_pair("bind_i", bind_i));
		// Function that gets num to bind_ii
		auto bind_ii = [&](int x_dub){
			if(n_bound_i_bindable_ > 0){
				return properties_->gsl.SampleBinomialDist(p_bind_ii_, 
						n_bound_i_bindable_, 1);
			}
			else return 0;
		};
		sampling_functs.insert(std::make_pair("bind_ii", bind_ii));
		// Function that gets num to unbind_i
		auto unbind_i = [&](int x_dub){
			if(n_bound_i_ > 0){
				return properties_->gsl.SampleBinomialDist(p_unbind_i_, 
						n_bound_i_, 2);
			}
			else return 0;
		};
		sampling_functs.insert(std::make_pair("unbind_i", unbind_i));
		// Function that gets num to unbind_ii
		auto unbind_ii = [&](int x_dub){
			if(n_bound_ii_ > 0){
				return properties_->gsl.SampleBinomialDist(p_unbind_ii_,
						n_bound_ii_, 3);
			}
			else return 0;
		};
		sampling_functs.insert(std::make_pair("unbind_ii", unbind_ii));
		// Function that gets num to step
		auto step = [&](int x_dub){
			if(n_stepable_ > 0){
			return properties_->gsl.SampleBinomialDist(p_step_, 
					n_stepable_, 4);
			}
			else return 0;
		};
		sampling_functs.insert(std::make_pair("step", step));
	}
	// Otherwise, serialize all populations of different tether extension
	else{
		// Function that gets num to bind_i_tethered
		auto bind_i_tethered = [&](int x_dub){
			return 0;
		};
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
    if(motors_[ID].front_site_ != nullptr
            || motors_[ID].rear_site_ != nullptr){
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
    if(motors_[ID].front_site_ == nullptr
            || motors_[ID].rear_site_ == nullptr){
        printf("Error: motor #%i is out of sync with motor_list\n", ID);
        exit(1);
    }
}

int KinesinManagement::GetNumBoundUntethered(){

	int n_untethered = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin* motor = &motors_[i_motor];
		if(motor->heads_active_ > 0
		&& motor->tethered_ == false){
			n_untethered++;
		}
	}
	return n_untethered;
}

Kinesin* KinesinManagement::GetFreeMotor(){

	// Randomly pick a motor from the reservoir
	int i_motor = properties_->gsl.GetRanInt(n_motors_);
	Kinesin *motor = &motors_[i_motor];
	int attempts = 0;
	while(motor->heads_active_ > 0 
	|| motor->tethered_ == true){
		i_motor++;
		if(i_motor == n_motors_)
			i_motor = 0;
		motor = &motors_[i_motor];
		attempts++;
		if(attempts > n_motors_){
			printf("error in get free motor\n");
			exit(1);
		}
	}
	UnboundCheck(motor);
	return motor;
}

Kinesin* KinesinManagement::GetBoundUntetheredMotor(){

	Kinesin* untethered_list[n_motors_];
	int i_entry = 0;
	int n_untethered = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ > 0
		&& motor->tethered_ == false){
			untethered_list[i_entry] = motor;
			i_entry++;
			n_untethered++; 
		}
	}
	int i_motor = properties_->gsl.GetRanInt(n_untethered); 
	Kinesin *motor = untethered_list[i_motor];
	return motor;
}

void KinesinManagement::UpdateAllLists(){

	#pragma omp parallel
	{
	#pragma omp single
	{
		#pragma omp task
		properties_->microtubules.UpdateUnoccupiedList();
		#pragma omp task
		UpdateBoundIList();
		#pragma omp task
		UpdateBoundIBindableList();
		#pragma omp task
		UpdateBoundIIList();
		#pragma omp task
		UpdateStepableList();
		if(parameters_->motors.tethers_active){
			#pragma omp task
			properties_->prc1.UpdateUntetheredList();
			#pragma omp task
			UpdateFreeTetheredList();
			#pragma omp task
			UpdateBoundUntethered();
			#pragma omp task
			UpdateBindableToTeth();
			#pragma omp task
			UpdateBindableFromTeth();
			#pragma omp task
			UpdateBoundITethered();
			#pragma omp task
			UpdateBoundIITetheredList();
			#pragma omp task
			UpdateBoundIITetheredTable();
			#pragma omp task
			UpdateBoundTethered();
			#pragma omp task
			UpdateStepableTetheredTables();
		}
	#pragma omp taskwait
	}
	}
}

void KinesinManagement::UpdateFreeTetheredList(){

	int i_entry = 0;
	int n_entries = 0; 
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
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

void KinesinManagement::UpdateBoundIList(){
	
	int i_entry = 0;
	int n_entries = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ == 1){
			if(motor->tethered_ == false){
				bound_i_list_[i_entry] = motor;
				i_entry++;
				n_entries++;
			}
			else if(motor->xlink_->heads_active_ == 0){
				bound_i_list_[i_entry] = motor;
				i_entry++;
				n_entries++;
			}
		}
	}
	if(n_entries != n_bound_i_){
		printf("something wrong in update_pseudo_bound (motors)\n");
		properties_->wallace.PrintMicrotubules(); 
		printf("local: %i, stats: %i\n", n_entries, n_bound_i_);
		exit(1);
	}
}

void KinesinManagement::UpdateBoundIBindableList(){

	n_bound_i_bindable_ = 0; 
	int i_entry = 0; 
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor]; 
		if(motor->heads_active_ == 1){
			if(motor->tethered_ == false){
				Microtubule *mt = motor->mt_;
				Tubulin *bound_site = motor->GetActiveHeadSite();
				int i_site = bound_site->index_;
				int i_plus = mt->plus_end_;
				int i_minus = mt->minus_end_;
				int dx = mt->delta_x_;
				// Don't access lattice sites that don't exist
				if(i_site == i_plus){
					if(mt->lattice_[i_site - dx].occupied_ == false){
						bound_i_bindable_list_[i_entry] = motor;
						i_entry++;
						n_bound_i_bindable_++;
					}
				}
				else if(i_site == i_minus){
					if(mt->lattice_[i_site + dx].occupied_ == false){
						bound_i_bindable_list_[i_entry] = motor;
						i_entry++;
						n_bound_i_bindable_++;
					}
				}
				else{
					// Add motor if it has an unoccupied site to either side
					if(mt->lattice_[i_site + dx].occupied_ == false
					|| mt->lattice_[i_site - dx].occupied_ == false){
						bound_i_bindable_list_[i_entry] = motor;
						i_entry++;
						n_bound_i_bindable_++; 
					}
				}
			}
			else if(motor->xlink_->heads_active_ == 0){
				Microtubule *mt = motor->mt_;
				Tubulin *bound_site = motor->GetActiveHeadSite();
				int i_site = bound_site->index_;
				int i_plus = mt->plus_end_;
				int i_minus = mt->minus_end_;
				int dx = mt->delta_x_;
				// Don't access lattice sites that don't exist
				if(i_site == i_plus){
					if(mt->lattice_[i_site - dx].occupied_ == false){
						bound_i_bindable_list_[i_entry] = motor;
						i_entry++;
						n_bound_i_bindable_++;
					}
				}
				else if(i_site == i_minus){
					if(mt->lattice_[i_site + dx].occupied_ == false){
						bound_i_bindable_list_[i_entry] = motor;
						i_entry++;
						n_bound_i_bindable_++;
					}
				}
				else{
					// Add motor if it has an unoccupied site to either side
					if(mt->lattice_[i_site + dx].occupied_ == false
					|| mt->lattice_[i_site - dx].occupied_ == false){
						bound_i_bindable_list_[i_entry] = motor;
						i_entry++;
						n_bound_i_bindable_++; 
					}
				}
			}
		}
	}
}

void KinesinManagement::UpdateBoundIIList(){

	int i_entry = 0;
	int n_entries = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ == 2){
			if(motor->tethered_ == false){
				bound_ii_list_[i_entry] = motor;
				i_entry++;
				n_entries++;
			}
			else if(motor->xlink_->heads_active_ == 0){
				bound_ii_list_[i_entry] = motor;
				i_entry++;
				n_entries++;
			}
		}
	}
	if(n_entries != n_bound_ii_){
		printf("something wrong in update_unteth_list (motors 1D)");
		printf(" %i in statistics, %i entries tho\n", n_bound_ii_, 
				n_entries);
		exit(1);
	}
}

void KinesinManagement::UpdateBoundUntethered(){

	n_bound_untethered_ = 0;
	int i_entry = 0;
	for(int i_motor(0); i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->tethered_ == false
		&& motor->heads_active_ > 0){
			bound_untethered_[i_entry] = motor;
			i_entry++;
			n_bound_untethered_++; 
		}
	}
}

void KinesinManagement::UpdateStepableList(){

	n_stepable_ = 0; 
	int i_entry = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
        if(motor->heads_active_ == 2){
			if(motor->tethered_ == false){
				Microtubule *mt = motor->mt_;
				int plus_end = mt->plus_end_;
				int dx = mt->delta_x_;
				int i_front_site = motor->front_site_->index_;
				// Exclude plus_end (can't step off of MTs)
				if(i_front_site != plus_end){
					if(mt->lattice_[i_front_site + dx].occupied_ == false){
						stepable_list_[i_entry] = motor;
						i_entry++; 
						n_stepable_++;
					}
				}
			}
			else if(motor->xlink_->heads_active_ == 0){
				Microtubule *mt = motor->mt_;
				int plus_end = mt->plus_end_;
				int dx = mt->delta_x_;
				int i_front_site = motor->front_site_->index_;
				// Exclude plus_end (can't step off of MTs)
				if(i_front_site != plus_end){
					if(mt->lattice_[i_front_site + dx].occupied_ == false){
						stepable_list_[i_entry] = motor;
						i_entry++; 
						n_stepable_++;
					}
				}
			}
		}
	}
}

void KinesinManagement::UpdateBindableToTeth(){

	int i_entry[2*dist_cutoff_ + 1];
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		i_entry[x_dist_dub] = 0;
		n_bindable_to_teth_[x_dist_dub] = 0;
	}
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ == 1
		&& motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				motor->UpdateExtension();
				if(motor->tethered_ == true){
					Microtubule *mt = motor->mt_;
					int i_plus = mt->plus_end_;
					int i_minus = mt->minus_end_;
					int dx = mt->delta_x_; 
					Tubulin *site = motor->GetActiveHeadSite();
					int i_site = site->index_; 
					int rest_dx = motor->GetDirectionTowardRest(); 
					// no seg faults plz
					if(!(i_site == i_plus && rest_dx == dx)
					&& !(i_site == i_minus && rest_dx == -dx)){
						if(mt->lattice_[i_site+rest_dx].occupied_ == false){
							int x_dist_dub = motor->x_dist_doubled_; 
							int index = i_entry[x_dist_dub]; 
							bindable_to_teth_[x_dist_dub][index] = motor;
							i_entry[x_dist_dub]++;
							n_bindable_to_teth_[x_dist_dub]++;
						}
					}
				}
			}
		}
	}
}

void KinesinManagement::UpdateBindableFromTeth(){

	int i_entry[2*dist_cutoff_ + 1];
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		i_entry[x_dist_dub] = 0;
		n_bindable_from_teth_[x_dist_dub] = 0;
	}
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ == 1
		&& motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				motor->UpdateExtension();
				if(motor->tethered_ == true){
					Microtubule *mt = motor->mt_;
					int i_plus = mt->plus_end_;
					int i_minus = mt->minus_end_;
					int dx = mt->delta_x_; 
					Tubulin *site = motor->GetActiveHeadSite();
					int i_site = site->index_; 
					int rest_dx = motor->GetDirectionTowardRest(); 
					// no seg faults plz
					if(!(i_site == i_plus && rest_dx == -dx)
					&& !(i_site == i_minus && rest_dx == dx)){
						if(mt->lattice_[i_site-rest_dx].occupied_ == false){
							int x_dist_dub = motor->x_dist_doubled_; 
							int index = i_entry[x_dist_dub]; 
							bindable_from_teth_[x_dist_dub][index] = motor;
							i_entry[x_dist_dub]++;
							n_bindable_from_teth_[x_dist_dub]++;
						}
					}
				}
			}
		}
	}
}

void KinesinManagement::UpdateBoundITethered(){

	int i_entry[2*dist_cutoff_ + 1];
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		i_entry[x_dist_dub] = 0;
		n_bound_i_tethered_[x_dist_dub] = 0;
	}
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ == 1
		&& motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				motor->UpdateExtension();
				if(motor->tethered_ == true){
					int x_dist_dub = motor->x_dist_doubled_; 
					int index = i_entry[x_dist_dub]; 
					bound_i_tethered_[x_dist_dub][index] = motor;
					i_entry[x_dist_dub]++;
					n_bound_i_tethered_[x_dist_dub]++;
				}
			}
		}
	}
}

void KinesinManagement::UpdateBoundIITetheredList(){

	int i_entry = 0;
	int n_entries = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ == 2
		&& motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				bound_ii_tethered_list_[i_entry] = motor;
				i_entry++;
				n_entries++;
			}
		}
	}
	if(n_entries != n_bound_ii_tethered_tot_){
		printf("something wrong in update_teth_list (motors 1D)");
		printf(" %in statistics, %i entries tho\n", 
				n_bound_ii_tethered_tot_, n_entries);
		exit(1);
	}
}

void KinesinManagement::UpdateBoundIITetheredTable(){

	int i_entry[2*dist_cutoff_ + 1];
	int n_entries[2*dist_cutoff_ + 1]; 
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		i_entry[x_dist_dub] = 0;
		n_entries[x_dist_dub] = 0;
	}
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ == 2
		&& motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				motor->UpdateExtension();
				if(motor->tethered_ == true){
					int x_dist_dub = motor->x_dist_doubled_; 
					int index = i_entry[x_dist_dub]; 
					bound_ii_tethered_table_[x_dist_dub][index] = motor;
					i_entry[x_dist_dub]++;
					n_entries[x_dist_dub]++;
				}
			}
		}
	}
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		if(n_entries[x_dist_dub] != n_bound_ii_tethered_[x_dist_dub]){
			printf("something wrong in update_bound_tethered (motors)");
			printf("for ext %i, \n%i in stats but %i entries counted\n", 
					x_dist_dub, n_bound_ii_tethered_[x_dist_dub], 
					n_entries[x_dist_dub]);
			exit(1);
		}
	}
}

void KinesinManagement::UpdateBoundTethered(){

	int i_entry[2*dist_cutoff_ + 1];
	for(int i_ext(0); i_ext <= 2*dist_cutoff_; i_ext++){
		n_bound_tethered_[i_ext] = 0;
		i_entry[i_ext] = 0;
	}
	for(int i_motor(0); i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				motor->UpdateExtension(); 
				if(motor->tethered_ == true){
					int x_dub = motor->x_dist_doubled_; 
					int index = i_entry[x_dub];
					bound_tethered_[x_dub][index] = motor;
					i_entry[x_dub]++;
					n_bound_tethered_[x_dub]++; 
				}
			}
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
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ == 2
		&& motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				motor->UpdateExtension();
				Microtubule *mt = motor->mt_;
				int plus_end = mt->plus_end_;
				int delta_x = mt->delta_x_;
				int dx_rest = motor->GetDirectionTowardRest();
				int i_front = motor->front_site_->index_;
				// end pausing
				if(i_front != plus_end
				&& motor->tethered_ == true){
					if(mt->lattice_[i_front + delta_x].occupied_ == false
					&& motor->AtCutoff() == false){
						int x_dist_dub = motor->x_dist_doubled_;
						// if MT's dx is towards rest, add to to_rest list
						if(delta_x == dx_rest){
							int index = i_entry_to[x_dist_dub];
							stepable_to_rest_table_[x_dist_dub][index] 
								= motor;	
							i_entry_to[x_dist_dub]++;
							n_stepable_to_teth_rest_[x_dist_dub]++;
						}
						// otherwise, add to from_rest list
						else if(delta_x == -1 * dx_rest){
							int index = i_entry_from[x_dist_dub];
							stepable_from_rest_table_[x_dist_dub][index] 
								= motor;
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
}

void KinesinManagement::GenerateDiffusionList(){

	int n_events = 0; 
	// Untethered statistics
	UpdateBoundIIList();
	int n_fwd_unteth = GetNumToStepForward();
	int n_bck_unteth = GetNumToStepBackward();
	while(n_fwd_unteth + n_bck_unteth > n_bound_ii_){
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
	UpdateBoundIITetheredTable();
	int n_toward_rest[2*dist_cutoff_ + 1];
	int n_from_rest[2*dist_cutoff_ + 1];
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		n_toward_rest[x_dist_dub] = GetNumToStepToTethRest(x_dist_dub);
		n_from_rest[x_dist_dub] = GetNumToStepFromTethRest(x_dist_dub);
		int n_to = n_toward_rest[x_dist_dub];
		int n_from = n_from_rest[x_dist_dub];
		int n_tethered = n_bound_ii_tethered_[x_dist_dub];
		while(n_to + n_from > n_tethered){
			double ran = properties_->gsl.GetRanProb();
			double p_to = p_diffuse_to_teth_rest_[x_dist_dub];
			double p_from = p_diffuse_from_teth_rest_[x_dist_dub];
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

int KinesinManagement::GetNumToStepForward(){

	int n_bound = n_bound_ii_;
	double p_step = p_diffuse_forward_;
	if(n_bound > 0){
		int n_to_step = 
			properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToStepBackward(){

	int n_bound = n_bound_ii_;
	double p_step = p_diffuse_backward_;
	if(n_bound > 0){
		int n_to_step = 
			properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToStepToTethRest(int x_dist_doubled){

	int n_bound = n_bound_ii_tethered_[x_dist_doubled];
	double p_step = p_diffuse_to_teth_rest_[x_dist_doubled];
	if(n_bound > 0){
		int n_to_step = 
			properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToStepFromTethRest(int x_dist_doubled){

	int n_bound = n_bound_ii_tethered_[x_dist_doubled];
	double p_step = p_diffuse_from_teth_rest_[x_dist_doubled];
	if(n_bound > 0){
		int n_to_step = 
			properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
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
					RunDiffusion_Forward();
					break;	
				case 20:
//					printf("unteth step bck\n");
					RunDiffusion_Backward();
					break;
				case 30:
//					printf("teth step to (%i)[%i avail]\n", x_dist_dub, 
//							n_bound_ii_tethered_[x_dist_dub]);
					RunDiffusion_To_Teth_Rest(x_dist_dub);
					break;
				case 40:
//					printf("teth step from (%i)[%i avail]\n", x_dist_dub, 
//							n_bound_ii_tethered_[x_dist_dub]);
					RunDiffusion_From_Teth_Rest(x_dist_dub);
					break;
			}
		}
	}
}

void KinesinManagement::RunDiffusion_Forward(){

	UpdateBoundIIList();
	int n_bound = n_bound_ii_;
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_ii_list_[i_entry];
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

void KinesinManagement::RunDiffusion_Backward(){

	UpdateBoundIIList();
	int n_bound = n_bound_ii_;
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_ii_list_[i_entry];
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

void KinesinManagement::RunDiffusion_To_Teth_Rest(int x_dist_doubled){

	int mt_length = parameters_->microtubules.length;
	int mt_array_length = mt_length - 1;	
	UpdateBoundIITetheredTable();
	int n_bound = n_bound_ii_tethered_[x_dist_doubled];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_ii_tethered_table_[x_dist_doubled][i_entry];
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
					n_bound_ii_tethered_[x_dub_pre]--;
					n_bound_ii_tethered_[x_dub_post]++;
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

void KinesinManagement::RunDiffusion_From_Teth_Rest(int x_dist_doubled){

	int mt_length = parameters_->microtubules.length;
	int mt_array_length = mt_length - 1;	
	UpdateBoundIITetheredTable();
	int n_bound = n_bound_ii_tethered_[x_dist_doubled];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_ii_tethered_table_[x_dist_doubled][i_entry];
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
					n_bound_ii_tethered_[x_dub_pre]--;
					n_bound_ii_tethered_[x_dub_post]++;
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

	bool tethers_on = parameters_->motors.tethers_active;
	UpdateSerializedPopulations();
	UpdateSerializedEvents();

	int n_events(0);

	int n_bind_i,
		n_bind_ii, 
		n_unbind_i, 
		n_unbind_ii, 
		n_step;
	
	int n_bind_i_tethered, 
		n_bind_ii_to_teth[2*dist_cutoff_ + 1],
		n_bind_ii_fr_teth[2*dist_cutoff_ + 1],
		n_unbind_i_tethered[2*dist_cutoff_ + 1],
		n_unbind_ii_to_teth[2*dist_cutoff_ + 1],
		n_unbind_ii_fr_teth[2*dist_cutoff_ + 1],
		n_tether_free, 
		n_tether_bound,
		n_untether_free,
		n_untether_bound[2*dist_cutoff_ + 1], 
		n_step_to_teth[2*dist_cutoff_ + 1], 
		n_step_fr_teth[2*dist_cutoff_ + 1];

	n_bind_i = serial_kmc_[0].n_entries_; 
	n_events += n_bind_i;

	if(tethers_on){
		/*
		int n_bind_i_teth = GetNumToBind_I_Tethered();
		n_events += n_bind_i_teth;
		UpdateBoundIBindableList();
		int n_bind_ii = GetNumToBind_II();
		n_events += n_bind_ii;
		UpdateStepableList();
		int n_unbind= GetNumToUnbind_II();
		n_events += n_unbind;
		int n_step = GetNumToStep(); 
		n_events += n_step;
		while(n_unbind + n_step > n_stepable_){
			if(n_step > 0){
				n_step--;
				n_events--;
			}
			else if(n_unbind > 0){
				n_unbind--;
				n_events--;
			}
		}
		UpdateBoundIITetheredList();
		UpdateBoundIList();
		int n_unbind_i = GetNumToUnbind_I();
		n_events += n_unbind_i;
		while(n_unbind_i + n_bind_ii > n_bound_i_){
			double p_bind = p_bind_ii_; 
			double p_unbind = p_unbind_i_;
			double p_tot = p_bind + p_unbind;
			double ran = properties_->gsl.GetRanProb();
			if(n_bind_ii > 0
			&& ran < p_bind / p_tot){
				n_bind_ii--;
				n_events--;
			}
			else if(n_unbind_i > 0){
				n_unbind_i--;
				n_events--;
			}
			else if(n_bind_ii > 0){
				n_bind_ii--;
				n_events--;
			}
		}
		int n_tether_free = GetNumToTether_Free();
		n_events += n_tether_free;
		UpdateBoundIIList();
		UpdateBoundUntethered();
		int n_tether_bound = GetNumToTether_Bound();
		n_events += n_tether_bound;
		int n_untether_free = GetNumToUntether_Free();
		n_events += n_untether_free;
		while(n_untether_free + n_bind_i_teth > n_free_tethered_){
			double p_unteth = p_untether_free_;
			double p_bind = p_bind_i_tethered_;
			double p_tot = p_unteth + p_bind; 
			double ran = properties_->gsl.GetRanProb();
			if(n_bind_i_teth > 0
					&& ran < p_unteth / p_tot){
				n_bind_i_teth--;
				n_events--;
			}
			else if(n_untether_free > 0){
				n_untether_free--;
				n_events--;
			}
			else if(n_bind_i_teth > 0){
				n_bind_i_teth--;
				n_events--;
			}
		}
		// Handle the statistics of differnt tether extensions separately 
		UpdateBoundITethered();
		UpdateBindableToTeth();
		UpdateBindableFromTeth();
		UpdateStepableTetheredTables();
		UpdateBoundIITetheredTable();
		UpdateBoundTethered();
		int n_bind_to_teth[2*dist_cutoff_ + 1];
		int n_bind_from_teth[2*dist_cutoff_ + 1]; 
		int n_unbind_tethered[2*dist_cutoff_ + 1];
		int n_unbind_to_teth[2*dist_cutoff_ + 1];
		int n_unbind_from_teth[2*dist_cutoff_ + 1];

		int n_untether_bound[2*dist_cutoff_ + 1];
		int n_step_to_rest[2*dist_cutoff_ + 1];
		int n_step_from_rest[2*dist_cutoff_ + 1];
		for(int x_dub = 0; x_dub <= 2*dist_cutoff_; x_dub++){

			n_bind_to_teth[x_dub] = GetNumToBind_II_To_Teth(x_dub);
			n_bind_from_teth[x_dub] = GetNumToBind_II_From_Teth(x_dub);
			n_unbind_tethered[x_dub] = GetNumToUnbind_I_Tethered(x_dub);
			n_unbind_to_teth[x_dub] = GetNumToUnbind_II_To_Teth(x_dub);
			n_unbind_from_teth[x_dub] = GetNumToUnbind_II_From_Teth(x_dub);

			n_untether_bound[x_dub] = GetNumToUntether_Bound(x_dub);
			n_step_to_rest[x_dub] = GetNumToStep_ToTethRest(x_dub);
			n_step_from_rest[x_dub] = GetNumToStep_FromTethRest(x_dub);

			int n_bind_to = n_bind_to_teth[x_dub];
			int n_bind_from = n_bind_from_teth[x_dub];
			int n_unbind_teth = n_unbind_tethered[x_dub];
			int n_bound_i_teth = n_bound_i_tethered_[x_dub]; 
			while(n_bind_to + n_bind_from + n_unbind_teth > n_bound_i_teth){
				double p_to = p_bind_ii_to_teth_[x_dub];
				double p_from = p_bind_ii_from_teth_[x_dub];
				double p_unbind = p_unbind_i_tethered_[x_dub];
				double p_tot = p_to + p_from + p_unbind; 
				double ran = properties_->gsl.GetRanProb();
				if(ran < p_to / p_tot
						&& n_bind_to > 0){
					n_bind_to_teth[x_dub]--;
					n_bind_to = n_bind_to_teth[x_dub]; 
				}
				else if(ran < p_from / p_tot
						&& n_bind_from > 0){
					n_bind_from_teth[x_dub]--;
					n_bind_from = n_bind_from_teth[x_dub];
				}
				else{
					n_unbind_tethered[x_dub]--;
					n_unbind_teth = n_unbind_tethered[x_dub];
				}

			}

			int n_unteth = n_untether_bound[x_dub];

			int n_unbind_to = n_unbind_to_teth[x_dub];
			int n_unbind_from = n_unbind_from_teth[x_dub];

			int n_step_to = n_step_to_rest[x_dub];
			int n_step_from = n_step_from_rest[x_dub];
			int n_tethered = n_bound_ii_tethered_[x_dub];
			// Prevent too many steps from occuring (imagine if 2 motors
			// are bound and we roll for 1 unbind but 2 steps)
			while(n_unteth + n_step_to + n_step_from +
			n_unbind_to + n_unbind_from > n_tethered){
				double p_to = p_step_to_teth_rest_[x_dub];
				double p_from = p_step_from_teth_rest_[x_dub];
				double p_tot = p_to + p_from;
				double ran = properties_->gsl.GetRanProb();
				if(ran < p_to / p_tot
						&& n_step_to > 0){
					n_step_to_rest[x_dub]--;
					n_step_to = n_step_to_rest[x_dub];
				}
				else if(n_step_from > 0){
					n_step_from_rest[x_dub]--;
					n_step_from = n_step_from_rest[x_dub];
				}
				else if(n_step_to > 0){
					n_step_to_rest[x_dub]--;
					n_step_to = n_step_to_rest[x_dub];
				}
				if(n_unteth + n_step_to + n_step_from + 
						n_unbind_to + n_unbind_from < n_tethered){
					break;
				}
				double p_unbind_to = p_unbind_ii_to_teth_[x_dub];
				double p_unbind_from = p_unbind_ii_from_teth_[x_dub];
				double p_unteth = p_untether_bound_[x_dub];
				double p_tot2 = p_unbind_to + p_unbind_from + p_unteth;
				double ran2 = properties_->gsl.GetRanProb();
				if(ran2 < p_unbind_to / p_tot2
						&& n_unbind_to > 0){
					n_unbind_to_teth[x_dub]--;
					n_unbind_to = n_unbind_to_teth[x_dub];
				}
				else if(ran2 < p_unbind_from / p_tot2
						&& n_unbind_from > 0){
					n_unbind_from_teth[x_dub]--;
					n_unbind_from = n_unbind_from_teth[x_dub];
				}
				else if (n_unteth > 0){
					n_untether_bound[x_dub]--;
					n_unteth = n_untether_bound[x_dub];
				}
			} 

			n_events += n_bind_to_teth[x_dub];
			n_events += n_bind_from_teth[x_dub];
			n_events += n_unbind_tethered[x_dub];
			n_events += n_unbind_to_teth[x_dub];
			n_events += n_unbind_from_teth[x_dub];

			n_events += n_untether_bound[x_dub];
			n_events += n_step_to_rest[x_dub];
			n_events += n_step_from_rest[x_dub];
		}
		*/
	}
	else{
		n_bind_ii = serial_kmc_[1].n_entries_;
		n_unbind_i = serial_kmc_[2].n_entries_;
		// Ensure there aren't more KMC events than singly-bound motors
		while(n_unbind_i + n_bind_ii > n_bound_i_){
			double p_tot = p_bind_ii_ + p_unbind_i_;
			double ran = properties_->gsl.GetRanProb();
			if(ran < p_bind_ii_/p_tot
			&& n_bind_ii > 0)
				n_bind_ii--;
			else if(n_unbind_i > 0)
				n_unbind_i--;
			else if(n_bind_ii > 0)
				n_bind_ii--;
		}

		n_unbind_ii = serial_kmc_[3].n_entries_;
		n_step = serial_kmc_[4].n_entries_;
		// Ensure there aren't more KMC events than doubly-bound motors
		while(n_unbind_ii + n_step > n_stepable_){
			double p_tot = p_step_ + p_unbind_ii_;
			double ran = properties_->gsl.GetRanProb();
			if(ran < p_step_/p_tot
			&& n_step > 0)
				n_step--;
			else if(n_unbind_ii > 0)
				n_unbind_ii--;
			else if(n_step > 0)
				n_bind_ii--;
		}
		
		n_events += n_bind_ii;
		n_events += n_unbind_i;
		n_events += n_unbind_ii;
		n_events += n_step; 
	}
	if(n_events > 0){
//		printf("n_events: %i \n", n_events);
		int pre_list[n_events];
		int kmc_index = 0;
		if(tethers_on){
		/*
		  	for(int i_event = 0; i_event < n_bind_i; i_event++){
//				printf("10\n");
				pre_list[kmc_index] = 10;
				kmc_index++;
			}
			for(int i_event = 0; i_event < n_bind_i_teth; i_event++){
//				printf("11\n");
				pre_list[kmc_index] = 11;
				kmc_index++;
			}
			for(int i_event = 0; i_event < n_bind_ii; i_event++){
//				printf("12\n");
				pre_list[kmc_index] = 12;
				kmc_index++;
			}
			for(int i_event = 0; i_event < n_unbind_i; i_event++){
//				printf("20\n");
				pre_list[kmc_index] = 20;
				kmc_index++;
			}
			for(int i_event = 0; i_event < n_unbind_ii; i_event++){
//				printf("22\n");
				pre_list[kmc_index] = 22;
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
			for(int i_event = 0; i_event < n_untether_free; i_event++){
				//			printf("50\n");
				pre_list[kmc_index] = 50;
				kmc_index++;
			}
			for(int i_event = 0; i_event < n_step; i_event++){
//				printf("60\n");
				pre_list[kmc_index] = 60;
				kmc_index++;
			}
			for(int x_dub = 0; x_dub <= 2*dist_cutoff_; x_dub++){
				int n_bind_to = n_bind_to_teth[x_dist_dub];
				for(int i_event = 0; i_event < n_bind_to; i_event++){
//					printf("%i\n", 1300 + x_dist_dub);
					pre_list[kmc_index] = 1300 + x_dist_dub; 
					kmc_index++;
				}
				int n_bind_from = n_bind_from_teth[x_dist_dub];
				for(int i_event = 0; i_event < n_bind_from; i_event++){
//					printf("%i\n", 1400 + x_dist_dub);
					pre_list[kmc_index] = 1400 + x_dist_dub; 
					kmc_index++;
				}
				int n_unbind_teth = n_unbind_tethered[x_dist_dub];
				for(int i_event = 0; i_event < n_unbind_teth; i_event++){
					//				printf("%i\n", 2100 + x_dist_dub);
					pre_list[kmc_index] = 2100 + x_dist_dub; 
					kmc_index++;
				}
				int n_unbind_to = n_unbind_to_teth[x_dist_dub];
				for(int i_event = 0; i_event < n_unbind_to; i_event++){
					//				printf("%i\n", 2300 + x_dist_dub);
					pre_list[kmc_index] = 2300 + x_dist_dub; 
					kmc_index++;
				}
				int n_unbind_from = n_unbind_from_teth[x_dist_dub];
				for(int i_event = 0; i_event < n_unbind_from; i_event++){
					//				printf("%i\n", 2400 + x_dist_dub);
					pre_list[kmc_index] = 2400 + x_dist_dub; 
					kmc_index++;
				}
				int n_unteth_b = n_untether_bound[x_dist_dub];
				for(int i_event = 0; i_event < n_unteth_b; i_event++){
					pre_list[kmc_index] = 500 + x_dist_dub; 	
					//				printf("%i\n", 500 + x_dist_dub);
					kmc_index++;
				}
				int n_step_to = n_step_to_rest[x_dist_dub];
				for(int i_event = 0; i_event < n_step_to; i_event++){
					pre_list[kmc_index] = 600 + x_dist_dub;
					//				printf("%i\n", 600 + x_dist_dub);
					kmc_index++;
				}
				int n_step_from = n_step_from_rest[x_dist_dub];
				for(int i_event = 0; i_event < n_step_from; i_event++){
					pre_list[kmc_index] = 700 + x_dist_dub;
					//				printf("%i\n", 700 + x_dist_dub);
					kmc_index++;
				}
			}
		*/
		}
		else{
			for(int i_event = 0; i_event < n_bind_i; i_event++){
//				printf("10\n");
				pre_list[kmc_index] = 10;
				kmc_index++;
			}
			for(int i_event = 0; i_event < n_bind_ii; i_event++){
//				printf("12\n");
				pre_list[kmc_index] = 12;
				kmc_index++;
			}
			for(int i_event = 0; i_event < n_unbind_i; i_event++){
//				printf("20\n");
				pre_list[kmc_index] = 20;
				kmc_index++;
			}
			for(int i_event = 0; i_event < n_unbind_ii; i_event++){
//				printf("22\n");
				pre_list[kmc_index] = 22;
				kmc_index++;
			}
			for(int i_event = 0; i_event < n_step; i_event++){
//				printf("60\n");
				pre_list[kmc_index] = 60;
				kmc_index++;
			}
		}
		if(n_events > 1){
			RandomNumberManagement *gsl = &properties_->gsl;
			gsl_ran_shuffle(gsl->rng, pre_list, n_events, sizeof(int));
		}
		// Transfer pre_list to class' kmc_list_ structure
        kmc_list_.resize(n_events);
        for(int i_event = 0; i_event < n_events; i_event++){
            kmc_list_[i_event] = pre_list[i_event];
        }
    }
    else{
        kmc_list_.clear();
    }
}

void KinesinManagement::UpdateSerializedPopulations(){

	bool tethers_on = parameters_->motors.tethers_active;

	// Update all pop. lists
	UpdateAllLists();
	#pragma omp parallel
	{
	#pragma omp single
	{
		// for bind_i
//		#pragma omp task
		serial_pop_[0].n_entries_ = properties_->microtubules.n_unoccupied_;

		if(tethers_on){
			/*
			// for bind_i_tethered
			if(n_free_tethered_ > 0){
				serial_pop_[1].n_entries_ = n_free_tethered_;
				active_serial_pop[n_active] = serial_pop_[1];
				n_active++;
			}

			// for bind_ii
			if(n_bound_i_bindable_ > 0){
				serial_pop_[2].n_entries_ = n_bound_i_bindable_; 
				active_serial_pop[n_active] = serial_pop_[2];
				n_active++;
			}

			// for bind_ii_to/fr_teth
			for(int x_dub(0); x_dub <= 2*dist_cutoff_; x_dub++){
				int n_to = n_bindable_to_teth_[x_dub];
				int n_fr = n_bindable_from_teth_[x_dub];
				if(n_to > 0){
					serial_pop_[3 + x_dub].n_entries_ = n_to; 
					active_serial_pop[n_active] = serial_pop_[3 + x_dub];
					n_active++;
				}
				if(n_fr > 0){
					serial_pop_[4+2*dist_cutoff_+x_dub].n_entries_ = n_fr; 
					active_serial_pop[n_active] = 
						serial_pop_[4 + 2 * dist_cutoff_ + x_dub];
					n_active++;
				}
			}

			int offset = 4 + 4*dist_cutoff_;

			// for unbind_i 
			if(n_bound_i_ > 0){
				serial_pop_[offset + 1].n_entries_ = n_bound_i_;
				active_serial_pop[n_active] = serial_pop_[offset + 1];
				n_active++;
			}

			pop_t unbind_ii = {n_bound_ii_, "unbind_ii", -1, -1};
			serial_pop_[offset + 2] = unbind_ii;

			for(int x_dub(0); x_dub <= 2*dist_cutoff_; x_dub++){
				int n = n_bound_ii_tethered_[x_dub];
				pop_t unbind_to_teth = {n, "unbind_ii_to_teth", -1, x_dub};
				pop_t unbind_fr_teth = {n, "unbind_ii_fr_teth", -1, x_dub};
				serial_pop_[3 + offset + x_dub] = unbind_to_teth;
				serial_pop_[4+2*dist_cutoff_+offset+x_dub] = unbind_fr_teth;
			}

			int offset2 = 8 + 8*dist_cutoff_; 

			int n_unteth = properties_->prc1.n_untethered_;
			pop_t tether_free = {n_unteth, "tether_free", -1, -1};
			serial_pop_[offset2 + 1] = tether_free;

			pop_t tether_bound = {0, "tether_bound", -1, -1};
			serial_pop_[offset2 + 2] = tether_bound;

			for(int x_dub(0); x_dub <= 2*dist_cutoff_; x_dub++){
				int n = n_bound_tethered_[x_dub];
				pop_t untether_bound = {n, "untether_bound", -1, x_dub};
				serial_pop_[offset2 + 3 + x_dub] = untether_bound; 
			}

			int offset3 = 11 + 10*dist_cutoff_; 

			pop_t untether_free = {n_free_tethered_,"untether_free", -1, -1};
			serial_pop_[offset3 + 1] = untether_free;

			pop_t to_step = {n_stepable_, "to_step", -1, -1};
			serial_pop_[offset3 + 2] = to_step;

			for(int x_dub(0); x_dub <= 2*dist_cutoff_; x_dub++){
				int n_to = n_stepable_to_teth_rest_[x_dub];
				int n_fr = n_stepable_from_teth_rest_[x_dub]; 
				pop_t step_to_teth = {n_to, "step_to_teth", -1, x_dub};
				pop_t step_fr_teth = {n_fr, "step_fr_teth", -1, x_dub};
				serial_pop_[3 + offset3 + x_dub] = step_to_teth;
				serial_pop_[4+2*dist_cutoff_+offset3+x_dub] = step_fr_teth; 
			}
			*/
		}
		else{
			// for bind_ii
//			#pragma omp task
			serial_pop_[1].n_entries_ = n_bound_i_bindable_;

			// for unbind_i
//			#pragma omp task
			serial_pop_[2].n_entries_ = n_bound_i_;

			// for unbind_ii
//			#pragma omp task
			serial_pop_[3].n_entries_ = n_bound_ii_;

			// for to_step
//			#pragma omp task
			serial_pop_[4].n_entries_ = n_stepable_;
		}
	}
	}
}

void KinesinManagement::UpdateSerializedEvents(){

	// Run through serialized population types
	#pragma omp parallel for
	for(int i_pop = 0; i_pop < serial_pop_.size(); i_pop++){
		// If population is greater than 0, sample for possible KMC event
		if(serial_pop_[i_pop].n_entries_ > 0){
			// Get iterator pointing to appropriate sampling function
			auto itr = sampling_functs.find(serial_pop_[i_pop].type_);
			// Ensure sampling function was properly found
			if(itr != sampling_functs.end()){
				// Function itself is 'second' in map's pair entry
				auto funct = itr->second;
				// Get number of KMC events expected for this population
				int n_events = funct(serial_pop_[i_pop].x_dist_dub_);
				serial_kmc_[i_pop].n_entries_ = n_events;
			}
			// If sampling function was not properly found, exit
			else{
				printf("error in cereal events\n");
				std::cout << serial_pop_[i_pop].type_ << std::endl;
				exit(1);
			}
		}
		// If population is 0, no. of kmc events is automatically 0
		else{
			serial_kmc_[i_pop].n_entries_ = 0;
		}
	}
}

int KinesinManagement::GetNumToBind_I(){

    int n_unocc = properties_->microtubules.n_unoccupied_;
	double p_bind = p_bind_i_; 
	if(n_unocc > 0){
		return properties_->gsl.SampleBinomialDist(p_bind, n_unocc, 0);
	}
	else{
		return 0;
	}
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
	double p_bind = p_bind_i_tethered_;
	double n_avg = p_bind * weights_summed;
	if(n_avg > 0){
		int n_to_bind = properties_->gsl.SamplePoissonDist(n_avg);
		return n_to_bind;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToBind_II(){

	int n_able = n_bound_i_bindable_;
	double p_bind = p_bind_ii_; 
	if(n_able > 0
	&& p_bind > 0){
		return properties_->gsl.SampleBinomialDist(p_bind, n_able, 1);
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToBind_II_To_Teth(int x_dist_doubled){

	int n_able = n_bindable_to_teth_[x_dist_doubled];
	double p_bind = p_bind_ii_to_teth_[x_dist_doubled];
	if(n_able > 0
	&& p_bind > 0){
		int n_to_bind = properties_->gsl.SampleBinomialDist(p_bind, n_able);
		return n_to_bind;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToBind_II_From_Teth(int x_dist_doubled){

	int n_able = n_bindable_from_teth_[x_dist_doubled];
	double p_bind = p_bind_ii_from_teth_[x_dist_doubled];
	if(n_able > 0
	&& p_bind > 0){
		int n_to_bind = properties_->gsl.SampleBinomialDist(p_bind, n_able);
		return n_to_bind;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToUnbind_I(){

	int n_bound = n_bound_i_; 
	double p_unbind = p_unbind_i_;
	if(n_bound > 0
	&& p_unbind > 0){
		return properties_->gsl.SampleBinomialDist(p_unbind, n_bound, 2);
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToUnbind_I_Tethered(int x_dist_doubled){
	
	// XXX what is this stat actually tracking?? where is it coming from?
	int n_bound = n_bound_i_tethered_[x_dist_doubled];
	double p_unbind = p_unbind_i_tethered_[x_dist_doubled];
	if(n_bound > 0
	&& p_unbind > 0)
		return properties_->gsl.SampleBinomialDist(p_unbind, n_bound);
	else
		return 0;
}

int KinesinManagement::GetNumToUnbind_II(){

	int n_bound = n_bound_ii_;
	double p_unbind = p_unbind_ii_;
	if(n_bound > 0
	&& p_unbind > 0){
		return properties_->gsl.SampleBinomialDist(p_unbind, n_bound, 3);
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToUnbind_II_To_Teth(int x_dist_doubled){

	int n_bound = n_bound_ii_tethered_[x_dist_doubled];
	double p_unbind = p_unbind_ii_to_teth_[x_dist_doubled];
	if(n_bound > 0
	&& p_unbind > 0)
		return properties_->gsl.SampleBinomialDist(p_unbind, n_bound);
	else
		return 0;
}

int KinesinManagement::GetNumToUnbind_II_From_Teth(int x_dist_doubled){

	int n_bound = n_bound_ii_tethered_[x_dist_doubled];
	double p_unbind = p_unbind_ii_from_teth_[x_dist_doubled];
	if(n_bound > 0
	&& p_unbind > 0)
		return properties_->gsl.SampleBinomialDist(p_unbind, n_bound);
	else
		return 0;
}

int KinesinManagement::GetNumToTether_Free(){

    int n_unteth = properties_->prc1.n_untethered_;
	// Calculate how many free motors tether within delta_t on avg
	double p_teth = p_tether_free_;
	if(n_unteth > 0
	&& p_teth > 0){
		int n_to_tether = 
			properties_->gsl.SampleBinomialDist(p_teth, n_unteth);
		return n_to_tether;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToTether_Bound(){

	double weights_summed = 0;
	for(int i_motor = 0; i_motor < n_bound_untethered_; i_motor++){
		Kinesin *motor = bound_untethered_[i_motor];
		motor->UpdateNeighborXlinks();
		int n_neighbs = motor->n_neighbor_xlinks_; 
		for(int i_neighb = 0; i_neighb < n_neighbs; i_neighb++){
			AssociatedProtein *xlink = motor->neighbor_xlinks_[i_neighb];
			double weight = motor->GetTetheringWeight(xlink);
			weights_summed += weight;
		}
	}
	double n_avg = p_tether_bound_ * weights_summed;
	if(n_avg > 0){
		int n_to_teth = properties_->gsl.SamplePoissonDist(n_avg);
		return n_to_teth;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToUntether_Bound(int x_dist_doubled){

	int n_teth = n_bound_tethered_[x_dist_doubled]; 
	double p_unteth = p_untether_bound_[x_dist_doubled];
	if(n_teth > 0
	&& p_unteth > 0){
		int n_to_unteth = 
			properties_->gsl.SampleBinomialDist(p_unteth, n_teth);
		return n_to_unteth;
	}
	else{
		return 0;
	}
}


int KinesinManagement::GetNumToUntether_Free(){

	int n_teth = n_free_tethered_;
	double p_unteth = p_untether_free_;
	if(n_teth > 0
	&& p_unteth > 0){
		int n_to_unteth = 
			properties_->gsl.SampleBinomialDist(p_unteth, n_teth);
		return n_to_unteth; 
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToStep(){

	int n_stepable = n_stepable_;
	double p_step = p_step_;
	if(n_stepable > 0
	&& p_step > 0)
		return properties_->gsl.SampleBinomialDist(p_step, n_stepable, 4);
	else
		return 0;
}

int KinesinManagement::GetNumToStep_ToTethRest(int x_dist_doubled){

	int n_stepable = n_stepable_to_teth_rest_[x_dist_doubled];
	double p_step = p_step_to_teth_rest_[x_dist_doubled];
	if(n_stepable > 0
	&& p_step > 0){
		int n_to_step = 
			properties_->gsl.SampleBinomialDist(p_step, n_stepable);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToStep_FromTethRest(int x_dist_doubled){

	int n_stepable = n_stepable_from_teth_rest_[x_dist_doubled];
	double p_step = p_step_from_teth_rest_[x_dist_doubled];
	if(n_stepable > 0
	&& p_step > 0){
		int n_to_step = 
			properties_->gsl.SampleBinomialDist(p_step, n_stepable);
		return n_to_step;
	}
	else{
		return 0;
	}
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

			if(kmc_event >= 1300 && kmc_event < 1400){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 13; 
			}
			if(kmc_event >= 1400 && kmc_event < 1500){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 14;
			}
			if(kmc_event >= 2100 && kmc_event < 2200){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 21;
			}
			if(kmc_event >= 2300 && kmc_event < 2400){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 23;
			}
			if(kmc_event >= 2400 && kmc_event < 2500){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 24;
			}
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

            switch(kmc_event){
				case 10:
//						printf("free motor pseudo-bound\n");
						KMC_Bind_I();
                        break;
                case 11:
//						printf("tethered motor pseudo-bound\n");
						KMC_Bind_I_Tethered();
                        break;
				case 12:
//						printf("pseudo motor bound\n");
						KMC_Bind_II();
						break;
				case 13:
//						printf("bind ii to teth\n");
						KMC_Bind_II_To_Teth(x_dist_doubled); 
						break;
				case 14: 
//						printf("bind ii from teth\n");
						KMC_Bind_II_From_Teth(x_dist_doubled);
						break;
				case 20:
//						printf("pseudo-bound motor unbound\n");
						KMC_Unbind_I();
						break;
				case 21:
//						printf("unbind i tethered\n");
						KMC_Unbind_I_Tethered(x_dist_doubled);
						break;
				case 22:
//						printf("untethered stepable motor unbound\n");
						KMC_Unbind_II();
                        break;
				case 23:
///						printf("unbind ii to teth\n");
						KMC_Unbind_II_To_Teth(x_dist_doubled);
						break;
				case 24:
//						printf("unbind ii from teth\n");
						KMC_Unbind_II_From_Teth(x_dist_doubled);
						break;
                case 30:
//						printf("free motor tethered\n");
						KMC_Tether_Free();
                        break;
				case 31:
//						printf("bound motor tethered\n");
						KMC_Tether_Bound();
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
						KMC_Step();
						break;
						
				case 62:
//						printf("tethered motor (ext %i) stepped\n", 
//								x_dist_doubled);
						KMC_Step_ToTethRest(x_dist_doubled);
						break;
				case 63:
//						printf("tethered motor (ext %i) stepped\n", 
//								x_dist_doubled);
						KMC_Step_FromTethRest(x_dist_doubled);
						break;
            }
        }
    }
}

void KinesinManagement::KMC_Bind_I(){

    // Make sure that at least one unbound motor exists
	properties_->microtubules.UpdateUnoccupiedList();
	if(properties_->microtubules.n_unoccupied_ > 0){	
        Kinesin *motor = GetFreeMotor();
		MicrotubuleManagement *mts = &properties_->microtubules;
		Tubulin *site = mts->GetUnoccupiedSite();
		Microtubule *mt = site->mt_;
		// Update site details
		site->motor_ = motor;
		site->occupied_ = true;
		// Update motor details
		motor->mt_ = mt;
		motor->front_site_ = site;
		motor->heads_active_ = 1;
		// Update statistics
		n_bound_i_++;
	}
	else{
		printf("Error in Bind_I: no unoccupied sites.\n");
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
		Tubulin *site = motor->GetWeightedNeighborSite();
		int attempts = 0;
		while(site == nullptr){
			if(attempts > 10*n_free_tethered_){
				break;
			}
			i_motor = properties_->gsl.GetRanInt(n_free_tethered_);
			motor = free_tethered_list_[i_motor];
			site = motor->GetWeightedNeighborSite();
			attempts++;	
		}
		if(site != nullptr){
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
//			n_bound_i_++; 
			// update sites for prc1 management
			motor->UpdateExtension();
			if(motor->tethered_ == true){
				if(motor->xlink_->heads_active_ > 0){
					AssociatedProtein *xlink = motor->xlink_;
					AssociatedProteinManagement *prc1 = &properties_->prc1;
					int x_dist_dub = motor->x_dist_doubled_;
					if(xlink->heads_active_ == 1){
						prc1->n_single_bound_--;
						prc1->n_sites_i_tethered_[x_dist_dub]++;	
						prc1->n_sites_i_untethered_--;
					}
					else if(xlink->heads_active_ == 2){
						int x_dist = xlink->x_dist_;
						prc1->n_double_bound_[x_dist]--;
						prc1->n_sites_ii_tethered_[x_dist_dub][x_dist] += 2;
						prc1->n_sites_ii_untethered_[x_dist] -= 2;
					}
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
		else{
			printf("Failed to bind_i_teth\n");
//			properties_->wallace.PrintMicrotubules(2);
		}
	}
	else if(properties_->microtubules.n_unoccupied_ == 0){
		printf("Error in Bind_I_Tethered: no unoccupied sites. ");
		//		properties_->wallace.PrintMicrotubules(0.5);
//		exit(1);
	}
	else{
		printf("Error in Bind_I_Tethered: no tethered free motors ");
//		properties_->wallace.PrintMicrotubules(0.5);
	}
}

void KinesinManagement::KMC_Bind_II(){

	UpdateBoundIBindableList();
	if(n_bound_i_bindable_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound_i_bindable_);
		Kinesin *motor = bound_i_bindable_list_[i_entry];
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
		// Choose where to bind second head
		if(i_site == i_plus){
			if(mt->lattice_[i_site - dx].occupied_ == false){
				front_site = bound_site;
				rear_site = &mt->lattice_[i_site - dx];	
			}
			else{
				printf("something wrong in bind i list A\n");
				exit(1);
			}
		}
		else if(i_site == i_minus){
			if(mt->lattice_[i_site + dx].occupied_ == false){
				front_site = &mt->lattice_[i_site + dx];
				rear_site = bound_site;
			}
			else{
				printf("something wrong in bind i list B\n");
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
			printf("Error with bind i bound motors list \n");
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
		n_bound_i_--;
		n_bound_i_bindable_--;
		n_bound_ii_++;
		// If motor is tethered, we gotta do a whole lot of bullshit
		if(motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				printf("umm wat\n");
				exit(1);
			}
		}
	}
	else{
//		printf("Error in Bind_II: no singly-bound motors. \n");
//		exit(1);
	}
}	

void KinesinManagement::KMC_Bind_II_To_Teth(int x_dist_doubled){

	UpdateBindableToTeth();
	int x_dub = x_dist_doubled; 
	int n_bound = n_bindable_to_teth_[x_dub];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bindable_to_teth_[x_dub][i_entry];
		Microtubule *mt = motor->mt_;
		Tubulin *bound_site = motor->GetActiveHeadSite();
		int i_site = bound_site->index_;
		int mt_dx = bound_site->mt_->delta_x_;
		int dx = motor->GetDirectionTowardRest();
		int i_plus = mt->plus_end_;
		int i_minus = mt->minus_end_;
		Tubulin *front_site,
				*rear_site; 
		if(motor->heads_active_ != 1){
			printf("nope in pseudos\n");
			exit(1);
		}
		// Choose where to bind second head
		if(i_site == i_plus){
			if(dx != mt_dx){
				if(mt->lattice_[i_site + dx].occupied_ == false){
					if(dx == mt_dx){
						rear_site = bound_site;
						front_site = &mt->lattice_[i_site + dx];	
					}
					else{
						front_site = bound_site;
						rear_site = &mt->lattice_[i_site + dx];	
					}
				}
			}
			else{
				printf("something wrong in bind i list ONE\n");
				exit(1);
			}
		}
		else if(i_site == i_minus){
			if(dx != -1*mt_dx){
				if(mt->lattice_[i_site + dx].occupied_ == false){
					if(dx == mt_dx){
						front_site = &mt->lattice_[i_site + dx];
						rear_site = bound_site;
					}
					else{
						rear_site = &mt->lattice_[i_site + dx];
						front_site = bound_site;
					}
				}
			}
			else{
				printf("something wrong in bind i list TWO\n");
				exit(1);
			}
		}
		else if(mt->lattice_[i_site + dx].occupied_ == false){
			if(dx == mt_dx){
				front_site = &mt->lattice_[i_site + dx];
				rear_site = bound_site;
			}
			else{
				rear_site = &mt->lattice_[i_site + dx];
				front_site = bound_site;
			}
		}
		else{
			printf("bruh. check bind_to_teth\n");
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
		// If motor is tethered, we gotta do a whole lot of bullshit
		if(motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				// Update statistics
				n_bound_i_tethered_[x_dub]--;
				n_bindable_to_teth_[x_dub]--;
				AssociatedProtein* xlink = motor->xlink_;
				int x_dub_pre = motor->x_dist_doubled_;
				if(x_dub_pre != x_dub){
					printf("pretty whack in bind ii to eth\n");
					exit(1);
				}
				motor->UpdateExtension();
				// weird exception if we triggered a force untether
				if(motor->tethered_ == false){
					// Cancel out the subtraction in the force untether since
					// this particular xlink never contributed to stats, 
					// the subtraction just fucks with stuff
					n_bound_ii_tethered_tot_++;
					n_bound_ii_tethered_[x_dub_pre]++;
				}
				// Otherwise, routine statistic stuff
				else{
					// KMC stuff
					int x_dist_dub = motor->x_dist_doubled_; 
					n_bound_ii_tethered_[x_dist_dub]++;
					n_bound_ii_tethered_tot_++;
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
			else{
				printf("why tho?? motor bind TO teth\n");
				exit(1);
			}
		}
		// Otherwise, error bruh
		else{
			printf("Error in Bind_II_TO_Teth\n");
			exit(1);
		}
	}
	else{
		//		printf("Error in Bind_II: no singly-bound motors. \n");
		//		exit(1);
	}
}

void KinesinManagement::KMC_Bind_II_From_Teth(int x_dist_doubled){

	UpdateBindableFromTeth();
	int x_dub = x_dist_doubled; 
	int n_bound = n_bindable_from_teth_[x_dub];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bindable_from_teth_[x_dub][i_entry];
		Microtubule *mt = motor->mt_;
		Tubulin *bound_site = motor->GetActiveHeadSite();
		int i_site = bound_site->index_;
		int mt_dx = bound_site->mt_->delta_x_;
		int dx = motor->GetDirectionTowardRest();
		int i_plus = mt->plus_end_;
		int i_minus = mt->minus_end_;
		Tubulin *front_site,
				*rear_site; 
		if(motor->heads_active_ != 1){
			printf("nope in pseudos\n");
			exit(1);
		}
		// Choose where to bind second head
		if(i_site == i_plus){
			if(dx != -1*mt_dx){
				if(mt->lattice_[i_site - dx].occupied_ == false){
					if(dx == mt_dx){
						front_site = bound_site;
						rear_site = &mt->lattice_[i_site - dx];	
					}
					else{
						rear_site = bound_site;
						front_site = &mt->lattice_[i_site - dx];	
					}
				}
			}
			else{
				printf("something wrong in bind i list THREE\n");
				exit(1);
			}
		}
		else if(i_site == i_minus){
			if(dx != mt_dx){
				if(mt->lattice_[i_site - dx].occupied_ == false){
					if(dx == mt_dx){
						front_site = bound_site;
						rear_site = &mt->lattice_[i_site - dx];
					}
					else{
						rear_site = bound_site;
						front_site = &mt->lattice_[i_site - dx];
					}
				}
			}
			else{
				printf("something wrong in bind i list FOUR\n");
				exit(1);
			}
		}
		else if(mt->lattice_[i_site - dx].occupied_ == false){
			if(dx == mt_dx){
				rear_site = &mt->lattice_[i_site - dx];
				front_site = bound_site;
			}
			else{
				front_site = &mt->lattice_[i_site - dx];
				rear_site = bound_site;
			}
		}
		else{
			printf("bruh. check bind_to_teth\n");
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
		// If motor is tethered, we gotta do a whole lot of bullshit
		if(motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				// Update statistics
				n_bound_i_tethered_[x_dub]--;
				n_bindable_to_teth_[x_dub]--;
				AssociatedProtein* xlink = motor->xlink_;
				int x_dub_pre = motor->x_dist_doubled_;
				if(x_dub_pre != x_dub){
					printf("pretty whack in bind ii to eth\n");
					exit(1);
				}
				motor->UpdateExtension();
				// weird exception if we triggered a force untether
				if(motor->tethered_ == false){
					// Cancel out the subtraction in the force untether since
					// this particular xlink never contributed to stats, 
					// the subtraction just fucks with stuff
					n_bound_ii_tethered_tot_++;
					n_bound_ii_tethered_[x_dub_pre]++;
				}
				// Otherwise, routine statistic stuff
				else{
					// KMC stuff
					int x_dist_dub = motor->x_dist_doubled_; 
					n_bound_ii_tethered_[x_dist_dub]++;
					n_bound_ii_tethered_tot_++;
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
			else{
				printf("why tho?? motor bind TO teth\n");
				exit(1);
			}
		}
		// Otherwise, error bruh
		else{
			printf("Error in Bind_II_TO_Teth\n");
			exit(1);
		}
	}
	else{
//		printf("Error in Bind_II: no singly-bound motors. \n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Unbind_I(){

	UpdateBoundIList();
	if(n_bound_i_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound_i_);
		Kinesin *motor = bound_i_list_[i_entry];
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
			if(motor->xlink_->heads_active_ == 0)
				motor->UntetherSatellite();
			else{
				printf("ummm whywouldudothis??? mot_unbind_i\n");
				exit(1);
			}
		}
		n_bound_i_--;
	}
	else{
		printf("Error in Unbind_I: no pseudo bound motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Unbind_I_Tethered(int x_dist_doubled){

	UpdateBoundITethered();
	int n_bound = n_bound_i_tethered_[x_dist_doubled]; 
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_i_tethered_[x_dist_doubled][i_entry];
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
			if(motor->xlink_->heads_active_ > 0){
				n_bound_i_tethered_[x_dist_doubled]--;
				n_free_tethered_++;
				int x_dist_dub = motor->x_dist_doubled_;
				if(x_dist_dub != x_dist_doubled){
					printf("\"not a fan\" -- unbind_i_teth\n");
					exit(1);
				}
				// Update sites for prc1 management
				AssociatedProtein *xlink = motor->xlink_;
				AssociatedProteinManagement *prc1 = &properties_->prc1;	
				if(xlink->heads_active_ == 1){
					prc1->n_single_bound_++;
					prc1->n_sites_i_tethered_[x_dist_dub]--;
					prc1->n_sites_i_untethered_++;
				}
				else{
					int x_dist = xlink->x_dist_;
					prc1->n_double_bound_[x_dist]++;
					prc1->n_sites_ii_tethered_[x_dist_dub][x_dist] -= 2;
					prc1->n_sites_ii_untethered_[x_dist] += 2;
				}
			}
			else{
				printf("ummm whywouldudothis??? mot_unbind_i_TETH TWOOO\n");
				exit(1);
			}
		}
		else{
			printf("ummm whywouldudothis??? mot_unbind_i_TETH\n");
			exit(1);
		}
	}
	else{
		printf("Error in Unbind_I: no pseudo bound motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Unbind_II(){

    // Make sure that at least one bound motor exists
	UpdateBoundIIList();
    if(n_bound_ii_ > 0){
        // Randomly pick a bound motor (that's not on the boundary)
		int i_entry = 0;
		i_entry = properties_->gsl.GetRanInt(n_bound_ii_);
		Kinesin *motor = bound_ii_list_[i_entry];
		// Roll and randomly pick a head to unbind
		double ran = properties_->gsl.GetRanProb();
		// Update site details
		if(ran < 0.5){
			motor->front_site_->motor_ = nullptr;
			motor->front_site_->occupied_ = false;
			motor->front_site_ = nullptr;
		}
		else{
			motor->rear_site_->motor_ = nullptr;
			motor->rear_site_->occupied_ = false;
			motor->rear_site_ = nullptr;
		}
		// Update motor details
		motor->heads_active_--;
		// Update statistics
		n_bound_ii_--; 
		n_bound_i_++; 
	}
	else{
		printf("Error in Unbind_Untethered MOBILE: ");
		printf("no bound untethered motors!\n");
//      exit(1);
    }
}

void KinesinManagement::KMC_Unbind_II_To_Teth(int x_dist_doubled){

	UpdateBoundIITetheredTable();
	int x_dub = x_dist_doubled;
	int n_bound = n_bound_ii_tethered_[x_dub];
    if(n_bound > 0){
        // Randomly pick a bound motor (that's not on the boundary)
		int i_entry = 0;
		i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_ii_tethered_table_[x_dub][i_entry];
		// Unbind site farther from rest to bring stalk closer to rest
		Tubulin *unbind_site = motor->GetSiteFartherFromRest();
		// Update site details
		if(unbind_site == motor->front_site_){
			motor->front_site_->motor_ = nullptr;
			motor->front_site_->occupied_ = false;
			motor->front_site_ = nullptr;
		}
		else{
			motor->rear_site_->motor_ = nullptr;
			motor->rear_site_->occupied_ = false;
			motor->rear_site_ = nullptr;
		}
		motor->heads_active_--;
		if(motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				// Update statistics
				AssociatedProtein* xlink = motor->xlink_;
				int x_dub_pre = motor->x_dist_doubled_;
				if(x_dub_pre != x_dub){
					printf("pretty whack in bind ii to eth\n");
					exit(1);
				}
				motor->UpdateExtension();
				// weird exception if we triggered a force untether
				if(motor->tethered_ == false){
					printf("HEY CHECK UNBIND II TO TETH IF ERROR THROWS\n");
					// Cancel out the subtraction in the force untether since
					// this particular xlink never contributed to stats, 
					// the subtraction just fucks with stuff
					//	XXX i dont think this is necessary
//					n_bound_ii_tethered_tot_++;
//					n_bound_ii_tethered_[x_dub_pre]++;
				}
				// Otherwise, routine statistic stuff
				else{
					// KMC stuff
					int x_dist_dub = motor->x_dist_doubled_; 
					n_bound_i_tethered_[x_dist_dub]++;
					n_bound_ii_tethered_[x_dub_pre]--;
					n_bound_ii_tethered_tot_--;
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
			else{
				printf("WUT in unbind_ii to teth TYPE 2\n");
				exit(1);
			}
		}
		else{
			printf("WUT in unbind_ii to teth\n");
			exit(1);
		}
	}
	else{
		printf("WUT in unbind_ii to teth: no teth bound motors\n");
		exit(1);
	}
}

void KinesinManagement::KMC_Unbind_II_From_Teth(int x_dist_doubled){

	UpdateBoundIITetheredTable();
	int x_dub = x_dist_doubled;
	int n_bound = n_bound_ii_tethered_[x_dub];
    if(n_bound > 0){
        // Randomly pick a bound motor (that's not on the boundary)
		int i_entry = 0;
		i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_ii_tethered_table_[x_dub][i_entry];
		// Unbind site closer to rest to bring stalk further from rest
		Tubulin *unbind_site = motor->GetSiteCloserToRest();
		// Update site details
		if(unbind_site == motor->front_site_){
			motor->front_site_->motor_ = nullptr;
			motor->front_site_->occupied_ = false;
			motor->front_site_ = nullptr;
		}
		else{
			motor->rear_site_->motor_ = nullptr;
			motor->rear_site_->occupied_ = false;
			motor->rear_site_ = nullptr;
		}
		motor->heads_active_--;
		if(motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				// Update statistics
				AssociatedProtein* xlink = motor->xlink_;
				int x_dub_pre = motor->x_dist_doubled_;
				if(x_dub_pre != x_dub){
					printf("pretty whack in bind ii to eth\n");
					exit(1);
				}
				motor->UpdateExtension();
				// weird exception if we triggered a force untether
				if(motor->tethered_ == false){
					printf("HEY CHECK UNBIND II TO TETH IF ERROR THROWS\n");
					// Cancel out the subtraction in the force untether since
					// this particular xlink never contributed to stats, 
					// the subtraction just fucks with stuff
					//	XXX i dont think this is necessary
//					n_bound_ii_tethered_tot_++;
//					n_bound_ii_tethered_[x_dub_pre]++;
				}
				// Otherwise, routine statistic stuff
				else{
					// KMC stuff
					int x_dist_dub = motor->x_dist_doubled_; 
					n_bound_i_tethered_[x_dist_dub]++;
					n_bound_ii_tethered_[x_dub_pre]--;
					n_bound_ii_tethered_tot_--;
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
			else{
				printf("WUT in unbind_ii to teth TYPE 2B\n");
				exit(1);
			}
		}
		else{
			printf("WUT in unbind_ii to teth\n");
			exit(1);
		}
	}
	else{
		printf("WUT in unbind_ii to teth: no teth bound motors\n");
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
		AssociatedProteinManagement *prc1 = &properties_->prc1;
		prc1->n_untethered_--;
	}
	else{
		printf("Error in Tether_Free: no bound untethered xlinks\n");
	}
}

void KinesinManagement::KMC_Tether_Bound(){

	UpdateBoundUntethered();
	if(n_bound_untethered_ > 0){
		int i_motor = properties_->gsl.GetRanInt(n_bound_untethered_);
		Kinesin *motor = bound_untethered_[i_motor];
		AssociatedProtein* xlink = nullptr;
		if(motor->tethered_ == false){
			xlink = motor->GetWeightedNeighborXlink();
		}
		int attempts = 0; 
		while(xlink == nullptr){
			if(attempts > 10*n_bound_untethered_){
				break;
			}
			i_motor = properties_->gsl.GetRanInt(n_bound_untethered_);
			motor = bound_untethered_[i_motor];
			if(motor->tethered_ == false){
				xlink = motor->GetWeightedNeighborXlink();
			}
			attempts++;
		}
		if(xlink != nullptr){
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
			// Update local statistics
			int x_dist_dub = motor->x_dist_doubled_; 
			n_bound_tethered_[x_dist_dub]++;
			n_bound_untethered_--; 
			if(motor->heads_active_ == 1){
				n_bound_i_--;
			}
			else if(motor->heads_active_ == 2){
				n_bound_ii_tethered_[x_dist_dub]++;
				n_bound_ii_tethered_tot_++;
				n_bound_ii_--;
			}
			else{
				printf("Error in tether bound motor >:(\n");
				exit(1);
			}
			// Update prc1 statistics 
			AssociatedProteinManagement *prc1 = &properties_->prc1;
			prc1->n_untethered_--;
			if(xlink->heads_active_ == 1){
				prc1->n_single_bound_--;
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
				prc1->n_double_bound_[x_dist]--;
				prc1->n_sites_ii_untethered_[x_dist] -= 2;
				prc1->n_sites_ii_tethered_[x_dist_dub][x_dist] += 2;
			}
			else{
				printf("error in tether_bound...\n");
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

void KinesinManagement::KMC_Untether_Bound(int x_dist_doubled){

	UpdateBoundTethered();
	int n_bound_tethered = n_bound_tethered_[x_dist_doubled];
	if(n_bound_tethered  > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound_tethered);
		Kinesin *motor = bound_tethered_[x_dist_doubled][i_entry];
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
		n_bound_untethered_--; 
		if(motor->heads_active_ == 2){
			n_bound_ii_tethered_[x_dub_pre]--;
			n_bound_ii_tethered_tot_--;
			n_bound_ii_++;
		}
		else{
			n_bound_i_++;
		}
		properties_->prc1.n_untethered_++;
		// Update sites for prc1_management
		AssociatedProteinManagement *prc1 = &properties_->prc1;
		if(xlink->heads_active_ == 1){
			prc1->n_single_bound_++;
			prc1->n_sites_i_tethered_[x_dub_pre]--;
			prc1->n_sites_i_untethered_++;
		}
		else if(xlink->heads_active_ ==2){
			int x_dist = xlink->x_dist_;
			prc1->n_double_bound_[x_dist]++;
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
		properties_->prc1.n_untethered_++;
	}
	else{
		printf("Error in Untether_Free: no free tethered motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Step(){

    // Make sure there is at least one stepable untethered motor
	UpdateStepableList();
	if(n_stepable_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_stepable_);
		Kinesin *motor = stepable_list_[i_entry];
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
//		printf("Error in Step_Untethered: no stepable untethered motors\n");
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
			exit(1);
		}
		n_bound_ii_tethered_[x_dub_pre]--;
		n_bound_ii_tethered_[x_dub_post]++;
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
//		printf("Error in Step_Tethered: no stepable motors!\n");
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
			exit(1);
		}
		n_bound_ii_tethered_[x_dub_pre]--;
		n_bound_ii_tethered_[x_dub_post]++;
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
//		printf("Error in Step_Tethered: no stepable motors!\n");
//		exit(1);
	}
}
