#include "master_header.h"
#include "associated_protein_management.h"

AssociatedProteinManagement::AssociatedProteinManagement(){
}

void AssociatedProteinManagement::Initialize(system_parameters *parameters, 
		system_properties *properties){

	parameters_ = parameters;
	properties_ = properties;

	GenerateXLinks();
	SetParameters();
	InitializeLists();
	InitializeDiffusionEvents();
	InitializeKMCEvents();
}

void AssociatedProteinManagement::GenerateXLinks(){

	int n_mts = parameters_->microtubules.count;
	n_xlinks_ = 0;
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		n_xlinks_ += parameters_->microtubules.length[i_mt];
	}
	// Since only one head has to be bound, the sim will at most
	// as many xlinks as sites in the bulk (all single-bound)
	xlinks_.resize(n_xlinks_);
	for(int ID = 0; ID < n_xlinks_; ID++){
		xlinks_[ID].Initialize(parameters_, properties_, ID);
	}
}

void AssociatedProteinManagement::SetParameters(){

	double delta_t = parameters_->delta_t;
	double site_size = parameters_->microtubules.site_size;
	// DIFFUSION STATISTICS FOR SELF BELOW
	double D_const_i = parameters_->xlinks.diffusion_const_i;
	double D_const_ii = parameters_->xlinks.diffusion_const_ii;
	double x_squared = (site_size/1000)*(site_size/1000); // in um^2
	double tau_i = x_squared / (2 * D_const_i);
	double tau_ii = x_squared / (2 * D_const_ii);
	p_diffuse_i_fwd_ = delta_t / tau_i;
	p_diffuse_i_bck_ = delta_t / tau_i;
	// Generate different stepping rates based on changes in
	// potential energy (dU) associated with that step
	teth_cutoff_ = properties_->kinesin4.dist_cutoff_; 
	dist_cutoff_ = xlinks_[0].dist_cutoff_;
	rest_dist_ = xlinks_[0].rest_dist_;
	if(parameters_->motors.tethers_active){
		properties_->wallace.Log("\nFor crosslinkers:\n");
		properties_->wallace.Log("  rest_dist is %i\n", rest_dist_);
		properties_->wallace.Log("  dist_cutoff is %i\n\n", dist_cutoff_);
	}
	p_diffuse_ii_to_rest_.resize(dist_cutoff_ + 1);
	p_diffuse_ii_from_rest_.resize(dist_cutoff_ + 1);
	double kbT = parameters_->kbT;
	double r_0 = xlinks_[0].r_0_;
	double k_spring = xlinks_[0].k_spring_;
	double r_y = parameters_->microtubules.y_dist;
	for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
		double r_x = x_dist * site_size;
		double r_x_to = (x_dist - 1) * site_size;
		double r_x_from = (x_dist + 1) * site_size;
		double r = sqrt(r_x*r_x + r_y*r_y);
		double r_to = sqrt(r_x_to*r_x_to + r_y*r_y);
		double r_from = sqrt(r_x_from*r_x_from + r_y*r_y);
		// Get extension for current dist and steps to/from spring rest
		double dr = r - r_0;
		double dr_to = r_to -  r_0;
		double dr_from = r_from - r_0;
		if(dr >= 0){
			// Get corresponding changes in potential energy
			double dU_to_rest = (k_spring/2)*(dr_to*dr_to - dr*dr);
			double dU_from_rest = (k_spring/2)*(dr_from*dr_from - dr*dr);
			// Weights according to Lanksy et al.
			double weight_to = exp(-dU_to_rest/(2*kbT));
			double weight_from = exp(-dU_from_rest/(2*kbT));
			double p_to = weight_to * delta_t / tau_ii;
			double p_from = weight_from * delta_t / tau_ii;
			if(x_dist == rest_dist_){
				p_diffuse_ii_to_rest_[x_dist] = 0;
				p_diffuse_ii_from_rest_[x_dist] = 2*p_from;
			}
			else if(x_dist == dist_cutoff_){
				p_diffuse_ii_to_rest_[x_dist] = p_to;
				p_diffuse_ii_from_rest_[x_dist] = 0;
			}
			else{
				p_diffuse_ii_to_rest_[x_dist] = p_to;
				p_diffuse_ii_from_rest_[x_dist] = p_from;
			}

			if(p_to > 1)
				printf("WARNING: p_diffuse_to_rest=%g for x=%i\n", 
						p_to, x_dist);
			if(2*p_from > 1)
				printf("WARNING: 2*p_diffuse_from_rest=%g for x=%i\n", 
						2*p_from, x_dist);
		}
		else{
			printf("woah mayne. xlink set parameters \n");
			exit(1);
		}
	}
	// DIFFUSION STATISTICS INVOLVING TETHER BELOW
	int teth_dist_cutoff = properties_->kinesin4.motors_[0].dist_cutoff_;
	int teth_comp_cutoff = properties_->kinesin4.motors_[0].comp_cutoff_;
	p_diffuse_i_to_teth_rest_.resize(2*teth_dist_cutoff + 1);
	p_diffuse_i_from_teth_rest_.resize(2*teth_dist_cutoff + 1);
	p_diffuse_ii_to_both_rest_.resize(2*teth_dist_cutoff + 1);
	p_diffuse_ii_to_self_from_teth_.resize(2*teth_dist_cutoff + 1);
	p_diffuse_ii_from_self_to_teth_.resize(2*teth_dist_cutoff + 1);
	p_diffuse_ii_from_both_rest_.resize(2*teth_dist_cutoff + 1);
	double k_teth_spring = properties_->kinesin4.motors_[0].k_spring_;
	double k_teth_slack = properties_->kinesin4.motors_[0].k_slack_;
	double r_0_teth = properties_->kinesin4.motors_[0].r_0_;
	double r_y_teth = parameters_->microtubules.y_dist / 2;
	double rest_dist_teth = properties_->kinesin4.motors_[0].rest_dist_;
	double r_rest_teth = 
		sqrt(site_size*rest_dist_teth*site_size*rest_dist_teth
		+ r_y_teth*r_y_teth);
	for(int x_dub = 0; x_dub <= 2*teth_dist_cutoff; x_dub++){
		p_diffuse_ii_to_both_rest_[x_dub].resize(dist_cutoff_ + 1);
		p_diffuse_ii_to_self_from_teth_[x_dub].resize(dist_cutoff_ + 1);
		p_diffuse_ii_from_self_to_teth_[x_dub].resize(dist_cutoff_ + 1);
		p_diffuse_ii_from_both_rest_[x_dub].resize(dist_cutoff_ + 1);
		// Calc x-distances (in nm) for tether
		double r_x_teth = x_dub * site_size / 2;
		double r_x_teth_bck = (x_dub - 1) * site_size / 2;
		double r_x_teth_fwd = (x_dub + 1) * site_size / 2;
		// Calc total r values 
		double r_teth = sqrt(r_x_teth*r_x_teth + r_y_teth*r_y_teth);
		double r_teth_bck = sqrt(r_x_teth_bck*r_x_teth_bck 
								 + r_y_teth*r_y_teth);
		double r_teth_fwd = sqrt(r_x_teth_fwd*r_x_teth_fwd 
								 + r_y_teth*r_y_teth);
		// Calc tether extensions for current dist and stepping to/from rest
		double dr_teth = r_teth - r_0_teth;	
		double dr_teth_to, 
			   dr_teth_from;
		if(dr_teth >= 0){
			dr_teth_to = r_teth_bck - r_0_teth;
			dr_teth_from = r_teth_fwd - r_0_teth;	
		}
		else{
			dr_teth_to = r_teth_fwd - r_0_teth;
			dr_teth_from = r_teth_bck - r_0_teth;
		}
		double dU_from_teth, 
			   dU_to_teth;
		if(x_dub == 2*rest_dist_teth){
			if(r_0_teth > r_rest_teth){
				dU_from_teth = (k_teth_slack/2)
					* (dr_teth_from*dr_teth_from - dr_teth*dr_teth);
				dU_to_teth = (0.5)*(k_teth_spring*dr_teth_to*dr_teth_to 
					- k_teth_slack*dr_teth*dr_teth);
			}
			else{
				dU_from_teth = (k_teth_spring/2)
					* (dr_teth_from*dr_teth_from - dr_teth*dr_teth);
				dU_to_teth = (0.5)*(k_teth_slack*dr_teth_to*dr_teth_to
					- k_teth_spring*dr_teth*dr_teth);
			}
		}
		else if(dr_teth > 0){
			dU_from_teth = (k_teth_spring/2)
						 * (dr_teth_from*dr_teth_from - dr_teth*dr_teth);
			dU_to_teth = (k_teth_spring/2)
					   * (dr_teth_to*dr_teth_to - dr_teth*dr_teth);
		}
		else{
			dU_from_teth = (k_teth_slack/2)
						 * (dr_teth_from*dr_teth_from - dr_teth*dr_teth);
			dU_to_teth = (k_teth_slack/2)
					   * (dr_teth_to*dr_teth_to - dr_teth*dr_teth);
		}
		double weight_to_teth;
		if(x_dub < 2*teth_comp_cutoff)
			weight_to_teth = 0;
		else
			weight_to_teth = exp(-dU_to_teth/(2*kbT));	
		double weight_from_teth;
		if(x_dub < 2*teth_dist_cutoff - 1
		&& x_dub > 2*teth_comp_cutoff + 1)
			weight_from_teth = exp(-dU_from_teth/(2*kbT));
		else
			weight_from_teth = 0; 
		if(!parameters_->motors.tethers_active){
			weight_to_teth = 0;
			weight_from_teth = 0;
		}
		double p_to_teth_i = weight_to_teth * delta_t / tau_i;
		double p_from_teth_i = weight_from_teth * delta_t / tau_i;
		if(p_to_teth_i > 1)
			printf("WARNING: p_diffuse_to_teth_i=%g for 2x=%i\n", 
					p_to_teth_i, x_dub);
		if(p_from_teth_i > 1)
			printf("WARNING: p_diffuse_from_teth_i=%g for 2x=%i\n", 
					p_from_teth_i, x_dub);
		// Input probabilities for stage_i / tethered xlinks
		p_diffuse_i_to_teth_rest_[x_dub] = p_to_teth_i;
		p_diffuse_i_from_teth_rest_[x_dub] = p_from_teth_i;
		// Run through x_dists to get probs for stage_ii / tethered xlinks
		for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
			double r_x = x_dist * site_size;
			double r_x_to = (x_dist - 1) * site_size;
			double r_x_from = (x_dist + 1) * site_size;
			double r = sqrt(r_x*r_x + r_y*r_y);
			double r_to = sqrt(r_x_to*r_x_to + r_y*r_y);
			double r_from = sqrt(r_x_from*r_x_from + r_y*r_y);
			// Get extension for current dist and steps to/from rest
			double dr = r - r_0;
			double dr_to = r_to -  r_0;
			double dr_from = r_from - r_0;
			if(dr >= 0){
				// Get corresponding changes in potential energy
				double dU_to_rest = (k_spring/2)*(dr_to*dr_to - dr*dr);
				double dU_fr_rest = (k_spring/2)*(dr_from*dr_from - dr*dr);
				// Weights according to Lanksy et al.
				double weight_to;
				double weight_from;
			   	if(x_dist == rest_dist_){
					weight_to = 0;
					weight_from = 2*exp(-dU_fr_rest/(2*kbT));
				}
				else if(x_dist == dist_cutoff_){
					weight_to = exp(-dU_to_rest/(2*kbT));
					weight_from = 0;
				}
				else{
					weight_to = exp(-dU_to_rest/(2*kbT));
					weight_from = exp(-dU_fr_rest/(2*kbT));
				}
				if(!parameters_->motors.tethers_active){
					weight_to = 0;
					weight_from = 0;
				}
				// Convolve these bitches
				double p_to_both = weight_to * weight_to_teth 
					 * delta_t / tau_ii;
				double p_to_self_from_teth = weight_to * weight_from_teth
					 * delta_t / tau_ii;
				double p_from_self_to_teth = weight_from * weight_to_teth 
					 * delta_t / tau_ii;
				double p_from_both = weight_from * weight_from_teth
					* delta_t / tau_ii;
				p_diffuse_ii_to_both_rest_[x_dub][x_dist]
					= p_to_both;
				p_diffuse_ii_to_self_from_teth_[x_dub][x_dist]
					= p_to_self_from_teth;
				p_diffuse_ii_from_self_to_teth_[x_dub][x_dist]
					= p_from_self_to_teth; 
				p_diffuse_ii_from_both_rest_[x_dub][x_dist]
					= p_from_both; 

				if(p_to_both > 1){
					printf("WARNING: p_diff_to_both=%g", 
							p_to_both);	
					printf(" for 2x=%i, x=%i\n", 
							x_dub, x_dist);
				}
				if(p_to_self_from_teth > 1){
					printf("WARNING: p_diff_to_self_fr_teth=%g", 
							p_to_self_from_teth);	
					printf(" for 2x=%i, x=%i\n", 
							x_dub, x_dist);
				}
				if(p_from_self_to_teth > 1){
					printf("WARNING: p_diff_fr_self_to_teth=%g", 
							p_from_self_to_teth);	
					printf(" for 2x=%i, x=%i\n", 
							x_dub, x_dist);
				}
				if(p_from_both > 1){
					printf("WARNING: p_diff_fr_both=%g", 
							p_from_both);	
					printf(" for 2x=%i, x=%i\n", 
							x_dub, x_dist);
				}
			}
			else{
				printf("woah mayne. xlink set parameters TWOO \n");
				exit(1);
			}
		}
	}
	// KMC STATISTICS BELOW
	double k_on = parameters_->xlinks.k_on; 
	double c_xlink = parameters_->xlinks.concentration;
	p_bind_i_.resize(3);
	p_bind_i_[0] = k_on * c_xlink * delta_t;
	double c_eff_teth = parameters_->motors.c_eff_tether;
	if(!parameters_->motors.tethers_active){
		c_eff_teth = 0;
	}
	p_bind_i_tethered_ = k_on * c_eff_teth * delta_t; 
	double c_eff_bind = parameters_->xlinks.conc_eff_bind;
	p_bind_ii_ = k_on * c_eff_bind * delta_t;
	double k_off_i = parameters_->xlinks.k_off_i;
	p_unbind_i_.resize(3);
	p_unbind_i_[0] = k_off_i * delta_t;
	// XXX janky PRC1 coop
	for(int n_neighbs = 1; n_neighbs < 3; n_neighbs++){
		double interaction_energy = 3 * n_neighbs;		// in kBT
		double weight = exp(interaction_energy / 2);
		p_bind_i_[n_neighbs] = p_bind_i_[0] * weight;
		p_unbind_i_[n_neighbs] = p_unbind_i_[0] / weight;
	}
	// Generate unbinding rates based on discretized spring extension
	double k_off_ii = parameters_->xlinks.k_off_ii;
	p_unbind_ii_.resize(dist_cutoff_ + 1);
	for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
		double r_x = x_dist*site_size;
		double r = sqrt(r_y*r_y + r_x*r_x);
		double dr = r - r_0;
		double U_xlink = (k_spring/2)*dr*dr;
		double unbind_weight = exp(U_xlink/(2*kbT));
//		printf("%i:  %g\n", distance, k_off_ii);
		p_unbind_ii_[x_dist] = k_off_ii * unbind_weight * delta_t;
	}
	p_unbind_i_tethered_.resize(2*teth_dist_cutoff + 1);
	p_unbind_ii_to_teth_.resize(2*teth_dist_cutoff + 1);
	p_unbind_ii_from_teth_.resize(2*teth_dist_cutoff + 1);
	for(int x_dub = 0; x_dub <= 2*teth_dist_cutoff; x_dub++){
		p_unbind_ii_to_teth_[x_dub].resize(dist_cutoff_ + 1); 
		p_unbind_ii_from_teth_[x_dub].resize(dist_cutoff_ + 1);
		// Calc x-distances (in nm) for tether
		double r_x_teth = x_dub * site_size / 2;
		// Calc total r values 
		double r_teth = sqrt(r_x_teth*r_x_teth + r_y_teth*r_y_teth);
		// Calc tether extensions for current dist and stepping to/from rest
		double dr_teth = r_teth - r_0_teth;	
		double U_at_teth;
		if(dr_teth > 0)
			U_at_teth = (k_teth_spring/2)*dr_teth*dr_teth;
		else
			U_at_teth = (k_teth_slack/2)*dr_teth*dr_teth;
		double weight_at_teth = exp(U_at_teth/(2*kbT));
		if(x_dub < 2*teth_comp_cutoff){
			weight_at_teth = 0;
		}
		if(!parameters_->motors.tethers_active){
			weight_at_teth = 0;
		}
		double p_unbind_teth = weight_at_teth * p_unbind_i_[0];
		if(p_unbind_teth > 1)
			printf("WARNING: p_unbind_teth (XLINK)=%g for 2x=%i\n", 
					p_unbind_teth, x_dub);
		p_unbind_i_tethered_[x_dub] = p_unbind_teth; 
		// XXX variable shift in delta-E based on x_dist XXX
		// Run through x_dists to get probs for stage_ii / tethered xlinks
		for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
			// change in teth extension if 2nd xlink head were to unbind
			double dx_teth = (double)x_dist / 2; 
			double dr_x_teth = dx_teth * site_size;
			double r_x_teth_to, 
				   r_x_teth_from; 
			if(dr_teth > 0){
				r_x_teth_to = r_x_teth - dr_x_teth;
				r_x_teth_from = r_x_teth + dr_x_teth; 	
			}
			else{
				r_x_teth_to = r_x_teth + dr_x_teth;
				r_x_teth_from = r_x_teth - dr_x_teth; 
			}
			double r_teth_to = 
				sqrt(r_x_teth_to*r_x_teth_to + r_y_teth*r_y_teth);
			double r_teth_from = 
				sqrt(r_x_teth_from*r_x_teth_from + r_y_teth*r_y_teth);
			double dr_teth_to = r_teth_to - r_0_teth; 
			double dr_teth_from = r_teth_from - r_0_teth; 
			double dU_to_teth, 
				   dU_from_teth;
			int x_from_rest_dub = abs(2*rest_dist_teth - x_dub);
			int dx_teth_dub = 2 * dx_teth;
			// Check if we're crossing over equil. point of tether 
			if(dx_teth_dub > x_from_rest_dub){
				if(r_teth < r_0_teth){
					dU_from_teth = (k_teth_slack/2)
						* (dr_teth_from*dr_teth_from - dr_teth*dr_teth);
					dU_to_teth = (0.5)*(k_teth_spring*dr_teth_to*dr_teth_to 
							- k_teth_slack*dr_teth*dr_teth);
				}
				else{
					dU_from_teth = (k_teth_spring/2)
						* (dr_teth_from*dr_teth_from - dr_teth*dr_teth);
					dU_to_teth = (0.5)*(k_teth_slack*dr_teth_to*dr_teth_to
							- k_teth_spring*dr_teth*dr_teth);
				}
			}
			else if(dr_teth > 0){
				dU_from_teth = (k_teth_spring/2)
					* (dr_teth_from*dr_teth_from - dr_teth*dr_teth);
				dU_to_teth = (k_teth_spring/2)
					* (dr_teth_to*dr_teth_to - dr_teth*dr_teth);
			}
			else{
				dU_from_teth = (k_teth_slack/2)
					* (dr_teth_from*dr_teth_from - dr_teth*dr_teth);
				dU_to_teth = (k_teth_slack/2)
					* (dr_teth_to*dr_teth_to - dr_teth*dr_teth);
			}
			double weight_to_teth = exp(-dU_to_teth/(2*kbT));
			double weight_from_teth = exp(-dU_from_teth/(2*kbT));
			if(x_dub < 2*teth_comp_cutoff){
				weight_to_teth = 0;
				weight_from_teth = 0;
			}
			if(!parameters_->motors.tethers_active){
				weight_to_teth = 0;
				weight_from_teth = 0;
			}
			double p_unbind_to = weight_to_teth * p_unbind_ii_[x_dist];
			double p_unbind_from = weight_from_teth * p_unbind_ii_[x_dist];
			if(p_unbind_to > 1)
				printf("WARNING: p_unbind_to = %g for 2x=%ix, x=%i\n", 
						p_unbind_to, x_dub, x_dist);
			if(p_unbind_from > 1)
				printf("WARNING: p_unbind_from = %g for 2x=%i, x=%i\n", 
						p_unbind_from, x_dub, x_dist);
			p_unbind_ii_to_teth_[x_dub][x_dist] = p_unbind_to;
			p_unbind_ii_from_teth_[x_dub][x_dist] = p_unbind_from; 
		}
	}
	double k_teth = parameters_->motors.k_tether;
	double k_unteth = parameters_->motors.k_untether;
	if(!parameters_->motors.tethers_active){
		k_teth = 0;
		k_unteth = 0;
	}
	p_tether_free_ = k_teth * c_xlink * delta_t; 
	p_untether_free_ = k_unteth * delta_t;
}

void AssociatedProteinManagement::InitializeLists(){

	// Stats (not a list ok bite me)
	n_bound_i_.resize(3); 				// Can have 0, 1, or 2 neighbs
	n_bound_ii_.resize(dist_cutoff_ + 1);
	n_sites_ii_.resize(dist_cutoff_ + 1);
	for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
		n_bound_ii_[x_dist] = 0;
		n_sites_ii_[x_dist] = 0;
	}
	n_bound_i_tethered_.resize(2*teth_cutoff_ + 1);
	n_bound_i_tethered_bindable_.resize(2*teth_cutoff_ + 1);
	n_bound_ii_tethered_.resize(2*teth_cutoff_ + 1);
	n_sites_i_tethered_.resize(2*teth_cutoff_ + 1);
	n_sites_ii_tethered_same_.resize(2*teth_cutoff_ + 1);
	n_sites_ii_tethered_oppo_.resize(2*teth_cutoff_ + 1);
	for(int x_dub = 0; x_dub <= 2*teth_cutoff_; x_dub++){
		n_bound_i_tethered_[x_dub] = 0;
		n_sites_i_tethered_[x_dub] = 0;
		n_bound_i_tethered_bindable_[x_dub].resize(dist_cutoff_ + 1);
		n_bound_ii_tethered_[x_dub].resize(dist_cutoff_ + 1);
		n_sites_ii_tethered_same_[x_dub].resize(dist_cutoff_ + 1);
		n_sites_ii_tethered_oppo_[x_dub].resize(dist_cutoff_ + 1);
		for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
			n_bound_i_tethered_bindable_[x_dub][x_dist] = 0;
			n_bound_ii_tethered_[x_dub][x_dist] = 0;
			n_sites_ii_tethered_same_[x_dub][x_dist] = 0;
			n_sites_ii_tethered_oppo_[x_dub][x_dist] = 0; 
		}
	}
	// Lists
	active_.resize(n_xlinks_);
	free_tethered_.resize(n_xlinks_);
	bound_i_all_.resize(n_xlinks_);
	bound_untethered_.resize(n_xlinks_);
	sites_i_untethered_.resize(n_xlinks_);
	bound_ii_.resize(dist_cutoff_ + 1);
	sites_ii_untethered_.resize(dist_cutoff_ + 1);
	for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
		bound_ii_[x_dist].resize(n_xlinks_);
		sites_ii_untethered_[x_dist].resize(n_xlinks_);
	}
	// XXX janky PRC1 coop
	bound_i_.resize(3);
	for(int n_neighbs(0); n_neighbs < 3; n_neighbs++){
		bound_i_[n_neighbs].resize(n_xlinks_);
	}
	sites_i_tethered_.resize(2*teth_cutoff_ + 1);
	bound_i_tethered_.resize(2*teth_cutoff_ + 1);
	bound_ii_tethered_.resize(2*teth_cutoff_ + 1);
	sites_ii_tethered_oppo_.resize(2*teth_cutoff_ + 1);
	sites_ii_tethered_same_.resize(2*teth_cutoff_ + 1);
	for(int x_dub = 0; x_dub <= 2*teth_cutoff_; x_dub++){
		sites_i_tethered_[x_dub].resize(n_xlinks_);
		bound_i_tethered_[x_dub].resize(n_xlinks_);
		bound_ii_tethered_[x_dub].resize(dist_cutoff_ + 1);
		sites_ii_tethered_oppo_[x_dub].resize(dist_cutoff_ + 1);
		sites_ii_tethered_same_[x_dub].resize(dist_cutoff_ + 1);
		for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
			bound_ii_tethered_[x_dub][x_dist].resize(n_xlinks_);
			sites_ii_tethered_oppo_[x_dub][x_dist].resize(n_xlinks_);
			sites_ii_tethered_same_[x_dub][x_dist].resize(n_xlinks_);
		}
	}
}

void AssociatedProteinManagement::InitializeDiffusionEvents(){

	auto binomial = [&](double p, int n){
		if(n > 0) return properties_->gsl.SampleBinomialDist(p, n);
		else return 0;
	};
	int index = 0; 
	
	diffu_events_.emplace_back(index++, 10, "i_fwd", "i_unteth", 
			binomial, &n_sites_i_, p_diffuse_i_fwd_);
	diffu_events_.emplace_back(index++, 20, "i_bck", "i_unteth", 
			binomial, &n_sites_i_, p_diffuse_i_bck_);
	for(int x(0); x <= dist_cutoff_; x++){
		diffu_events_.emplace_back(index++, 300 + x, "ii_to", 
				"ii_unteth_" + std::to_string(x), binomial, 
				&n_sites_ii_[x], p_diffuse_ii_to_rest_[x]);
		diffu_events_.emplace_back(index++, 400 + x, "ii_fr", 
				"ii_unteth_" + std::to_string(x), binomial, 
				&n_sites_ii_[x], p_diffuse_ii_from_rest_[x]);
	}
	// If tethers are enabled, initialize those events as well
	if(parameters_->motors.tethers_active){
		int comp_cutoff = properties_->kinesin4.comp_cutoff_;
		for(int x_dub(2*comp_cutoff); x_dub <= 2*teth_cutoff_; x_dub++){
			diffu_events_.emplace_back(index++, 500 + x_dub, "i_to_teth", 
					"i_teth_" + std::to_string(x_dub), binomial, 
					&n_sites_i_tethered_[x_dub], 
					p_diffuse_i_to_teth_rest_[x_dub]);
			diffu_events_.emplace_back(index++, 600 + x_dub, "i_fr_teth", 
					"i_teth_" + std::to_string(x_dub), binomial, 
					&n_sites_i_tethered_[x_dub], 
					p_diffuse_i_from_teth_rest_[x_dub]);
			for(int x(0); x <= dist_cutoff_; x++){
				diffu_events_.emplace_back(index++, 7000 + 10*x_dub + x, 
					"ii_to_both", "ii_teth_same_" + std::to_string(x_dub)
					+ "_" + std::to_string(x), binomial, 
					&n_sites_ii_tethered_same_[x_dub][x], 
					p_diffuse_ii_to_both_rest_[x_dub][x]);
				diffu_events_.emplace_back(index++, 8000 + 10*x_dub + x, 
					"ii_fr_both", "ii_teth_same_" + std::to_string(x_dub)
					+ "_" + std::to_string(x), binomial, 
					&n_sites_ii_tethered_same_[x_dub][x], 
					p_diffuse_ii_from_both_rest_[x_dub][x]);
				diffu_events_.emplace_back(index++, 9000 + 10*x_dub + x, 
					"ii_to_fr", "ii_teth_oppo_" + std::to_string(x_dub)
					+ "_" + std::to_string(x), binomial, 
					&n_sites_ii_tethered_oppo_[x_dub][x], 
					p_diffuse_ii_to_self_from_teth_[x_dub][x]);
				diffu_events_.emplace_back(index++, 10000 + 10*x_dub + x, 
					"ii_fr_to", "ii_teth_oppo_" + std::to_string(x_dub)
					+ "_" + std::to_string(x), binomial, 
					&n_sites_ii_tethered_oppo_[x_dub][x], 
					p_diffuse_ii_from_self_to_teth_[x_dub][x]);
			}
		}
	}
	int n_pops = 0;
	int n_entries[diffu_events_.size()];
	int indices[diffu_events_.size()]
		[2*(2*teth_cutoff_+1)*(dist_cutoff_+1)];
	for(int i_entry(0); i_entry < diffu_events_.size(); i_entry++){
		bool new_pop(true);
		int pop_index(0); 
		for(int i_pop(0); i_pop < n_pops; i_pop++){
			std::string pop = diffu_events_[indices[i_pop][0]].target_pop_;
			if(diffu_events_[i_entry].target_pop_ == pop){
				new_pop = false;
				pop_index = i_pop;
			}
		}
		if(new_pop){
			indices[n_pops][0] = diffu_events_[i_entry].index_;
			n_entries[n_pops] = 1; 
			n_pops++;
		}
		else{
			int entry_no = n_entries[pop_index];
			indices[pop_index][entry_no] = diffu_events_[i_entry].index_;
			n_entries[pop_index]++;
		}
	}
	diffu_by_pop_.resize(n_pops);
	for(int i_pop(0); i_pop < n_pops; i_pop++){
		printf("DIFF Pop #%i:", i_pop);
		std::cout << diffu_events_[indices[i_pop][0]].target_pop_;
		diffu_by_pop_[i_pop].resize(n_entries[i_pop]);
		for(int i_entry(0); i_entry < n_entries[i_pop]; i_entry++){
			diffu_by_pop_[i_pop][i_entry] = indices[i_pop][i_entry];
			printf(" %i,", indices[i_pop][i_entry]);
		}
		printf("\n");
	}
}

void AssociatedProteinManagement::InitializeKMCEvents(){

	int index(0);
	auto binomial = [&](double p, int n){
		if(n > 0) return properties_->gsl.SampleBinomialDist(p, n);
		else return 0;
	};
	auto poisson_bind = [&](double p, int n){
		if(n > 0){ 
			double n_wt = GetWeightBindII();
			if(n_wt > 0) 
				return properties_->gsl.SamplePoissonDist(p*n_wt);
			else return 0;
		}
		else return 0;
	};
	auto poisson_bind_i_teth = [&](double p, int n){
		if(n > 0){
			double n_wt = GetWeightBindITethered();
			if(n_wt > 0) 
				return properties_->gsl.SamplePoissonDist(p*n_wt);
			else return 0;
		}
		else return 0;
	};
	auto poisson_bind_ii_teth = [&](double p, int n){
		if(n > 0){
			double n_wt = GetWeightBindIITethered();
			if(n_wt > 0)
				return properties_->gsl.SamplePoissonDist(p*n_wt);
			else return 0;
		}
		else return 0;
	};
	for(int n_neighbs(0); n_neighbs < 3; n_neighbs++){
		kmc_events_.emplace_back(index++, 100 + n_neighbs, "bind_i", 
				"unocc_" + std::to_string(n_neighbs), binomial, 
				&properties_->microtubules.n_unoccupied_xl_[n_neighbs], 
				p_bind_i_[n_neighbs]);
	}
	kmc_events_.emplace_back(index++, 20, "bind_ii", "bound_i", 
			poisson_bind, &n_bound_i_tot_, p_bind_ii_);
	for(int n_neighbs(0); n_neighbs < 3; n_neighbs++){
		kmc_events_.emplace_back(index++, 200 + n_neighbs, "unbind_i", 
				"bound_i", binomial, &n_bound_i_[n_neighbs], 
				p_unbind_i_[n_neighbs]);	
	}
	for(int x(0); x <= dist_cutoff_; x++){
		kmc_events_.emplace_back(index++, 40 + x, "unbind_ii", 
				"bound_ii_" + std::to_string(x), binomial, 
				&n_bound_ii_[x], p_unbind_ii_[x]);
	}
	// If tethering is enabled, add those event populations as well
	if(parameters_->motors.tethers_active){
		kmc_events_.emplace_back(index++, 11, "bind_i_teth", "free_teth", 
				poisson_bind_i_teth, &n_free_tethered_, p_bind_i_tethered_);
		kmc_events_.emplace_back(index++, 21, "bind_ii_teth", 
				"bound_i_teth", poisson_bind_ii_teth,
				&n_bound_i_tethered_tot_, p_bind_ii_);

		int comp_cutoff = properties_->kinesin4.comp_cutoff_;
		for(int x_dub(2*comp_cutoff); x_dub <= 2*teth_cutoff_; x_dub++){
			kmc_events_.emplace_back(index++, 300 + x_dub, "unbind_i_teth",
					"bound_i_teth", binomial, &n_bound_i_tethered_[x_dub], 
					p_unbind_i_tethered_[x_dub]);
		}
		for(int x_dub(2*comp_cutoff); x_dub <= 2*teth_cutoff_; x_dub++){
			for(int x(0); x <= dist_cutoff_; x++){
				kmc_events_.emplace_back(index++, 4000 + 10*x_dub + x, 
						"unbind_ii_to_teth", "bound_ii_teth_" + 
						std::to_string(x_dub) + "_" + std::to_string(x), 
						binomial, &n_bound_ii_tethered_[x_dub][x], 
						p_unbind_ii_to_teth_[x_dub][x]);
			}
		}
		for(int x_dub(2*comp_cutoff); x_dub <= 2*teth_cutoff_; x_dub++){
			for(int x(0); x <= dist_cutoff_; x++){
				kmc_events_.emplace_back(index++, 5000 + 10*x_dub + x,
						"unbind_ii_fr_teth", "bound_ii_teth_" + 
						std::to_string(x_dub) + "_" + std::to_string(x), 
						binomial, &n_bound_ii_tethered_[x_dub][x], 
						p_unbind_ii_from_teth_[x_dub][x]);
			}
		}
		kmc_events_.emplace_back(index++, 50, "tether_free", "unteth_mots",
				binomial, &properties_->kinesin4.n_bound_untethered_, 
				p_tether_free_);
		kmc_events_.emplace_back(index++, 60, "untether_free", "free_teth", 
				binomial, &n_free_tethered_, p_untether_free_);
	}
	int n_pops = 0;
	int n_entries[kmc_events_.size()];
	int indices[kmc_events_.size()][2*(2*teth_cutoff_+1)*(dist_cutoff_+1)];
	for(int i_entry(0); i_entry < kmc_events_.size(); i_entry++){
		bool new_pop(true);
		int pop_index(0); 
		for(int i_pop(0); i_pop < n_pops; i_pop++){
			std::string pop = kmc_events_[indices[i_pop][0]].target_pop_;
			if(kmc_events_[i_entry].target_pop_ == pop){
				new_pop = false;
				pop_index = i_pop;
			}
		}
		if(new_pop){
			indices[n_pops][0] = kmc_events_[i_entry].index_;
			n_entries[n_pops] = 1; 
			n_pops++;
		}
		else{
			int entry_no = n_entries[pop_index];
			indices[pop_index][entry_no] = kmc_events_[i_entry].index_;
			n_entries[pop_index]++;
		}
	}
	kmc_by_pop_.resize(n_pops);
	for(int i_pop(0); i_pop < n_pops; i_pop++){
		printf("KMC Pop #%i:", i_pop);
		std::cout << kmc_events_[indices[i_pop][0]].target_pop_;
		kmc_by_pop_[i_pop].resize(n_entries[i_pop]);
		for(int i_entry(0); i_entry < n_entries[i_pop]; i_entry++){
			kmc_by_pop_[i_pop][i_entry] = indices[i_pop][i_entry];
			printf(" %i,", indices[i_pop][i_entry]);
		}
		printf("\n");
	}
}

void AssociatedProteinManagement::UpdateAllLists(){
	
	UpdateSingleBoundList();
	UpdateDoubleBoundList();
	properties_->microtubules.UpdateUnoccupied();
	if(parameters_->motors.tethers_active){
		UpdateBoundITethered();
		UpdateBoundIITethered();
		UpdateFreeTetheredList();
		UpdateBoundUntethered();
		properties_->kinesin4.UpdateBoundUntethered();
	}
}

void AssociatedProteinManagement::UpdateSingleBoundList(){
	
	n_bound_i_tot_ = 0;
	for(int n_neighbs = 0; n_neighbs < 3; n_neighbs++){
		n_bound_i_[n_neighbs] = 0;
	}
	for(int i_xlink = 0; i_xlink < n_active_; i_xlink++){
		AssociatedProtein *xlink = active_[i_xlink];
		if(xlink->heads_active_ == 1
		&& xlink->tethered_ == false){
			int n_neighbs = 0;
			Tubulin *site = xlink->GetActiveHeadSite(); 
			int i_site = site->index_;
			int i_plus = site->mt_->plus_end_;
			int i_minus = site->mt_->minus_end_;
			int dx = site->mt_->delta_x_; 
			if(i_site == i_plus){
				if(site->mt_->lattice_[i_site-dx].xlink_ != nullptr)
					n_neighbs++;
			}
			else if(i_site == i_minus){
				if(site->mt_->lattice_[i_site+dx].xlink_ != nullptr)
					n_neighbs++;
			}
			else{
				if(site->mt_->lattice_[i_site-dx].xlink_ != nullptr)
					n_neighbs++;
				if(site->mt_->lattice_[i_site+dx].xlink_ != nullptr)
					n_neighbs++;
			}
			bound_i_[n_neighbs][n_bound_i_[n_neighbs]] = xlink;
			n_bound_i_[n_neighbs]++;
			bound_i_all_[n_bound_i_tot_] = xlink;
			n_bound_i_tot_++;
		}
		else if(xlink->tethered_ == true){
			if(xlink->heads_active_ == 1
			&& xlink->motor_->heads_active_ == 0){
				int n_neighbs = 0;
				Tubulin *site = xlink->GetActiveHeadSite(); 
				int i_site = site->index_;
				int i_plus = site->mt_->plus_end_;
				int i_minus = site->mt_->minus_end_;
				int dx = site->mt_->delta_x_; 
				if(i_site == i_plus){
					if(site->mt_->lattice_[i_site-dx].xlink_ != nullptr)
						n_neighbs++;
				}
				else if(i_site == i_minus){
					if(site->mt_->lattice_[i_site+dx].xlink_ != nullptr)
						n_neighbs++;
				}
				else{
					if(site->mt_->lattice_[i_site-dx].xlink_ != nullptr)
						n_neighbs++;
					if(site->mt_->lattice_[i_site+dx].xlink_ != nullptr)
						n_neighbs++;
				}
				bound_i_[n_neighbs][n_bound_i_[n_neighbs]] = xlink;
				n_bound_i_[n_neighbs]++;
				bound_i_all_[n_bound_i_tot_] = xlink;
				n_bound_i_tot_++;
			}
		}
	}
}

void AssociatedProteinManagement::UpdateBoundITethered(){
	
	int comp_cutoff = properties_->kinesin4.comp_cutoff_;
	for(int x_dub(2*comp_cutoff); x_dub <= 2*teth_cutoff_; x_dub++){
		n_bound_i_tethered_[x_dub] = 0;
	}
	n_bound_i_tethered_tot_ = 0;
	for(int i_xlink = 0; i_xlink < n_active_; i_xlink++){
		AssociatedProtein *xlink = active_[i_xlink];
		if(xlink->heads_active_ == 1
		&& xlink->tethered_ == true){
			if(xlink->motor_->heads_active_ > 0){
				xlink->motor_->UpdateExtension(); 
				// Make sure we didn't force an untether event
				if(xlink->tethered_){
					int x_dub = xlink->motor_->x_dist_doubled_;
					int index = n_bound_i_tethered_[x_dub];
					bound_i_tethered_[x_dub][index] = xlink;
					n_bound_i_tethered_[x_dub]++; 
					n_bound_i_tethered_tot_++;
				}
			}
		}
	}
}

void AssociatedProteinManagement::UpdateDoubleBoundList(){

	for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
		n_bound_ii_[x_dist] = 0;
	}
	for(int i_xlink = 0; i_xlink < n_active_; i_xlink++){
		AssociatedProtein *xlink = active_[i_xlink];
		if(xlink->heads_active_ == 2
		&& xlink->tethered_ == false){
			int x_dist = xlink->x_dist_;
			int index = n_bound_ii_[x_dist];
			bound_ii_[x_dist][index] = xlink;
			n_bound_ii_[x_dist]++;
		}
		else if(xlink->tethered_ == true){
			if(xlink->heads_active_ == 2
			&& xlink->motor_->heads_active_ == 0){
				int x_dist = xlink->x_dist_;
				int index = n_bound_ii_[x_dist];
				bound_ii_[x_dist][index] = xlink;
				n_bound_ii_[x_dist]++;
			}
		}
	}
}

void AssociatedProteinManagement::UpdateBoundIITethered(){

	int comp_cutoff = properties_->kinesin4.comp_cutoff_;
	for(int x_dub(2*comp_cutoff); x_dub <= 2*teth_cutoff_; x_dub++){
		for(int x_dist(0); x_dist <= dist_cutoff_; x_dist++){
			n_bound_ii_tethered_[x_dub][x_dist] = 0;
		}
	}
	for(int i_xlink = 0; i_xlink < n_active_; i_xlink++){
		AssociatedProtein *xlink = active_[i_xlink];
		if(xlink->heads_active_ == 2
		&& xlink->tethered_ == true){
			if(xlink->motor_->heads_active_ > 0){
				xlink->UpdateExtension();
				xlink->motor_->UpdateExtension(); 
				// Make sure we didn't force an untether or unbind event
				if(xlink->heads_active_ == 2
				&& xlink->tethered_){
					int x_dist = xlink->x_dist_; 
					int x_dub = xlink->motor_->x_dist_doubled_;
					int index = n_bound_ii_tethered_[x_dub][x_dist]; 
					bound_ii_tethered_[x_dub][x_dist][index] = xlink;
					n_bound_ii_tethered_[x_dub][x_dist]++; 
				}
			}
		}
	}
}

void AssociatedProteinManagement::UpdateFreeTetheredList(){

	n_free_tethered_ = 0;
	for(int i_xlink = 0; i_xlink < n_active_; i_xlink++){
		AssociatedProtein *xlink = active_[i_xlink];
		if(xlink->heads_active_ == 0
		&& xlink->tethered_ == true){
			if(xlink->motor_->heads_active_ > 0){
				free_tethered_[n_free_tethered_] = xlink;
				n_free_tethered_++;
			}
			else{
				printf("woah. error in update_free_teth_list (XLINKS)\n");
				exit(1);
			}
		}
	}
}

void AssociatedProteinManagement::UpdateBoundUntethered(){
	
	n_bound_untethered_ = 0;
	for(int i_xlink = 0; i_xlink < n_active_; i_xlink++){
		AssociatedProtein *xlink = active_[i_xlink];
		if(xlink->heads_active_ > 0){
			if(xlink->tethered_ == false){
				bound_untethered_[n_bound_untethered_] = xlink;
				n_bound_untethered_++;
			}
		}
	}
}

void AssociatedProteinManagement::UpdateAllSiteLists(){

	UpdateSingleUntetheredSites();
	UpdateDoubleUntetheredSites();
	if(parameters_->motors.tethers_active){
		UpdateSingleTetheredSites();
		UpdateDoubleTetheredSites();
	}
}

void AssociatedProteinManagement::UpdateSingleUntetheredSites(){

	n_sites_i_ = 0;
	for(int i_xlink = 0; i_xlink < n_active_; i_xlink++){
		AssociatedProtein *xlink = active_[i_xlink];
		if(xlink->heads_active_ == 1
		&& xlink->tethered_ == false){
			Tubulin *site = xlink->GetActiveHeadSite();
			sites_i_untethered_[n_sites_i_] = site;
			n_sites_i_++;
		}
		// xlinks tethered to satellite motors diffuse as untethered xlinks
		else if(xlink->tethered_ == true){
			if(xlink->heads_active_ == 1
			&& xlink->motor_->heads_active_ == 0){
				Tubulin *site = xlink->GetActiveHeadSite();
				sites_i_untethered_[n_sites_i_] = site;
				n_sites_i_++;
			}
		}
	}
}

void AssociatedProteinManagement::UpdateDoubleUntetheredSites(){

	for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
		n_sites_ii_[x_dist] = 0;
	}
	for(int i_xlink = 0; i_xlink < n_active_; i_xlink++){
		AssociatedProtein *xlink = active_[i_xlink];
		if(xlink->heads_active_ == 2
		&& xlink->tethered_ == false){
			xlink->UpdateExtension();
			// Make sure we didn't force an unbind event
			if(xlink->heads_active_ == 2){
				int x_dist = xlink->x_dist_;
				Tubulin *site_one = xlink->site_one_;
				int index_one = n_sites_ii_[x_dist];
				sites_ii_untethered_[x_dist][index_one] = site_one;
				n_sites_ii_[x_dist]++;
				Tubulin *site_two = xlink->site_two_;
				int index_two = n_sites_ii_[x_dist];
				sites_ii_untethered_[x_dist][index_two] = site_two;
				n_sites_ii_[x_dist]++;
			}
			else{
				printf("wat in update_dub_unteth sites xlink\n");
			}
		}
		// xlinks tethered to satellite motors diffuse as untethered xlinks
		else if(xlink->tethered_ == true){
			if(xlink->heads_active_ == 2
			&& xlink->motor_->heads_active_ == 0){
				xlink->UpdateExtension();
				// Make sure we didn't force an unbind event
				if(xlink->heads_active_ == 2){
					int x_dist = xlink->x_dist_;
					Tubulin *site_one = xlink->site_one_;
					int index_one = n_sites_ii_[x_dist];
					sites_ii_untethered_[x_dist][index_one] = site_one;
					n_sites_ii_[x_dist]++;
					Tubulin *site_two = xlink->site_two_;
					int index_two = n_sites_ii_[x_dist];
					sites_ii_untethered_[x_dist][index_two] = site_two;
					n_sites_ii_[x_dist]++;
				}
				else{
					printf("wat in update_dub_unteth sites xlink\n");
				}
			}
		}
	}
}

void AssociatedProteinManagement::UpdateSingleTetheredSites(){

	int comp_cutoff = properties_->kinesin4.comp_cutoff_;
	for(int x_dub(2*comp_cutoff); x_dub <= 2*teth_cutoff_; x_dub++){
		n_sites_i_tethered_[x_dub] = 0;
	}
	for(int i_xlink = 0; i_xlink < n_active_; i_xlink++){
		AssociatedProtein *xlink = active_[i_xlink];
        if(xlink->heads_active_ == 1
		&& xlink->tethered_ == true){
			// Dont count xlinks tethered to satellite motors
			if(xlink->motor_->heads_active_ > 0){
				xlink->motor_->UpdateExtension();
				// Make sure we didn't force an untether event
				if(xlink->tethered_ == true){
					int x_dub = xlink->motor_->x_dist_doubled_;
					Tubulin *site = xlink->GetActiveHeadSite();
					int index = n_sites_i_tethered_[x_dub]; 
					sites_i_tethered_[x_dub][index] = site;
					n_sites_i_tethered_[x_dub];
				}
			}
		}
	}
}

void AssociatedProteinManagement::UpdateDoubleTetheredSites(){

	int comp_cutoff = properties_->kinesin4.comp_cutoff_;
	for(int x_dub(2*comp_cutoff); x_dub <= 2*teth_cutoff_; x_dub++){
		for(int x = 0; x <= dist_cutoff_; x++){
			n_sites_ii_tethered_same_[x_dub][x] = 0;
			n_sites_ii_tethered_oppo_[x_dub][x] = 0;
		}
	}
	for(int i_xlink = 0; i_xlink < n_active_; i_xlink++){
		AssociatedProtein *xlink = active_[i_xlink];
		if(xlink->heads_active_ == 2
		&& xlink->tethered_ == true){
			// Dont count xlinks tethered to satellite motors
			if(xlink->motor_->heads_active_ > 0){
				xlink->UpdateExtension();
				xlink->motor_->UpdateExtension();
				// Make sure an untether or unbind event wasn't forced
				if(xlink->heads_active_ == 2
				&& xlink->tethered_ == true){
					int x = xlink->x_dist_;
					int x_dub = xlink->motor_->x_dist_doubled_;
					// Site one
					Tubulin *site_one = xlink->site_one_;
					if(x == rest_dist_){
						int i_one = n_sites_ii_tethered_same_[x_dub][x];
						sites_ii_tethered_same_[x_dub][x][i_one] = 
							site_one;
						int i_two = n_sites_ii_tethered_oppo_[x_dub][x];
						sites_ii_tethered_oppo_[x_dub][x][i_two] =
							site_one;
						n_sites_ii_tethered_same_[x_dub][x]++;
						n_sites_ii_tethered_oppo_[x_dub][x]++;
					}
					else if(site_one->EquilibriumInSameDirection()){
						int i = n_sites_ii_tethered_same_[x_dub][x];
						sites_ii_tethered_same_[x_dub][x][i] = site_one;
						n_sites_ii_tethered_same_[x_dub][x]++;
					}
					else{
						int i = n_sites_ii_tethered_oppo_[x_dub][x]; 
						sites_ii_tethered_oppo_[x_dub][x][i] = site_one;
						n_sites_ii_tethered_oppo_[x_dub][x]++;
					}
					// Site two
					Tubulin *site_two = xlink->site_two_;
					if(x == rest_dist_){
						int i_one = n_sites_ii_tethered_same_[x_dub][x];
						sites_ii_tethered_same_[x_dub][x][i_one] = 
							site_two;
						int i_two = n_sites_ii_tethered_oppo_[x_dub][x];
						sites_ii_tethered_oppo_[x_dub][x][i_two] = 
							site_two;
						n_sites_ii_tethered_same_[x_dub][x]++;
						n_sites_ii_tethered_oppo_[x_dub][x]++;
					}
					else if(site_two->EquilibriumInSameDirection()){
						int i = n_sites_ii_tethered_same_[x_dub][x];
						sites_ii_tethered_same_[x_dub][x][i] = site_two;
						n_sites_ii_tethered_same_[x_dub][x]++;
					}
					else{
						int i = n_sites_ii_tethered_oppo_[x_dub][x];
						sites_ii_tethered_oppo_[x_dub][x][i] = site_two;
						n_sites_ii_tethered_oppo_[x_dub][x]++;
					}
				}
			}
		}
	}
}

AssociatedProtein* AssociatedProteinManagement::GetFreeXlink(){


	// Randomly choose an unbound xlink
	int i_xlink = properties_->gsl.GetRanInt(n_xlinks_);
	AssociatedProtein *xlink = &xlinks_[i_xlink];
	int attempts = 0;
	while(xlink->heads_active_ > 0 || xlink->tethered_){
		i_xlink++;
		if(i_xlink == n_xlinks_)
			i_xlink = 0;
		xlink = &xlinks_[i_xlink];
		attempts++;
		if(attempts > n_xlinks_){
			printf("error in get_free_xlink\n");
			exit(1);
		}
	}
	return xlink;
}

AssociatedProtein* AssociatedProteinManagement::GetBoundUntetheredXlink(){

	UpdateBoundUntethered();
	if(n_bound_untethered_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound_untethered_);
		AssociatedProtein* xlink = bound_untethered_[i_entry];
		return xlink;
	}
	else{
		printf("Error in GetUnTetheredXlink: no untethered xlinks!\n");
		exit(1);
	}
}

void AssociatedProteinManagement::GenerateDiffusionList(){

	sys_time start = sys_clock::now();

	// Update all different population lists
	UpdateAllSiteLists();
	// Run through serialized population types
	int n_events_tot = 0;
	for(int i_entry = 0; i_entry < diffu_events_.size(); i_entry++){
		diffu_events_[i_entry].SampleStatistics();
		int n_expected = diffu_events_[i_entry].n_expected_;
		for(int i_event(0); i_event < n_expected; i_event++) 
			n_events_tot++;
	}

	sys_time finish = sys_clock::now();
	auto elapsed = std::chrono::duration_cast<t_unit>(finish - start);
	properties_->wallace.t_xlinks_dif_[1] += elapsed.count();

	start = sys_clock::now();

	// Scan through all target populations & ensure none become negative 
	start = sys_clock::now();
	for(int i_pop(0); i_pop < diffu_by_pop_.size(); i_pop++){
		int n_entries = diffu_by_pop_[i_pop].size();
		if(n_entries > 1){
			int n_events = 0;
			double p_tot = 0;
			for(int i_entry(0); i_entry < n_entries; i_entry++){
				int index = diffu_by_pop_[i_pop][i_entry];
				n_events += diffu_events_[index].n_expected_;
				p_tot += diffu_events_[index].p_occur_; 
			}
			int n_available = *diffu_events_
				[diffu_by_pop_[i_pop][0]].pop_ptr_;
			while(n_events > n_available){
				double p_cum = 0;
				double ran = properties_->gsl.GetRanProb();
				for(int i_entry(0); i_entry < n_entries; i_entry++){
					int index = diffu_by_pop_[i_pop][i_entry];
					p_cum += diffu_events_[index].p_occur_ / p_tot;
					if(p_cum >= ran
					&& diffu_events_[index].n_expected_ > 0){
						diffu_events_[index].n_expected_--;
						n_events_tot--;
						n_events--;
						break;
					}
				}
			}
		}
	}
	// Construct KMC list if any events are predicted, otherwise clear it
	if(n_events_tot > 0){
		int pre_array[n_events_tot];
		int i_event = 0;
		for(int i_entry(0); i_entry < diffu_events_.size(); i_entry++){
			int n_expected = diffu_events_[i_entry].n_expected_;
			for(int i_event(0); i_event < n_expected; i_event++){
				if(i_event >= n_events_tot){
					printf("WOAH BUDDY. check XLINK genKMC\n");
					exit(1);
				}
				pre_array[i_event] = diffu_events_[i_entry].kmc_code_;
				i_event++;
			}
		}
		if(n_events_tot > 1) gsl_ran_shuffle(properties_->gsl.rng_, 
				pre_array, n_events_tot, sizeof(int));
		dif_list_.resize(n_events_tot);
		for(int i_entry(0); i_entry < n_events_tot; i_entry++){
			dif_list_[i_entry] = pre_array[i_entry];
		}
	}
	else dif_list_.clear();

	finish = sys_clock::now();
	elapsed = std::chrono::duration_cast<t_unit>(finish - start);
	properties_->wallace.t_xlinks_dif_[2] += elapsed.count();
}

void AssociatedProteinManagement::RunDiffusion(){

//	printf("start of xlink diffusion cycle\n");
	sys_time start1 = sys_clock::now();

	if(parameters_->xlinks.concentration > 0) GenerateDiffusionList();
	else return;
	
	sys_time start2 = sys_clock::now();

	if(dif_list_.empty() == false){
		int n_events = dif_list_.size();
//		printf("%i XLINK DIFFUSION EVENTS\n", n_events);
		int x_dist, 			// refers to dist (in sites) of xlink ext
			x_dist_dub;			// refers to 2*dist ('') of motor ext
		for(int i_event = 0; i_event < n_events; i_event++){
			int diff_event = dif_list_[i_event];
			if(diff_event >= 300 && diff_event < 400){
				x_dist = diff_event % 100; 
				diff_event = 30;
			}
			else if(diff_event >= 400 && diff_event < 500){
				x_dist = diff_event % 100;
				diff_event = 40;
			}	
			else if(diff_event >= 500 && diff_event < 600){
				x_dist_dub = diff_event % 100;
				diff_event = 50;
			}
			else if(diff_event >= 600 && diff_event < 700){
				x_dist_dub = diff_event % 100;
				diff_event = 60;
			}
			else if(diff_event >= 7000 && diff_event < 8000){
				x_dist = diff_event % 10;
				x_dist_dub = (diff_event % 1000 - x_dist)/10; 
				diff_event = 70;
			}
			else if(diff_event >= 8000 && diff_event < 9000){
				x_dist = diff_event % 10;
				x_dist_dub = (diff_event % 1000 - x_dist)/10; 
				diff_event = 80;
			}
			else if(diff_event >= 9000 && diff_event < 10000){
				x_dist = diff_event % 10;
				x_dist_dub = (diff_event % 1000 - x_dist)/10; 
				diff_event = 90;
			}
			else if(diff_event >= 10000 && diff_event < 11000){
				x_dist = diff_event % 10;
				x_dist_dub = (diff_event % 1000 - x_dist)/10; 
				diff_event = 100;
			}
//			printf("ole:\n");
//			properties_->wallace.PrintMicrotubules(0.000);
			switch(diff_event){
				case 10:
//					printf("unteth_i step fwd\n");
					Diffuse_I_Forward();
					break;
				case 20:
//					printf("unteth_i step bck\n");
					Diffuse_I_Backward();
					break;
				case 30:
//					printf("unteth_ii step to rest (%i)[%i avail]\n", 
//						x_dist, n_sites_ii_[x_dist]);
					Diffuse_II_ToRest(x_dist);
					break;
				case 40:
//					printf("unteth_ii step from rest (%i)[%i avail]\n", 
//						x_dist, n_sites_ii_[x_dist]);
					Diffuse_II_FromRest(x_dist);
					break;
				case 50:
///					printf("teth_i step to rest (%i)\n", x_dist_dub);
					Diffuse_I_ToTethRest(x_dist_dub);
					break;
				case 60:
//					printf("teth_i step from rest (%i)\n", x_dist_dub);
					Diffuse_I_FromTethRest(x_dist_dub);
					break;
				case 70:
//					printf("teth_ii step to both (2x: %i, x: %i)\n", 
//						x_dist_dub, x_dist);
					Diffuse_II_ToBothRest(x_dist_dub, x_dist);
					break;
				case 80:
//					printf("teth_ii step from both (2x: %i, x: %i)\n", 
//							x_dist_dub, x_dist);
					Diffuse_II_FromBothRest(x_dist_dub, x_dist);
					break;
				case 90:
//					printf("teth_ii step to self from teth 
//						(2x: %i, x: %i)\n" ,x_dist_dub, x_dist);
					Diffuse_II_ToSelf_FromTeth(x_dist_dub, x_dist);
					break;
				case 100:
//					printf("teth_ii step from self to teth 
//						(2x: %i, x: %i)\n" ,x_dist_dub, x_dist);
					Diffuse_II_FromSelf_ToTeth(x_dist_dub, x_dist);
					break;
			}
		}
	}
	sys_time finish = sys_clock::now();
	auto elapsed = std::chrono::duration_cast<t_unit>(finish - start1);
	properties_->wallace.t_xlinks_dif_[0] += elapsed.count();
	elapsed = std::chrono::duration_cast<t_unit>(finish - start2);
	properties_->wallace.t_xlinks_dif_[3] += elapsed.count();
}

void AssociatedProteinManagement::Diffuse_I_Forward(){

	UpdateSingleUntetheredSites();
	if(n_sites_i_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_sites_i_);
		Tubulin* site = sites_i_untethered_[i_entry]; 
		AssociatedProtein* xlink = site->xlink_;
		int i_site = site->index_;
		// cant step off them MTs
		if(i_site != site->mt_->n_sites_ - 1){
			Tubulin *old_site = site;
			Tubulin *new_site = &site->mt_->lattice_[i_site + 1];
			if(new_site->occupied_ == false){
				old_site->xlink_ = nullptr;
				old_site->occupied_ = false;
				if(old_site == xlink->site_one_)
					xlink->site_one_ = new_site;
				else
					xlink->site_two_ = new_site;
				new_site->xlink_ = xlink;
				new_site->occupied_ = true;
			}
		}
	}
	else{
		printf("come on kid we cant step_i_fwd\n");
		exit(1);
	}
}

void AssociatedProteinManagement::Diffuse_I_Backward(){

	UpdateSingleUntetheredSites();
	if(n_sites_i_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_sites_i_);
		Tubulin* site = sites_i_untethered_[i_entry]; 
		AssociatedProtein* xlink = site->xlink_;
		int i_site = site->index_;
		// can't step off them MTs
		if(i_site != 0){
			Tubulin *old_site = site;
			Tubulin *new_site = &site->mt_->lattice_[i_site - 1];
			if(new_site->occupied_ == false){
				old_site->xlink_ = nullptr;
				old_site->occupied_ = false;
				if(old_site == xlink->site_one_)
					xlink->site_one_ = new_site;
				else
					xlink->site_two_ = new_site;
				new_site->xlink_ = xlink;
				new_site->occupied_ = true;
			}
		}
	}
	else{
		printf("come on kid we cant step_i_bck\n");
		exit(1);
	}
}

void AssociatedProteinManagement::Diffuse_II_ToRest(int x){

	UpdateDoubleUntetheredSites();
	int n_bound = n_sites_ii_[x];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site = sites_ii_untethered_[x][i_entry];
		int i_site = site->index_;
		AssociatedProtein* xlink = site->xlink_;
		int dx = xlink->GetDirectionTowardRest(site);
		// cant step off them MTs
		if(!(i_site == site->mt_->n_sites_ - 1 && dx == 1)
		&& !(i_site == 0 && dx == -1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &site->mt_->lattice_[i_site + dx];
			if(new_site->occupied_ == false){
				old_site->xlink_ = nullptr;
				old_site->occupied_ = false;
				if(old_site == xlink->site_one_)
					xlink->site_one_ = new_site;
				else
					xlink->site_two_ = new_site;
				new_site->xlink_ = xlink;
				new_site->occupied_ = true;
				// Update statistics
				xlink->UpdateExtension();
			}
		}
	}
	else{
		printf("we cannot step_ii_toward (xlink UT)\n");
		exit(1);
	}
}

void AssociatedProteinManagement::Diffuse_II_FromRest(int x){

	UpdateDoubleUntetheredSites();
	int n_bound = n_sites_ii_[x];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site = sites_ii_untethered_[x][i_entry];
		int i_site = site->index_;
		AssociatedProtein* xlink = site->xlink_;
		int dx = -1 * xlink->GetDirectionTowardRest(site);
		// cant step off them MTs
		if(!(i_site == site->mt_->n_sites_ - 1 && dx == 1)
		&& !(i_site == 0 && dx == -1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &site->mt_->lattice_[i_site + dx];
			if(new_site->occupied_ == false){
				old_site->xlink_ = nullptr;
				old_site->occupied_ = false;
				if(old_site == xlink->site_one_)
					xlink->site_one_ = new_site;
				else
					xlink->site_two_ = new_site;
				new_site->xlink_ = xlink;
				new_site->occupied_ = true;
				// Update statistics
				xlink->UpdateExtension();
			}
		}
	}
	else{
		printf("we cannot step_ii_from (xlink UT)\n");
		exit(1);
	}
}

void AssociatedProteinManagement::Diffuse_I_ToTethRest(int x_dub){

	UpdateSingleTetheredSites();
	int n_bound = n_sites_i_tethered_[x_dub];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site = sites_i_tethered_[x_dub][i_entry];
		int i_site = site->index_;
		AssociatedProtein *xlink = site->xlink_;
		// direction xlink needs to move to rest is opposite of motor's
		int dx = -1 * xlink->motor_->GetDirectionTowardRest(); 
		// cant step off them MTs
		if(!(i_site == site->mt_->n_sites_ - 1 && dx == 1)
		&& !(i_site == 0 && dx == -1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &site->mt_->lattice_[i_site + dx];
			if(new_site->occupied_ == false){
				old_site->xlink_ = nullptr;
				old_site->occupied_ = false;
				if(old_site == xlink->site_one_)
					xlink->site_one_ = new_site;
				else
					xlink->site_two_ = new_site;
				new_site->xlink_ = xlink;
				new_site->occupied_ = true;
				// Update statistics
				xlink->motor_->UpdateExtension();
			}
		}
	}
	else{
		printf("we cannot step_i_to_teth_rest (xlink)\n");
		exit(1);
	}
}

void AssociatedProteinManagement::Diffuse_I_FromTethRest(int x_dub){

	UpdateSingleTetheredSites();
	int n_bound = n_sites_i_tethered_[x_dub];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site = sites_i_tethered_[x_dub][i_entry];
		int i_site = site->index_;
		AssociatedProtein *xlink = site->xlink_;
		// direction xlink needs to move to rest is opposite of motor's
		int dx = xlink->motor_->GetDirectionTowardRest(); 
		// cant step off them MTs
		if(!(i_site == site->mt_->n_sites_ - 1 && dx == 1)
		&& !(i_site == 0 && dx == -1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &site->mt_->lattice_[i_site + dx];
			if(new_site->occupied_ == false){
				old_site->xlink_ = nullptr;
				old_site->occupied_ = false;
				if(old_site == xlink->site_one_)
					xlink->site_one_ = new_site;
				else
					xlink->site_two_ = new_site;
				new_site->xlink_ = xlink;
				new_site->occupied_ = true;
				// Update statistics
				xlink->motor_->UpdateExtension();
			}
		}
	}
	else{
		printf("we cannot step_i_fr_teth_rest (xlink)\n");
		exit(1);
	}
}

void AssociatedProteinManagement::Diffuse_II_ToBothRest(
		int x_dub, int x){

	UpdateDoubleTetheredSites();
	int n_bound = n_sites_ii_tethered_same_[x_dub][x];
	if(n_bound > 0){
		int i = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site = sites_ii_tethered_same_[x_dub][x][i];
		int i_site = site->index_;
		AssociatedProtein *xlink = site->xlink_;
		// xlink->GetDirToRest is ill-defined for x = 0, so use motor's 
		// function instead (multiply by -1 since xlink is the one moving)
		int dx = -1 * xlink->motor_->GetDirectionTowardRest();
		if(!(i_site == site->mt_->n_sites_ - 1 && dx == 1)
		&& !(i_site == 0 && dx == -1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &site->mt_->lattice_[i_site + dx];
			if(new_site->occupied_ == false){
				old_site->xlink_ = nullptr;
				old_site->occupied_ = false;
				if(old_site == xlink->site_one_)
					xlink->site_one_ = new_site;
				else
					xlink->site_two_ = new_site;
				new_site->xlink_ = xlink;
				new_site->occupied_ = true;
				// Update statistics
				xlink->UpdateExtension();
				xlink->motor_->UpdateExtension();
			}
		}
	}
	else{
		printf("issues w/ run_diffusion_ii_to_both_rest\n");
		exit(1);
	}
}

void AssociatedProteinManagement::Diffuse_II_FromBothRest(
		int x_dub, int x){

	UpdateDoubleTetheredSites();
	int n_bound = n_sites_ii_tethered_same_[x_dub][x];
	if(n_bound > 0){
		int i = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site = sites_ii_tethered_same_[x_dub][x][i];
		int i_site = site->index_;
		AssociatedProtein *xlink = site->xlink_;
		// xlink->GetDirToRest is ill-defined for x = 0, so use motor's 
		// function instead (multiply by -1 since xlink is the one moving)
		int dx = xlink->motor_->GetDirectionTowardRest();
		if(!(i_site == site->mt_->n_sites_ - 1 && dx == 1)
		&& !(i_site == 0 && dx == -1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &site->mt_->lattice_[i_site + dx];
			if(new_site->occupied_ == false){
				old_site->xlink_ = nullptr;
				old_site->occupied_ = false;
				if(old_site == xlink->site_one_)
					xlink->site_one_ = new_site;
				else
					xlink->site_two_ = new_site;
				new_site->xlink_ = xlink;
				new_site->occupied_ = true;
				// Update statistics
				xlink->UpdateExtension();
				xlink->motor_->UpdateExtension();
			}
		}
	}
	else{
		printf("issues w/ run_diffusion_ii_from_both_rest\n");
		exit(1);
	}
}

void AssociatedProteinManagement::Diffuse_II_ToSelf_FromTeth(
		int x_dub, int x){

	UpdateDoubleTetheredSites();
	int	n_bound = n_sites_ii_tethered_oppo_[x_dub][x];
	if(n_bound > 0){
		int i = properties_->gsl.GetRanInt(n_bound);
	   	Tubulin *site = sites_ii_tethered_oppo_[x_dub][x][i];
		int i_site = site->index_;
		AssociatedProtein *xlink = site->xlink_;
		// xlink->GetDirToRest is ill-defined for x = 0, so use motor's 
		// function instead (multiply by -1 since xlink is the one moving)
		int dx = xlink->motor_->GetDirectionTowardRest();
		if(!(i_site == site->mt_->n_sites_ - 1 && dx == 1)
		&& !(i_site == 0 && dx == -1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &site->mt_->lattice_[i_site + dx];
			if(new_site->occupied_ == false){
				old_site->xlink_ = nullptr;
				old_site->occupied_ = false;
				if(old_site == xlink->site_one_)
					xlink->site_one_ = new_site;
				else
					xlink->site_two_ = new_site;
				new_site->xlink_ = xlink;
				new_site->occupied_ = true;
				int x_pre = xlink->x_dist_;
				xlink->UpdateExtension();
				xlink->motor_->UpdateExtension();
			}
		}
	}
	else{
		printf("issues w/ run_diffusion_ii_toself_fromteth_rest\n");
		printf("x=%i, 2x=%i\n", x, x_dub);
		exit(1);
	}
}

void AssociatedProteinManagement::Diffuse_II_FromSelf_ToTeth(
		int x_dub, int x){

	UpdateDoubleTetheredSites();
	int n_bound = n_sites_ii_tethered_oppo_[x_dub][x];
	if(n_bound > 0){
		int i = properties_->gsl.GetRanInt(n_bound);
	   	Tubulin *site = sites_ii_tethered_oppo_[x_dub][x][i];
		int i_site = site->index_;
		AssociatedProtein *xlink = site->xlink_;
		// xlink->GetDirToRest is ill-defined for x = 0, so use motor's 
		// function instead (multiply by -1 since xlink is the one moving)
		int dx = -1 * xlink->motor_->GetDirectionTowardRest();
		if(!(i_site == site->mt_->n_sites_ - 1 && dx == 1)
		&& !(i_site == 0 && dx == -1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &site->mt_->lattice_[i_site + dx];
			if(new_site->occupied_ == false){
				old_site->xlink_ = nullptr;
				old_site->occupied_ = false;
				if(old_site == xlink->site_one_)
					xlink->site_one_ = new_site;
				else
					xlink->site_two_ = new_site;
				new_site->xlink_ = xlink;
				new_site->occupied_ = true;
				// Update statistics
				xlink->UpdateExtension();
				xlink->motor_->UpdateExtension();
			}
		}
	}
	else{
		printf("issues w/ run_diffusion_ii_fromself_toteth_rest\n");
		printf("x=%i, 2x=%i\n", x, x_dub);
		exit(1);
	}
}

void AssociatedProteinManagement::GenerateKMCList(){

	// Scan through all events & calculate number expected to occur
	sys_time start = sys_clock::now();
	UpdateAllLists();
	int n_events_tot = 0;
	for(int i_entry(0); i_entry < kmc_events_.size(); i_entry++){
		kmc_events_[i_entry].SampleStatistics();
		int n_expected = kmc_events_[i_entry].n_expected_;
		for(int i_event(0); i_event < n_expected; i_event++) 
			n_events_tot++;
	}
	sys_time finish = sys_clock::now();
	auto elapsed = std::chrono::duration_cast<t_unit>(finish - start);
	properties_->wallace.t_xlinks_kmc_[1] += elapsed.count();

	// Scan through all target populations & ensure none become negative 
	start = sys_clock::now();
	for(int i_pop(0); i_pop < kmc_by_pop_.size(); i_pop++){
		int n_entries = kmc_by_pop_[i_pop].size();
		if(n_entries > 1){
			int n_events = 0;
			double p_tot = 0;
			for(int i_entry(0); i_entry < n_entries; i_entry++){
				int index = kmc_by_pop_[i_pop][i_entry];
				n_events += kmc_events_[index].n_expected_;
				p_tot += kmc_events_[index].p_occur_; 
			}
			int n_available = *kmc_events_[kmc_by_pop_[i_pop][0]].pop_ptr_;
			while(n_events > n_available){
				double p_cum = 0;
				double ran = properties_->gsl.GetRanProb();
				for(int i_entry(0); i_entry < n_entries; i_entry++){
					int index = kmc_by_pop_[i_pop][i_entry];
					p_cum += kmc_events_[index].p_occur_ / p_tot;
					if(p_cum >= ran
					&& kmc_events_[index].n_expected_ > 0){
						kmc_events_[index].n_expected_--;
						n_events_tot--;
						n_events--;
						break;
					}
				}
			}
		}
	}
	// Construct KMC list if any events are predicted, otherwise clear it
	if(n_events_tot > 0){
		int pre_array[n_events_tot];
		int i_event = 0;
		for(int i_entry(0); i_entry < kmc_events_.size(); i_entry++){
			int n_expected = kmc_events_[i_entry].n_expected_;
			for(int i_event(0); i_event < n_expected; i_event++){
				if(i_event >= n_events_tot){
					printf("WOAH BUDDY. check XLINK genKMC\n");
					exit(1);
				}
				pre_array[i_event] = kmc_events_[i_entry].kmc_code_;
				i_event++;
			}
		}
		if(n_events_tot > 1) gsl_ran_shuffle(properties_->gsl.rng_, 
				pre_array, n_events_tot, sizeof(int));
		kmc_list_.resize(n_events_tot);
		for(int i_entry(0); i_entry < n_events_tot; i_entry++){
			kmc_list_[i_entry] = pre_array[i_entry];
		}
	}
	else kmc_list_.clear();

	finish = sys_clock::now();
	elapsed = std::chrono::duration_cast<t_unit>(finish - start);
	properties_->wallace.t_xlinks_kmc_[2] += elapsed.count();
}

double AssociatedProteinManagement::GetWeightBindII(){

	double weights_summed = 0;
	// Sum over all single-bound xlinks
	for(int n_neighbs(0); n_neighbs < 3; n_neighbs++){
		for(int i_xlink = 0; i_xlink < n_bound_i_[n_neighbs]; i_xlink++){
			AssociatedProtein *xlink = bound_i_[n_neighbs][i_xlink];
			xlink->UpdateNeighborSites();
			// Get weight of every possible orientation w/ neighbors
			int n_neighbors = xlink->n_neighbor_sites_;
			for(int i_neighb = 0; i_neighb < n_neighbors; i_neighb++){
				Tubulin *site = xlink->neighbor_sites_[i_neighb];
				double weight = xlink->GetBindingWeight(site);
				weights_summed += weight;
			}
		}
	}
	return weights_summed;
}

double AssociatedProteinManagement::GetWeightBindITethered(){

	double weights_summed = 0;
	// Sum over all free_tethered xlinks
//	printf("%i free tethered\n", n_free_tethered_);
	for(int i_xlink = 0; i_xlink < n_free_tethered_; i_xlink++){
		AssociatedProtein *xlink = free_tethered_[i_xlink];
		xlink->UpdateTethNeighborSites(); 
		// Get weight of every possible orientation with neighbors
		int n_neighbs = xlink->n_teth_neighbor_sites_; 
		for(int i_neighb = 0; i_neighb < n_neighbs; i_neighb++){
			Tubulin *site = xlink->teth_neighbor_sites_[i_neighb];
			double weight = xlink->GetTethBindingWeight(site); 
			weights_summed += weight;
		}
	}
	return weights_summed;
}

double AssociatedProteinManagement::GetWeightBindIITethered(){

	int comp_cutoff = properties_->kinesin4.comp_cutoff_;
	double weights_summed = 0;
	// Sum over all single-bound tethered xlink extensions
	for(int x_dub(2*comp_cutoff); x_dub <= 2*teth_cutoff_; x_dub++){
		int n_bound = n_bound_i_tethered_[x_dub];
		// Sum over xlinks at this specific extension
		for(int i_xlink(0); i_xlink < n_bound; i_xlink++){
			AssociatedProtein *xlink = bound_i_tethered_[x_dub][i_xlink];
			xlink->UpdateTethNeighborSitesII(); 
			int n_neighbs = xlink->n_teth_neighbor_sites_ii_;
			// Get weight of every possible orientation with neighbors
			for(int i_neighb(0); i_neighb < n_neighbs; i_neighb++){
				Tubulin *site = xlink->teth_neighbor_sites_ii_[i_neighb];
				double weight = xlink->GetTethBindingWeightII(site); 
				weights_summed += weight;
			}
		}
	}
	return weights_summed;
}

void AssociatedProteinManagement::RunKMC(){
	
//	printf("Start of xlink KMC cycle\n");
	
	sys_time start1 = sys_clock::now();

	if(parameters_->xlinks.concentration > 0) GenerateKMCList();
	else return;

	sys_time start2 = sys_clock::now();

	int x_dist;		// extension of xlink for stage2 unbinding
	int x_dub;		// twice the extension of tether for stage1&2 binding
	int n_neighbs;
	if(kmc_list_.empty() == false){
		int n_events = kmc_list_.size();
//		printf("%i XLINK KMC EVENTS\n", n_events);
		for(int i_event = 0; i_event < n_events; i_event++){
			int kmc_event = kmc_list_[i_event];
			// Allows us to encode xlink extension as second digit
			if(kmc_event >= 40 && kmc_event < 50){
				x_dist = kmc_event % 10;
				kmc_event = 40;
			}
			else if(kmc_event >= 100 && kmc_event < 200){
				n_neighbs = kmc_event % 100;
				kmc_event = 10;
			}	
			else if(kmc_event >= 200 && kmc_event < 300){
				n_neighbs = kmc_event % 100;
				kmc_event = 30;
			}
			else if(kmc_event >= 300 && kmc_event < 400){
				x_dub = kmc_event % 100;
				kmc_event = 31;	
			}
			else if(kmc_event >= 4000 && kmc_event < 5000){
				x_dist = kmc_event % 10;
				x_dub = (kmc_event % 1000 - x_dist) / 10;
				kmc_event = 41;
			}
			else if(kmc_event >= 5000 && kmc_event < 6000){
				x_dist = kmc_event % 10;
				x_dub = (kmc_event % 1000 - x_dist) / 10;
				kmc_event = 42;
			}
//			properties_->wallace.PrintMicrotubules(0.000);
			switch(kmc_event){
				case 10: 
//						printf("xlink bound stage 1\n");
						Bind_I(n_neighbs);
						break;
				case 11:
//						printf("xlink TETH bound stage 1\n");
						Bind_I_Tethered();
						break;
				case 20: 
//						printf("xlink bound stage 2\n");
						Bind_II();
						break;
				case 21:
//						printf("bind ii teth xlink\n");
						Bind_II_Tethered();
						break;
				case 30: 
//						printf("xlink fully unbound\n");
						Unbind_I(n_neighbs);
						break;
				case 31:
//						printf("unbind i teth xlink\n");
						Unbind_I_Tethered(x_dub);
						break;
				case 40: 
//						printf("xlink unbound 2nd head (ext was %i)\n", 
//								x_dist);
						Unbind_II(x_dist);
						break;
				case 41:
//						printf("unbind ii to teth XLINK\n");
						Unbind_II_To_Teth(x_dub, x_dist);
						break;
				case 42:
//						printf("unbind ii from teth XLINK\n");
						Unbind_II_From_Teth(x_dub, x_dist);
						break;
				case 50:
//						printf("xlink tethered free\n");
						Tether_Free();
						break;
				case 60:
//						printf("xlink untethered free\n");
						Untether_Free();
						break;
			}
		}
	}
	sys_time finish = sys_clock::now();
	auto elapsed = std::chrono::duration_cast<t_unit>(finish - start1);
	properties_->wallace.t_xlinks_kmc_[0] += elapsed.count();
//	printf("added %lu nanoseconds\n", elapsed.count());
	elapsed = std::chrono::duration_cast<t_unit>(finish - start2);
	properties_->wallace.t_xlinks_kmc_[3] += elapsed.count();
}

void AssociatedProteinManagement::Bind_I(int n_neighbs){

	// Make sure unoccupied sites are available
	properties_->microtubules.UpdateUnoccupied();
	if(properties_->microtubules.n_unoccupied_xl_[n_neighbs]> 0){
		// Randomly choose an unbound xlink
		AssociatedProtein *xlink = GetFreeXlink();
		if(xlink->tethered_ == true){
			printf("error in xlink bind_i\n");
			exit(1);
		}
		// Get random unoccupied site
		Tubulin *site = 
			properties_->microtubules.GetUnoccupiedSite(n_neighbs);
		// Place xlink onto site
		site->xlink_ = xlink;
		site->occupied_ = true;
		// Update xlink details
		xlink->heads_active_++;
		xlink->site_one_ = site; 
		// Update active_ list
		active_[n_active_] = xlink;
		xlink->active_index_ = n_active_;
		n_active_++;
	}
	else{
		printf("Error in xlink Bind_I: no unoccupied sites\n");
//		exit(1);
	}
}

void AssociatedProteinManagement::Bind_II(){

	// Make sure stage 1 xlinks and unoccupied sites are available
	UpdateSingleBoundList();
	properties_->microtubules.UpdateUnoccupied();
	if(n_bound_i_tot_ > 0
	&& properties_->microtubules.n_unoccupied_> 0){
		// Randomly pick single-bound xlink
		int i_xlink = properties_->gsl.GetRanInt(n_bound_i_tot_);
		AssociatedProtein *xlink = bound_i_all_[i_xlink];
		if(xlink->tethered_ == true){
			if(xlink->motor_->heads_active_ != 0){
				printf("error in xlink bind_ii_\n");
				exit(1);
			}
		}
		Tubulin *site = xlink->GetWeightedNeighborSite();
		int attempts = 0;
		while(site == nullptr){
			// roll for a new distance after a certain # of tries
			if(attempts > 10*n_bound_i_tot_){
				break;
			}
			i_xlink = properties_->gsl.GetRanInt(n_bound_i_tot_);
			xlink = bound_i_all_[i_xlink];	
			site = xlink->GetWeightedNeighborSite();
			attempts++;
		}
		if(site != nullptr){
			// Place  xlink onto site
			site->xlink_ = xlink;
			site->occupied_ = true;
			// Update xlink details
			xlink->heads_active_++;
			if(xlink->site_one_ == nullptr)
				xlink->site_one_ = site;
			else if(xlink->site_two_ == nullptr)
				xlink->site_two_ = site;
			// Update statistics
			xlink->UpdateExtension();
		}
		else{
			printf("failed to xlink bind_ii\n");
		}
	}
	else{
		printf("Error in xlink Bind_II: no unoccupied sites\n");
//		exit(1);
	}
}

void AssociatedProteinManagement::Unbind_I(int n_neighbs){

	UpdateSingleBoundList();
	if(n_bound_i_[n_neighbs] > 0){
		// Randomly pick a singly-bound xlink
		int i_entry = properties_->gsl.GetRanInt(n_bound_i_[n_neighbs]);
		AssociatedProtein *xlink = bound_i_[n_neighbs][i_entry];
		Tubulin *site = xlink->GetActiveHeadSite();
		// Remove xlink from site
		site->xlink_ = nullptr;
		site->occupied_ = false;
		// Update xlink details
		xlink->site_one_ = nullptr;
		xlink->site_two_ = nullptr;
		xlink->heads_active_--;
		// If this xlink has a satellite motor, untether it
		if(xlink->tethered_ == true){ 
			if(xlink->motor_->heads_active_ == 0){
				xlink->UntetherSatellite();
			}
			else{
				printf("error in xlink unbind_i\n");
				exit(1);
			}
		}
		// Remove this xlink from active_, replace with last entry
		AssociatedProtein *last_entry = active_[n_active_ - 1];
		int this_index = xlink->active_index_; 
		if(this_index != n_active_ - 1){
			active_[this_index] = last_entry; 
			last_entry->active_index_ = this_index; 
		}
		n_active_--;
	}
	else{
		printf("Error in Unbind: no bound xlinks\n");
//		exit(1);
	}
}

void AssociatedProteinManagement::Unbind_II(int x_dist){

	UpdateDoubleBoundList();
	if(n_bound_ii_[x_dist] > 0){
		// Randomly pick a double-bound xlink
		int i_entry = properties_->gsl.GetRanInt(n_bound_ii_[x_dist]);
		AssociatedProtein* xlink = bound_ii_[x_dist][i_entry];
		if(x_dist != xlink->x_dist_){
			printf("error in xink unbind_ii \n");
			exit(1);
		}
		if(xlink->tethered_ == true){
			if(xlink->motor_->heads_active_ != 0){
				printf("error in xlink unbind_ii NON TETH\n");
				exit(1);
			}
		}
		Tubulin* site;
		double ran = properties_->gsl.GetRanProb();
		if(ran < 0.5)
			site = xlink->site_one_;
		else
			site = xlink->site_two_;
		// Remove xlink from site
		site->xlink_ = nullptr;
		site->occupied_ = false;
		// Update xlink details
		if(ran < 0.5)
			xlink->site_one_ = nullptr;
		else
			xlink->site_two_ = nullptr;
		xlink->heads_active_--;
	}
	else{
		printf("Error in Unbind_II:no double bound xlinks\n");
//		exit(1);
	}
}

void AssociatedProteinManagement::Bind_I_Tethered(){

	UpdateFreeTetheredList();
	properties_->microtubules.UpdateUnoccupied();
	if(n_free_tethered_ > 0
	&& properties_->microtubules.n_unoccupied_> 0){
		int i_xlink = properties_->gsl.GetRanInt(n_free_tethered_);
		AssociatedProtein* xlink = free_tethered_[i_xlink];
		if(xlink->motor_->heads_active_ == 0){
			printf("error in xlink bind_i_teth\n");
			exit(1);
		}
		Tubulin* site = xlink->GetWeightedTethNeighborSite();
		if(site != nullptr){
			// Update site details
			site->xlink_ = xlink; 
			site->occupied_ = true;
			// Update xlink details
			xlink->site_one_ = site; 
			xlink->heads_active_++; 
		}
		else{
			printf("failed to XLINK Bind_I_Free_Tethered\n");
		}
	}
	else{
		printf("Error in XLINK bind_i_free_teth; no avail xlinks\n");
	}
}

void AssociatedProteinManagement::Bind_II_Tethered(){

	// Make sure stage 1 xlinks and unoccupied sites are available
	UpdateBoundITethered();
	properties_->microtubules.UpdateUnoccupied();
	if(n_bound_i_tethered_tot_ > 0
	&& properties_->microtubules.n_unoccupied_> 0){
		/* Get a weighted teth extension */
		int x_dist_dub = -1;
		// total weight summed over all extensions
		double weight_tot(0);
		// cumulative weight up to and including index extension
		double weight_cum[2*teth_cutoff_ + 1];
		// scan over all all extensions, add up weights appropriately
		for(int x_dub = 0; x_dub <= 2*teth_cutoff_; x_dub++){
			double weight = xlinks_[0].teth_binding_weight_lookup_[x_dub];
			weight_cum[x_dub] = weight * n_bound_i_tethered_[x_dub];
			weight_tot += weight_cum[x_dub];
			if(weight_cum[x_dub] > 0){
				for(int x_dub_dub = 0; x_dub_dub < x_dub; x_dub_dub++){
					weight_cum[x_dub] += weight_cum[x_dub_dub];
				}
			}
		}
		// now that we have weight_tot, we can use weight_cum[i] to
		// map the available x_dub populations onto different ranges 
		// in [0, 1), which allows us to use GetRanProb() to select one
		double ran = properties_->gsl.GetRanProb();
		for(int x_dub = 0; x_dub <= 2*teth_cutoff_; x_dub++){
			if(ran <= weight_cum[x_dub] / weight_tot){
				x_dist_dub = x_dub;
				break;
			}
		}
		if(x_dist_dub == -1){
			printf("error in xlink bind_ii_teth ZERO\n");
			exit(1);
		}
		int n_bound = n_bound_i_tethered_[x_dist_dub];
		if(n_bound > 0){
			// Randomly pick single-bound xlink
			int i_xl = properties_->gsl.GetRanInt(n_bound);
			AssociatedProtein *xlink = bound_i_tethered_[x_dist_dub][i_xl];
			if(xlink->tethered_ == false){
				printf("error in xlink bind_ii_teth ONE\n");
				exit(1);
			}
			else if(xlink->motor_->heads_active_ == 0){
				printf("error in xlink bind_ii_teth TWO\n");
				exit(1);
			}
			Tubulin *site = xlink->GetWeightedTethNeighborSiteII();
			int attempts = 0;
			int switches = 0;
			bool searching = true;
			while(site == nullptr && searching){
				// roll for a new distance after a certain # of tries
				if(attempts > 50*n_bound){
					switches++;
					// stop the search after a certain # of dist switches
					if(switches > 50*n_bound_i_tethered_tot_){
						searching = false;
						break;
					}
					double ran = properties_->gsl.GetRanProb();
					for(int x_dub = 0; x_dub <= 2*teth_cutoff_; x_dub++){
						if(ran <= weight_cum[x_dub] / weight_tot){
							x_dist_dub = x_dub;
							break;
						}
					}
					n_bound = n_bound_i_tethered_[x_dist_dub];
					attempts = 0;
				}
				i_xl = properties_->gsl.GetRanInt(n_bound);
				xlink = bound_i_tethered_[x_dist_dub][i_xl];	
				site = xlink->GetWeightedTethNeighborSiteII();
				attempts++;
			}
			if(site != nullptr){
				// Place  xlink onto site
				site->xlink_ = xlink;
				site->occupied_ = true;
				// Update xlink details
				xlink->heads_active_++;
				if(xlink->site_one_ == nullptr)
					xlink->site_one_ = site;
				else if(xlink->site_two_ == nullptr)
					xlink->site_two_ = site;
				// Update statistics
				xlink->UpdateExtension();
				xlink->motor_->UpdateExtension();
			}
			else{
				printf("failed to bind_ii_teth xlink\n");
			}
		}
		else{
			printf("error in xlink bind_ii_teth \n");
			exit(1);
		}
	}
	else{
		printf("Error in xlink Bind_II_Tethered: no unoccupied sites\n");
//		exit(1);
	}
}

void AssociatedProteinManagement::Unbind_I_Tethered(int x_dist_dub){

	UpdateBoundITethered();
	int n_bound = n_bound_i_tethered_[x_dist_dub];
	if(n_bound > 0){
		// Randomly pick a single-bound xlink
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		AssociatedProtein *xlink = bound_i_tethered_[x_dist_dub][i_entry];
		if(xlink->tethered_ == false){ 
			printf("error in xlink unbind_i_teth TWO\n");
			exit(1);
		}
		else if(xlink->motor_->heads_active_ == 0){
			printf("error in xlink unbind_i_tethhh\n");
			exit(1);
		}
		Tubulin *site = xlink->GetActiveHeadSite();
		// Remove xlink from site
		site->xlink_ = nullptr;
		site->occupied_ = false;
		// Update xlink details
		xlink->site_one_ = nullptr;
		xlink->site_two_ = nullptr;
		xlink->heads_active_--;
		// Update statistics
		xlink->motor_->UpdateExtension();
	}
	else{
		printf("Error in Unbind: no bound xlinks\n");
//		exit(1);
	}
}

void AssociatedProteinManagement::Unbind_II_To_Teth(int x_dist_dub, 
		int x_dist){

	UpdateBoundIITethered();
	int n_bound = n_bound_ii_tethered_[x_dist_dub][x_dist];
	if(n_bound > 0){
		// Randomly pick a double-bound xlink
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		AssociatedProtein* xlink = 
			bound_ii_tethered_[x_dist_dub][x_dist][i_entry];
		if(x_dist != xlink->x_dist_){
			printf("error in xlink unbind_ii_to_teth ZERO \n");
			exit(1);
		}
		if(xlink->tethered_ == false){
			printf("error in xlink unbind_ii_to_teth ONE\n");
			exit(1);
		}
		else if(xlink->motor_->heads_active_ == 0){
			printf("error in xlink unbind_ii_to_teth TWO\n");
			exit(1);
		}
		Tubulin* site = xlink->GetSiteFartherFromTethRest();;
		// Remove xlink from site
		site->xlink_ = nullptr;
		site->occupied_ = false;
		// Update xlink details
		if(xlink->site_one_ == site)
			xlink->site_one_ = nullptr;
		else
			xlink->site_two_ = nullptr;
		xlink->heads_active_--;
		// Update statistics
		xlink->UpdateExtension();
		xlink->motor_->UpdateExtension();
	}
	else{
		printf("Error in Unbind_II_To_Teth: no doubly-bound xlinks\n");
//		exit(1);
	}
}

void AssociatedProteinManagement::Unbind_II_From_Teth(
		int x_dist_dub, int x_dist){

	UpdateBoundIITethered();
	int n_bound = n_bound_ii_tethered_[x_dist_dub][x_dist];
	if(n_bound > 0){
		// Randomly pick a double-bound xlink
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		AssociatedProtein* xlink = 
			bound_ii_tethered_[x_dist_dub][x_dist][i_entry];
		if(x_dist != xlink->x_dist_){
			printf("error in xlink unbind_ii_to_teth ZERO \n");
			exit(1);
		}
		if(xlink->tethered_ == false){
			printf("error in xlink unbind_ii_to_teth ONE\n");
			exit(1);
		}
		else if(xlink->motor_->heads_active_ == 0){
			printf("error in xlink unbind_ii_to_teth TWO\n");
			exit(1);
		}
		Tubulin* site = xlink->GetSiteCloserToTethRest();;
		// Remove xlink from site
		site->xlink_ = nullptr;
		site->occupied_ = false;
		// Update xlink details
		if(xlink->site_one_ == site)
			xlink->site_one_ = nullptr;
		else
			xlink->site_two_ = nullptr;
		xlink->heads_active_--;
		// Update statistics
		xlink->UpdateExtension();
		xlink->motor_->UpdateExtension();
	}
	else{
		printf("Error in Unbind_II_Fr_Teth: no doubly-bound xlinks\n");
//		exit(1);
	}
}

void AssociatedProteinManagement::Tether_Free(){

	int n_motors_unteth = properties_->kinesin4.GetNumBoundUntethered();
	if(n_motors_unteth > 0){
		AssociatedProtein *xlink = GetFreeXlink();
		Kinesin *motor = properties_->kinesin4.GetBoundUntetheredMotor();
		// Update xlink and motor details
		xlink->motor_ = motor;
		xlink->tethered_ = true;
		motor->xlink_ = xlink;
		motor->tethered_ = true;
		// Add xlink to active_ list
		active_[n_active_] = xlink;
		xlink->active_index_ = n_active_;
		n_active_++;
	}
	else{
		printf("Error in XLINK Tether_free: no untethered bound motors\n");
	}
}

void AssociatedProteinManagement::Untether_Free(){

	UpdateFreeTetheredList();
	if(n_free_tethered_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_free_tethered_);
		AssociatedProtein* xlink = free_tethered_[i_entry];
		Kinesin* motor = xlink->motor_;
		// Update xlink and motor details
		xlink->motor_ = nullptr;
		xlink->tethered_ = false;
		motor->xlink_ = nullptr;
		motor->tethered_ = false;
		// Remove this xlink from active_, replace with last entry
		AssociatedProtein *last_entry = active_[n_active_ - 1];
		int this_index = xlink->active_index_; 
		if(this_index != n_active_ - 1){
			active_[this_index] = last_entry; 
			last_entry->active_index_ = this_index; 
		}
		n_active_--;
	}
	else{
		printf("Error in XLINK Untether_free: no free tethered xlinks\n");
	}
}
