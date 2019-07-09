#include "master_header.h"
#include "associated_protein_management.h"

AssociatedProteinManagement::AssociatedProteinManagement(){
}

void AssociatedProteinManagement::Initialize(
		system_parameters *parameters, system_properties *properties){

	parameters_ = parameters;
	properties_ = properties;
	GenerateXLinks();
	SetParameters();
	InitializeLists();
	InitializeEvents();
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

	/* 		Assume lambda = 0.5 for binding &unbinding,		  *
	 *  	whereas lambda = 1/0 for diffusing away/towards	  */
	max_neighbs_ = 2; 			// can have one in front & one behind
	interaction_energy_ = -1 * parameters_->xlinks.interaction_energy;
	double delta_t = parameters_->delta_t;
	double site_size = parameters_->microtubules.site_size;
	// DIFFUSION STATISTICS FOR SELF BELOW
	double D_coeff = parameters_->xlinks.diffusion_coeff;
	double x_squared = (site_size/1000)*(site_size/1000); // in um^2
	double tau = x_squared / (2 * D_coeff);
	p_diffuse_i_fwd_.resize(max_neighbs_ + 1);
	p_diffuse_i_bck_.resize(max_neighbs_ + 1); 
	for(int n_neighbs(0); n_neighbs <= max_neighbs_; n_neighbs++){
		double tot_E = n_neighbs * interaction_energy_;
		double weight = exp(tot_E);
		if(n_neighbs == 2) weight = 0;
		p_diffuse_i_fwd_[n_neighbs] = (delta_t / tau) * weight;
		p_diffuse_i_bck_[n_neighbs] = (delta_t / tau) * weight;
	}
	// Generate different stepping rates based on changes in
	// potential energy (dU) associated with that step
	dist_cutoff_ = xlinks_[0].dist_cutoff_;
	rest_dist_ = xlinks_[0].rest_dist_;
	teth_cutoff_ = properties_->kinesin4.dist_cutoff_; 
	comp_cutoff_ = properties_->kinesin4.comp_cutoff_;
	if(parameters_->motors.tethers_active){
		properties_->wallace.Log("\nFor crosslinkers:\n");
		properties_->wallace.Log("  rest_dist is %i\n", rest_dist_);
		properties_->wallace.Log("  dist_cutoff is %i\n\n", dist_cutoff_);
	}
	double kbT = parameters_->kbT;
	double r_0 = xlinks_[0].r_0_;
	double k_spring = xlinks_[0].k_spring_;
	double r_y = parameters_->microtubules.y_dist;
	p_diffuse_ii_to_rest_.resize(max_neighbs_ + 1);
	p_diffuse_ii_fr_rest_.resize(max_neighbs_ + 1);
	for(int n_neighbs(0); n_neighbs <= max_neighbs_; n_neighbs++){
		p_diffuse_ii_to_rest_[n_neighbs].resize(dist_cutoff_ + 1);
		p_diffuse_ii_fr_rest_[n_neighbs].resize(dist_cutoff_ + 1);
		double tot_E = n_neighbs * interaction_energy_;
		double weight_neighb = exp(tot_E);
		if(n_neighbs == max_neighbs_) weight_neighb = 0;
		for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
			double r_x = x_dist * site_size;
			double r_x_to = (x_dist - 1) * site_size;
			double r_x_fr = (x_dist + 1) * site_size;
			double r = sqrt(r_x*r_x + r_y*r_y);
			double r_to = sqrt(r_x_to*r_x_to + r_y*r_y);
			double r_fr = sqrt(r_x_fr*r_x_fr + r_y*r_y);
			// Get extension for current dist and steps to/from spring rest
			double dr = r - r_0;
			double dr_to = r_to - r_0;
			double dr_fr = r_fr - r_0;
			if(dr >= 0){
				// Get corresponding changes in potential energy
				double dU_to_rest = (k_spring/2)*(dr_to*dr_to - dr*dr);
				double dU_fr_rest = (k_spring/2)*(dr_fr*dr_fr - dr*dr);
				// Weights according to Lanksy et al.
				double weight_to = exp(-dU_to_rest/(2*kbT));
				double weight_fr = exp(-dU_fr_rest/(2*kbT));
				double p_to = weight_neighb * weight_to * delta_t / tau;
				double p_fr = weight_neighb * weight_fr * delta_t / tau;
				if(x_dist == rest_dist_){
					p_diffuse_ii_to_rest_[n_neighbs][x_dist] = 0;
					p_diffuse_ii_fr_rest_[n_neighbs][x_dist] = 2 * p_fr;
				}
				else if(x_dist == dist_cutoff_){
					p_diffuse_ii_to_rest_[n_neighbs][x_dist] = p_to;
					p_diffuse_ii_fr_rest_[n_neighbs][x_dist] = 0;
				}
				else{
					p_diffuse_ii_to_rest_[n_neighbs][x_dist] = p_to;
					p_diffuse_ii_fr_rest_[n_neighbs][x_dist] = p_fr;
				}
				if(p_to > 1)
					printf("WARNING: p_diffuse_to_rest=%g for x=%i\n", 
							p_to, x_dist);
				if(2*p_fr > 1)
					printf("WARNING: 2*p_diffuse_fr_rest=%g for x=%i\n", 
							2*p_fr, x_dist);
			}
			else{
				printf("woah mayne. xlink set parameters \n");
				exit(1);
			}
		}
	}
	// DIFFUSION STATISTICS INVOLVING TETHER BELOW
	int teth_cutoff = properties_->kinesin4.motors_[0].dist_cutoff_;
	double k_teth_spring = properties_->kinesin4.motors_[0].k_spring_;
	double k_teth_slack = properties_->kinesin4.motors_[0].k_slack_;
	double r_0_teth = properties_->kinesin4.motors_[0].r_0_;
	double r_y_teth = parameters_->microtubules.y_dist / 2;
	double rest_teth = properties_->kinesin4.motors_[0].rest_dist_;
	double r_rest_teth = 
		sqrt(site_size*rest_teth*site_size*rest_teth + r_y_teth*r_y_teth);
	p_diffuse_i_to_teth_rest_.resize(max_neighbs_ + 1);
	p_diffuse_i_fr_teth_rest_.resize(max_neighbs_ + 1);
	p_diffuse_ii_to_both_.resize(max_neighbs_ + 1);
	p_diffuse_ii_to_self_fr_teth_.resize(max_neighbs_ + 1);
	p_diffuse_ii_fr_self_to_teth_.resize(max_neighbs_ + 1);
	p_diffuse_ii_fr_both_.resize(max_neighbs_ + 1);
	for(int n_neighbs = 0; n_neighbs < max_neighbs_; n_neighbs++){
		double tot_E = n_neighbs * interaction_energy_;
		double weight_neighb = exp(tot_E);
		if(n_neighbs == 2) weight_neighb = 0;
		p_diffuse_i_to_teth_rest_[n_neighbs].resize(2*teth_cutoff + 1);
		p_diffuse_i_fr_teth_rest_[n_neighbs].resize(2*teth_cutoff + 1);
		p_diffuse_ii_to_both_[n_neighbs].resize(2*teth_cutoff + 1);
		p_diffuse_ii_to_self_fr_teth_[n_neighbs].resize(2*teth_cutoff + 1);
		p_diffuse_ii_fr_self_to_teth_[n_neighbs].resize(2*teth_cutoff + 1);
		p_diffuse_ii_fr_both_[n_neighbs].resize(2*teth_cutoff + 1);
		for(int x_dub = 0; x_dub <= 2*teth_cutoff; x_dub++){
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
			// Calc tether exts for current dist and stepping to/from rest
			double dr_teth = r_teth - r_0_teth,
				   dr_teth_to, 
				   dr_teth_fr;
			if(dr_teth >= 0){
				dr_teth_to = r_teth_bck - r_0_teth;
				dr_teth_fr = r_teth_fwd - r_0_teth;	
			}
			else{
				dr_teth_to = r_teth_fwd - r_0_teth;
				dr_teth_fr = r_teth_bck - r_0_teth;
			}
			double dU_fr_teth, 
				   dU_to_teth;
			if(x_dub == 2*rest_teth){
				if(r_0_teth > r_rest_teth){
					dU_fr_teth = (k_teth_slack/2)
						* (dr_teth_fr*dr_teth_fr - dr_teth*dr_teth);
					dU_to_teth = (0.5)*(k_teth_spring*dr_teth_to*dr_teth_to 
						- k_teth_slack*dr_teth*dr_teth);
				}
				else{
					dU_fr_teth = (k_teth_spring/2)
						* (dr_teth_fr*dr_teth_fr - dr_teth*dr_teth);
					dU_to_teth = (0.5)*(k_teth_slack*dr_teth_to*dr_teth_to
						- k_teth_spring*dr_teth*dr_teth);
				}
			}
			else if(dr_teth > 0){
				dU_fr_teth = (k_teth_spring/2)
					* (dr_teth_fr*dr_teth_fr - dr_teth*dr_teth);
				dU_to_teth = (k_teth_spring/2)
					* (dr_teth_to*dr_teth_to - dr_teth*dr_teth);
			}
			else{
				dU_fr_teth = (k_teth_slack/2)
					* (dr_teth_fr*dr_teth_fr - dr_teth*dr_teth);
				dU_to_teth = (k_teth_slack/2)
					* (dr_teth_to*dr_teth_to - dr_teth*dr_teth);
			}
			double weight_to_teth = exp(-dU_to_teth/(2*kbT));
			double weight_fr_teth = exp(-dU_fr_teth/(2*kbT));
			if(x_dub < 2*comp_cutoff_
			|| !parameters_->motors.tethers_active){
				weight_to_teth = 0;
				weight_fr_teth = 0;
			}
			if(x_dub < 2*(comp_cutoff_ + 1)
			|| x_dub > 2*(teth_cutoff - 1)){
				weight_fr_teth = 0;
			}
			double p_to_teth_i = weight_neighb*weight_to_teth*delta_t/tau;
			double p_fr_teth_i = weight_neighb*weight_fr_teth*delta_t/tau;
			// Input probabilities for stage_i / tethered xlinks
			p_diffuse_i_to_teth_rest_[n_neighbs][x_dub] = p_to_teth_i;
			p_diffuse_i_fr_teth_rest_[n_neighbs][x_dub] = p_fr_teth_i;
			if(p_to_teth_i > 1)
				printf("WARNING: p_diffuse_to_teth_i=%g for 2x=%i\n", 
						p_to_teth_i, x_dub);
			if(p_fr_teth_i > 1)
				printf("WARNING: p_diffuse_fr_teth_i=%g for 2x=%i\n", 
						p_fr_teth_i, x_dub);
			// Run through x to get probs for stage_ii tethered xlinks
			p_diffuse_ii_to_both_[n_neighbs][x_dub]
				.resize(dist_cutoff_ + 1);
			p_diffuse_ii_to_self_fr_teth_[n_neighbs][x_dub]
				.resize(dist_cutoff_ + 1);
			p_diffuse_ii_fr_self_to_teth_[n_neighbs][x_dub]
				.resize(dist_cutoff_ + 1);
			p_diffuse_ii_fr_both_[n_neighbs][x_dub]
				.resize(dist_cutoff_ + 1);
			for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
				double r_x = x_dist * site_size;
				double r_x_to = (x_dist - 1) * site_size;
				double r_x_fr = (x_dist + 1) * site_size;
				double r = sqrt(r_x*r_x + r_y*r_y);
				double r_to = sqrt(r_x_to*r_x_to + r_y*r_y);
				double r_fr = sqrt(r_x_fr*r_x_fr + r_y*r_y);
				// Get extension for current dist and steps to/from rest
				double dr = r - r_0;
				double dr_to = r_to -  r_0;
				double dr_fr = r_fr - r_0;
				if(dr >= 0){
					// Get corresponding changes in potential energy
					double dU_to_rest = (k_spring/2)*(dr_to*dr_to - dr*dr);
					double dU_fr_rest = (k_spring/2)*(dr_fr*dr_fr - dr*dr);
					// Weights according to Lanksy et al.
					double weight_to;
					double weight_fr;
					if(x_dist == rest_dist_){
						weight_to = 0;
						weight_fr = 2*exp(-dU_fr_rest/(2*kbT));
					}
					else if(x_dist == dist_cutoff_){
						weight_to = exp(-dU_to_rest/(2*kbT));
						weight_fr = 0;
					}
					else{
						weight_to = exp(-dU_to_rest/(2*kbT));
						weight_fr = exp(-dU_fr_rest/(2*kbT));
					}
					// Convolve these bitches
					double p_to_both = weight_neighb * weight_to 
						* weight_to_teth * delta_t / tau;
					double p_to_self_fr_teth = weight_neighb * weight_to 
						* weight_fr_teth * delta_t / tau;
					double p_fr_self_to_teth = weight_neighb * weight_fr 
						* weight_to_teth * delta_t / tau;
					double p_fr_both = weight_neighb * weight_fr 
						* weight_fr_teth * delta_t / tau;
					p_diffuse_ii_to_both_[n_neighbs][x_dub][x_dist]
						= p_to_both;
					p_diffuse_ii_to_self_fr_teth_[n_neighbs][x_dub][x_dist]
						= p_to_self_fr_teth;
					p_diffuse_ii_fr_self_to_teth_[n_neighbs][x_dub][x_dist]
						= p_fr_self_to_teth; 
					p_diffuse_ii_fr_both_[n_neighbs][x_dub][x_dist]
						= p_fr_both; 

					if(p_to_both > 1){
						printf("WARNING: p_diff_to_both=%g", p_to_both);	
						printf(" for 2x=%i, x=%i\n", x_dub, x_dist);
					}
					if(p_to_self_fr_teth > 1){
						printf("WARNING: p_diff_to_self_fr_teth=%g", 
								p_to_self_fr_teth);	
						printf(" for 2x=%i, x=%i\n", x_dub, x_dist);
					}
					if(p_fr_self_to_teth > 1){
						printf("WARNING: p_diff_fr_self_to_teth=%g", 
								p_fr_self_to_teth);	
						printf(" for 2x=%i, x=%i\n", x_dub, x_dist);
					}
					if(p_fr_both > 1){
						printf("WARNING: p_diff_fr_both=%g", p_fr_both);	
						printf(" for 2x=%i, x=%i\n", x_dub, x_dist);
					}
				}
				else{
					printf("woah mayne. xlink set parameters TWOO \n");
					exit(1);
				}
			}
		}
	}
	// KMC STATISTICS BELOW
	double k_on = parameters_->xlinks.k_on; 
	double c_xlink = parameters_->xlinks.c_bulk;
	double c_eff_teth = parameters_->motors.c_eff_tether;
	if(!parameters_->motors.tethers_active){
		c_eff_teth = 0;
	}
	p_bind_i_teth_base_ = k_on * c_eff_teth * delta_t; 
	double c_eff_bind = parameters_->xlinks.c_eff_bind;
	p_bind_ii_base_ = k_on * c_eff_bind * delta_t;
	double k_off = parameters_->xlinks.k_off;
	p_bind_i_.resize(max_neighbs_ + 1);
	p_unbind_i_.resize(max_neighbs_ + 1);
	p_unbind_ii_.resize(max_neighbs_ + 1);
	p_unbind_i_teth_.resize(max_neighbs_ + 1);
	p_unbind_ii_to_teth_.resize(max_neighbs_ + 1);
	p_unbind_ii_fr_teth_.resize(max_neighbs_ + 1);
	for(int n_neighbs(0); n_neighbs <= max_neighbs_; n_neighbs++){
		double tot_E = n_neighbs * interaction_energy_;
		double weight_neighb = exp(tot_E/2);
		p_bind_i_[n_neighbs] = (k_on * c_xlink * delta_t) / weight_neighb;
		p_unbind_i_[n_neighbs] = weight_neighb * k_off * delta_t;
		// Rates involving xlink spring only
		p_unbind_ii_[n_neighbs].resize(dist_cutoff_ + 1);
		for(int x = 0; x <= dist_cutoff_; x++){
			double r_x = x*site_size;
			double r = sqrt(r_y*r_y + r_x*r_x);
			double dr = r - r_0;
			double U_xlink = (k_spring/2)*dr*dr;
			double unbind_weight = exp(U_xlink/(2*kbT));
			if(x == rest_dist_) unbind_weight = 1;
			p_unbind_ii_[n_neighbs][x] = weight_neighb * unbind_weight 
				* k_off * delta_t;
		}
		// Rates involving both xlink & tether spring
		p_unbind_i_teth_[n_neighbs].resize(2*teth_cutoff + 1);
		p_unbind_ii_to_teth_[n_neighbs].resize(2*teth_cutoff + 1);
		p_unbind_ii_fr_teth_[n_neighbs].resize(2*teth_cutoff + 1);
		for(int x_dub = 0; x_dub <= 2*teth_cutoff; x_dub++){
			// Calc x-distances (in nm) for tether
			double r_x_teth = x_dub * site_size / 2;
			// Calc total r values 
			double r_teth = sqrt(r_x_teth*r_x_teth + r_y_teth*r_y_teth);
			// Calc tether exts for current dist and stepping to/from rest
			double dr_teth = r_teth - r_0_teth;	
			double U_at_teth;
			if(dr_teth > 0)
				U_at_teth = (k_teth_spring/2)*dr_teth*dr_teth;
			else
				U_at_teth = (k_teth_slack/2)*dr_teth*dr_teth;
			double weight_at_teth = exp(U_at_teth/(2*kbT));
			if(x_dub < 2*comp_cutoff_
			|| !parameters_->motors.tethers_active){
				weight_at_teth = 0;
			}
			p_unbind_i_teth_[n_neighbs][x_dub] = weight_neighb 
				* weight_at_teth * k_off * delta_t; 
			if(p_unbind_i_teth_[n_neighbs][x_dub] > 1)
				printf("WARNING: p_unbind_teth (XLINK)=%g for 2x=%i\n", 
						p_unbind_i_teth_[n_neighbs][x_dub], x_dub);
			// Run through x to get probs for stage_ii / tethered xlinks
			p_unbind_ii_to_teth_[n_neighbs][x_dub].resize(dist_cutoff_+1);
			p_unbind_ii_fr_teth_[n_neighbs][x_dub].resize(dist_cutoff_+1);
			for(int x = 0; x <= dist_cutoff_; x++){
				// change in teth extension if 2nd xlink head were to unbind
				double dx_teth = (double)x/2; 
				double dr_x_teth = dx_teth * site_size;
				double r_x_teth_to, 
					   r_x_teth_fr; 
				if(dr_teth >= 0){
					r_x_teth_to = r_x_teth - dr_x_teth;
					r_x_teth_fr = r_x_teth + dr_x_teth; 	
				}
				else{
					r_x_teth_to = r_x_teth + dr_x_teth;
					r_x_teth_fr = r_x_teth - dr_x_teth; 
				}
				double r_teth_to = 
					sqrt(r_x_teth_to*r_x_teth_to + r_y_teth*r_y_teth);
				double r_teth_fr = 
					sqrt(r_x_teth_fr*r_x_teth_fr + r_y_teth*r_y_teth);
				double dr_teth_to = r_teth_to - r_0_teth; 
				double dr_teth_fr = r_teth_fr - r_0_teth; 
				double dU_to_teth, 
					   dU_fr_teth;
				int x_fr_rest_dub = abs(2*rest_teth - x_dub);
				int dx_teth_dub = 2 * dx_teth;
				// Check if we're crossing over equil. point of tether 
				if(dx_teth_dub > x_fr_rest_dub){
					if(r_teth < r_0_teth){
						dU_fr_teth = (k_teth_slack/2)
							* (dr_teth_fr*dr_teth_fr - dr_teth*dr_teth);
						dU_to_teth = (0.5)
							*(k_teth_spring*dr_teth_to*dr_teth_to 
							- k_teth_slack*dr_teth*dr_teth);
					}
					else{
						dU_fr_teth = (k_teth_spring/2)
							* (dr_teth_fr*dr_teth_fr - dr_teth*dr_teth);
						dU_to_teth = (0.5)
							*(k_teth_slack*dr_teth_to*dr_teth_to
							- k_teth_spring*dr_teth*dr_teth);
					}
				}
				else if(dr_teth > 0){
					dU_fr_teth = (k_teth_spring/2)
						* (dr_teth_fr*dr_teth_fr - dr_teth*dr_teth);
					dU_to_teth = (k_teth_spring/2)
						* (dr_teth_to*dr_teth_to - dr_teth*dr_teth);
				}
				else{
					dU_fr_teth = (k_teth_slack/2)
						* (dr_teth_fr*dr_teth_fr - dr_teth*dr_teth);
					dU_to_teth = (k_teth_slack/2)
						* (dr_teth_to*dr_teth_to - dr_teth*dr_teth);
				}
				double weight_to_teth = exp(-dU_to_teth/(2*kbT));
				double weight_fr_teth = exp(-dU_fr_teth/(2*kbT));
				// To ensure we don't double-count sites for x = 0, 
				// only use unbind_ii_fr_teth and ignore to_teth
				if(x == rest_dist_){
					weight_to_teth = 0;
					weight_fr_teth = 1;
				}
				if(x_dub < 2*comp_cutoff_
				|| !parameters_->motors.tethers_active){
					weight_to_teth = 0;
					weight_fr_teth = 0;
				}
				if(x_dub < 2*(comp_cutoff_ + 1)
				|| x_dub > 2*(teth_cutoff - 1)){
					weight_fr_teth = 0;
				}
				p_unbind_ii_to_teth_[n_neighbs][x_dub][x] = weight_neighb
					* weight_to_teth * k_off * delta_t;
				p_unbind_ii_fr_teth_[n_neighbs][x_dub][x] = weight_neighb
					* weight_fr_teth * k_off * delta_t;
				if(p_unbind_ii_to_teth_[n_neighbs][x_dub][x] > 1)
					printf("WARNING: p_unbind_to = %g for 2x=%ix, x=%i\n", 
							p_unbind_ii_to_teth_[n_neighbs][x_dub][x], 
							x_dub, x);
				if(p_unbind_ii_fr_teth_[n_neighbs][x_dub][x] > 1)
					printf("WARNING: p_unbind_fr = %g for 2x=%ix, x=%i\n", 
							p_unbind_ii_fr_teth_[n_neighbs][x_dub][x], 
							x_dub, x);
			}
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

	n_scratch_.resize(max_neighbs_ + 1);
	scratch_.resize(max_neighbs_ + 1);
	for(int n_neighbs(0); n_neighbs <= max_neighbs_; n_neighbs++){
		n_scratch_[n_neighbs] = 0;
		scratch_[n_neighbs].resize(n_xlinks_);
	}
	// Stats (not a list ok bite me)
	n_heads_i_.resize(max_neighbs_ + 2);
	n_heads_ii_.resize(max_neighbs_ + 2);
	n_heads_i_teth_.resize(max_neighbs_ + 2);
	n_heads_ii_teth_same_.resize(max_neighbs_ + 2);
	n_heads_ii_teth_oppo_.resize(max_neighbs_ + 2);
	// Final entry, [max_neighbs_+1], holds ALL stats regardless of neighbs
	for(int n_neighbs(0); n_neighbs <= max_neighbs_ + 1; n_neighbs++){
		n_heads_i_[n_neighbs] = 0;
		n_heads_ii_[n_neighbs].resize(dist_cutoff_ + 1);
		for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++)
			n_heads_ii_[n_neighbs][x_dist] = 0;
		n_heads_i_teth_[n_neighbs].resize(2*teth_cutoff_ + 1);
		n_heads_ii_teth_same_[n_neighbs].resize(2*teth_cutoff_ + 1);
		n_heads_ii_teth_oppo_[n_neighbs].resize(2*teth_cutoff_ + 1);
		for(int x_dub = 0; x_dub <= 2*teth_cutoff_; x_dub++){
			n_heads_i_teth_[n_neighbs][x_dub] = 0;
			n_heads_ii_teth_same_[n_neighbs][x_dub].resize(dist_cutoff_+1);
			n_heads_ii_teth_oppo_[n_neighbs][x_dub].resize(dist_cutoff_+1);
			for(int x = 0; x <= dist_cutoff_; x++){
				n_heads_ii_teth_same_[n_neighbs][x_dub][x] = 0;
				n_heads_ii_teth_oppo_[n_neighbs][x_dub][x] = 0; 
			}
		}
	}
	// Lists
	active_.resize(n_xlinks_);
	free_teth_.resize(n_xlinks_);
	bound_unteth_.resize(n_xlinks_);
	heads_i_.resize(max_neighbs_ + 2);
	heads_ii_.resize(max_neighbs_ + 2);
	heads_i_teth_.resize(max_neighbs_ + 2);
	heads_ii_teth_oppo_.resize(max_neighbs_ + 2);
	heads_ii_teth_same_.resize(max_neighbs_ + 2);
	// Final entry, [max_neighbs_+1], holds ALL stats regardless of neighbs
	for(int n_neighbs(0); n_neighbs <= max_neighbs_ + 1; n_neighbs++){
		heads_i_[n_neighbs].resize(n_xlinks_);
		heads_ii_[n_neighbs].resize(dist_cutoff_ + 1);
		for(int x = 0; x <= dist_cutoff_; x++)
			heads_ii_[n_neighbs][x].resize(n_xlinks_);
		heads_i_teth_[n_neighbs].resize(2*teth_cutoff_ + 1);
		heads_ii_teth_oppo_[n_neighbs].resize(2*teth_cutoff_ + 1);
		heads_ii_teth_same_[n_neighbs].resize(2*teth_cutoff_ + 1);
		for(int x_dub = 0; x_dub <= 2*teth_cutoff_; x_dub++){
			heads_i_teth_[n_neighbs][x_dub].resize(n_xlinks_);
			heads_ii_teth_oppo_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
			heads_ii_teth_same_[n_neighbs][x_dub].resize(dist_cutoff_ + 1);
			for(int x = 0; x <= dist_cutoff_; x++){
				heads_ii_teth_oppo_[n_neighbs][x_dub][x].resize(n_xlinks_);
				heads_ii_teth_same_[n_neighbs][x_dub][x].resize(n_xlinks_);
			}
		}
	}
}

void AssociatedProteinManagement::InitializeEvents(){

	/* *** Serialized & unique index of each KMC event *** */
	int ID(0);
	/* *** Basic random integer generator ** */
	auto ran_int = [&](int n){
		if(n > 0) return properties_->gsl.GetRanInt(n);
		else return 0;
	};
	/* *** Probability distributions *** */
	// Baseline probability dist. is Binomial since sim is discretized
	auto binomial = [&](double p, int n){
		if(n > 0) return properties_->gsl.SampleBinomialDist(p, n);
		else return 0;
	};
	// Poisson dist. is used w/ partition function for E-dependent binding
	auto poisson_ii = [&](double p, int n){
		if(n > 0){ 
			double n_wt = GetWeight_Bind_II();
			if(n_wt > 0) 
				return properties_->gsl.SamplePoissonDist(p*n_wt);
			else return 0;
		}
		else return 0;
	};
	auto poisson_i_teth = [&](double p, int n){
		if(n > 0){
			double n_wt = GetWeight_Bind_I_Teth();
			if(n_wt > 0) 
				return properties_->gsl.SamplePoissonDist(p*n_wt);
			else return 0;
		}
		else return 0;
	};
	auto poisson_ii_teth = [&](double p, int n){
		if(n > 0){
			double n_wt = GetWeight_Bind_II_Teth();
			if(n_wt > 0)
				return properties_->gsl.SamplePoissonDist(p*n_wt);
			else return 0;
		}
		else return 0;
	};
	/* ***    Actual instantiation of event objects; format is:    *** */
	//		  (*management, ID, kmc_code, 'event_name', 'target_pop', 
	//		    p_occur, *n_avail, *pop_pool, ran_int, p_dist)
	for(int n_neighbs(0); n_neighbs <= max_neighbs_; n_neighbs++){
		std::string N = std::to_string(n_neighbs);
		events_.emplace_back(this, ID++, -10, "diff_i_fwd", "bound_i_" + N, 
				p_diffuse_i_fwd_[n_neighbs], &n_heads_i_[n_neighbs], 
				heads_i_[n_neighbs], ran_int, binomial);
		events_.emplace_back(this, ID++, -11, "diff_i_bck", "bound_i_" + N, 
				p_diffuse_i_bck_[n_neighbs], &n_heads_i_[n_neighbs], 
				heads_i_[n_neighbs], ran_int, binomial);
		events_.emplace_back(this, ID++, 10, "bind_i", "unocc_" + N, 
				p_bind_i_[n_neighbs], 
				&properties_->microtubules.n_unoccupied_xl_[n_neighbs], 
				&properties_->microtubules.unoccupied_list_xl_[n_neighbs], 
				ran_int, binomial);
		if(n_neighbs == 0 && parameters_->microtubules.count > 1){
			events_.emplace_back(this, ID++, 20, "bind_ii", "bound_i_ALL", 
					p_bind_ii_base_, &n_heads_i_[max_neighbs_+1], 
					heads_i_[max_neighbs_+1], ran_int, poisson_ii);
		}
		events_.emplace_back(this, ID++, 30, "unbind_i", "bound_i_" + N, 
				p_unbind_i_[n_neighbs], &n_heads_i_[n_neighbs], 
				heads_i_[n_neighbs], ran_int, binomial);
		/*
		// Only create doubly-bound events if n_MTs > 1
		if(parameters_->microtubules.count > 1){
			for(int x(0); x <= dist_cutoff_; x++){
				events_.emplace_back(index++, -(20000 + x), n_neighbs,
						"diff_ii_to", "*heads_ii_" 
						+ std::to_string(n_neighbs) + "_" + 
						std::to_string(x), binomial, 
						&n_heads_ii_[n_neighbs][x], 
						p_diffuse_ii_to_rest_[n_neighbs][x]);
				events_.emplace_back(index++, -(21000 + x), n_neighbs, 
						"diff_ii_fr", "*heads_ii_" 
						+ std::to_string(n_neighbs) + "_" 
						+ std::to_string(x), binomial, 
						&n_heads_ii_[n_neighbs][x], 
						p_diffuse_ii_fr_rest_[n_neighbs][x]);
				events_.emplace_back(index++, 40000 + x, n_neighbs, 
						"unbind_ii", "*heads_ii_" 
						+ std::to_string(n_neighbs) + "_" 
						+ std::to_string(x), binomial, 
						&n_heads_ii_[n_neighbs][x], 
						p_unbind_ii_[n_neighbs][x]);
			}
		}
		*/
	}
	// If tethering is enabled, add those event populations as well
	if(parameters_->motors.tethers_active){
		/*
		events_.emplace_back(index++, 11, -1, 
				"bind_i_teth", "free_teth", poisson_bind_i_teth, 
				&n_free_teth_, p_bind_i_teth_base_);
		events_.emplace_back(index++, 21, -1, 
				"bind_ii_teth", "heads_i_teth_ALL", poisson_bind_ii_teth, 
				&n_heads_i_teth_tot_, p_bind_ii_base_);
		for(int n_neighbs(0); n_neighbs <= max_neighbs_; n_neighbs++){
			for(int x_dub(2*comp_cutoff_); x_dub<=2*teth_cutoff_; x_dub++){
				events_.emplace_back(index++, -(30000 + 10*x_dub),
						n_neighbs, "diff_i_to_teth", "heads_i_teth_" 
						+ std::to_string(x_dub) + "_"
						+ std::to_string(n_neighbs), binomial, 
						&n_heads_i_teth_[n_neighbs][x_dub], 
						p_diffuse_i_to_teth_rest_[n_neighbs][x_dub]);
				events_.emplace_back(index++, -(31000 + 10*x_dub), 
						n_neighbs, "diff_i_fr_teth", "heads_i_teth_" 
						+ std::to_string(x_dub) + "_"
						+ std::to_string(n_neighbs), binomial, 
						&n_heads_i_teth_[n_neighbs][x_dub], 
						p_diffuse_i_fr_teth_rest_[n_neighbs][x_dub]);
				events_.emplace_back(index++, 31000 + 10*x_dub, 
						n_neighbs, "unbind_i_teth", "heads_i_teth_" 
						+ std::to_string(x_dub) + "_" 
						+ std::to_string(n_neighbs), binomial, 
						&n_heads_i_teth_[n_neighbs][x_dub], 
						p_unbind_i_teth_[n_neighbs][x_dub]);
				// Only create doubly-bound events if n_MTs > 1
				if(parameters_->microtubules.count > 1){
					for(int x(0); x <= dist_cutoff_; x++){
						events_.emplace_back(index++, -(40000+10*x_dub+x), 
								n_neighbs, "diff_ii_to_both", 
								"heads_ii_teth_same_" 
								+ std::to_string(x_dub) + "_" 
								+ std::to_string(x) + "_"
								+ std::to_string(n_neighbs), binomial, 
								&n_heads_ii_teth_same_[n_neighbs][x_dub][x],
								p_diffuse_ii_to_both_[n_neighbs][x_dub][x]);
						events_.emplace_back(index++, -(41000+10*x_dub+x), 
								n_neighbs, "diff_ii_fr_both", 
								"heads_ii_teth_same_" 
								+ std::to_string(x_dub) + "_" 
								+ std::to_string(x) + "_"
								+ std::to_string(n_neighbs), binomial, 
								&n_heads_ii_teth_same_[n_neighbs][x_dub][x],
								p_diffuse_ii_fr_both_[n_neighbs][x_dub][x]);
						events_.emplace_back(index++, -(50000+10*x_dub+x), 
								n_neighbs, "diff_ii_to_fr", 
								"heads_ii_teth_oppo_" 
								+ std::to_string(x_dub) + "_" 
								+ std::to_string(x) + "_"
								+ std::to_string(n_neighbs), binomial, 
								&n_heads_ii_teth_oppo_[n_neighbs][x_dub][x],
								p_diffuse_ii_to_self_fr_teth_
								[n_neighbs][x_dub][x]);
						events_.emplace_back(index++, -(51000+10*x_dub+x), 
								n_neighbs, "diff_ii_fr_to", 
								"heads_ii_teth_oppo_" 
								+ std::to_string(x_dub) + "_" 
								+ std::to_string(x) + "_"
								+ std::to_string(n_neighbs), binomial, 
								&n_heads_ii_teth_oppo_[n_neighbs][x_dub][x],
								p_diffuse_ii_fr_self_to_teth_
								[n_neighbs][x_dub][x]);
						events_.emplace_back(index++, 41000+10*x_dub+x, 
								n_neighbs, "unbind_ii_to_teth", 
								"heads_ii_teth_same_" 
								+ std::to_string(x_dub) + "_" 
								+ std::to_string(x) + "_" 
								+ std::to_string(n_neighbs), binomial, 
								&n_heads_ii_teth_same_[n_neighbs][x_dub][x],
								p_unbind_ii_to_teth_[n_neighbs][x_dub][x]);
						events_.emplace_back(index++, 42000 + 10*x_dub + x,
								n_neighbs, "unbind_ii_fr_teth", 
								"heads_ii_teth_oppo_"
							   	+ std::to_string(x_dub) + "_" 
								+ std::to_string(x) + "_"
								+ std::to_string(n_neighbs), binomial, 
								&n_heads_ii_teth_oppo_[n_neighbs][x_dub][x],
								p_unbind_ii_fr_teth_[n_neighbs][x_dub][x]);
					}
				}
			}
		}
		events_.emplace_back(index++, 50, -1, "tether_free", "unteth_mots",
				binomial, &properties_->kinesin4.n_bound_untethered_, 
				p_tether_free_);
		events_.emplace_back(index++, 60, -1, "untether_free", "free_teth", 
				binomial, &n_free_teth_, p_untether_free_);
		*/
	}
	/* ** Segregate events_ into IDs_by_pop_ ** */
	int n_pops = 0;					// Total no. of distinct populations
	int n_entries[events_.size()];	// No. of entries for each pop.
	// Index of each entry for each population:
	int entries[events_.size()][2*(2*teth_cutoff_+1)*(dist_cutoff_+1)];
	// Certain events have their stats corrected SECONDARY to others,
	// e.g., events that affect all extensions will be corrected after
	// events that affect a specific extension; 'sec' tracks these 
	// 'root' refers to the root pop., e.g., teth_i w/o extension info
	int n_sec_roots = 0;
	int i_sec_roots[events_.size()];
	std::string root_pops[events_.size()];
	for(int i_entry(0); i_entry < events_.size(); i_entry++){
		std::string tar = events_[i_entry].target_pop_;
		// '*' flag at beginning -> both primary and secondary scan
		// (first over x & n_neighbs, then over x for all n_neighbs)
		if(tar.substr(0, 1) == "*"){
			std::cout << tar;
			printf(" is secondary!\n");
			std::string base = tar.substr(0, 10);
			std::string x = tar.substr(tar.length() - 1);
			std::string root = base + x;
			std::cout << root;
			printf(" is its root!\n");
			bool counted = false;
			for(int i_root(0); i_root < n_sec_roots; i_root++){
				std::string pop = events_[i_sec_roots[i_root]].target_pop_;
				if(pop.substr(0, 1) == "*" 
				&& pop.substr(pop.length() - 1) == x)
					counted = true; 
			}    
			if(!counted){
				root_pops[n_sec_roots] = root; 
				i_sec_roots[n_sec_roots] = events_[i_entry].ID_;
				n_sec_roots++;
			}
		}
		// 'ALL' flag at end -> secondary scan only
		if(tar.substr(tar.length() - 3) == "ALL"){
			std::cout << tar;
			printf(" is secondary!\n");
			std::string root = tar.substr(0, tar.length() - 4);
			std::cout << root;
			printf(" is its root!\n");
			root_pops[n_sec_roots] = root; 
			i_sec_roots[n_sec_roots] = events_[i_entry].ID_;
			n_sec_roots++;
		}
		// No 'ALL' flag -> primary scan; compare to recorded pops.
		else{
			bool new_pop(true); 	// Assume population type is new
			int pop_index(0); 		// Index of population in n_entries
			for(int i_pop(0); i_pop < n_pops; i_pop++){
				std::string rec = events_[entries[i_pop][0]].target_pop_;
				// If type matches a recorded type; pop. isn't new
				if(tar == rec){
					new_pop = false;
					pop_index = i_pop;
				}
			}
			// If indeed a new pop., record it in a new row
			if(new_pop){
				entries[n_pops][0] = events_[i_entry].ID_;
				n_entries[n_pops] = 1; 
				n_pops++;
			}
			// Otherwise, add entry to row of already-seen pop.
			else{
				int entry_no = n_entries[pop_index];
				entries[pop_index][entry_no] = events_[i_entry].ID_;
				n_entries[pop_index]++;
			}
		}
	}
	// Scan through primary pops. & place them into IDs_by_pop_
	IDs_by_pop_.resize(n_pops + n_sec_roots);
	// As we transfer primary pops; tally up all secondary events
	int n_sec_entries[n_sec_roots];
	for(int i_root(0); i_root < n_sec_roots; i_root++)
		n_sec_entries[i_root] = 0;
	int sec_entries[n_sec_roots][2*(2*teth_cutoff_+1)*(dist_cutoff_+1)];
	// 1st index of IDs_by_pop_ corresponds to population type
	for(int i_pop(0); i_pop < n_pops; i_pop++){
		std::string target = events_[entries[i_pop][0]].target_pop_;
		printf("KMC Pop #%i - ", i_pop);
		std::cout << target;
		printf(":");
		// 2nd index corresponds to event entries which target that pop.
		IDs_by_pop_[i_pop].resize(n_entries[i_pop]);
		for(int i_entry(0); i_entry < n_entries[i_pop]; i_entry++){
			IDs_by_pop_[i_pop][i_entry] = entries[i_pop][i_entry];
			printf(" %i,", entries[i_pop][i_entry]);
		}
		printf("\n");
		// Check to see if pop. root matches any secondary-scan pops. 
		for(int i_root(0); i_root < n_sec_roots; i_root++){
			std::string pop_root = root_pops[i_root];
			std::string tar_root = target.substr(0, pop_root.length());
			if(target.substr(0, 1) == "*"){
				std::string base = target.substr(0, 10);
				std::string x = target.substr(target.length() - 1);
				tar_root = base + x;
			}
			// If it does, add all entries of this pop. to sec_entries
			if(tar_root == pop_root){
				for(int i_entry(0); i_entry < n_entries[i_pop]; i_entry++){
//					printf("%i\n", n_sec_entries[i_root]);
					sec_entries[i_root][n_sec_entries[i_root]] 
						= entries[i_pop][i_entry];
					n_sec_entries[i_root]++;
				}
			}
		}
	}
	// Scan through secondary pops. & place them into IDs_by_pop_
	for(int i_root(0); i_root < n_sec_roots; i_root++){
		std::string label = events_[i_sec_roots[i_root]].target_pop_;
		printf("KMC SEC ROOT #%i - ", i_root);
		std::cout << label;
		printf(":");
		bool coupled = false;
		if(label.substr(0, 1) == "*") coupled = true;
		int i_pop = n_pops + i_root;
		IDs_by_pop_[i_pop].resize(n_sec_entries[i_root]);
		for(int i_entry(0); i_entry < n_sec_entries[i_root]; i_entry++){
			IDs_by_pop_[i_pop][i_entry] = sec_entries[i_root][i_entry];
			// Make coupled pops. (e.g., sites_ii) neg. to distinguish
			if(coupled) IDs_by_pop_[i_pop][i_entry] *= -1;
			printf(" %i,", IDs_by_pop_[i_pop][i_entry]);
		}
		printf("\n");
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

	Update_Bound_Unteth();
	if(n_bound_unteth_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound_unteth_);
		AssociatedProtein* xlink = bound_unteth_[i_entry];
		return xlink;
	}
	else{
		printf("Error in GetUnTetheredXlink: no untethered xlinks!\n");
		exit(1);
	}
}

double AssociatedProteinManagement::GetWeight_Bind_II(){

	double weights_summed = 0;
	// Sum over all single-bound xlinks
	for(int i_xlink = 0; i_xlink < n_heads_i_[max_neighbs_+1]; i_xlink++){
		AssociatedProtein *xlink 
			= heads_i_[max_neighbs_+1][i_xlink]->xlink_;
		xlink->UpdateNeighborSites();
		// Get weight of every possible orientation w/ neighbors
		int n_neighbors = xlink->n_neighbor_sites_;
		for(int i_neighb = 0; i_neighb < n_neighbors; i_neighb++){
			Tubulin *site = xlink->neighbor_sites_[i_neighb];
			double weight = xlink->GetBindingWeight(site);
			weights_summed += weight;
		}
	}
	return weights_summed;
}

double AssociatedProteinManagement::GetWeight_Bind_I_Teth(){

	double weights_summed = 0;
	// Sum over all free_tethered xlinks
//	printf("%i free tethered\n", n_free_tethered_);
	for(int i_xlink = 0; i_xlink < n_free_teth_; i_xlink++){
		AssociatedProtein *xlink = free_teth_[i_xlink];
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

double AssociatedProteinManagement::GetWeight_Bind_II_Teth(){

	double weights_summed = 0;
	// Sum over all single-bound tethered xlink extensions
	for(int x_dub(2*comp_cutoff_); x_dub <= 2*teth_cutoff_; x_dub++){
		int n_bound = n_heads_i_teth_[max_neighbs_+1][x_dub];
		// Sum over xlinks at this specific extension
		for(int i_xlink(0); i_xlink < n_bound; i_xlink++){
			AssociatedProtein *xlink = heads_i_teth_
				[max_neighbs_+1][x_dub][i_xlink]->xlink_;
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

void AssociatedProteinManagement::Update_All(){
	
	Update_Heads_I();
	Update_Heads_II();
	properties_->microtubules.UpdateUnoccupied();
	if(parameters_->motors.tethers_active){
		Update_Free_Teth();
		Update_Bound_Unteth();
		Update_Heads_I_Teth();
		Update_Heads_II_Teth();
		properties_->kinesin4.UpdateBoundUntethered();
	}
}

void AssociatedProteinManagement::Update_Free_Teth(){

	n_free_teth_ = 0;
	for(int i_xlink = 0; i_xlink < n_active_; i_xlink++){
		AssociatedProtein *xlink = active_[i_xlink];
		if(xlink->heads_active_ == 0
		&& xlink->tethered_ == true){
			if(xlink->motor_->heads_active_ > 0){
				free_teth_[n_free_teth_] = xlink;
				n_free_teth_++;
			}
			else{
				printf("woah. error in update_free_teth_list (XLINKS)\n");
				exit(1);
			}
		}
	}
}

void AssociatedProteinManagement::Update_Heads_I(){
	
	for(int n_neighbs = 0; n_neighbs <= max_neighbs_ + 1; n_neighbs++){
		n_heads_i_[n_neighbs] = 0;
	}
	for(int i_xlink = 0; i_xlink < n_active_; i_xlink++){
		AssociatedProtein *xlink = active_[i_xlink];
		if(xlink->heads_active_ == 1
		&& xlink->tethered_ == false){
			AssociatedProtein::Monomer* head = xlink->GetActiveHead();
			int n_neighbs = head->site_->GetPRC1NeighborCount();
			heads_i_[n_neighbs][n_heads_i_[n_neighbs]] = head;
			n_heads_i_[n_neighbs]++;
			heads_i_[max_neighbs_+1][n_heads_i_[max_neighbs_+1]] = head;
			n_heads_i_[max_neighbs_+1]++;
		}
		else if(xlink->tethered_ == true){
			if(xlink->heads_active_ == 1
			&& xlink->motor_->heads_active_ == 0){
				AssociatedProtein::Monomer* head = xlink->GetActiveHead();
				int n_neighbs = head->site_->GetPRC1NeighborCount();
				heads_i_[n_neighbs][n_heads_i_[n_neighbs]] = head;
				n_heads_i_[n_neighbs]++;
				heads_i_[max_neighbs_+1][n_heads_i_[max_neighbs_+1]] 
					= head;
				n_heads_i_[max_neighbs_+1]++;
			}
		}
	}
}

void AssociatedProteinManagement::Update_Heads_I_Teth(){
	
	n_heads_i_teth_tot_ = 0;
	for(int n_neighbs(0); n_neighbs <= max_neighbs_ + 1; n_neighbs++){
		for(int x_dub(2*comp_cutoff_); x_dub <= 2*teth_cutoff_; x_dub++){
			n_heads_i_teth_[n_neighbs][x_dub] = 0;
		}
	}
	for(int i_xlink = 0; i_xlink < n_active_; i_xlink++){
		AssociatedProtein *xlink = active_[i_xlink];
		if(xlink->heads_active_ == 1
		&& xlink->tethered_){
			if(xlink->motor_->heads_active_ > 0){
				xlink->motor_->UpdateExtension(); 
				// Make sure we didn't force an untether event
				if(xlink->tethered_){
					int x_dub = xlink->motor_->x_dist_doubled_;
					AssociatedProtein::Monomer* head 
						= xlink->GetActiveHead();
					int n_neighbs = head->site_->GetPRC1NeighborCount();
//					printf("FOUND %i NEBS\n", n_neighbs);
					int index = n_heads_i_teth_[n_neighbs][x_dub];
					heads_i_teth_[n_neighbs][x_dub][index] = head;
					n_heads_i_teth_[n_neighbs][x_dub]++; 
					heads_i_teth_[max_neighbs_+1][x_dub][index] = head;
					n_heads_i_teth_[max_neighbs_+1][x_dub]++; 
					n_heads_i_teth_tot_++;
				}
			}
		}
	}
}

void AssociatedProteinManagement::Update_Heads_II(){

	for(int n_neighbs(0); n_neighbs <= max_neighbs_ + 1; n_neighbs++){
		for(int x = 0; x <= dist_cutoff_; x++){
			n_heads_ii_[n_neighbs][x] = 0;
		}
	}
	for(int i_xlink = 0; i_xlink < n_active_; i_xlink++){
		AssociatedProtein *xlink = active_[i_xlink];
		if(xlink->heads_active_ == 2
		&& xlink->tethered_ == false){
			xlink->UpdateExtension();
			// Make sure we didn't force an unbind event
			if(xlink->heads_active_ == 2){
				int x = xlink->x_dist_;
				int neighbs_one 
					= xlink->head_one_.site_->GetPRC1NeighborCount();
				int i_one = n_heads_ii_[neighbs_one][x];
				heads_ii_[neighbs_one][x][i_one] = &xlink->head_one_;
				n_heads_ii_[neighbs_one][x]++;
				int neighbs_two 
					= xlink->head_two_.site_->GetPRC1NeighborCount();
				int i_two = n_heads_ii_[neighbs_two][x];
				heads_ii_[neighbs_two][x][i_two] = &xlink->head_two_;
				n_heads_ii_[neighbs_two][x]++;
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
					int x = xlink->x_dist_;
					int neighbs_one = 
						xlink->head_one_.site_->GetPRC1NeighborCount();
					int i_one = n_heads_ii_[neighbs_one][x];
					heads_ii_[neighbs_one][x][i_one] = &xlink->head_one_;
					n_heads_ii_[neighbs_one][x]++;
					int neighbs_two = 
						xlink->head_two_.site_->GetPRC1NeighborCount();
					int i_two = n_heads_ii_[neighbs_two][x];
					heads_ii_[neighbs_two][x][i_two] = &xlink->head_two_;
					n_heads_ii_[neighbs_two][x]++;
				}
				else{
					printf("wat in update_dub_unteth sites xlink\n");
				}
			}
		}
	}
}

void AssociatedProteinManagement::Update_Heads_II_Teth(){

	for(int n_neighbs(0); n_neighbs <= max_neighbs_ + 1; n_neighbs++){
		for(int x_dub(2*comp_cutoff_); x_dub <= 2*teth_cutoff_; x_dub++){
			for(int x(0); x <= dist_cutoff_; x++){
				n_heads_ii_teth_same_[n_neighbs][x_dub][x] = 0;
				n_heads_ii_teth_oppo_[n_neighbs][x_dub][x] = 0;
			}
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
					Tubulin *site_one = xlink->head_one_.site_;
					int n_one = site_one->GetPRC1NeighborCount();
					if(x == rest_dist_){
						int i_one = n_heads_ii_teth_same_[n_one][x_dub][x];
						heads_ii_teth_same_[n_one][x_dub][x][i_one] 
							= &xlink->head_one_;
						n_heads_ii_teth_same_[n_one][x_dub][x]++;
						int i_two = n_heads_ii_teth_oppo_[n_one][x_dub][x];
						heads_ii_teth_oppo_[n_one][x_dub][x][i_two] 
							= &xlink->head_one_;
						n_heads_ii_teth_oppo_[n_one][x_dub][x]++;
					}
					else if(site_one->EquilibriumInSameDirection()){
						int i = n_heads_ii_teth_same_[n_one][x_dub][x];
						heads_ii_teth_same_[n_one][x_dub][x][i]
							= &xlink->head_one_;
						n_heads_ii_teth_same_[n_one][x_dub][x]++;
					}
					else{
						int i = n_heads_ii_teth_oppo_[n_one][x_dub][x]; 
						heads_ii_teth_oppo_[n_one][x_dub][x][i]
							= &xlink->head_one_;
						n_heads_ii_teth_oppo_[n_one][x_dub][x]++;
					}
					// Site two
					Tubulin *site_two = xlink->head_two_.site_;
					int n_two = site_two->GetPRC1NeighborCount();
					if(x == rest_dist_){
						int i_one = n_heads_ii_teth_same_[n_two][x_dub][x];
						heads_ii_teth_same_[n_two][x_dub][x][i_one]
							= &xlink->head_two_;
						n_heads_ii_teth_same_[n_two][x_dub][x]++;
						int i_two = n_heads_ii_teth_oppo_[n_two][x_dub][x];
						heads_ii_teth_oppo_[n_two][x_dub][x][i_two]
							= &xlink->head_two_;
						n_heads_ii_teth_oppo_[n_two][x_dub][x]++;
					}
					else if(site_two->EquilibriumInSameDirection()){
						int i = n_heads_ii_teth_same_[n_two][x_dub][x];
						heads_ii_teth_same_[n_two][x_dub][x][i]
							= &xlink->head_two_;
						n_heads_ii_teth_same_[n_two][x_dub][x]++;
					}
					else{
						int i = n_heads_ii_teth_oppo_[n_two][x_dub][x];
						heads_ii_teth_oppo_[n_two][x_dub][x][i]
							= &xlink->head_two_;
						n_heads_ii_teth_oppo_[n_two][x_dub][x]++;
					}
				}
			}
		}
	}
}

void AssociatedProteinManagement::Update_Bound_Unteth(){
	
	n_bound_unteth_ = 0;
	for(int i_xlink = 0; i_xlink < n_active_; i_xlink++){
		AssociatedProtein *xlink = active_[i_xlink];
		if(xlink->heads_active_ > 0
		&& xlink->tethered_ == false){
			bound_unteth_[n_bound_unteth_] = xlink;
			n_bound_unteth_++;
		}
	}
}

void AssociatedProteinManagement::Run_KMC(){
	
	
	sys_time start1 = sys_clock::now();
	if(parameters_->xlinks.c_bulk > 0) Generate_KMC_List();
	else return;
	sys_time start2 = sys_clock::now();
//	printf("Start of xlink KMC cycle\n");
	for(int i_event = 0; i_event < IDs_to_exe_.size(); i_event++)
		events_[IDs_to_exe_[i_event]].Execute();
	sys_time finish = sys_clock::now();
	auto elapsed = std::chrono::duration_cast<t_unit>(finish - start1);
	properties_->wallace.t_xlinks_kmc_[0] += elapsed.count();
	elapsed = std::chrono::duration_cast<t_unit>(finish - start2);
	properties_->wallace.t_xlinks_kmc_[3] += elapsed.count();
}

void AssociatedProteinManagement::Generate_KMC_List(){

	// Scan through all events & calculate number expected to occur
	sys_time start = sys_clock::now();
	Update_All();
	int n_events_tot = 0;
	for(int i_event(0); i_event < events_.size(); i_event++){
		events_[i_event].SampleStatistics();
		n_events_tot += events_[i_event].n_expected_;
	}
	sys_time finish = sys_clock::now();
	auto elapsed = std::chrono::duration_cast<t_unit>(finish - start);
	properties_->wallace.t_xlinks_kmc_[1] += elapsed.count();
	// Scan through all target populations & ensure none become negative 
	start = sys_clock::now();
	for(int i_pop(0); i_pop < IDs_by_pop_.size(); i_pop++){
		int n_competitors = IDs_by_pop_[i_pop].size();
		if(n_competitors > 1){
			bool coupled(false);
			int n_events_loc = 0;
			double p_tot = 0;
			for(int i_entry(0); i_entry < n_competitors; i_entry++){
				int i_event = IDs_by_pop_[i_pop][i_entry];
				if(i_event < 0){
					coupled = true;
					i_event = abs(i_event);
				}
				n_events_loc += events_[i_event].n_expected_;
				p_tot += (events_[i_event].n_expected_ 
						* events_[i_event].p_occur_); 
			}
			int n_avail_loc 
				= *events_[abs(IDs_by_pop_[i_pop][0])].n_avail_;
//			if(coupled)
//			if(events_[abs(IDs_by_pop_[i_pop][0])].n_neighbs_ == 1) 
//				n_avail_loc /= 2; 
			while(n_events_loc > n_avail_loc){
				double p_cum = 0;
				double ran = properties_->gsl.GetRanProb();
				for(int i_entry(0); i_entry < n_competitors; i_entry++){
					int i_event = abs(IDs_by_pop_[i_pop][i_entry]);
					p_cum += (events_[i_event].n_expected_
							* events_[i_event].p_occur_) / p_tot;
					if(p_cum >= ran
					&& events_[i_event].n_expected_ > 0){
						events_[i_event].n_expected_--;
						n_events_tot--;
						n_events_loc--;
						p_tot -= events_[i_event].p_occur_;
						break;
					}
				}
			}
		}
	}
	// Construct KMC list if any events are predicted, otherwise clear it
	if(n_events_tot > 0){
		int pre_array[n_events_tot];
		int i_array(0);
		for(int i_event(0); i_event < events_.size(); i_event++){
			int n_expected = events_[i_event].n_expected_;
			for(int i_entry(0); i_entry < n_expected; i_entry++){
				if(i_array >= n_events_tot){
					printf("WOAH BUDDY. check XLINK genKMC\n");
					exit(1);
				}
				// Store & pass unique index of this event
				pre_array[i_array] = events_[i_event].ID_;
				i_array++;
			}
		}
		if(i_array != n_events_tot){
			printf("NOT SURE in xlink GEN KMC \n");
			exit(1);
		}
		if(n_events_tot > 1) gsl_ran_shuffle(properties_->gsl.rng_, 
				pre_array, n_events_tot, sizeof(int));
		IDs_to_exe_.resize(n_events_tot);
		for(int i_entry(0); i_entry < n_events_tot; i_entry++){
			IDs_to_exe_[i_entry] = pre_array[i_entry];
		}
	}
	else IDs_to_exe_.clear();
	finish = sys_clock::now();
	elapsed = std::chrono::duration_cast<t_unit>(finish - start);
	properties_->wallace.t_xlinks_kmc_[2] += elapsed.count();
}

void AssociatedProteinManagement::KMC_Relay(monomer target, int code){

	int x, x_dub; 
	bool diffusion_event(false);
	// Diffusion events are encoded as negative values
	if(kmc_code < 0){
		diffusion_event = true;
		kmc_code = abs(kmc_code);
	}
	if(kmc_code > 100){
		x = kmc_code % 10; 
		x_dub = (kmc_code % 1000 - x) / 10;
		kmc_code = (int)(kmc_code / 1000);
	}
	if(diffusion_event){
		switch(kmc_code){
			case 10:
//				printf("unteth_i step fwd\n");
				Diffuse_I_Fwd(n_neighbs);
				break;
			case 11:
//				printf("unteth_i step bck\n");
				Diffuse_I_Bck(n_neighbs);
				break;
			case 20:
//				printf("unteth_ii step to rest (%i)[%i avail]\n", 
//					x_dist, n_heads_ii_[x_dist]);
				Diffuse_II_To_Rest(n_neighbs, x);
				break;
			case 21:
//				printf("unteth_ii step from rest (%i)[%i avail]\n", 
//					x_dist, n_heads_ii_[x_dist]);
				Diffuse_II_Fr_Rest(n_neighbs, x);
				break;
			case 30:
//				printf("teth_i step to rest (%i)\n", x_dist_dub);
				Diffuse_I_To_Teth(n_neighbs, x_dub);
				break;
			case 31:
//				printf("teth_i step from rest (%i)\n", x_dist_dub);
				Diffuse_I_Fr_Teth(n_neighbs, x_dub);
				break;
			case 40:
//				printf("teth_ii step to both (2x: %i, x: %i)\n", x_dub, x);
				Diffuse_II_To_Both(n_neighbs, x_dub, x);
				break;
			case 41:
//				printf("teth_ii step fr both (2x: %i, x: %i)\n", x_dub, x);
				Diffuse_II_Fr_Both(n_neighbs, x_dub, x);
				break;
			case 50:
//				printf("teth_ii step to self from teth (2x: %i, x: %i)\n", 
//					x_dub, x);
				Diffuse_II_To_Self_Fr_Teth(n_neighbs, x_dub, x);
				break;
			case 51:
//				printf("teth_ii step from self to teth (2x: %i, x: %i)\n",
//					x_dub, x);
				Diffuse_II_Fr_Self_To_Teth(n_neighbs, x_dub, x);
				break;
		}
	}
	else{
		switch(kmc_code){
			case 10: 
//				printf("xlink bound stage 1\n");
				Bind_I(n_neighbs);
				break;
			case 11:
//				printf("xlink TETH bound stage 1\n");
				Bind_I_Teth();
				break;
			case 20: 
//				printf("xlink bound stage 2\n");
				Bind_II();
				break;
			case 21:
//				printf("bind ii teth xlink\n");
				Bind_II_Teth();
				break;
			case 30: 
//				printf("xlink fully unbound\n");
				Unbind_I(n_neighbs);
				break;
			case 31:
//				printf("unbind i teth xlink\n");
				Unbind_I_Teth(n_neighbs, x_dub);
				break;
			case 40: 
//				printf("xlink unbound 2nd head (ext was %i)\n", x);
				Unbind_II(n_neighbs, x);
				break;
			case 41:
//				printf("unbind ii to teth XLINK\n");
				Unbind_II_To_Teth(n_neighbs, x_dub, x);
				break;
			case 42:
//				printf("unbind ii from teth XLINK\n");
				Unbind_II_Fr_Teth(n_neighbs, x_dub, x);
				break;
			case 50:
//				printf("xlink tethered free\n");
				Tether_Free();
				break;
			case 60:
//				printf("xlink untethered free\n");
				Untether_Free();
				break;
		}
	}
}

void AssociatedProteinManagement::Diffuse_I_Fwd(int n_neighbs){

	Update_Heads_I();
	int n_bound = n_heads_i_[n_neighbs];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		AssociatedProtein::Monomer *head = heads_i_[n_neighbs][i_entry];
		Tubulin *site = head->site_;
		// cant step off them MTs
		if(site->index_ != site->mt_->n_sites_ - 1){
			Tubulin *old_site = site;
			Tubulin *new_site = &site->mt_->lattice_[site->index_ + 1];
			if(!new_site->occupied_){
				old_site->xlink_head_ = nullptr;
				old_site->occupied_ = false;
				head->site_ = new_site;
				new_site->xlink_head_ = head;
				new_site->occupied_ = true;
				scratch_[n_neighbs][n_scratch_[n_neighbs]] = head;
				n_scratch_[n_neighbs]++;
				if(n_neighbs == 1){
					scratch_[n_neighbs][n_scratch_[n_neighbs]] 
						= site->mt_->lattice_[site->index_-1].xlink_head_;
					n_scratch_[n_neighbs]++;
				}
			}
		}
	}
	else if(n_scratch_[n_neighbs] > 0){
		AssociatedProtein::Monomer *head(nullptr);
		for(int i_scratch(0); i_scratch<n_scratch_[n_neighbs]; i_scratch++){
			head = scratch_[n_neighbs][i_scratch];
			if(head->state_ == "bound_i") break;
			else head = nullptr;
		}
		if(head != nullptr){
			Tubulin *site = head->site_;
			// cant step off them MTs
			if(site->index_ != site->mt_->n_sites_ - 1){
				Tubulin *old_site = site;
				Tubulin *new_site = &site->mt_->lattice_[site->index_ + 1];
				if(!new_site->occupied_){
					old_site->xlink_head_ = nullptr;
					old_site->occupied_ = false;
					head->site_ = new_site;
					new_site->xlink_head_ = head;
					new_site->occupied_ = true;
				}
			}
		}
		else{
			printf("come on kid we cant step_i_fwd TWOO\n");
		}
	}
	else{
		printf("come on kid we cant step_i_fwd\n");
//		exit(1);
	}
}

void AssociatedProteinManagement::Diffuse_I_Bck(int n_neighbs){

	xlink_monomer head = GetActiveHead("bound_i", n_neighbs);
	Update_Heads_I();
	int n_bound = n_heads_i_[n_neighbs];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		AssociatedProtein::Monomer* head = heads_i_[n_neighbs][i_entry];
		Tubulin* site = head->site_;
		// can't step off them MTs
		if(site->index_ != 0){
			Tubulin *old_site = site;
			Tubulin *new_site = &site->mt_->lattice_[site->index_ - 1];
			if(new_site->occupied_ == false){
				old_site->xlink_head_ = nullptr;
				old_site->occupied_ = false;
				head->site_ = new_site;
				new_site->xlink_head_ = head;
				new_site->occupied_ = true;
				scratch_[n_neighbs][n_scratch_[n_neighbs]] = head;
				n_scratch_[n_neighbs]++;
				if(n_neighbs == 1){
					scratch_[n_neighbs][n_scratch_[n_neighbs]] 
						= site->mt_->lattice_[site->index_+1].xlink_head_;
					n_scratch_[n_neighbs]++;
				}
			}
		}
	}
	else if(n_scratch_[n_neighbs] > 0){
		AssociatedProtein::Monomer *head(nullptr);
		for(int i_scratch(0); i_scratch<n_scratch_[n_neighbs]; i_scratch++){
			head = scratch_[n_neighbs][i_scratch];
			if(head->state_ == "bound_i") break;
			else head = nullptr;
		}
		if(head != nullptr){
			Tubulin* site = head->site_;
			// can't step off them MTs
			if(site->index_ != 0){
				Tubulin *old_site = site;
				Tubulin *new_site = &site->mt_->lattice_[site->index_ - 1];
				if(new_site->occupied_ == false){
					old_site->xlink_head_ = nullptr;
					old_site->occupied_ = false;
					head->site_ = new_site;
					new_site->xlink_head_ = head;
					new_site->occupied_ = true;
				}
			}
		}
		else{
			printf("come on kid we cant step_i_bck TWOO\n");
		}
	}
	else{
		printf("come on kid we cant step_i_bck\n");
//		exit(1);
	}
}

void AssociatedProteinManagement::Diffuse_II_To_Rest(int n_neighbs, int x){

	/*
	Update_Bound_II_Sites();
	int n_bound = n_heads_ii_[n_neighbs][x];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site = heads_ii_[n_neighbs][x][i_entry];
		if(site == nullptr){
			printf("n: %i, x: %i, i: %i\n", n_neighbs, x, i_entry);
			printf("WHOO in DIFF II TO REST \n");
			exit(1);
		}
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
//		exit(1);
	}
	*/
}

void AssociatedProteinManagement::Diffuse_II_Fr_Rest(int n_neighbs, int x){

	/*
	Update_Bound_II_Sites();
	int n_bound = n_heads_ii_[n_neighbs][x];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site = heads_ii_[n_neighbs][x][i_entry];
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
//		exit(1);
	}
	*/
}

void AssociatedProteinManagement::Diffuse_I_To_Teth(
		int n_neighbs, int x_dub){

	/*
	Update_Bound_I_Teth();
	int n_bound = n_heads_i_teth_[n_neighbs][x_dub];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		AssociatedProtein *xlink = heads_i_teth_[n_neighbs][x_dub][i_entry];
		Tubulin *site = xlink->GetActiveHeadSite();
		int i_site = site->index_;
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
//		exit(1);
	}
	*/
}

void AssociatedProteinManagement::Diffuse_I_Fr_Teth(
		int n_neighbs, int x_dub){

	/*
	Update_Bound_I_Teth();
	int n_bound = n_heads_i_teth_[n_neighbs][x_dub];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		AssociatedProtein *xlink = heads_i_teth_[n_neighbs][x_dub][i_entry];
		Tubulin *site = xlink->GetActiveHeadSite();
		int i_site = site->index_;
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
//		exit(1);
	}
	*/
}

void AssociatedProteinManagement::Diffuse_II_To_Both(
		int n_neighbs, int x_dub, int x){

	/*
	Update_Bound_II_Sites_Teth();
	int n_bound = n_heads_ii_teth_same_[n_neighbs][x_dub][x];
	if(n_bound > 0){
		int i = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site = heads_ii_teth_same_[n_neighbs][x_dub][x][i];
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
//		exit(1);
	}
	*/
}

void AssociatedProteinManagement::Diffuse_II_Fr_Both(
		int n_neighbs, int x_dub, int x){

	/*
	Update_Bound_II_Sites_Teth();
	int n_bound = n_heads_ii_teth_same_[n_neighbs][x_dub][x];
	if(n_bound > 0){
		int i = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site = heads_ii_teth_same_[n_neighbs][x_dub][x][i];
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
//		exit(1);
	}
	*/
}

void AssociatedProteinManagement::Diffuse_II_To_Self_Fr_Teth(
		int n_neighbs, int x_dub, int x){

	/*
	Update_Bound_II_Sites_Teth();
	int	n_bound = n_heads_ii_teth_oppo_[n_neighbs][x_dub][x];
	if(n_bound > 0){
		int i = properties_->gsl.GetRanInt(n_bound);
	   	Tubulin *site = heads_ii_teth_oppo_[n_neighbs][x_dub][x][i];
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
//		exit(1);
	}
	*/
}

void AssociatedProteinManagement::Diffuse_II_Fr_Self_To_Teth(
		int n_neighbs, int x_dub, int x){

	/*
	Update_Bound_II_Sites_Teth();
	int n_bound = n_heads_ii_teth_oppo_[n_neighbs][x_dub][x];
	if(n_bound > 0){
		int i = properties_->gsl.GetRanInt(n_bound);
	   	Tubulin *site = heads_ii_teth_oppo_[n_neighbs][x_dub][x][i];
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
//		exit(1);
	}
	*/
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
		Tubulin *site
			= properties_->microtubules.GetUnoccupiedSite(n_neighbs);
		// Place xlink onto site
		site->xlink_head_ = &xlink->head_one_;
		site->occupied_ = true;
		// Update xlink details
		xlink->heads_active_++;
		xlink->head_one_.site_ = site; 
		xlink->head_one_.state_ = std::string("bound_i");
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

	 /*
	// Make sure stage 1 xlinks and unoccupied sites are available
	Update_Bound_I();
	properties_->microtubules.UpdateUnoccupied();
	int n_bound = n_heads_i_[max_neighbs_+1];
	if(n_bound > 0 && properties_->microtubules.n_unoccupied_> 0){
		// Randomly pick single-bound xlink
		int i_xlink = properties_->gsl.GetRanInt(n_bound);
		AssociatedProtein *xlink = heads_i_[max_neighbs_+1][i_xlink];
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
			if(attempts > 10*n_bound){
				break;
			}
			i_xlink = properties_->gsl.GetRanInt(n_bound);
			xlink = heads_i_[max_neighbs_+1][i_xlink];	
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
	*/
}

void AssociatedProteinManagement::Unbind_I(int n_neighbs){

	Update_Heads_I();
	int n_bound = n_heads_i_[n_neighbs];
	if(n_bound > 0){
		// Randomly pick a singly-bound xlink
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		AssociatedProtein::Monomer *head = heads_i_[n_neighbs][i_entry];
		Tubulin *site = head->site_;
		// Remove xlink from site
		site->xlink_head_ = nullptr;
		site->occupied_ = false;
		// Update xlink details
		head->xlink_->heads_active_--;
		head->state_ = std::string("unbound");
		head->site_ = nullptr;
		// If this xlink has a satellite motor, untether it
		if(head->xlink_->tethered_ == true){ 
			if(head->xlink_->motor_->heads_active_ == 0){
				head->xlink_->UntetherSatellite();
			}
			else{
				printf("error in xlink unbind_i\n");
				exit(1);
			}
		}
		// Remove this xlink from active_, replace with last entry
		AssociatedProtein *last_entry = active_[n_active_ - 1];
		int this_index = head->xlink_->active_index_; 
		if(this_index != n_active_ - 1){
			active_[this_index] = last_entry; 
			last_entry->active_index_ = this_index; 
		}
		n_active_--;
	}
	else if(n_scratch_[n_neighbs] > 0){
		AssociatedProtein::Monomer *head(nullptr);
		for(int i_scratch(0); i_scratch<n_scratch_[n_neighbs]; i_scratch++){
			head = scratch_[n_neighbs][i_scratch];
			if(head->state_ == "bound_i") break;
			else head = nullptr;
		}
		if(head != nullptr){
			Tubulin *site = head->site_;
			// Remove xlink from site
			site->xlink_head_ = nullptr;
			site->occupied_ = false;
			// Update xlink details
			head->xlink_->heads_active_--;
			head->state_ = std::string("unbound");
			head->site_ = nullptr;
			// If this xlink has a satellite motor, untether it
			if(head->xlink_->tethered_ == true){ 
				if(head->xlink_->motor_->heads_active_ == 0){
					head->xlink_->UntetherSatellite();
				}
				else{
					printf("error in xlink unbind_i\n");
					exit(1);
				}
			}
			// Remove this xlink from active_, replace with last entry
			AssociatedProtein *last_entry = active_[n_active_ - 1];
			int this_index = head->xlink_->active_index_; 
			if(this_index != n_active_ - 1){
				active_[this_index] = last_entry; 
				last_entry->active_index_ = this_index; 
			}
			n_active_--;
		}
		else{
			printf("Error in Unbind: no bound xlinks TWO\n");
			//		exit(1);
		}
	}
	else{
		printf("Error in Unbind: no bound xlinks\n");
//		exit(1);
	}
}

void AssociatedProteinManagement::Unbind_II(int n_neighbs, int x){

	/*
	Update_Bound_II_Sites();
	int n_bound = n_heads_ii_[n_neighbs][x];
	if(n_bound > 0){
		// Randomly pick a double-bound xlink
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site = heads_ii_[n_neighbs][x][i_entry];
		AssociatedProtein* xlink = site->xlink_;
		if(x != xlink->x_dist_){
			printf("error in xink unbind_ii \n");
			exit(1);
		}
		if(xlink->tethered_ == true){
			if(xlink->motor_->heads_active_ != 0){
				printf("error in xlink unbind_ii NON TETH\n");
				exit(1);
			}
		}
		// Remove xlink from site
		site->xlink_ = nullptr;
		site->occupied_ = false;
		// Update xlink details
		if(site == xlink->site_one_)
			xlink->site_one_ = nullptr;
		else
			xlink->site_two_ = nullptr;
		xlink->heads_active_--;
	}
	else{
		printf("Error in Unbind_II:no double bound xlinks\n");
//		exit(1);
	}
	*/
}

void AssociatedProteinManagement::Bind_I_Teth(){

	/*
	Update_Free_Teth();
	properties_->microtubules.UpdateUnoccupied();
	int n_bound = n_free_teth_;
	if(n_bound > 0 && properties_->microtubules.n_unoccupied_> 0){
		int i_xlink = properties_->gsl.GetRanInt(n_bound);
		AssociatedProtein* xlink = free_teth_[i_xlink];
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
	*/
}

void AssociatedProteinManagement::Bind_II_Teth(){

	/*
	Update_Bound_I_Teth();
	properties_->microtubules.UpdateUnoccupied();
	int n_bound_tot = n_heads_i_teth_tot_;
	if(n_bound_tot > 0 && properties_->microtubules.n_unoccupied_> 0){
		// Get a weighted teth extension
		int x_dist_dub = -1;
		// total weight summed over all extensions
		double weight_tot(0);
		// cumulative weight up to and including index extension
		double weight_cum[2*teth_cutoff_ + 1];
		// scan over all all extensions, add up weights appropriately
		for(int x_dub(2*comp_cutoff_); x_dub <= 2*teth_cutoff_; x_dub++){
			double wt = xlinks_[0].teth_binding_weight_lookup_[x_dub];
			weight_cum[x_dub] = wt * n_heads_i_teth_[max_neighbs_+1][x_dub];
			weight_tot += weight_cum[x_dub];
			if(weight_cum[x_dub] > 0){
				for(int X_DUB(2*comp_cutoff_); X_DUB < x_dub; X_DUB++){
					weight_cum[x_dub] += weight_cum[X_DUB];
				}
			}
		}
		// now that we have weight_tot, we can use weight_cum[i] to
		// map the available x_dub populations onto different ranges 
		// in [0, 1), which allows us to use GetRanProb() to select one
		double ran = properties_->gsl.GetRanProb();
		for(int x_dub(2*comp_cutoff_); x_dub <= 2*teth_cutoff_; x_dub++){
			if(ran <= weight_cum[x_dub] / weight_tot){
				x_dist_dub = x_dub;
				break;
			}
		}
		if(x_dist_dub == -1){
			printf("error in xlink bind_ii_teth ZERO\n");
			exit(1);
		}
		int n_bound = n_heads_i_teth_[max_neighbs_+1][x_dist_dub];
		if(n_bound > 0){
			// Randomly pick single-bound xlink
			int i_xl = properties_->gsl.GetRanInt(n_bound);
			AssociatedProtein *xlink 
				= heads_i_teth_[max_neighbs_+1][x_dist_dub][i_xl];
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
					if(switches > 50*n_bound_tot){
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
					n_bound = n_heads_i_teth_[max_neighbs_+1][x_dist_dub];
					attempts = 0;
				}
				i_xl = properties_->gsl.GetRanInt(n_bound);
				xlink = heads_i_teth_[max_neighbs_+1][x_dist_dub][i_xl];	
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
	*/
}

void AssociatedProteinManagement::Unbind_I_Teth(int n_neighbs, int x_dub){

	/*
	Update_Bound_I_Teth();
	int n_bound = n_heads_i_teth_[n_neighbs][x_dub];
	if(n_bound > 0){
		// Randomly pick a single-bound xlink
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		AssociatedProtein *xlink = heads_i_teth_[n_neighbs][x_dub][i_entry];
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
	*/
}

void AssociatedProteinManagement::Unbind_II_To_Teth(
		int n_neighbs, int x_dub, int x){

	/*
	Update_Bound_II_Sites_Teth();
	int n_bound = n_heads_ii_teth_same_[n_neighbs][x_dub][x];
	if(n_bound > 0){
		// Randomly pick a double-bound xlink
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Tubulin* site = heads_ii_teth_same_[n_neighbs][x_dub][x][i_entry];
		AssociatedProtein* xlink = site->xlink_;
		if(x != xlink->x_dist_){
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
	*/
}

void AssociatedProteinManagement::Unbind_II_Fr_Teth(
		int n_neighbs, int x_dub, int x){

	/*
	Update_Bound_II_Sites_Teth();
	int n_bound = n_heads_ii_teth_oppo_[n_neighbs][x_dub][x];
	if(n_bound > 0){
		// Randomly pick a double-bound xlink
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Tubulin* site = heads_ii_teth_oppo_[n_neighbs][x_dub][x][i_entry];
		AssociatedProtein* xlink = site->xlink_;
		if(x != xlink->x_dist_){
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
	*/
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

	Update_Free_Teth();
	if(n_free_teth_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_free_teth_);
		AssociatedProtein* xlink = free_teth_[i_entry];
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
