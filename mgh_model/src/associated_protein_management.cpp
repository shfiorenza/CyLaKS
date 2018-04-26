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
	InitiateLists();
}

void AssociatedProteinManagement::GenerateXLinks(){

	int n_mts = parameters_->n_microtubules;
	int n_sites = parameters_->length_of_microtubule;
	// Since only one head has to be bound, the sim will at most
	// as many xlinks as sites in the bulk (all single-bound)
	n_xlinks_ = n_mts*n_sites;
	xlink_list_.resize(n_xlinks_);
	for(int ID = 0; ID < n_xlinks_; ID++){
		xlink_list_[ID].Initialize(parameters_, properties_, ID);
	}
}

void AssociatedProteinManagement::SetParameters(){

	double delta_t = parameters_->delta_t;
	double site_size = parameters_->site_size;
	// DIFFUSION STATISTICS FOR SELF BELOW
	double D_const_i = parameters_->D_xlink_i;
	double D_const_ii = parameters_->D_xlink_ii;
	double x_squared = (site_size/1000)*(site_size/1000); // in um^2
	tau_i_ = x_squared / (2 * D_const_i);
	tau_ii_ = x_squared / (2 * D_const_ii);
	p_diffuse_i_fwd_ = delta_t / tau_i_;
	p_diffuse_i_bck_ = delta_t / tau_i_;
	// Generate different stepping rates based on changes in
	// potential energy (dU) associated with that step
	dist_cutoff_ = xlink_list_[0].dist_cutoff_;
	rest_dist_ = xlink_list_[0].rest_dist_;
	p_diffuse_ii_to_rest_.resize(dist_cutoff_ + 1);
	p_diffuse_ii_from_rest_.resize(dist_cutoff_ + 1);
	double r_0 = xlink_list_[0].r_0_;
	double kbT = xlink_list_[0].kbT_;
	double k_spring = xlink_list_[0].k_spring_;
	double r_y = 35;	// in nm
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
			double p_to = weight_to * delta_t / tau_ii_; //XXX
			double p_from = weight_from * delta_t / tau_ii_;
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
			printf("for x of %i (to): %g\n", x_dist, p_diffuse_ii_to_rest_[x_dist]);
			printf("for x of %i (from): %g\n", x_dist, p_diffuse_ii_from_rest_[x_dist]);
		}
		else{
			printf("woah mayne. xlink set parameters \n");
			exit(1);
		}
	}
	// DIFFUSION STATISTICS INVOLVING TETHER BELOW
	int teth_dist_cutoff = properties_->kinesin4.motor_list_[0].dist_cutoff_;
	int teth_comp_cutoff = properties_->kinesin4.motor_list_[0].comp_cutoff_;
	p_diffuse_i_to_teth_rest_.resize(2*teth_dist_cutoff + 1);
	p_diffuse_i_from_teth_rest_.resize(2*teth_dist_cutoff + 1);
	p_diffuse_ii_to_both_rest_.resize(2*teth_dist_cutoff + 1);
	p_diffuse_ii_to_self_from_teth_.resize(2*teth_dist_cutoff + 1);
	p_diffuse_ii_from_self_to_teth_.resize(2*teth_dist_cutoff + 1);
	p_diffuse_ii_from_both_rest_.resize(2*teth_dist_cutoff + 1);
	double k_teth_spring = properties_->kinesin4.motor_list_[0].k_spring_;
	double k_teth_slack = properties_->kinesin4.motor_list_[0].k_eff_slack_;
	double r_0_teth = properties_->kinesin4.motor_list_[0].r_0_;
	double r_y_teth = 17.5;
	double rest_dist_teth = properties_->kinesin4.motor_list_[0].rest_dist_;
	for(int x_dist_dub = 0; x_dist_dub <= 2*teth_dist_cutoff; x_dist_dub++){
		p_diffuse_ii_to_both_rest_[x_dist_dub].resize(dist_cutoff_ + 1);
		p_diffuse_ii_to_self_from_teth_[x_dist_dub].resize(dist_cutoff_ + 1);
		p_diffuse_ii_from_self_to_teth_[x_dist_dub].resize(dist_cutoff_ + 1);
		p_diffuse_ii_from_both_rest_[x_dist_dub].resize(dist_cutoff_ + 1);
		// Calc x-distances (in nm) for tether
		double r_x_teth = x_dist_dub * site_size / 2;
		double r_x_teth_bck = (x_dist_dub - 1) * site_size / 2;
		double r_x_teth_fwd = (x_dist_dub + 1) * site_size / 2;
		// Calc total r values 
		double r_teth = sqrt(r_x_teth*r_x_teth + r_y_teth*r_y_teth);
		double r_teth_bck = sqrt(r_x_teth_bck*r_x_teth_bck + r_y_teth*r_y_teth);
		double r_teth_fwd = sqrt(r_x_teth_fwd*r_x_teth_fwd+ r_y_teth*r_y_teth);
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
		if(x_dist_dub == 2*rest_dist_teth){
			dU_from_teth = (k_teth_slack/2)
						 * (dr_teth_from*dr_teth_from - dr_teth*dr_teth);
			dU_to_teth = (1/2)*(k_teth_spring*dr_teth_to*dr_teth_to 
							 - k_teth_slack*dr_teth*dr_teth);
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
		if(x_dist_dub < 2*teth_comp_cutoff)
			weight_to_teth = 0;
		else
			weight_to_teth = exp(-dU_to_teth/(2*kbT));	
		double weight_from_teth;
		if(x_dist_dub < 2*teth_dist_cutoff
		&& x_dist_dub > 2*teth_comp_cutoff)
			weight_from_teth = exp(-dU_from_teth/(2*kbT));
		else
			weight_from_teth = 0; 
		double p_to_teth_i = weight_to_teth * delta_t / tau_i_;
		double p_from_teth_i = weight_from_teth * delta_t / tau_i_;
		// Input probabilities for stage_i / tethered xlinks
		p_diffuse_i_to_teth_rest_[x_dist_dub] = p_to_teth_i;
		p_diffuse_i_from_teth_rest_[x_dist_dub] = p_from_teth_i;
		// Run through x_dists to get probs for stage_ii / tethered xlinks
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
				double weight_to;
				double weight_from;
			   	if(x_dist == rest_dist_){
					weight_to = 0;
					weight_from = 2*exp(-dU_from_rest/(2*kbT));
				}
				else if(x_dist == dist_cutoff_){
					weight_to = exp(-dU_to_rest/(2*kbT));
					weight_from = 0;
				}
				else{
					weight_to = exp(-dU_to_rest/(2*kbT));
					weight_from = exp(-dU_from_rest/(2*kbT));
				}
				// Convolve these bitches
				p_diffuse_ii_to_both_rest_[x_dist_dub][x_dist]
					= weight_to * weight_to_teth * delta_t / tau_ii_;
				p_diffuse_ii_to_self_from_teth_[x_dist_dub][x_dist]
					= weight_to * weight_from_teth * delta_t / tau_ii_;
				p_diffuse_ii_from_self_to_teth_[x_dist_dub][x_dist]
					= weight_from * weight_to_teth * delta_t / tau_ii_;
				p_diffuse_ii_from_both_rest_[x_dist_dub][x_dist]
					= weight_from * weight_from_teth * delta_t / tau_ii_;
			}
			else{
				printf("woah mayne. xlink set parameters TWOO \n");
				exit(1);
			}
		}
	}
	// KMC STATISTICS BELOW
	c_eff_ = parameters_->c_eff_xlink;
	double k_on = parameters_->k_on_xlink; 
	double c_xlink = parameters_->c_xlink;
	p_bind_i_ = k_on * c_xlink * delta_t;
	double k_off_i = parameters_->k_off_xlink_i;
	p_unbind_i_ = k_off_i * delta_t;
	// Generate unbinding rates based on discretized spring extension
	double k_off_ii = parameters_->k_off_xlink_ii;
	p_unbind_ii_.resize(dist_cutoff_ + 1);
	for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
		double r_x = x_dist*site_size;
		double r = sqrt(r_y*r_y + r_x*r_x);
		double dr = r - r_0;
		double unbind_weight = exp(dr*dr*k_spring/(2*kbT));
//		printf("%i:  %g\n", distance, k_off_ii);
		p_unbind_ii_[x_dist] = k_off_ii * unbind_weight * delta_t;
	}
}

void AssociatedProteinManagement::InitiateLists(){

	// Stats (not a list ok bite me)
	n_double_bound_.resize(dist_cutoff_ + 1);
	n_sites_ii_untethered_.resize(dist_cutoff_ + 1);
	for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
		n_double_bound_[x_dist] = 0;
		n_sites_ii_untethered_[x_dist] = 0;
	}
	int teth_cutoff = properties_->kinesin4.motor_list_[0].dist_cutoff_;
	n_sites_i_tethered_.resize(2*teth_cutoff + 1);
	n_sites_ii_tethered_.resize(2*teth_cutoff + 1);
	n_sites_ii_tethered_same_.resize(2*teth_cutoff + 1);
	n_sites_ii_tethered_oppo_.resize(2*teth_cutoff + 1);
	for(int x_dist_dub = 0; x_dist_dub <= 2*teth_cutoff; x_dist_dub++){
		n_sites_i_tethered_[x_dist_dub] = 0;
		n_sites_ii_tethered_[x_dist_dub].resize(dist_cutoff_ + 1);
		n_sites_ii_tethered_same_[x_dist_dub].resize(dist_cutoff_ + 1);
		n_sites_ii_tethered_oppo_[x_dist_dub].resize(dist_cutoff_ + 1);
		for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
			n_sites_ii_tethered_[x_dist_dub][x_dist] = 0;
			n_sites_ii_tethered_same_[x_dist_dub][x_dist] = 0;
			n_sites_ii_tethered_oppo_[x_dist_dub][x_dist] = 0; 
		}
	}
	// Lists
	single_bound_list_.resize(n_xlinks_);
	untethered_list_.resize(n_xlinks_);
	single_untethered_sites_.resize(n_xlinks_);
	double_bound_list_.resize(dist_cutoff_ + 1);
	double_untethered_sites_.resize(dist_cutoff_ + 1);
	for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
		double_bound_list_[x_dist].resize(n_xlinks_);
		double_untethered_sites_[x_dist].resize(n_xlinks_);
	}
	single_tethered_sites_.resize(2*teth_cutoff + 1);
	double_tethered_sites_oppo_.resize(2*teth_cutoff + 1);
	double_tethered_sites_same_.resize(2*teth_cutoff + 1);
	for(int x_dist_dub = 0; x_dist_dub <= 2*teth_cutoff; x_dist_dub++){
		single_tethered_sites_[x_dist_dub].resize(n_xlinks_);
		double_tethered_sites_oppo_[x_dist_dub].resize(dist_cutoff_ + 1);
		double_tethered_sites_same_[x_dist_dub].resize(dist_cutoff_ + 1);
		for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
			double_tethered_sites_oppo_[x_dist_dub][x_dist]
				.resize(n_xlinks_);
			double_tethered_sites_same_[x_dist_dub][x_dist]
				.resize(n_xlinks_);
		}
	}
}

void AssociatedProteinManagement::BoundCheck(AssociatedProtein *xlink){

	// Check xlink
	if(xlink->site_one_ == nullptr
	&& xlink->site_two_ == nullptr){
		printf("ERROR WITH XLINK BINDING\n");
		exit(1);
	}
	// Check xlink_list
	int ID = xlink->ID_;
	if(xlink_list_[ID].site_one_ == nullptr
	&& xlink_list_[ID].site_two_ == nullptr){
		printf("ERROR WITH XLINK BINDING\n");
		exit(1);
	}
}
void AssociatedProteinManagement::UntetheredCheck(AssociatedProtein *xlink){

	// Check xlink
	if(xlink->motor_ != nullptr){
		printf("ERROR WITH XLINK TETHERING\n");
		exit(1);
	}
	// Check xlink_list 
	int ID = xlink->ID_;
	if(xlink_list_[ID].motor_ != nullptr){
		printf("ERROR WITH XLINK TETHERING 2\n");
		exit(1);
	}
}

void AssociatedProteinManagement::UpdateSingleBoundList(){
	
	int i_entry = 0;
	int n_entries = 0;
	for(int i_xlink = 0; i_xlink < n_xlinks_; i_xlink++){
		AssociatedProtein *xlink = &xlink_list_[i_xlink];
		if(xlink->heads_active_ == 1){
			single_bound_list_[i_entry] = xlink;
			i_entry++;
			n_entries++;
		}
	}
	if(n_entries != n_single_bound_){
		printf("something bad in update_single_bound_list (xlink)\n");
		printf("%i counted but %i in stats\n", n_entries, n_single_bound_);
		exit(1);
	}
}

void AssociatedProteinManagement::UpdateDoubleBoundList(){

	int i_entry[dist_cutoff_ + 1]; 
	int n_entries[dist_cutoff_ + 1];
	for(int dist = 0; dist <= dist_cutoff_; dist++){
		i_entry[dist] = 0;
		n_entries[dist] = 0;
	}
	for(int i_xlink = 0; i_xlink < n_xlinks_; i_xlink++){
		AssociatedProtein *xlink = &xlink_list_[i_xlink]; 
		if(xlink->heads_active_ == 2){
			int x_dist = xlink->x_dist_;
			int index = i_entry[x_dist];
			double_bound_list_[x_dist][index] = xlink;
			i_entry[x_dist]++;
			n_entries[x_dist]++;	
		}
	}
	for(int dist = 0; dist <= dist_cutoff_; dist++){
		if(n_entries[dist] != n_double_bound_[dist]){
			printf("something wrong in update double bound (xlink)\n");
			exit(1);
		}
	}
}

void AssociatedProteinManagement::UpdateUntetheredList(){
	
	int i_entry = 0;
	int n_entries = 0;
	for(int i_xlink = 0; i_xlink < n_xlinks_; i_xlink++){
		AssociatedProtein *xlink = &xlink_list_[i_xlink];
		if(xlink->heads_active_ > 0){
			if(xlink->tethered_ == false){
				untethered_list_[i_entry] = xlink;
				i_entry++;
				n_entries++;
			}
		}
	}
	if(n_entries != n_untethered_){
			printf("ugh, updateuntetheredlist:");
			printf(" %i in stats, %i in list\n", n_untethered_, n_entries);
	}
}

void AssociatedProteinManagement::UpdateSingleUntetheredSites(){

	int i_site = 0;
	int n_sites = 0;
	for(int i_xlink = 0; i_xlink < n_xlinks_; i_xlink++){
		AssociatedProtein *xlink = &xlink_list_[i_xlink];
		if(xlink->heads_active_ == 1
		&& xlink->tethered_ == false){
			Tubulin *site = xlink->GetActiveHeadSite();
			single_untethered_sites_[i_site] = site;
			i_site++;
			n_sites++;
		}
		// Treat xlinks tethered to free motors as untethered xlinks
		// (for diffusion at least)
		if(xlink->heads_active_ == 1
		&& xlink->tethered_ == true){
			if(xlink->motor_->heads_active_ == 0){
				Tubulin *site = xlink->GetActiveHeadSite();
				single_untethered_sites_[i_site] = site;
				i_site++;
				n_sites++;
			}
		}
	}
	if(n_sites != n_sites_i_untethered_){
		printf("bad mojo single bound untethered sites\n");
		printf("%i in stats but %i counted\n", n_sites_i_untethered_, 
				n_sites);
		exit(1);
	}
}

void AssociatedProteinManagement::UpdateDoubleUntetheredSites(){

	int i_site[dist_cutoff_ + 1];
	int n_sites[dist_cutoff_ + 1];
	for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
		i_site[x_dist] = 0;
		n_sites[x_dist] = 0;
	}
	for(int i_xlink = 0; i_xlink < n_xlinks_; i_xlink++){
		AssociatedProtein *xlink = &xlink_list_[i_xlink];
		if(xlink->heads_active_ == 2
		&& xlink->tethered_ == false){
			xlink->UpdateExtension();
			// Make sure we didn't force an unbind event
			if(xlink->heads_active_ == 2){
				int x_dist = xlink->x_dist_;
				Tubulin *site_one = xlink->site_one_;
				int index_one = i_site[x_dist];
				double_untethered_sites_[x_dist][index_one] = site_one;
			   	i_site[x_dist]++;	
				n_sites[x_dist]++;
				Tubulin *site_two = xlink->site_two_;
				int index_two = i_site[x_dist];
				double_untethered_sites_[x_dist][index_two] = site_two;
				i_site[x_dist]++;
				n_sites[x_dist]++;
			}
			else{
				printf("wat in update_dub_unteth sites xlink\n");
			}
		}
		if(xlink->heads_active_ == 2
		&& xlink->tethered_ == true){
			// Treat xlinks tethered to free motors as untethered xlinks
			// (for diffusion at least)
			if(xlink->motor_->heads_active_ == 0){
				xlink->UpdateExtension();
				// Make sure we didn't force an unbind event
				if(xlink->heads_active_ == 2){
					int x_dist = xlink->x_dist_;
					Tubulin *site_one = xlink->site_one_;
					int index_one = i_site[x_dist];
					double_untethered_sites_[x_dist][index_one] = site_one;
					i_site[x_dist]++;	
					n_sites[x_dist]++;
					Tubulin *site_two = xlink->site_two_;
					int index_two = i_site[x_dist];
					double_untethered_sites_[x_dist][index_two] = site_two;
					i_site[x_dist]++;
					n_sites[x_dist]++;
				}
				else{
					printf("wat in update_dub_unteth sites xlink\n");
				}
			}
		}
	}
	for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
		if(n_sites[x_dist] != n_sites_ii_untethered_[x_dist]){
			printf("bad mojo double bound untethered sites\n");
			printf("%i in stats but %i counted (x: %i)\n", 
					n_sites_ii_untethered_[x_dist], n_sites[x_dist], 
					x_dist);
			exit(1);
		}
	}
}

void AssociatedProteinManagement::UpdateSingleTetheredSites(){

	int teth_cutoff = properties_->kinesin4.motor_list_[0].dist_cutoff_;
    int i_site[2*teth_cutoff + 1];
	int n_sites[2*teth_cutoff + 1];
	for(int x_dist_dub = 0; x_dist_dub <= 2*teth_cutoff; x_dist_dub++){
		i_site[x_dist_dub] = 0;
		n_sites[x_dist_dub] = 0;
	}
    for(int i_xlink = 0; i_xlink < n_xlinks_; i_xlink++){
        AssociatedProtein *xlink = &xlink_list_[i_xlink];
        if(xlink->heads_active_ == 1
		&& xlink->tethered_ == true){
			// Dont count xlinks tethered to free motors
			if(xlink->motor_->heads_active_ > 0){
				xlink->motor_->UpdateExtension();
				// Make sure we didn't force an untether event
				if(xlink->tethered_ == true){
					int x_dist_dub = xlink->motor_->x_dist_doubled_;
					Tubulin *site = xlink->GetActiveHeadSite();
					int index = i_site[x_dist_dub];
					single_tethered_sites_[x_dist_dub][index] = site;
					i_site[x_dist_dub]++;
					n_sites[x_dist_dub]++;
				}
			}
		}
	}
	for(int x_dist_dub = 0; x_dist_dub <= 2*teth_cutoff; x_dist_dub++){
		if(n_sites[x_dist_dub] != n_sites_i_tethered_[x_dist_dub]){
			printf("bad mojo single bound tethered sites\n");
			printf("%i in stats but %i counted (for 2x: %i)\n", 
					n_sites_i_tethered_[x_dist_dub], n_sites[x_dist_dub], 
					x_dist_dub);
			exit(1);
		}
	}
}

void AssociatedProteinManagement::UpdateDoubleTetheredSites(){

	int teth_cutoff = properties_->kinesin4.motor_list_[0].dist_cutoff_;
	int i_site_oppo[2*teth_cutoff + 1][dist_cutoff_ + 1];
	int i_site_same[2*teth_cutoff + 1][dist_cutoff_ + 1];
	int n_sites[2*teth_cutoff + 1][dist_cutoff_ + 1];
	int n_sites_ii_unstretched[2*teth_cutoff + 1];
	for(int x_dist_dub = 0; x_dist_dub <= 2*teth_cutoff; x_dist_dub++){
		n_sites_ii_unstretched[x_dist_dub] = 0;
		for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
			i_site_oppo[x_dist_dub][x_dist] = 0;
			i_site_same[x_dist_dub][x_dist] = 0;
			n_sites[x_dist_dub][x_dist] = 0;
			// Reset internal statistics
			n_sites_ii_tethered_oppo_[x_dist_dub][x_dist] = 0;
			n_sites_ii_tethered_same_[x_dist_dub][x_dist] = 0;
		}
	}
	for(int i_xlink = 0; i_xlink < n_xlinks_; i_xlink++){
		AssociatedProtein *xlink = &xlink_list_[i_xlink];
		if(xlink->heads_active_ == 2
		&& xlink->tethered_ == true){
			// Dont count xlinks tethered to free motors
			if(xlink->motor_->heads_active_ > 0){
				xlink->motor_->UpdateExtension();
				xlink->UpdateExtension();
				// Make sure an untether or unbind event wasn't forced
				if(xlink->heads_active_ == 2
				&& xlink->tethered_ == true){
					int x_dist_dub = xlink->motor_->x_dist_doubled_;
					int x_dist = xlink->x_dist_;
					// Site one
					Tubulin *site_one = xlink->site_one_;
					n_sites[x_dist_dub][x_dist]++;
					if(x_dist == 0){
						int i = i_site_same[x_dist_dub][x_dist];
						int j = i_site_oppo[x_dist_dub][x_dist];
						double_tethered_sites_same_[x_dist_dub][x_dist][i]
							= site_one;
						double_tethered_sites_oppo_[x_dist_dub][x_dist][j]
							= site_one;
						i_site_same[x_dist_dub][x_dist]++;
						i_site_oppo[x_dist_dub][x_dist]++;
						n_sites_ii_tethered_same_[x_dist_dub][x_dist]++;
						n_sites_ii_tethered_oppo_[x_dist_dub][x_dist]++;
						n_sites_ii_unstretched[x_dist_dub]++;
					}
					else if(site_one->SpringEquilOnSameSide() == true){
						int i = i_site_same[x_dist_dub][x_dist];
						double_tethered_sites_same_[x_dist_dub][x_dist][i]
							= site_one;
						i_site_same[x_dist_dub][x_dist]++;
						n_sites_ii_tethered_same_[x_dist_dub][x_dist]++;
					}
					else{
						int i = i_site_oppo[x_dist_dub][x_dist];
						double_tethered_sites_oppo_[x_dist_dub][x_dist][i]
							= site_one;
						i_site_oppo[x_dist_dub][x_dist]++;
						n_sites_ii_tethered_oppo_[x_dist_dub][x_dist]++;
					}
					// Site two
					Tubulin *site_two = xlink->site_two_;
					n_sites[x_dist_dub][x_dist]++;
					if(x_dist == 0){
						int i = i_site_same[x_dist_dub][x_dist];
						int j = i_site_oppo[x_dist_dub][x_dist];
						double_tethered_sites_same_[x_dist_dub][x_dist][i]
							= site_two;
						double_tethered_sites_oppo_[x_dist_dub][x_dist][j]
							= site_two;
						i_site_same[x_dist_dub][x_dist]++;
						i_site_oppo[x_dist_dub][x_dist]++;
						n_sites_ii_tethered_same_[x_dist_dub][x_dist]++;
						n_sites_ii_tethered_oppo_[x_dist_dub][x_dist]++;
						n_sites_ii_unstretched[x_dist_dub]++;
					}
					else if(site_two->SpringEquilOnSameSide() == true){
						int i = i_site_same[x_dist_dub][x_dist];
						double_tethered_sites_same_[x_dist_dub][x_dist][i]
							= site_two;
						i_site_same[x_dist_dub][x_dist]++;
						n_sites_ii_tethered_same_[x_dist_dub][x_dist]++;
					}
					else{
						int i = i_site_oppo[x_dist_dub][x_dist];
						double_tethered_sites_oppo_[x_dist_dub][x_dist][i]
							= site_two;
						i_site_oppo[x_dist_dub][x_dist]++;
						n_sites_ii_tethered_oppo_[x_dist_dub][x_dist]++;
					}
				}
			}
		}
	}
	for(int x_dist_dub = 0; x_dist_dub <= 2*teth_cutoff; x_dist_dub++){
		for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
			if(n_sites[x_dist_dub][x_dist] 
					!= n_sites_ii_tethered_[x_dist_dub][x_dist]){
				printf("bad mojo double bound tethered sites\n");
				printf("%i in stats but %i counted (2x: %i, x: %i)\n", 
						n_sites_ii_tethered_[x_dist_dub][x_dist], 
						n_sites[x_dist_dub][x_dist], x_dist_dub, x_dist);
				exit(1);
			}
			if(x_dist == 0){
				if(n_sites[x_dist_dub][x_dist]
					!= n_sites_ii_tethered_same_[x_dist_dub][x_dist]
					 + n_sites_ii_tethered_oppo_[x_dist_dub][x_dist]
					 - n_sites_ii_unstretched[x_dist_dub]){
				printf("literally the worst error in ii_teth_sites\n");
				exit(1);
				}
			}
			else if(n_sites[x_dist_dub][x_dist]
					!= n_sites_ii_tethered_same_[x_dist_dub][x_dist]
					 + n_sites_ii_tethered_oppo_[x_dist_dub][x_dist]){
				printf("dont u like, believe in karma? ii_teth sites\n");
				printf("ticks: %i, same: %i, oppo: %i\n", 
						n_sites[x_dist_dub][x_dist], 
						n_sites_ii_tethered_same_[x_dist_dub][x_dist],
						n_sites_ii_tethered_oppo_[x_dist_dub][x_dist]);
				printf("for 2x: %i and x: %i\n", x_dist_dub, x_dist);
				exit(1);
			}
		}
	}
}

AssociatedProtein* AssociatedProteinManagement::GetUntetheredXlink(){

	UpdateUntetheredList();
	if(n_untethered_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_untethered_);
		AssociatedProtein* xlink = untethered_list_[i_entry];
		BoundCheck(xlink);
		UntetheredCheck(xlink);
		return xlink;
	}
	else{
		printf("Error in GetUnTetheredXlink: no untethered xlinks!\n");
		exit(1);
	}
}

void AssociatedProteinManagement::GenerateDiffusionList(){

	int n_events = 0;
	// Untethered statistics
	UpdateSingleUntetheredSites();
	int n_i_fwd = GetNumToStepI_Forward();
	int n_i_bck = GetNumToStepI_Backward();
	while(n_i_fwd + n_i_bck > n_sites_i_untethered_){
//		printf("yah for n_i_unteth\n");
		double ran = properties_->gsl.GetRanProb();
		if(ran < 0.5)
			n_i_fwd--;
		else
			n_i_bck--;
	}
	n_events += n_i_fwd;
	n_events += n_i_bck;
	UpdateDoubleUntetheredSites();
	int n_ii_to[dist_cutoff_ + 1];
	int n_ii_from[dist_cutoff_ + 1];
	// Handle statistics for each xlink extension separately
	for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
		n_ii_to[x_dist] = 0; //GetNumToStepII_ToRest(x_dist);
		n_ii_from[x_dist] = 0; //GetNumToStepII_FromRest(x_dist);
		int n_to = n_ii_to[x_dist];
		int n_from = n_ii_from[x_dist];
		int n_avail = n_sites_ii_untethered_[x_dist];
//		printf("%i out of %i\n", x_dist, dist_cutoff_);
		while(2*(n_to + n_from) > n_avail && n_avail != 0){
//		printf("yah for n_ii_unteth (x: %i)\n", x_dist);
//		printf("n_to: %i n_from: %i n_avail: %i\n", n_to, n_from, n_avail);
			double ran = properties_->gsl.GetRanProb();
			double p_to = p_diffuse_ii_to_rest_[x_dist];
			double p_from = p_diffuse_ii_from_rest_[x_dist];
			double p_tot = p_to + p_from;
			if(ran < p_to/p_tot
			&& n_to > 0){
				n_ii_to[x_dist]--;
				n_to = n_ii_to[x_dist];
			}
			else if(n_from > 0){
				n_ii_from[x_dist]--;
				n_from = n_ii_from[x_dist];
			}
			else if(n_to > 0){
				n_ii_to[x_dist]--;
				n_to = n_ii_to[x_dist];
			}
		}
		n_events += n_ii_to[x_dist];
		n_events += n_ii_from[x_dist];
	}
	// Handle statistics for each tether extension separately
	int teth_cutoff = properties_->kinesin4.motor_list_[0].dist_cutoff_;
	UpdateSingleTetheredSites();
	UpdateDoubleTetheredSites();
	int n_i_to_teth[2*teth_cutoff + 1];
	int n_i_from_teth[2*teth_cutoff + 1];
	int n_ii_to_both[2*teth_cutoff + 1][dist_cutoff_ + 1];
	int n_ii_to_self_from_teth[2*teth_cutoff + 1][dist_cutoff_ + 1];
	int n_ii_from_self_to_teth[2*teth_cutoff + 1][dist_cutoff_ + 1];
	int n_ii_from_both[2*teth_cutoff + 1][dist_cutoff_ + 1];
	int n_to, n_from, n_avail,
		n_to_one, n_from_one, 
		n_to_two, n_from_two,
		n_avail_one, n_avail_two;
	double p_to, p_from, p_tot,
		   p_to_one, p_from_one, 
		   p_to_two, p_from_two; 
	for(int x_dist_dub = 0; x_dist_dub <= 2*teth_cutoff; x_dist_dub++){
		n_i_to_teth[x_dist_dub] = GetNumToStepI_ToTethRest(x_dist_dub);
		n_i_from_teth[x_dist_dub] = GetNumToStepI_FromTethRest(x_dist_dub);
		// Make sure we don't get more events than sites that exist
		n_to = n_i_to_teth[x_dist_dub];
		n_from = n_i_from_teth[x_dist_dub];
		n_avail = n_sites_i_tethered_[x_dist_dub];
		while(n_to + n_from > n_avail){
			double ran = properties_->gsl.GetRanProb();
//			printf("yah for n_i_teth (2x: %i)\n", x_dist_dub);
//			printf("n_to: %i n_from: %i n_avail: %i\n", 
//					n_to, n_from, n_avail);
			p_to = p_diffuse_i_to_teth_rest_[x_dist_dub];
			p_from = p_diffuse_i_from_teth_rest_[x_dist_dub];
			p_tot = p_to + p_from;
			if(ran < p_to/p_tot
			&& n_to > 0){
				n_i_to_teth[x_dist_dub]--;
				n_to = n_i_to_teth[x_dist_dub];
			}
			else if(n_from > 0){
				n_i_from_teth[x_dist_dub]--;
				n_from = n_i_from_teth[x_dist_dub];
			}
			else if(n_to > 0){
				n_i_to_teth[x_dist_dub]--;
				n_to = n_i_to_teth[x_dist_dub];
			}
		}
		n_events += n_i_to_teth[x_dist_dub];
		n_events += n_i_from_teth[x_dist_dub];
		KinesinManagement *kinesin4 = &properties_->kinesin4;
		int rest_dist_dub = 2*kinesin4->motor_list_[0].rest_dist_;
		for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
			// First do stats for sites with xlink/teth rest on SAME side
			n_ii_to_both[x_dist_dub][x_dist] 
				= GetNumToStepII_ToBothRest(x_dist_dub, x_dist);
			n_ii_from_both[x_dist_dub][x_dist]
				= GetNumToStepII_FromBothRest(x_dist_dub, x_dist);
			// Next do stats for sites with xlink/teth on OPPOSITE sides
			n_ii_to_self_from_teth[x_dist_dub][x_dist]
				= GetNumToStepII_ToSelf_FromTeth(x_dist_dub, x_dist);
			n_ii_from_self_to_teth[x_dist_dub][x_dist]
				= GetNumToStepII_FromSelf_ToTeth(x_dist_dub, x_dist);
			// Make sure we don't get more events than available sites
			if(x_dist == 0){
				n_to = n_ii_from_self_to_teth[x_dist_dub][x_dist];
				n_from = n_ii_from_both[x_dist_dub][x_dist];
				n_avail = n_sites_ii_tethered_oppo_[x_dist_dub][x_dist];
				while(2*(n_to + n_from) > n_avail){
//					printf("yah for ii_teth_serr (2x: %i, x: %i)\n", 
//							x_dist_dub, x_dist);
//					printf("n_to: %i n_from: %i n_avail: %i\n", 
//							n_to, n_from, n_avail);
					double ran = properties_->gsl.GetRanProb();
					p_to = p_diffuse_ii_from_self_to_teth_
									[x_dist_dub][x_dist];
					p_from = p_diffuse_ii_from_both_rest_
									[x_dist_dub][x_dist];
					p_tot = p_to + p_from;
					if(ran < p_to/p_tot
							&& n_to > 0){
						n_ii_from_self_to_teth[x_dist_dub][x_dist]--;
						n_to = n_ii_from_self_to_teth[x_dist_dub][x_dist];
					}
					else if(n_from > 0){
						n_ii_from_both[x_dist_dub][x_dist]--;
						n_from = n_ii_from_both[x_dist_dub][x_dist];
					}
					else if (n_to > 0){
						n_ii_from_self_to_teth[x_dist_dub][x_dist]--;
						n_to = n_ii_from_self_to_teth[x_dist_dub][x_dist];
					}
				}
			}
			// For compression, sides are switched
			else if(x_dist_dub <= rest_dist_dub){
				n_to_one = n_ii_to_both[x_dist_dub][x_dist];
				n_from_one = n_ii_from_both[x_dist_dub][x_dist];
				n_avail_one = n_sites_ii_tethered_oppo_[x_dist_dub][x_dist];
				// Make sure there aren't more events than sites w/ oppo equil
				while((n_to_one + n_from_one) > n_avail_one){
//					printf("yah for to/from both_ext (2x: %i, x: %i)\n", 
//							x_dist_dub, x_dist);
//					printf("n_to: %i n_from: %i n_avail: %i\n", 
//							n_to_one, n_from_one, n_avail_one);
					double ran = properties_->gsl.GetRanProb();
					p_to_one = p_diffuse_ii_to_both_rest_
								[x_dist_dub][x_dist];
					p_from_one = p_diffuse_ii_from_both_rest_
								[x_dist_dub][x_dist];
					p_tot = p_to_one + p_from_one;
					if(ran < p_to_one/p_tot
							&& n_to_one > 0){
						n_ii_to_both[x_dist_dub][x_dist]--;
						n_to_one = n_ii_to_both[x_dist_dub][x_dist];
					}
					else if(n_from_one > 0){
						n_ii_from_both[x_dist_dub][x_dist]--;
						n_from_one = n_ii_from_both[x_dist_dub][x_dist];
					}
					else if (n_to_one > 0){
						n_ii_to_both[x_dist_dub][x_dist]--;
						n_to_one = n_ii_to_both[x_dist_dub][x_dist];
					}

				}
				n_to_two = n_ii_to_self_from_teth[x_dist_dub][x_dist];
				n_from_two = n_ii_from_self_to_teth[x_dist_dub][x_dist];
				n_avail_two = n_sites_ii_tethered_same_[x_dist_dub][x_dist];
				// Make sure there aren't more events than sites w/ same equil
				while((n_to_two + n_from_two) > n_avail_two){
//					printf("yah for to/from self/teth_ext (2x: %i, x: %i)\n",
//								x_dist_dub, x_dist);
//					printf("n_to: %i n_from: %i n_avail: %i\n", 
//								n_to_two, n_from_two, n_avail_two);
					double ran = properties_->gsl.GetRanProb();
					p_to_two = p_diffuse_ii_to_self_from_teth_
								[x_dist_dub][x_dist];
					p_from_two = p_diffuse_ii_from_self_to_teth_
								[x_dist_dub][x_dist];
					p_tot = p_to_two + p_from_two;
					if(ran < p_to_two/p_tot
					&& n_to_two > 0){
						n_ii_to_self_from_teth[x_dist_dub][x_dist]--;
						n_to_two = n_ii_to_self_from_teth
									[x_dist_dub][x_dist];
					}
					else if(n_from_two > 0){
						n_ii_from_self_to_teth[x_dist_dub][x_dist]--;
						n_from_two = n_ii_from_self_to_teth
									[x_dist_dub][x_dist];
					}
					else if(n_to_two > 0){
						n_ii_to_self_from_teth[x_dist_dub][x_dist]--;
						n_to_two = n_ii_to_self_from_teth
									[x_dist_dub][x_dist];
					}
				}
				n_to = n_to_one + n_to_two;
				n_from = n_from_one + n_from_two;
				n_avail = n_avail_one + n_avail_two;
				// Make sure there aren't more events than ALL sites
				while(2*(n_to + n_from) > n_avail){
//					printf("yah for teth_ii ext FINALE comp (2x: %i, x: %i)\n",
//								x_dist_dub, x_dist);
//					printf("n_to: %i n_from: %i n_avail: %i\n", 
//								n_to, n_from, n_avail);
					double ran = properties_->gsl.GetRanProb();
					p_to_one = p_diffuse_ii_to_both_rest_
								[x_dist_dub][x_dist];
					p_from_one = p_diffuse_ii_from_both_rest_
								[x_dist_dub][x_dist];
					p_to_two = p_diffuse_ii_to_self_from_teth_
								[x_dist_dub][x_dist];
					p_from_two = p_diffuse_ii_from_self_to_teth_
								[x_dist_dub][x_dist];
					p_to = p_to_one + p_to_two;
					p_from = p_from_one + p_from_two;
					p_tot = p_to + p_from;
//					printf("%g, %g, %g, %g\n", p_to_one, p_to_two, 
//							p_from_one, p_from_two);
					if(ran < p_to_one/p_tot
					&& n_to_one > 0){
						n_ii_to_both[x_dist_dub][x_dist]--;
						n_to_one = n_ii_to_both[x_dist_dub][x_dist];
						n_to = n_to_one + n_to_two;
					}
					else if(ran < p_to/p_tot
					&& n_to_two > 0){
						n_ii_to_self_from_teth[x_dist_dub][x_dist]--;
						n_to_two = n_ii_to_self_from_teth
									[x_dist_dub][x_dist];
						n_to = n_to_one + n_to_two;
					}
					else if(ran < (p_to + p_from_one)/p_tot
					&& n_from_one > 0){
						n_ii_from_both[x_dist_dub][x_dist]--;
						n_from_one = n_ii_from_both[x_dist_dub][x_dist];
						n_from = n_from_one + n_from_two;
					}
					else if(n_from_two > 0){
						n_ii_from_self_to_teth[x_dist_dub][x_dist]--;
						n_from_two = n_ii_from_self_to_teth
									[x_dist_dub][x_dist];
						n_from = n_from_one + n_from_two;
					}
				}
			}
			// Extension
			else{
				n_to_one = n_ii_to_both[x_dist_dub][x_dist];
				n_from_one = n_ii_from_both[x_dist_dub][x_dist];
				n_avail_one = n_sites_ii_tethered_same_[x_dist_dub][x_dist];
				// Make sure there aren't more events than sites w/ same equil
				while((n_to_one + n_from_one) > n_avail_one){
//					printf("yah for to/from both_ext (2x: %i, x: %i)\n", 
//							x_dist_dub, x_dist);
//					printf("n_to: %i n_from: %i n_avail: %i\n", 
//							n_to_one, n_from_one, n_avail_one);
					double ran = properties_->gsl.GetRanProb();
					p_to_one = p_diffuse_ii_to_both_rest_
								[x_dist_dub][x_dist];
					p_from_one = p_diffuse_ii_from_both_rest_
								[x_dist_dub][x_dist];
					p_tot = p_to_one + p_from_one;
					if(ran < p_to_one/p_tot
							&& n_to_one > 0){
						n_ii_to_both[x_dist_dub][x_dist]--;
						n_to_one = n_ii_to_both[x_dist_dub][x_dist];
					}
					else if(n_from_one > 0){
						n_ii_from_both[x_dist_dub][x_dist]--;
						n_from_one = n_ii_from_both[x_dist_dub][x_dist];
					}
					else if (n_to_one > 0){
						n_ii_to_both[x_dist_dub][x_dist]--;
						n_to_one = n_ii_to_both[x_dist_dub][x_dist];
					}

				}
				n_to_two = n_ii_to_self_from_teth[x_dist_dub][x_dist];
				n_from_two = n_ii_from_self_to_teth[x_dist_dub][x_dist];
				n_avail_two = n_sites_ii_tethered_oppo_[x_dist_dub][x_dist];
				// Make sure there aren't more events than sites w/ oppo equil
				while((n_to_two + n_from_two) > n_avail_two){
//					printf("yah for to/from self/teth_ext (2x: %i, x: %i)\n",
//								x_dist_dub, x_dist);
//					printf("n_to: %i n_from: %i n_avail: %i\n", 
//								n_to_two, n_from_two, n_avail_two);
					double ran = properties_->gsl.GetRanProb();
					p_to_two = p_diffuse_ii_to_self_from_teth_
								[x_dist_dub][x_dist];
					p_from_two = p_diffuse_ii_from_self_to_teth_
								[x_dist_dub][x_dist];
					p_tot = p_to_two + p_from_two;
					if(ran < p_to_two/p_tot
					&& n_to_two > 0){
						n_ii_to_self_from_teth[x_dist_dub][x_dist]--;
						n_to_two = n_ii_to_self_from_teth
									[x_dist_dub][x_dist];
					}
					else if(n_from_two > 0){
						n_ii_from_self_to_teth[x_dist_dub][x_dist]--;
						n_from_two = n_ii_from_self_to_teth
									[x_dist_dub][x_dist];
					}
					else if(n_to_two > 0){
						n_ii_to_self_from_teth[x_dist_dub][x_dist]--;
						n_to_two = n_ii_to_self_from_teth
									[x_dist_dub][x_dist];
					}
				}
				n_to = n_to_one + n_to_two;
				n_from = n_from_one + n_from_two;
				n_avail = n_avail_one + n_avail_two;
				// Make sure there aren't more events than ALL sites
				while(2*(n_to + n_from) > n_avail){
//					printf("yah for teth_ii ext FINALE ext (2x: %i, x: %i)\n",
//								x_dist_dub, x_dist);
//					printf("n_to: %i n_from: %i n_avail: %i\n", 
//								n_to, n_from, n_avail);
					double ran = properties_->gsl.GetRanProb();
					p_to_one = p_diffuse_ii_to_both_rest_
								[x_dist_dub][x_dist];
					p_from_one = p_diffuse_ii_from_both_rest_
								[x_dist_dub][x_dist];
					p_to_two = p_diffuse_ii_to_self_from_teth_
								[x_dist_dub][x_dist];
					p_from_two = p_diffuse_ii_from_self_to_teth_
								[x_dist_dub][x_dist];
					p_to = p_to_one + p_to_two;
					p_from = p_from_one + p_from_two;
					p_tot = p_to + p_from;
					if(ran < p_to_one/p_tot
					&& n_to_one > 0){
						n_ii_to_both[x_dist_dub][x_dist]--;
						n_to_one = n_ii_to_both[x_dist_dub][x_dist];
						n_to = n_to_one + n_to_two;
					}
					else if(ran < p_to/p_tot
					&& n_to_two > 0){
						n_ii_to_self_from_teth[x_dist_dub][x_dist]--;
						n_to_two = n_ii_to_self_from_teth
									[x_dist_dub][x_dist];
						n_to = n_to_one + n_to_two;
					}
					else if(ran < (p_to + p_from_one)/p_tot
					&& n_from_one > 0){
						n_ii_from_both[x_dist_dub][x_dist]--;
						n_from_one = n_ii_from_both[x_dist_dub][x_dist];
						n_from = n_from_one + n_from_two;
					}
					else if(n_from_two > 0){
						n_ii_from_self_to_teth[x_dist_dub][x_dist]--;
						n_from_two = n_ii_from_self_to_teth
									[x_dist_dub][x_dist];
						n_from = n_from_one + n_from_two;
					}
				}
			}
			n_events += n_ii_to_both[x_dist_dub][x_dist];
			n_events += n_ii_from_both[x_dist_dub][x_dist];
			n_events += n_ii_to_self_from_teth[x_dist_dub][x_dist];
			n_events += n_ii_from_self_to_teth[x_dist_dub][x_dist];
		}
	}
	// Only generate a list if more than zero events occured
	if(n_events > 0){
		int pre_list[n_events];
		int diff_index = 0;
		for(int i_event = 0; i_event < n_i_fwd; i_event++){
			pre_list[diff_index] = 10;
			diff_index++;
		}
		for(int i_event = 0; i_event < n_i_bck; i_event++){
			pre_list[diff_index] = 20;
			diff_index++;
		}
		for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
			int n_step_to = n_ii_to[x_dist];
			for(int i_event = 0; i_event < n_step_to; i_event++){
				pre_list[diff_index] = 300 + x_dist;
				diff_index++;
			}
			int n_step_from = n_ii_from[x_dist];
			for(int i_event = 0; i_event < n_step_from; i_event++){
				pre_list[diff_index] = 400 + x_dist;
				diff_index++;
			}
		}
		for(int x_dist_dub = 0; x_dist_dub <= 2*teth_cutoff; x_dist_dub++){
			int n_step_to_i = n_i_to_teth[x_dist_dub];
			for(int i_event = 0; i_event < n_step_to_i; i_event++){
				pre_list[diff_index] = 500 + x_dist_dub;
				diff_index++;
			}
			int n_step_from_i = n_i_from_teth[x_dist_dub];
			for(int i_event = 0; i_event < n_step_from_i; i_event++){
				pre_list[diff_index] = 600 + x_dist_dub;
				diff_index++;
			}
			for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
				int n_step_to_same = n_ii_to_both[x_dist_dub][x_dist];	
				for(int i_event = 0; i_event < n_step_to_same; i_event++){
					pre_list[diff_index] = 7000 + 10*x_dist_dub + x_dist;
					diff_index++;
				}
				int n_step_from_same = n_ii_from_both[x_dist_dub][x_dist];
				for(int i_event = 0; i_event < n_step_from_same; i_event++){
					pre_list[diff_index] = 8000 + 10*x_dist_dub + x_dist;
					diff_index++;
				}
				int n_step_to_oppo = n_ii_to_self_from_teth
						[x_dist_dub][x_dist];
				for(int i_event = 0; i_event < n_step_to_oppo; i_event++){
					pre_list[diff_index] = 9000 +10*x_dist_dub + x_dist;
					diff_index++;
				}
				int n_step_from_oppo = n_ii_from_self_to_teth
						[x_dist_dub][x_dist];
				for(int i_event = 0; i_event < n_step_from_oppo; i_event++){
					pre_list[diff_index] = 10000 + 10*x_dist_dub + x_dist;
					diff_index++;
				}
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

int AssociatedProteinManagement::GetNumToStepI_Forward(){

	int n_bound = n_sites_i_untethered_;
	if(n_bound > 0){
		double p_step = p_diffuse_i_fwd_;
		int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int AssociatedProteinManagement::GetNumToStepI_Backward(){

	int n_bound = n_sites_i_untethered_;
	if(n_bound > 0){
		double p_step = p_diffuse_i_bck_;
		int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int AssociatedProteinManagement::GetNumToStepII_ToRest(int x_dist){

	int n_bound = n_sites_ii_untethered_[x_dist];
	if(n_bound > 0){
		double p_step = p_diffuse_ii_to_rest_[x_dist];
		int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int AssociatedProteinManagement::GetNumToStepII_FromRest(int x_dist){

	int n_bound = n_sites_ii_untethered_[x_dist];
	if(n_bound > 0){
		double p_step = p_diffuse_ii_from_rest_[x_dist];
		int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int AssociatedProteinManagement::GetNumToStepI_ToTethRest(int x_dist_dub){

	int n_bound = n_sites_i_tethered_[x_dist_dub];
	if(n_bound > 0){
		double p_step = p_diffuse_i_to_teth_rest_[x_dist_dub];
		int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int AssociatedProteinManagement::GetNumToStepI_FromTethRest(int x_dist_dub){

	int n_bound = n_sites_i_tethered_[x_dist_dub];
	if(n_bound > 0){
		double p_step = p_diffuse_i_from_teth_rest_[x_dist_dub];
		int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int AssociatedProteinManagement::GetNumToStepII_ToBothRest(int x_dist_dub, 
		int x_dist){

	int rest_dist_dub = 2*properties_->kinesin4.motor_list_[0].rest_dist_;
	int n_bound;
	// this means we compressed, yo
	if(x_dist_dub <= rest_dist_dub){
		n_bound = n_sites_ii_tethered_oppo_[x_dist_dub][x_dist];
	}
	else{
		n_bound = n_sites_ii_tethered_same_[x_dist_dub][x_dist];
	}
	if(n_bound > 0){
		double p_step = p_diffuse_ii_to_both_rest_[x_dist_dub][x_dist];
		int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int AssociatedProteinManagement::GetNumToStepII_FromBothRest(int x_dist_dub, 
		int x_dist){

	int rest_dist_dub = 2*properties_->kinesin4.motor_list_[0].rest_dist_;
	int n_bound;
	// this means we compressed, yo
	if(x_dist_dub <= rest_dist_dub){
		n_bound = n_sites_ii_tethered_oppo_[x_dist_dub][x_dist];
	}
	else{
		n_bound = n_sites_ii_tethered_same_[x_dist_dub][x_dist];
	}
	if(n_bound > 0){
		double p_step = p_diffuse_ii_from_both_rest_[x_dist_dub][x_dist];
		int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int AssociatedProteinManagement::GetNumToStepII_ToSelf_FromTeth
(int x_dist_dub, int x_dist){

	int rest_dist_dub = 2*properties_->kinesin4.motor_list_[0].rest_dist_;
	int n_bound;
	// this means we compressed, yo
	if(x_dist_dub <= rest_dist_dub){
		n_bound = n_sites_ii_tethered_same_[x_dist_dub][x_dist];
	}
	else{
		n_bound = n_sites_ii_tethered_oppo_[x_dist_dub][x_dist];
	}
	if(n_bound > 0){
		double p_step = p_diffuse_ii_to_self_from_teth_[x_dist_dub][x_dist];
		int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int AssociatedProteinManagement::GetNumToStepII_FromSelf_ToTeth
									(int x_dist_dub, int x_dist){

	int rest_dist_dub = 2*properties_->kinesin4.motor_list_[0].rest_dist_;
	int n_bound;
	// this means we compressed, yo
	if(x_dist_dub <= rest_dist_dub){
		n_bound = n_sites_ii_tethered_same_[x_dist_dub][x_dist];
	}
	else{
		n_bound = n_sites_ii_tethered_oppo_[x_dist_dub][x_dist];
	}
	if(n_bound > 0){
		double p_step = p_diffuse_ii_from_self_to_teth_[x_dist_dub][x_dist];
		int n_to_step = properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
}

void AssociatedProteinManagement::RunDiffusion(){

	GenerateDiffusionList();
	if(diffusion_list_.empty() == false){
		int n_events = diffusion_list_.size();
//		printf("%i XLINK DIFFUSION EVENTS\n", n_events);
		int x_dist, 			// refers to dist (in sites) of xlink ext
			x_dist_dub;			// refers to 2*dist ('') of motor ext
		for(int i_event = 0; i_event < n_events; i_event++){
			int diff_event = diffusion_list_[i_event];
			if(diff_event >= 300 && diff_event < 400){
				x_dist = diff_event % 100; 
				diff_event = 30;
			}
			if(diff_event >= 400 && diff_event < 500){
				x_dist = diff_event % 100;
				diff_event = 40;
			}	
			if(diff_event >= 500 && diff_event < 600){
				x_dist_dub = diff_event % 100;
				diff_event = 50;
			}
			if(diff_event >= 600 && diff_event < 700){
				x_dist_dub = diff_event % 100;
				diff_event = 60;
			}
			if(diff_event >= 7000 && diff_event < 8000){
				x_dist = diff_event % 10;
				x_dist_dub = (diff_event % 1000 - x_dist)/10; 
				diff_event = 70;
			}
			if(diff_event >= 8000 && diff_event < 9000){
				x_dist = diff_event % 10;
				x_dist_dub = (diff_event % 1000 - x_dist)/10; 
				diff_event = 80;
			}
			if(diff_event >= 9000 && diff_event < 10000){
				x_dist = diff_event % 10;
				x_dist_dub = (diff_event % 1000 - x_dist)/10; 
				diff_event = 90;
			}
			if(diff_event >= 10000 && diff_event < 11000){
				x_dist = diff_event % 10;
				x_dist_dub = (diff_event % 1000 - x_dist)/10; 
				diff_event = 100;
			}
//			printf("ole:\n");
//			properties_->wallace.PrintMicrotubules(0.000);
			switch(diff_event){
				case 10:
//					printf("unteth_i step fwd\n");
					RunDiffusionI_Forward();
					break;
				case 20:
//					printf("unteth_i step bck\n");
					RunDiffusionI_Backward();
					break;
				case 30:
//						printf("unteth_ii step to rest (%i)[%i avail]\n", 
//							x_dist, n_sites_ii_untethered_[x_dist]);
					RunDiffusionII_ToRest(x_dist);
					break;
				case 40:
//					printf("unteth_ii step from rest (%i)[%i avail]\n", 
//							x_dist, n_sites_ii_untethered_[x_dist]);
					RunDiffusionII_FromRest(x_dist);
					break;
				case 50:
//					printf("teth_i step to rest (%i)\n", x_dist_dub);
					RunDiffusionI_ToTethRest(x_dist_dub);
					break;
				case 60:
//					printf("teth_i step from rest (%i)\n", x_dist_dub);
					RunDiffusionI_FromTethRest(x_dist_dub);
					break;
				case 70:
//					printf("teth_ii step to both (2x: %i, x: %i)\n", 
//						x_dist_dub, x_dist);
					RunDiffusionII_ToBothRest(x_dist_dub, x_dist);
					break;
				case 80:
//					printf("teth_ii step from both (2x: %i, x: %i)\n", 
//							x_dist_dub, x_dist);
					RunDiffusionII_FromBothRest(x_dist_dub, x_dist);
					break;
				case 90:
//					printf("teth_ii step to self from teth (2x: %i, x: %i)\n",
//							x_dist_dub, x_dist);
					RunDiffusionII_ToSelf_FromTeth(x_dist_dub, x_dist);
					break;
				case 100:
//					printf("teth_ii step from self to teth (2x: %i, x: %i)\n", 
//							x_dist_dub, x_dist);
					RunDiffusionII_FromSelf_ToTeth(x_dist_dub, x_dist);
					break;
			}
		}
	}
}

void AssociatedProteinManagement::RunDiffusionI_Forward(){

	UpdateSingleUntetheredSites();
	int n_bound = n_sites_i_untethered_;
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Tubulin* site = single_untethered_sites_[i_entry]; 
		AssociatedProtein* xlink = site->xlink_;
		Microtubule* mt = site->mt_;
		int i_site = site->index_;
		int dx = mt->delta_x_;
		int plus_end = mt->plus_end_;
		int minus_end = mt->minus_end_;
		// cant step off them MTs
		if(!(i_site == plus_end && dx == 1)
		&& !(i_site == minus_end && dx == -1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &mt->lattice_[i_site + 1];
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
			else{
	//			printf("oh well fwd xlink_i\n");
			}
		}
	}
	else{
		printf("come on kid we cant step_i_fwd\n");
		exit(1);
	}
}

void AssociatedProteinManagement::RunDiffusionI_Backward(){

	UpdateSingleUntetheredSites();
	int n_bound = n_sites_i_untethered_;
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Tubulin* site = single_untethered_sites_[i_entry]; 
		AssociatedProtein* xlink = site->xlink_;
		Microtubule* mt = site->mt_;
		int i_site = site->index_;
		int dx = mt->delta_x_;
		int plus_end = mt->plus_end_;
		int minus_end = mt->minus_end_;
		if(!(i_site == minus_end && dx == 1)
		&& !(i_site == plus_end && dx == -1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &mt->lattice_[i_site - 1];
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
			else{
	//			printf("oh well bck xlink_i\n");
			}
		}
	}
	else{
		printf("come on kid we cant step_i_bck\n");
		exit(1);
	}
}

void AssociatedProteinManagement::RunDiffusionII_ToRest(int x_dist){

	int mt_length = parameters_->length_of_microtubule;
	int mt_array_length = mt_length - 1;	
	UpdateDoubleUntetheredSites();
	int n_bound = n_sites_ii_untethered_[x_dist];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site = double_untethered_sites_[x_dist][i_entry];
		int i_site = site->index_;
		AssociatedProtein* xlink = site->xlink_;
		int dx = xlink->GetDirectionTowardRest(site);
		Microtubule* mt = site->mt_;
//		printf("i: %i, dx: %i\n", i_site, dx);
		if(!(i_site == mt_array_length && dx == 1)
		&& !(i_site == 0 && dx == -1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &mt->lattice_[i_site + dx];
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
				int x_dist_pre = xlink->x_dist_;
				xlink->UpdateExtension();
				// Make sure we didn't force an unbind event
				if(xlink->heads_active_ == 2){
					int x_dist_post = xlink->x_dist_;
					n_double_bound_[x_dist_pre]--;
					n_double_bound_[x_dist_post]++;
					n_sites_ii_untethered_[x_dist_pre] -= 2;
					n_sites_ii_untethered_[x_dist_post] += 2;
				}
			}
			else{
	//			printf("oh well towardrest xlink_ii\n");
			}
		}
	}
	else{
		printf("we cannot step_ii_toward (xlink UT)\n");
		exit(1);
	}
}

void AssociatedProteinManagement::RunDiffusionII_FromRest(int x_dist){

	int mt_length = parameters_->length_of_microtubule;
	int mt_array_length = mt_length - 1;	
	UpdateDoubleUntetheredSites();
	int n_bound = n_sites_ii_untethered_[x_dist];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site = double_untethered_sites_[x_dist][i_entry];
		int i_site = site->index_;
		AssociatedProtein* xlink = site->xlink_;
		int dx = xlink->GetDirectionTowardRest(site);
		Microtubule* mt = site->mt_;
		if(!(i_site == mt_array_length && dx == -1)
		&& !(i_site == 0 && dx == 1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &mt->lattice_[i_site - dx];
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
				int x_dist_pre = xlink->x_dist_;
				xlink->UpdateExtension();
				// Make sure we didn't force an unbind event
				if(xlink->heads_active_ == 2){
					int x_dist_post = xlink->x_dist_;
					n_double_bound_[x_dist_pre]--;
					n_double_bound_[x_dist_post]++;
					n_sites_ii_untethered_[x_dist_pre] -= 2;
					n_sites_ii_untethered_[x_dist_post] += 2;
				}
			}
			else{
	//			printf("oh well fromrest xlink_ii\n");
			}
		}
	}
	else{
		printf("we cannot step_ii_from (xlink UT)\n");
		exit(1);
	}
}

void AssociatedProteinManagement::RunDiffusionI_ToTethRest(int x_dist_dub){

	int mt_length = parameters_->length_of_microtubule; 
	int mt_array_length = mt_length - 1;
	UpdateSingleTetheredSites();
	int n_bound = n_sites_i_tethered_[x_dist_dub];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site = single_tethered_sites_[x_dist_dub][i_entry];
		int i_site = site->index_;
		AssociatedProtein *xlink = site->xlink_;
		Kinesin *motor = xlink->motor_;
		// direction xlink needs to move to rest is opposite of motor's
		int dx = -1 * motor->GetDirectionTowardRest(); 
		Microtubule *mt = site->mt_;
		if(!(i_site == mt_array_length && dx == 1)
		&& !(i_site == 0 && dx == -1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &mt->lattice_[i_site + dx];
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
				int x_dub_pre = motor->x_dist_doubled_;
				motor->UpdateExtension();
				// Make sure we didn't force an untether event
				if(motor->tethered_ == true){
					int x_dub_post = motor->x_dist_doubled_;
					n_sites_i_tethered_[x_dub_pre]--;
					n_sites_i_tethered_[x_dub_post]++;
					// Update kinesin4 stats
					KinesinManagement *kinesin4 = &properties_->kinesin4;
					if(motor->heads_active_ == 2){
						kinesin4->n_bound_tethered_[x_dub_pre]--;
						kinesin4->n_bound_tethered_[x_dub_post]++;
					}
				}
			}
			else{
	//			printf("oh well fromrest xlink_ii\n");
			}
		}
	}
	else{
		printf("cant work under these cndtns. xlink RD_i_to_tethrest\n");
		exit(1);
	}
}

void AssociatedProteinManagement::RunDiffusionI_FromTethRest(int x_dist_dub){

	int mt_length = parameters_->length_of_microtubule; 
	int mt_array_length = mt_length - 1;
	UpdateSingleTetheredSites();
	int n_bound = n_sites_i_tethered_[x_dist_dub];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site = single_tethered_sites_[x_dist_dub][i_entry];
		int i_site = site->index_;
		AssociatedProtein *xlink = site->xlink_;
		Kinesin *motor = xlink->motor_;
		if(motor == nullptr){
			printf("why\n");
			exit(1);
		}	
		// direction xlink needs to move to rest is opposite of motor's
		int dx = -1 * motor->GetDirectionTowardRest(); 
		Microtubule *mt = site->mt_;
		if(!(i_site == mt_array_length && dx == -1)
		&& !(i_site == 0 && dx == 1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &mt->lattice_[i_site - dx];
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
				int x_dub_pre = motor->x_dist_doubled_;
				motor->UpdateExtension();
				// Make sure we didn't force an untether event
				if(motor->tethered_ == true){
					int x_dub_post = motor->x_dist_doubled_;
					n_sites_i_tethered_[x_dub_pre]--;
					n_sites_i_tethered_[x_dub_post]++;
					// Update kinesin4 stats
					KinesinManagement *kinesin4 = &properties_->kinesin4;
					if(motor->heads_active_ == 2){
						kinesin4->n_bound_tethered_[x_dub_pre]--;
						kinesin4->n_bound_tethered_[x_dub_post]++;
					}
				}
			}
			else{
	//			printf("oh well fromrest xlink_ii\n");
			}
		}
	}
	else{
		printf("cant work under these cndtns. xlink RD_i_to_tethrest\n");
		exit(1);
	}
}

void AssociatedProteinManagement::RunDiffusionII_ToBothRest
										(int x_dist_dub, int x_dist){

	int mt_length = parameters_->length_of_microtubule;
	int mt_array_length = mt_length - 1;
	int rest_dist_dub = 2*properties_->kinesin4.motor_list_[0].rest_dist_;
	UpdateDoubleTetheredSites();
	int n_bound;
	// this means we compressed, yo
	if(x_dist_dub <= rest_dist_dub){
		n_bound = n_sites_ii_tethered_oppo_[x_dist_dub][x_dist];
	}
	else{
		n_bound = n_sites_ii_tethered_same_[x_dist_dub][x_dist];
	}
	if(n_bound > 0){
		int i = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site;
		if(x_dist_dub <= rest_dist_dub){
			site = double_tethered_sites_oppo_[x_dist_dub][x_dist][i];
		}
		else{
	   		site = double_tethered_sites_same_[x_dist_dub][x_dist][i];
		}
		Microtubule *mt = site->mt_;
		int i_site = site->index_;
		AssociatedProtein *xlink = site->xlink_;
		Kinesin *motor = xlink->motor_;
		int dx = -1 * motor->GetDirectionTowardRest();
		if(!(i_site == mt_array_length && dx == 1)
		&& !(i_site == 0 && dx == -1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &mt->lattice_[i_site + dx];
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
				int x_pre = xlink->x_dist_;
				xlink->UpdateExtension();
				// Make sure we didn't force an unbinding event
				if(xlink->heads_active_ == 2){
					int x_post = xlink->x_dist_;
					n_double_bound_[x_pre]--;
					n_double_bound_[x_post]++;
					int x_dub_pre = motor->x_dist_doubled_;
					motor->UpdateExtension();
					// Make sure we didn't force an untether event
					if(motor->tethered_ == true){
						int x_dub_post = motor->x_dist_doubled_;
						n_sites_ii_tethered_[x_dub_pre][x_pre] -= 2;
						n_sites_ii_tethered_[x_dub_post][x_post] += 2;
						KinesinManagement *kinesin4 = &properties_->kinesin4;
						if(motor->heads_active_ == 2){
							kinesin4->n_bound_tethered_[x_dub_pre]--;
							kinesin4->n_bound_tethered_[x_dub_post]++;
						}
					}
					// If we did, correct stats
					else{
						n_sites_ii_tethered_[x_dub_pre][x_post] += 2;
						n_sites_ii_tethered_[x_dub_pre][x_pre] -= 2;
					}
				}
			}
			else{
	//			printf("oh well to_bothrest xlink_ii\n");
			}
		}
	}
	else{
		printf("issues w/ run_diffusion_ii_to_both_rest\n");
		exit(1);
	}
}

void AssociatedProteinManagement::RunDiffusionII_FromBothRest
										(int x_dist_dub, int x_dist){

	int mt_length = parameters_->length_of_microtubule;
	int mt_array_length = mt_length - 1;
	int rest_dist_dub = 2*properties_->kinesin4.motor_list_[0].rest_dist_;
	UpdateDoubleTetheredSites();
	int n_bound;
	// this means we compressed, yo
	if(x_dist_dub <= rest_dist_dub){
		n_bound = n_sites_ii_tethered_oppo_[x_dist_dub][x_dist];
	}
	else{
		n_bound = n_sites_ii_tethered_same_[x_dist_dub][x_dist];
	}
	if(n_bound > 0){
		int i = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site;
		if(x_dist_dub <= rest_dist_dub){
			site = double_tethered_sites_oppo_[x_dist_dub][x_dist][i];
		}
		else{
	   		site = double_tethered_sites_same_[x_dist_dub][x_dist][i];
		}
		Microtubule *mt = site->mt_;
		int i_site = site->index_;
		AssociatedProtein *xlink = site->xlink_;
		Kinesin *motor = xlink->motor_;
		// xlink
		int dx = motor->GetDirectionTowardRest();
		if(!(i_site == mt_array_length && dx == 1)
		&& !(i_site == 0 && dx == -1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &mt->lattice_[i_site + dx];
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
				int x_pre = xlink->x_dist_;
				xlink->UpdateExtension();
				// Make sure we didn't force an unbinding event
				if(xlink->heads_active_ == 2){
					int x_post = xlink->x_dist_;
					n_double_bound_[x_pre]--;
					n_double_bound_[x_post]++;
					int x_dub_pre = motor->x_dist_doubled_;
					motor->UpdateExtension();
					// Make sure we didn't force an untether event
					if(motor->tethered_ == true){
						int x_dub_post = motor->x_dist_doubled_;
						n_sites_ii_tethered_[x_dub_pre][x_pre] -= 2;
						n_sites_ii_tethered_[x_dub_post][x_post] += 2;
						KinesinManagement *kinesin4 = &properties_->kinesin4;
						if(motor->heads_active_ == 2){
							kinesin4->n_bound_tethered_[x_dub_pre]--;
							kinesin4->n_bound_tethered_[x_dub_post]++;
						}
					}
					// If we did, correct stats
					else{
						n_sites_ii_tethered_[x_dub_pre][x_post] += 2;
						n_sites_ii_tethered_[x_dub_pre][x_pre] -= 2;
					}
				}
			}
			else{
	//			printf("oh well to_fromrest xlink_ii\n");
			}
		}
	}
	else{
		printf("issues w/ run_diffusion_ii_from_both_rest\n");
		exit(1);
	}
}

void AssociatedProteinManagement::RunDiffusionII_ToSelf_FromTeth
										(int x_dist_dub, int x_dist){

	int mt_length = parameters_->length_of_microtubule;
	int mt_array_length = mt_length - 1;
	int rest_dist_dub = 2*properties_->kinesin4.motor_list_[0].rest_dist_;
	UpdateDoubleTetheredSites();
	int n_bound;
	// this means we compressed, yo
	if(x_dist_dub <= rest_dist_dub){
		n_bound = n_sites_ii_tethered_same_[x_dist_dub][x_dist];
	}
	else{
		n_bound = n_sites_ii_tethered_oppo_[x_dist_dub][x_dist];
	}
	if(n_bound > 0){
		int i = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site;
		if(x_dist_dub <= rest_dist_dub){
			site = double_tethered_sites_same_[x_dist_dub][x_dist][i];
		}
		else{
	   		site = double_tethered_sites_oppo_[x_dist_dub][x_dist][i];
		}
		Microtubule *mt = site->mt_;
		int i_site = site->index_;
		AssociatedProtein *xlink = site->xlink_;
		Kinesin *motor = xlink->motor_;
		int dx = motor->GetDirectionTowardRest();
		if(!(i_site == mt_array_length && dx == 1)
		&& !(i_site == 0 && dx == -1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &mt->lattice_[i_site + dx];
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
				int x_pre = xlink->x_dist_;
				xlink->UpdateExtension();
				// Make sure we didn't force an unbinding event
				if(xlink->heads_active_ == 2){
					int x_post = xlink->x_dist_;
					n_double_bound_[x_pre]--;
					n_double_bound_[x_post]++;
					int x_dub_pre = motor->x_dist_doubled_;
					motor->UpdateExtension();
					// Make sure we didn't force an untether event
					if(motor->tethered_ == true){
						int x_dub_post = motor->x_dist_doubled_;
						n_sites_ii_tethered_[x_dub_pre][x_pre] -= 2;
						n_sites_ii_tethered_[x_dub_post][x_post] += 2;
						KinesinManagement *kinesin4 = &properties_->kinesin4;
						if(motor->heads_active_ == 2){
							kinesin4->n_bound_tethered_[x_dub_pre]--;
							kinesin4->n_bound_tethered_[x_dub_post]++;
						}
					}
					// If we did, correct stats
					else{
						n_sites_ii_tethered_[x_dub_pre][x_post] += 2;
						n_sites_ii_tethered_[x_dub_pre][x_pre] -= 2;
					}
				}
			}
			else{
	//			printf("oh well toself_fromteth xlink_ii\n");
			}
		}
	}
	else{
		printf("issues w/ run_diffusion_ii_toself_fromteth_rest\n");
		exit(1);
	}
}

void AssociatedProteinManagement::RunDiffusionII_FromSelf_ToTeth
										(int x_dist_dub, int x_dist){

	int mt_length = parameters_->length_of_microtubule;
	int mt_array_length = mt_length - 1;
	int rest_dist_dub = 2*properties_->kinesin4.motor_list_[0].rest_dist_;
	UpdateDoubleTetheredSites();
	int n_bound;
	// this means we compressed, yo
	if(x_dist_dub <= rest_dist_dub){
		n_bound = n_sites_ii_tethered_same_[x_dist_dub][x_dist];
	}
	else{
		n_bound = n_sites_ii_tethered_oppo_[x_dist_dub][x_dist];
	}
	if(n_bound > 0){
		int i = properties_->gsl.GetRanInt(n_bound);
		Tubulin *site;
		if(x_dist_dub <= rest_dist_dub){
			site = double_tethered_sites_same_[x_dist_dub][x_dist][i];
		}
		else{
	   		site = double_tethered_sites_oppo_[x_dist_dub][x_dist][i];
		}
		Microtubule *mt = site->mt_;
		int i_site = site->index_;
		AssociatedProtein *xlink = site->xlink_;
		Kinesin *motor = xlink->motor_;
		int dx = -1 * motor->GetDirectionTowardRest();
		if(!(i_site == mt_array_length && dx == 1)
		&& !(i_site == 0 && dx == -1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &mt->lattice_[i_site + dx];
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
				int x_pre = xlink->x_dist_;
				xlink->UpdateExtension();
				// Make sure we didn't force an unbinding event
				if(xlink->heads_active_ == 2){
					int x_post = xlink->x_dist_;
					n_double_bound_[x_pre]--;
					n_double_bound_[x_post]++;
					int x_dub_pre = motor->x_dist_doubled_;
					motor->UpdateExtension();
					// Make sure we didn't force an untether event
					if(motor->tethered_ == true){
						int x_dub_post = motor->x_dist_doubled_;
						n_sites_ii_tethered_[x_dub_pre][x_pre] -= 2;
						n_sites_ii_tethered_[x_dub_post][x_post] += 2;
						KinesinManagement *kinesin4 = &properties_->kinesin4;
						if(motor->heads_active_ == 2){
							kinesin4->n_bound_tethered_[x_dub_pre]--;
							kinesin4->n_bound_tethered_[x_dub_post]++;
						}
					}
					// If we did, correct stats
					else{
						n_sites_ii_tethered_[x_dub_pre][x_post] += 2;
						n_sites_ii_tethered_[x_dub_pre][x_pre] -= 2;
					}
				}
			}
			else{
	//			printf("oh well fromself_toteth xlink_ii\n");
			}
		}
	}
	else{
		printf("issues w/ run_diffusion_ii_fromself_toteth_rest\n");
		exit(1);
	}
}

void AssociatedProteinManagement::GenerateKMCList(){

	int n_events = 0;
	int n_to_bind_i = GetNumToBind_I();
	n_events += n_to_bind_i; 
	int n_to_bind_ii = GetNumToBind_II();
	n_events += n_to_bind_ii;
	int n_to_unbind_i = GetNumToUnbind_I();
	n_events += n_to_unbind_i; 
	while(n_to_bind_ii + n_to_unbind_i > n_single_bound_){
		if(n_to_bind_ii > 0){
			n_to_bind_ii--;
			n_events--;
		}
		else{
			n_to_unbind_i--;
			n_events--;
		}
	}
	int n_to_unbind_ii[dist_cutoff_ + 1]; 
	for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
		n_to_unbind_ii[x_dist] = GetNumToUnbind_II(x_dist);
		n_events += n_to_unbind_ii[x_dist];
	}
	int kmc_index = 0;
	if(n_events > 0){
		int pre_list[n_events];
		// "0" corresponds to a stage 1 binding event
		for(int i_event = 0; i_event < n_to_bind_i; i_event++){
			pre_list[kmc_index] = 0;
			kmc_index++;
		}
		// "1" corresponds to a stage 2 binding event
		for(int i_event = 0; i_event < n_to_bind_ii; i_event++){
			pre_list[kmc_index] = 1;
			kmc_index++;
		}
		// "2" corresponds to a stage 1 (single head) unbinding event
		for(int i_event = 0; i_event < n_to_unbind_i; i_event++){
			pre_list[kmc_index] = 2;
			kmc_index++;
		}
		// "30 + x" corresponds to a stage 2 unbinding event, where
		// x is the specific xlink x_dist we wish to unbind 
		for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
			// Adjust n_events to include all events up to unbinding
			// of this x_dist, but no distances past it
			int n_to_unbind = n_to_unbind_ii[x_dist];
			int kmc_event = 30 + x_dist;
			for(int i_event = 0; i_event < n_to_unbind; i_event++){
				pre_list[kmc_index] = kmc_event;
				kmc_index++; 
			}
		}
		// Shuffle using GSL (why an array is necessary in the 1st place)
		gsl_ran_shuffle(properties_->gsl.rng,pre_list,n_events,sizeof(int));
		kmc_list_.resize(n_events);
		// Transfer shuffled array into our class kmc vector 
		for(int j_event = 0; j_event < n_events; j_event++){
			kmc_list_[j_event] = pre_list[j_event];
		}
	}
	else{
		kmc_list_.clear();
	}
}

int AssociatedProteinManagement::GetNumToBind_I(){
	
	properties_->microtubules.UpdateNumUnoccupied();
	int n_unocc = properties_->microtubules.n_unoccupied_;
	double p_bind = p_bind_i_;
	int n_to_bind = properties_->gsl.SampleBinomialDist(p_bind, n_unocc);
	return n_to_bind;
}

int AssociatedProteinManagement::GetNumToBind_II(){
	
	UpdateSingleBoundList();
	double weights_summed = 0;
	// Sum over all single-bound xlinks
	for(int i_xlink = 0; i_xlink < n_single_bound_; i_xlink++){
		AssociatedProtein *xlink = single_bound_list_[i_xlink];
		xlink->UpdateNeighborSites();
		// Get weight of every possible orientation w/ neighbor
		int n_neighbors = xlink->n_neighbor_sites_;
		for(int i_neighb = 0; i_neighb < n_neighbors; i_neighb++){
			Tubulin *site = xlink->neighbor_sites_[i_neighb];
			double weight = xlink->GetBindingWeight(site);
			weights_summed += weight;
		}
	}
	// Scale summed weights by an effective conc. to get  n_bound at equil
	double k_on = parameters_->k_on_xlink;
	double delta_t = parameters_->delta_t;
	int n_equil = (int) weights_summed; 
//	printf("yo n_equil is %i\n", n_equil);
	double p_bind = k_on * c_eff_ * delta_t;
//	printf("yo p_bind is %g\n", c_eff_);
	int n_to_bind = properties_->gsl.SampleBinomialDist(p_bind, n_equil);
//	printf("yo n_to_bind is %i\n", n_to_bind);
	return n_to_bind;
}

int AssociatedProteinManagement::GetNumToUnbind_I(){

	int n_to_unbind = properties_->gsl.SampleBinomialDist(p_unbind_i_, 
			n_single_bound_);
	return n_to_unbind;
}

// Input extension of crosslinker (in # of sites) and get expected unbind #
int AssociatedProteinManagement::GetNumToUnbind_II(int x_dist){

	double p_unbind = p_unbind_ii_[x_dist];
	int n_bound = n_double_bound_[x_dist];
	int n_to_unbind = properties_->gsl.SampleBinomialDist(p_unbind, n_bound);
	return n_to_unbind;
}

void AssociatedProteinManagement::RunKMC(){
	
	int x_dist;		// extension of xlink for stage2 unbinding
//	printf("Start of xlink KMC cycle\n");
	GenerateKMCList();
	if(kmc_list_.empty() == false){
		int n_events = kmc_list_.size();
//		printf("%i XLINK KMC EVENTS\n", n_events);
		for(int i_event = 0; i_event < n_events; i_event++){
			int kmc_event = kmc_list_[i_event];
			// Allows us to encode xlink extension as second digit
			if(kmc_event >= 30 
			&& kmc_event < 40){
				x_dist = kmc_event % 10;
				kmc_event = 3;
			}
//			properties_->wallace.PrintMicrotubules(0.000);
			switch(kmc_event){
				case 0: 
//						printf("xlink bound stage 1\n");
						RunKMC_Bind_I();
						break;
				case 1: 
//						printf("xlink bound stage 2\n");
						RunKMC_Bind_II();
						break;
				case 2: 
//						printf("xlink fully unbound\n");
						RunKMC_Unbind_I();
						break;
				case 3: 
//						printf("xlink unbound 2nd head (ext was %i)\n", 
//								x_dist);
						RunKMC_Unbind_II(x_dist);
						break;
			}
		}
	}
}

void AssociatedProteinManagement::RunKMC_Bind_I(){

	// Make sure unoccupied sites are available
	properties_->microtubules.UpdateUnoccupiedList();
	if(properties_->microtubules.n_unoccupied_ > 0){
		// Randomly choose an unbound xlink
		int i_xlink = properties_->gsl.GetRanInt(n_xlinks_);
		AssociatedProtein *xlink = &xlink_list_[i_xlink];
		while(xlink->heads_active_ != 0){
			i_xlink++;
			if(i_xlink == n_xlinks_)
				i_xlink = 0;
			xlink = &xlink_list_[i_xlink];
		}
		// Get random unoccupied site
		Tubulin *site = properties_->microtubules.GetUnoccupiedSite();
		// Place xlink onto site
		site->xlink_ = xlink;
		site->occupied_ = true;
		// Update xlink details
		xlink->heads_active_++;
		xlink->site_one_ = site; 
		// Update statistics
		n_single_bound_++;
		if(xlink->tethered_ == true){
			n_tethered_++;
			xlink->motor_->UpdateExtension();
			// Tethers attached to free motors diffuse as untethered
			if(xlink->motor_->heads_active_ == 0){
				n_sites_i_untethered_++;
			}
			// Make sure we didn't force an unteter
			else if(xlink->tethered_ == true){
				int x_dist_dub = xlink->motor_->x_dist_doubled_;
				n_sites_i_tethered_[x_dist_dub]++;
			}
			else{
				printf("hummm in bind_i (XLINKS)\n");
				exit(1);
			}
		}
		else{
			n_untethered_++;
			n_sites_i_untethered_++;
		}
	}
	else{
		printf("Error in RunKMC_BindFirst: no unoccupied sites\n");
//		exit(1);
	}
}

void AssociatedProteinManagement::RunKMC_Bind_II(){

	// Make sure stage 1 xlinks and unoccupied sites are available
	UpdateSingleBoundList();
	properties_->microtubules.UpdateUnoccupiedList();
	if(n_single_bound_ > 0
	&& properties_->microtubules.n_unoccupied_ > 0){
		// Randomly pick single-bound xlink
		int i_xlink = properties_->gsl.GetRanInt(n_single_bound_);
		AssociatedProtein *xlink = single_bound_list_[i_xlink];
		xlink->UpdateNeighborSites();
		// Sample normal distribution for x-dist of xlink to insert 
		int x_dist = xlink->SampleSpringExtension();
//		printf("x_dist be %i\n", x_dist);
		double ran = properties_->gsl.GetRanProb();
//		x_dist = 0;
//		if(ran < 0.5)
//			x_dist = -1;
//		else //if(ran < 0.5)
//			x_dist = 1;
//		else
//			x_dist = 0;

		// Ensure xlink has a neighbor with desired distance
		int attempts = 0;
		int switches = 0;
		bool failure = false;
		while(xlink->NeighborExists(x_dist) == false){
			// roll for a new distance after a certain # of tries
			if(attempts > 2*n_single_bound_){
				if(switches > 50){
//					printf("failed to bind_ii xlink\n");
					failure = true;
					break;
				}
				int x_dist_pre = x_dist;
				while(x_dist == x_dist_pre){
					x_dist = xlink->SampleSpringExtension();
				}
				switches++;
				attempts = 0;
			}
			i_xlink = properties_->gsl.GetRanInt(n_single_bound_);
			xlink = single_bound_list_[i_xlink];	
			xlink->UpdateNeighborSites();
			attempts++;
		}
		Tubulin *site = xlink->GetNeighborSite(x_dist);
		if(failure == false){
			// Place  xlink onto site
			site->xlink_ = xlink;
			site->occupied_ = true;
			// Update xlink details
			xlink->heads_active_++;
			if(xlink->site_one_ == nullptr)
				xlink->site_one_ = site;
			else if(xlink->site_two_ == nullptr)
				xlink->site_two_ = site;
			else{
				printf("bruhhhhhhhh");
				exit(1);
			}
			xlink->UpdateExtension();
			// Take absolute value of x_dist for array access/etc
			x_dist = abs(x_dist);
			if(xlink->x_dist_ != x_dist){
				printf("why thoooooouaz: ");
				printf(" %i in xlink, %i in kmc\n", xlink->x_dist_, x_dist);
				exit(1);
			}
			// Update statistics
			n_double_bound_[x_dist]++;
			n_single_bound_--;
			if(xlink->tethered_ == true){
				Kinesin *motor = xlink->motor_;
				int x_dub_pre = motor->x_dist_doubled_;
				motor->UpdateExtension();
				if(motor->heads_active_ == 0){
					n_sites_ii_untethered_[x_dist] += 2;	
					n_sites_i_untethered_--;
				}
				// Make sure an untethering event wasn't forced
				else if(motor->tethered_ == true){
					int x_dub_post = xlink->motor_->x_dist_doubled_;
					if(motor->heads_active_ == 2){
						properties_->kinesin4.n_bound_tethered_[x_dub_pre]--;
						properties_->kinesin4.n_bound_tethered_[x_dub_post]++;
					}
					n_sites_ii_tethered_[x_dub_post][x_dist] += 2;
					n_sites_i_tethered_[x_dub_pre]--;
				}
				// If it was, counteract statistic change 
				else{
					n_sites_ii_tethered_[x_dub_pre][x_dist] += 2;
					n_sites_i_tethered_[x_dub_pre]--;
					printf("XLINK BIND_II ???? *** \n");
				}
			}
			else{
				n_sites_ii_untethered_[x_dist] += 2;	
				n_sites_i_untethered_--;
			}
		}
	}
	else{
		printf("Error in RunKMC_BindSecond: no unoccupied sites\n");
//		exit(1);
	}
}

void AssociatedProteinManagement::RunKMC_Unbind_I(){
	
	UpdateSingleBoundList();
	if(n_single_bound_ > 0){
		// Randomly pick a single-bound xlink
		int i_entry = properties_->gsl.GetRanInt(n_single_bound_);
		AssociatedProtein *xlink = single_bound_list_[i_entry];
		Tubulin *site = xlink->GetActiveHeadSite();
		if(xlink->heads_active_ != 1){
			printf("nope. xlink unbind_i\n");
			exit(1);
		}
		// Remove xlink from site
		site->xlink_ = nullptr;
		site->occupied_ = false;
		// Update statistics
		n_single_bound_--;
		if(xlink->tethered_ == true){ 
			n_tethered_--;
			Kinesin *motor = xlink->motor_;
			int x_dub_pre = motor->x_dist_doubled_;
			motor->UpdateExtension();
			if(motor->heads_active_ == 0){
				n_sites_i_untethered_--;
				properties_->kinesin4.n_free_tethered_--;
			}
			// Make sure an untethering event didn't occur
			else if(motor->tethered_ == true){
				if(motor->heads_active_ == 2){
					properties_->kinesin4.n_bound_untethered_++;
					properties_->kinesin4.n_bound_tethered_[x_dub_pre]--;
					properties_->kinesin4.n_bound_tethered_tot_--;
				}
				n_sites_i_tethered_[x_dub_pre]--;
			}
			xlink->tethered_ = false;
			motor->xlink_ = nullptr;
			motor->tethered_ = false;
			xlink->motor_ = nullptr;
			motor->UpdateExtension(); 		// Sets extension to 0
		}
		else{
			n_untethered_--;
			n_sites_i_untethered_--;
		}
		// Update xlink details
		xlink->site_one_ = nullptr;
		xlink->site_two_ = nullptr;
		xlink->heads_active_--;
	}
	else{
		printf("Error in RunKMC_Unbind: no bound xlinks\n");
		exit(1);
	}
}

void AssociatedProteinManagement::RunKMC_Unbind_II(int x_dist){

	UpdateDoubleBoundList();
	if(n_double_bound_[x_dist] > 0){
		// Randomly pick a double-bound xlink
		int i_entry = properties_->gsl.GetRanInt(n_double_bound_[x_dist]);
		AssociatedProtein* xlink = double_bound_list_[x_dist][i_entry];
		if(xlink->heads_active_ != 2){
			printf("what...error in xlink unbind ii");
			exit(1);
		}
		Tubulin* site;
		// FIXME Randomly choose a head to remove ... FIXME
		double ran = properties_->gsl.GetRanProb();
		if(ran < 0.5)
			site = xlink->site_one_;
		else
			site = xlink->site_two_;
		// Remove xlink from site
		site->xlink_ = nullptr;
		site->occupied_ = false;
//		printf("index is %i\n", site->index_);
		// Update xlink details
		if(ran < 0.5)
			xlink->site_one_ = nullptr;
		else
			xlink->site_two_ = nullptr;
		xlink->heads_active_--;
		if(x_dist != xlink->x_dist_){
			printf("why theuaroux (xink unbind_ii) \n");
			exit(1);
		}
		// Update statistics
		n_single_bound_++;
		n_double_bound_[x_dist]--;
		if(xlink->tethered_ == true){
			xlink->UpdateExtension();
			Kinesin *motor = xlink->motor_;
			int x_dub_pre = motor->x_dist_doubled_;
			motor->UpdateExtension();
			// xlinks tethered to free motors diffuse as untethered
			if(motor->heads_active_ == 0){
				n_sites_i_untethered_++;
				n_sites_ii_untethered_[x_dist] -= 2;
			}
			// Make sure we didn't force an untether event
			else if(motor->tethered_ == true){
				int x_dub_post = motor->x_dist_doubled_;
				if(motor->heads_active_ == 2){
					properties_->kinesin4.n_bound_tethered_[x_dub_pre]--;
					properties_->kinesin4.n_bound_tethered_[x_dub_post]++;
				}
				n_sites_i_tethered_[x_dub_post]++;
				n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
			}
			// If we did, counteract statistic skew
			else{
				n_sites_i_tethered_[x_dub_pre]++;
				n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
			}
		}
		else{
			xlink->UpdateExtension();
			n_sites_i_untethered_++;
			n_sites_ii_untethered_[x_dist] -= 2;
		}
	}
	else{
		printf("Error in RunKMC_Unbind_II:no double bound xlinks\n");
		exit(1);
	}
}
