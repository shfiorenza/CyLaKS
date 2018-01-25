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
	n_xlinks_ = n_mts*(n_sites - 2);
	xlink_list_.resize(n_xlinks_);
	for(int ID = 0; ID < n_xlinks_; ID++){
		xlink_list_[ID].Initialize(parameters_, properties_, ID);
	}
}

void AssociatedProteinManagement::SetParameters(){

	double c_xlink = 50;		// compare to 200 for motors; in nM
	double k_on = parameters_->k_on; 
	double k_off = parameters_->k_off;
	double delta_t = parameters_->delta_t;
	p_bind_i_ = k_on*c_xlink*delta_t;
	p_unbind_i_ = k_off*delta_t*2;
	// Generate unbinding rates based on discretized spring extension
	dist_cutoff_ = (int) xlink_list_[0].dist_cutoff_; 	//FIXME
	p_unbind_ii_.resize(dist_cutoff_ + 1);
	for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
		double site_size = xlink_list_[0].site_size_;
		double r_0 = xlink_list_[0].r_0_;
		double kbT = xlink_list_[0].kbT_;
		double k_spring = xlink_list_[0].k_spring_;
		double r_x = x_dist*site_size;
		double r_y = 35;									// in nm
		double r = sqrt(r_y*r_y + r_x*r_x);
		double dr = r - r_0;
		double k_off_ii = exp(dr*dr*k_spring/(2*kbT));
//		printf("%i:  %g\n", distance, k_off_ii);
		p_unbind_ii_[x_dist] = k_off_ii*delta_t;
	}
}

void AssociatedProteinManagement::InitiateLists(){

	single_bound_list_.resize(n_xlinks_);
	// Have a double-bound list for each spring extension up to cutoff
	int cutoff = (int) xlink_list_[0].dist_cutoff_;		//FIXME
	n_double_bound_.resize(cutoff + 1);
	double_bound_list_.resize(cutoff + 1);
	for(int i_entry = 0; i_entry <= cutoff; i_entry++){
		double_bound_list_[i_entry].resize(n_xlinks_);
		n_double_bound_[i_entry] = 0;
	}
	untethered_list_.resize(n_xlinks_);
	// The following lists contain the sites each xlink is bound to
	single_tethered_sites_.resize(n_xlinks_);
	double_tethered_sites_.resize(n_xlinks_);
	single_untethered_sites_.resize(n_xlinks_);
	double_untethered_sites_.resize(n_xlinks_);

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
	int n_unbound = 0;
	int n_single_bound = 0;
	int n_double_bound = 0;
	for(int i_xlink = 0; i_xlink < n_xlinks_; i_xlink++){
		AssociatedProtein *xlink = &xlink_list_[i_xlink];
		switch(xlink->heads_active_){
			case 0: n_unbound++;
					break;
			case 1: single_bound_list_[i_entry] = xlink;
					n_single_bound++;
					i_entry++;
					break;
			case 2: n_double_bound++;
					break;
		}
	}
	if(n_unbound + n_single_bound + n_double_bound != n_xlinks_)
		printf("idk bro 45\n");
	if(n_single_bound != n_single_bound_)
		printf("not good in single bound list %i_%i\n", 
				n_single_bound, n_single_bound_);
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
	
	int n_untethered = 0;
	int i_entry = 0;
	for(int i_xlink = 0; i_xlink < n_xlinks_; i_xlink++){
		AssociatedProtein *xlink = &xlink_list_[i_xlink];
		if(xlink->heads_active_ > 0){
			if(xlink->tethered_ == false){
				n_untethered++;
				untethered_list_[i_entry] = xlink;
				i_entry++;
			}
		}
	}
	if(n_untethered != n_untethered_){
			printf("ugh, updateuntetheredlist:");
			printf(" %i in stats, %i in list\n", n_untethered_, n_untethered);
	}
}

void AssociatedProteinManagement::UpdateSingleTetheredSites(){

    int n_single_tethered = 0;
    int i_site = 0;
    for(int i_xlink = 0; i_xlink < n_xlinks_; i_xlink++){
        AssociatedProtein *xlink = &xlink_list_[i_xlink];
        if(xlink->heads_active_ == 1
		&& xlink->tethered_ == true){
			n_single_tethered++;
			Tubulin *site = xlink->GetActiveHeadSite();
			single_tethered_sites_[i_site] = site;
			i_site++;
		}
    }
	if(n_single_tethered != n_sites_single_tethered_)
			printf("bad mojo single bound tethered sites\n");
}

void AssociatedProteinManagement::UpdateDoubleTetheredSites(){

	int n_double_tethered = 0;
	int i_site = 0;
	for(int i_xlink = 0; i_xlink < n_xlinks_; i_xlink++){
		AssociatedProtein *xlink = &xlink_list_[i_xlink];
		if(xlink->heads_active_ == 2
		&& xlink->tethered_ == true){
			n_double_tethered+=2;
			Tubulin *site_one = xlink->site_one_;
			double_tethered_sites_[i_site] = site_one;
			i_site++;
			Tubulin *site_two = xlink->site_two_;
			double_tethered_sites_[i_site] = site_two;
			i_site++;
		}
	}
	if(n_double_tethered != n_sites_double_tethered_)
		printf("bad mojo double bound tethered sites\n");
}

void AssociatedProteinManagement::UpdateSingleUntetheredSites(){

	int n_single_untethered = 0;
	int i_site = 0;
	for(int i_xlink = 0; i_xlink < n_xlinks_; i_xlink++){
		AssociatedProtein *xlink = &xlink_list_[i_xlink];
		if(xlink->heads_active_ == 1
		&& xlink->tethered_ == false){
			n_single_untethered++;
			Tubulin *site = xlink->GetActiveHeadSite();
			single_untethered_sites_[i_site] = site;
			i_site++;
		}
	}
	if(n_single_untethered != n_sites_single_untethered_)
		printf("bad mojo single bound untethered sites\n");
}

void AssociatedProteinManagement::UpdateDoubleUntetheredSites(){

	int n_double_untethered = 0;
	int i_site = 0;
	for(int i_xlink = 0; i_xlink < n_xlinks_; i_xlink++){
		AssociatedProtein *xlink = &xlink_list_[i_xlink];
		if(xlink->heads_active_ == 2
		&& xlink->tethered_ == false){
			n_double_untethered+=2;
			Tubulin *site_one = xlink->site_one_;
			double_untethered_sites_[i_site] = site_one;
			i_site++;
			Tubulin *site_two = xlink->site_two_;
			double_untethered_sites_[i_site] = site_two;
			i_site++;
		}
	}
	if(n_double_untethered != n_sites_double_untethered_)
		printf("bad mojo double bound untethered sites\n");
}

void AssociatedProteinManagement::UpdateExtensions(){
	
	for(int i_xlink = 0; i_xlink < n_xlinks_; i_xlink++){
		AssociatedProtein *xlink = &xlink_list_[i_xlink];
		if(xlink->heads_active_ == 2){
			int x_dist_pre = xlink->x_dist_;
			xlink->UpdateExtension();
			int x_dist_post = xlink->x_dist_;
			// If an unbinding event is forced, update statistics
			if(xlink->heads_active_ == 1){
				n_double_bound_[x_dist_pre]--; 
				n_single_bound_++;
				if(xlink->tethered_ == true){
					n_sites_double_tethered_ -= 2; 	
					n_sites_single_tethered_++;
				}
				else{
					n_sites_double_untethered_ -= 2;
					n_sites_single_untethered_++; 
				}
			}
			// If x_distance changed, update those statistics
			else if(x_dist_pre != x_dist_post){
				n_double_bound_[x_dist_pre]--;
				n_double_bound_[x_dist_post]++;

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

void AssociatedProteinManagement::RunDiffusion(){

	int i_step = properties_->current_step_;
	// FIXME the following is only valid for delta_t = 0.0005 FIXME
	int	i_tau_st = 4,  //(int)(tau_single_tethered_/parameters_->delta_t),
	    i_tau_dt = 58, //(int)(tau_double_tethered_/parameters_->delta_t),
		i_tau_su = 1,  //(int)(tau_single_untethered_/parameters_->delta_t),
		i_tau_du = 10; //(int)(tau_double_untethered_/parameters_->delta_t);
	if(i_step % i_tau_st == 0)
		RunDiffusion_Single_Tethered();
	if(i_step % i_tau_dt == 0)
		RunDiffusion_Double_Tethered();
	if(i_step % i_tau_su == 0)
		RunDiffusion_Single_Untethered();
	if(i_step % i_tau_du == 0)
		RunDiffusion_Double_Untethered();
	UpdateExtensions();
}

void AssociatedProteinManagement::RunDiffusion_Single_Tethered(){

	int mt_length = parameters_->length_of_microtubule;
	UpdateSingleTetheredSites();
	int n_sites = n_sites_single_tethered_;
	int step_dir[n_sites];
	if(n_sites > 1){
		int n_half = (int) (n_sites/2);
		// Half of all xlink heads step backwards
		for(int i_entry = 0; i_entry < n_half; i_entry++)
			step_dir[i_entry] = -1;
		// Other half step forward
		for(int i_entry = n_half; i_entry < n_sites; i_entry++)
			step_dir[i_entry] = 1;
		// Shuffle the order of forward/backwards steps
		gsl_ran_shuffle(properties_->gsl.rng, step_dir, n_sites, sizeof(int));
	}
	else if(n_sites > 0){
		double ran = properties_->gsl.GetRanProb();
		if(ran < 0.5)
			step_dir[0] = 1;
		else
			step_dir[0] = -1;
	}
	// Run through list of occupied sites and step xlink heads appriopriately
	for(int i_entry = 0; i_entry < n_sites; i_entry++){
		Tubulin *site = single_tethered_sites_[i_entry];
		AssociatedProtein *xlink = site->xlink_;
		int dx = step_dir[i_entry];
		int i_site = site->index_;
		// Forbid diffusion onto boundary sites
		if(!(i_site == 1 && dx == -1)
		&& !(i_site == mt_length - 2 && dx == 1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &site->mt_->lattice_[i_site + dx];
			// Only diffuse if not blocked
			if(new_site->occupied_ == false){
				old_site->xlink_ = nullptr;
				old_site->occupied_ = false;
				if(old_site == xlink->site_one_)
					xlink->site_one_ = new_site;
				else
					xlink->site_two_ = new_site;
				new_site->xlink_ = xlink;
				new_site->occupied_ = true;

				// Update motor statistics
				int x_dub_pre = xlink->motor_->x_dist_doubled_;
				xlink->motor_->UpdateExtension();
				int x_dub_post = xlink->motor_->x_dist_doubled_;
				properties_->kinesin4.n_bound_tethered_[x_dub_pre]--;
				properties_->kinesin4.n_bound_tethered_[x_dub_post]++;
			}
		}
	}
}

void AssociatedProteinManagement::RunDiffusion_Double_Tethered(){

	int mt_length = parameters_->length_of_microtubule;
	UpdateDoubleTetheredSites();
	int n_sites = n_sites_double_tethered_;
	int step_dir[n_sites];
	if(n_sites > 0){
		int n_half = (int) (n_sites/2);
		// Half of all xlink heads step backwards
		for(int i_entry = 0; i_entry < n_half; i_entry++)
			step_dir[i_entry] = -1;
		// Other half step forward
		for(int i_entry = n_half; i_entry < n_sites; i_entry++)
			step_dir[i_entry] = 1;
		// Shuffle the order of forward/backwards steps
		gsl_ran_shuffle(properties_->gsl.rng, step_dir, n_sites, sizeof(int));
	}	
	// Run through list of occupied sites and step xlink heads appriopriately
	for(int i_entry = 0; i_entry < n_sites; i_entry++){
		Tubulin *site = double_tethered_sites_[i_entry];
		AssociatedProtein *xlink = site->xlink_;
		int dx = step_dir[i_entry];
		int i_site = site->index_;
		// Forbid diffusion onto boundary sites
		if(!(i_site == 1 && dx == -1)
		&& !(i_site == mt_length - 2 && dx == 1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &site->mt_->lattice_[i_site + dx];
			// Only diffuse if not blocked
			if(new_site->occupied_ == false){
				old_site->xlink_ = nullptr;
				old_site->occupied_ = false;
				if(old_site == xlink->site_one_)
					xlink->site_one_ = new_site;
				else
					xlink->site_two_ = new_site;
				new_site->xlink_ = xlink;
				new_site->occupied_ = true;
				
				// Update motor statistics
				int x_dub_pre = xlink->motor_->x_dist_doubled_;
				xlink->motor_->UpdateExtension();
				int x_dub_post = xlink->motor_->x_dist_doubled_;
				properties_->kinesin4.n_bound_tethered_[x_dub_pre]--;
				properties_->kinesin4.n_bound_tethered_[x_dub_post]++;
			}
		}
	}
}

void AssociatedProteinManagement::RunDiffusion_Single_Untethered(){

	int mt_length = parameters_->length_of_microtubule;
	UpdateSingleUntetheredSites();
	int n_sites = n_sites_single_untethered_;
	int step_dir[n_sites];
	if(n_sites > 1){
		int n_half = (int) (n_sites/2);
		// Half of all xlink heads step backwards
		for(int i_entry = 0; i_entry < n_half; i_entry++)
			step_dir[i_entry] = -1;
		// Other half step forward
		for(int i_entry = n_half; i_entry < n_sites; i_entry++)
			step_dir[i_entry] = 1;
		// Shuffle the order of forward/backwards steps
		gsl_ran_shuffle(properties_->gsl.rng, step_dir, n_sites, sizeof(int));
	}	
	else if(n_sites > 0){
		double ran = properties_->gsl.GetRanProb();
		if(ran < 0.5)
			step_dir[0] = 1;
		else
			step_dir[0] = -1;
	}
	// Run through list of occupied sites and step xlink heads appriopriately
	for(int i_entry = 0; i_entry < n_sites; i_entry++){
		Tubulin *site = single_untethered_sites_[i_entry];
		AssociatedProtein *xlink = site->xlink_;
		int dx = step_dir[i_entry];
		int i_site = site->index_;
		// Forbid diffusion onto boundary sites
		if(!(i_site == 1 && dx == -1)
		&& !(i_site == mt_length - 2 && dx == 1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &site->mt_->lattice_[i_site + dx];
			// Only diffuse if not blocked
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
}

void AssociatedProteinManagement::RunDiffusion_Double_Untethered(){

	int mt_length = parameters_->length_of_microtubule;
	UpdateDoubleUntetheredSites();
	int n_sites = n_sites_double_untethered_;
	int step_dir[n_sites];
	if(n_sites > 0){
		int n_half = (int) (n_sites/2);
		// Half of all xlink heads step backwards
		for(int i_entry = 0; i_entry < n_half; i_entry++)
			step_dir[i_entry] = -1;
		// Other half step forward
		for(int i_entry = n_half; i_entry < n_sites; i_entry++)
			step_dir[i_entry] = 1;
		// Shuffle the order of forward/backwards steps
		gsl_ran_shuffle(properties_->gsl.rng, step_dir, n_sites, sizeof(int));
	}	
	// Run through list of occupied sites and step xlink heads appriopriately
	for(int i_entry = 0; i_entry < n_sites; i_entry++){
		Tubulin *site = double_untethered_sites_[i_entry];
		AssociatedProtein *xlink = site->xlink_;
		int dx = step_dir[i_entry];
		int i_site = site->index_;
		// Forbid diffusion onto boundary sites
		if(!(i_site == 1 && dx == -1)
		&& !(i_site == mt_length - 2 && dx == 1)){
			Tubulin *old_site = site;
			Tubulin *new_site = &site->mt_->lattice_[i_site + dx];
			// Only diffuse if not blocked
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
}

void AssociatedProteinManagement::GenerateKMCList(){

	UpdateExtensions();
	int n_events = 0;
	int n_to_bind_i = GetNumToBind_I();
	n_events += n_to_bind_i; 
	int n_to_bind_ii = GetNumToBind_II();
	n_events += n_to_bind_ii;
	int n_to_unbind_i = GetNumToUnbind_I();
	n_events += n_to_unbind_i; 
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
	int n_unoccupied = properties_->microtubules.n_unoccupied_;
	double n_avg = p_bind_i_*n_unoccupied;
	int n_to_bind = properties_->gsl.SamplePoissonDist(n_avg);
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
	double n_equil = c_eff_*weights_summed; 
	double n_avg = p_bind_i_*n_equil;
	int n_to_bind = properties_->gsl.SamplePoissonDist(n_avg);
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
	
	int x_dist = 0;		// extension of xlink for stage2 unbinding
	GenerateKMCList();
	if(kmc_list_.empty() == false){
		int n_events = kmc_list_.size();
		for(int i_event = 0; i_event < n_events; i_event++){
			int kmc_event = kmc_list_[i_event];
			// Allows us to encode xlink extension as second digit
			if(kmc_event >= 30 
			&& kmc_event < 40){
				x_dist = kmc_event % 10;
				kmc_event = 3;
			}
			switch(kmc_event){
				case 0: 
						printf("xlink bound stage 1\n");
						RunKMC_Bind_I();
						break;
				case 1: 
						printf("xlink bound stage 2\n");
						RunKMC_Bind_II();
						break;
				case 2: 
						printf("xlink fully unbound\n");
						RunKMC_Unbind_I();
						break;
				case 3: 
						printf("xlink unbound 2nd head\n");
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
			n_sites_single_tethered_++;
		}
		else{
			n_untethered_++;
			n_sites_single_untethered_++;
		}
	}
	else{
		printf("Error in RunKMC_BindFirst: no unoccupied sites\n");
		exit(1);
	}
}

void AssociatedProteinManagement::RunKMC_Bind_II(){

	// Make sure stage 1 xlinks and unoccupied sites are available
	UpdateSingleBoundList();
	properties_->microtubules.UpdateUnoccupiedList();
	if(n_single_bound_ > 0
	&& properties_->microtubules.n_unoccupied_ > 0){
		// Sample normal distribution for x-dist of xlink to insert 
		double sigma = ((double)dist_cutoff_)/3; 
		int x_dist = properties_->gsl.SampleNormalDist(sigma);
		if(x_dist > dist_cutoff_)
			x_dist = dist_cutoff_;
		// Randomly pick single-bound xlink
		int i_xlink = properties_->gsl.GetRanInt(n_single_bound_);
		AssociatedProtein *xlink = single_bound_list_[i_xlink];
		xlink->UpdateNeighborSites();
		// Ensure xlink has a neighbor with desired distance
		// roll for a new distance after a certain # of tries
		int attempts = 0;
		while(xlink->NeighborExists(x_dist) == false){
			if(attempts > n_single_bound_){
				attempts = 0;
				x_dist = properties_->gsl.SampleNormalDist(sigma);
				if(x_dist > dist_cutoff_)
					x_dist = dist_cutoff_;
			}
			i_xlink = properties_->gsl.GetRanInt(n_single_bound_);
			xlink = single_bound_list_[i_xlink];	
			xlink->UpdateNeighborSites();
			attempts++;
		}
		Tubulin *site = xlink->GetNeighborSite(x_dist);
		if(site == nullptr){
			printf("issue with xlink bind 2\n");
			exit(1);
		}
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
		if(xlink->x_dist_ != x_dist){
			printf("why thoooooouaz: ");
			printf(" %i in xlink, %i in kmc\n", xlink->x_dist_, x_dist);
			exit(1);
		}
		// Update statistics
		n_double_bound_[x_dist]++;
		n_single_bound_--;
		if(xlink->tethered_ == true){
			n_sites_double_tethered_ += 2;
			n_sites_single_tethered_--;
		}
		else{
			n_sites_double_untethered_ += 2;	
			n_sites_single_untethered_--;
		}
		// Check to see if motor tether extension was affected
		if(xlink->tethered_ == true){
			Kinesin *motor = xlink->motor_;
			if(motor->heads_active_ != 1){
				int x_dub_pre = motor->x_dist_doubled_;
				motor->UpdateExtension();
				if(motor->tethered_ == true){
					int x_dub_post = motor->x_dist_doubled_;
					if(x_dub_pre != x_dub_post){
						properties_->kinesin4.n_bound_tethered_[x_dub_pre]--;
						properties_->kinesin4.n_bound_tethered_[x_dub_post]++;
					}
				}
			}
		}
	}
	else{
		printf("Error in RunKMC_BindSecond: no unoccupied sites\n");
		exit(1);
	}
}

void AssociatedProteinManagement::RunKMC_Unbind_I(){
	
	UpdateSingleBoundList();
	if(n_single_bound_ > 0){
		// Randomly pick a single-bound xlink
		int i_entry = properties_->gsl.GetRanInt(n_single_bound_);
		AssociatedProtein *xlink = single_bound_list_[i_entry];
		Tubulin *site = xlink->GetActiveHeadSite();
		Kinesin *motor;
		int x_dist_dub = 0;
		if(xlink->tethered_ == true){
			motor = xlink->motor_;
			if(xlink->motor_->heads_active_ != 1){
				motor->UpdateExtension();
				x_dist_dub = motor->x_dist_doubled_;
			}
		}
		// Remove xlink from site
		site->xlink_ = nullptr;
		site->occupied_ = false;
		// Update xlink details
		xlink->site_one_ = nullptr;
		xlink->site_two_ = nullptr;
		xlink->heads_active_--;
		// Update statistics
		n_single_bound_--;
		if(xlink->tethered_ == true){ 
			n_tethered_--;
			n_sites_single_tethered_--;
			if(motor->heads_active_ == 0){
				properties_->kinesin4.n_free_tethered_--;
			}
			else if(motor->heads_active_ == 2){
				properties_->kinesin4.n_bound_untethered_++;
				properties_->kinesin4.n_bound_tethered_[x_dist_dub]--;
				properties_->kinesin4.n_bound_tethered_tot_--;
			}
			xlink->tethered_ = false;
			motor->xlink_ = nullptr;
			motor->tethered_ = false;
			xlink->motor_ = nullptr;
		}
		else{
			n_untethered_--;
			n_sites_single_untethered_--;
		}
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
		xlink->UpdateExtension();
		// Update statistics
		n_single_bound_++;
		n_double_bound_[x_dist]--;
		if(xlink->tethered_ == true){
			n_sites_double_tethered_ -= 2;
			n_sites_single_tethered_++; 
		}
		else{
			n_sites_double_untethered_ -= 2;
			n_sites_single_untethered_++;
		}
		// Check to see if motor tether extension was affected
		if(xlink->tethered_ == true){
			Kinesin *motor = xlink->motor_;
			if(motor->heads_active_ != 1){
				int x_dub_pre = motor->x_dist_doubled_;
				motor->UpdateExtension();
				if(motor->tethered_ == true){
					int x_dub_post = motor->x_dist_doubled_;
					if(x_dub_pre != x_dub_post){
						properties_->kinesin4.n_bound_tethered_[x_dub_pre]--;
						properties_->kinesin4.n_bound_tethered_[x_dub_post]++;
					}
				}
			}
		}
	}
	else{
		printf("Error in RunKMC_Unbind_II:no double bound xlinks\n");
		exit(1);
	}
}
