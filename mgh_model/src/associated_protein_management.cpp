#include "master_header.h"
#include "associated_protein_management.h"

AssociatedProteinManagement::AssociatedProteinManagement(){
}

void AssociatedProteinManagement::Initialize(system_parameters *parameters, 
		system_properties *properties){

	parameters_ = parameters;
	properties_ = properties;

	SetParameters();
	GenerateXLinks();
}

void AssociatedProteinManagement::SetParameters(){

	// FIXME xlinks need their own stats bruh FIXME
	double k_on = parameters_->k_on/2;
	double c_motor = parameters_->c_motor;
	double k_off = parameters_->k_off*2;
	double delta_t = parameters_->delta_t;
	p_bind_first_ = k_on*c_motor*delta_t;
	p_unbind_ = k_off*delta_t;

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
	tethered_list_.resize(n_xlinks_);
	untethered_list_.resize(n_xlinks_);	
	occupied_site_list_.resize(n_xlinks_);
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

void AssociatedProteinManagement::UpdateTetherLists(){

	int n_tethered = 0;
	int i_entry_tethered = 0;
	int i_entry_untethered = 0;
	int n_untethered = 0;
	int n_unbound = 0;
	for(int i_xlink = 0; i_xlink < n_xlinks_; i_xlink++){
		AssociatedProtein *xlink = &xlink_list_[i_xlink];
		if(xlink->heads_active_ > 0){
			if(xlink->tethered_ == true){
				n_tethered++;
				tethered_list_[i_entry_tethered] = xlink;
				i_entry_tethered++;
			}
			else{
				n_untethered++;
				untethered_list_[i_entry_untethered] = xlink;
				i_entry_untethered++;
			}
		}
		else{
			n_unbound++;
		}
	}
	if(n_tethered != n_tethered_ 
	|| n_untethered != n_untethered_)
		printf("ugh (%i, %i | %i, %i)\n", n_tethered, n_tethered_, 
			n_untethered, n_untethered_);
	if(n_tethered + n_untethered != n_single_bound_ + n_double_bound_)
		printf("why (%i, %i | %i, %i)\n", n_tethered_, n_untethered, 
				n_single_bound_, n_double_bound_);
}


void AssociatedProteinManagement::UpdateOccupiedSiteList(){

	int n_single_bound = 0, 
		n_double_bound = 0, 
		n_sites_occupied = 0;
	int i_site = 0;
	int n_unactive = 0;
	for(int i_xlink = 0; i_xlink < n_xlinks_; i_xlink++){
		AssociatedProtein *xlink = &xlink_list_[i_xlink];
		switch(xlink->heads_active_){
			case 0: n_unactive+=2;
					break;
			case 1: n_unactive++;
					if(xlink->site_one_ != nullptr)
						occupied_site_list_[i_site] = xlink->site_one_;
					else if(xlink->site_two_ != nullptr)
						occupied_site_list_[i_site] = xlink->site_two_;
					i_site++;
					n_single_bound++;
					n_sites_occupied++;
					break;
			case 2: occupied_site_list_[i_site] = xlink->site_one_;
					i_site++;
					occupied_site_list_[i_site] = xlink->site_two_;
					i_site++;
					n_double_bound++;
					n_sites_occupied+=2;
					break;
		}
	}	
	if(i_site != n_sites_occupied_)
		printf("something awful in update_bound_list for xlinks bruh\n");
	if(n_unactive != 2*n_xlinks_ - n_sites_occupied_)
		printf("something awful THE SEQUEL in update_bound_list for xlinks\n");
	if(n_single_bound != n_single_bound_ 
	|| n_double_bound != n_double_bound_
	|| n_sites_occupied != n_sites_occupied)
		printf("something awful THA THIRDD  (%i|%i)(%i|%i)(%i|%i)\n", 
			n_single_bound, n_single_bound_, n_double_bound, n_double_bound_, 
			n_sites_occupied, n_sites_occupied_);
}

AssociatedProtein* AssociatedProteinManagement::GetUntetheredXlink(){

	if(n_untethered_ > 0){
		UpdateTetherLists();
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

	UpdateOccupiedSiteList();
	int mt_length = parameters_->length_of_microtubule;	
	int step_dir [n_sites_occupied_];
	if(n_sites_occupied_ > 0){
		int n_half = (int) (n_sites_occupied_/2);
		// Half of all xlink heads step backwards
		for(int i_entry = 0; i_entry < n_half; i_entry ++)
			step_dir[i_entry] = -1;
		// The other half step forwards
		for(int i_entry = n_half; i_entry < n_sites_occupied_; i_entry++)
			step_dir[i_entry] = 1;
		// Shuffle the order to randomize which heads step forward or backwards 
		gsl_ran_shuffle(properties_->gsl.rng, step_dir, n_sites_occupied_, sizeof(int));
	}
	// Run through list of occupied sites and step the respective xlink heads appropriately
	for(int i_entry = 0; i_entry < n_sites_occupied_; i_entry++){
		Tubulin *site = occupied_site_list_[i_entry];
		AssociatedProtein *xlink = site->xlink_;
		int dx = step_dir[i_entry] * site->mt_->delta_x_;
		int i_site = site->index_;
		// Do not let xlinks diffuse onto boundary sites
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
				else if(old_site == xlink->site_two_)
					xlink->site_two_ = new_site;
				new_site->xlink_ = xlink;
				new_site->occupied_ = true;
			}
		}
	}
}

void AssociatedProteinManagement::GenerateKMCList(){

	int n_to_bind_first = GetNumToBindFirst();
	int n_to_bind_second = 0;
	int n_to_unbind = GetNumToUnbind();
	int n_events = n_to_bind_first + n_to_bind_second + n_to_unbind;
	if(n_events > 0){
		int pre_list[n_events];
		int i_event = 0;
		for(i_event; i_event < n_to_bind_first; i_event++){
			pre_list[i_event] = 0;
		}
		for(i_event; i_event < n_to_bind_first + n_to_bind_second; i_event++){
			pre_list[i_event] = 1;
		}
		for(i_event; i_event < n_events; i_event++){
			pre_list[i_event] = 2;
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

int AssociatedProteinManagement::GetNumToBindFirst(){
	
	properties_->microtubules.UpdateNumUnoccupied();
	int n_unoccupied = properties_->microtubules.n_unoccupied_;
	double n_avg = p_bind_first_*n_unoccupied;
	int n_to_bind = properties_->gsl.SamplePoissonDist(n_avg);
	return n_to_bind;
}

int AssociatedProteinManagement::GetNumToUnbind(){

	int n_to_unbind = properties_->gsl.SampleBinomialDist(p_unbind_, n_sites_occupied_);
	return n_to_unbind;
}

void AssociatedProteinManagement::RunKMC(){
	
	GenerateKMCList();
	if(kmc_list_.empty() == false){
		int n_events = kmc_list_.size();
		for(int i_event; i_event < n_events; i_event++){
			int kmc_event = kmc_list_[i_event];
			switch(kmc_event){
				case 0: RunKMC_BindFirst();
						break;
				case 1: RunKMC_BindSecond();
						break;
				case 2: RunKMC_Unbind();
						break;
			}
		}
	}
}

void AssociatedProteinManagement::RunKMC_BindFirst(){

	if(n_sites_occupied_ != n_xlinks_){
		int i_xlink = properties_->gsl.GetRanInt(n_xlinks_);
		AssociatedProtein *xlink = &xlink_list_[i_xlink];
		while(xlink->heads_active_ != 0){
			i_xlink++;
			if(i_xlink == n_xlinks_)
				i_xlink = 0;
			xlink = &xlink_list_[i_xlink];
		}
		properties_->microtubules.UpdateUnoccupiedList();
		if(properties_->microtubules.n_unoccupied_ > 0){
			Tubulin *site = properties_->microtubules.GetUnoccupiedSite();
			Microtubule *mt = site->mt_;
			// Place xlink onto site
			site->xlink_ = xlink;
			site->occupied_ = true;
			// Update xlink details
			xlink->heads_active_++;
			xlink->site_one_ = site; 
			xlink->mt_one_ = mt;
			// Update statistics
			n_sites_occupied_++;
			n_single_bound_++;
		}
		else{
			printf("Failed to bind xlink\n");
		}
	}
	else{
		printf("Error in RunKMC_BindFirst: no unoccupied sites\n");
		exit(1);
	}
}

void AssociatedProteinManagement::RunKMC_BindSecond(){

}

void AssociatedProteinManagement::RunKMC_Unbind(){
	
	if(n_sites_occupied_ > 0){
		UpdateOccupiedSiteList();
		// Randomly pick an xlink
		int i_entry = properties_->gsl.GetRanInt(n_sites_occupied_);
		Tubulin* site = occupied_site_list_[i_entry];
		AssociatedProtein* xlink = site->xlink_;
		// Rem_eve(ve xlink from site
		site->xlink_ = nullptr;
		site->occupied_ = false;
		// Update xlink details
		if(site == xlink->site_one_){
			xlink->site_one_ = nullptr;
			xlink->mt_one_ = nullptr;
		}
		else if(site == xlink->site_two_){
			xlink->site_two_ = nullptr;
			xlink->mt_two_ = nullptr;
		}
		xlink->heads_active_--;
		// Update statistics
		if(xlink->heads_active_ == 1){
			n_double_bound_--;
			n_single_bound_++;
		}
		else if(xlink->heads_active_ == 0){
			n_single_bound_--;
		}
		n_sites_occupied_--;
	}
	else{
		printf("Error in RunKMC_Unbind: no bound xlinks\n");
		exit(1);
	}
}
