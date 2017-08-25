#include "master_header.h"
#include "microtubule_management.h"

MicrotubuleManagement::MicrotubuleManagement(){

}

void MicrotubuleManagement::Initialize(system_parameters *parameters, 
									   system_properties *properties){

	parameters_ = parameters;
	properties_ = properties;

	SetParameters();
	GenerateMicrotubules();
	UpdateUnoccupiedList();
	UpdateUnoccupiedPairsList();
}

void MicrotubuleManagement::SetParameters(){

	int n_mts = parameters_->n_microtubules;
	int mt_length = parameters_->length_of_microtubule;
	n_sites_tot_ = n_mts*mt_length;
	int n_sites_bulk = n_sites_tot_ - 2*n_mts;
	unoccupied_list_.resize(n_sites_bulk);
	unoccupied_pairs_list_.resize(n_sites_bulk);
}

void MicrotubuleManagement::GenerateMicrotubules(){

	int n_mts = parameters_->n_microtubules;
	mt_list_.resize(n_mts);
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		mt_list_[i_mt].Initialize(parameters_, properties_, i_mt);
	}
}

void MicrotubuleManagement::UnoccupiedCheck(Tubulin *site){

	if(site->motor_ != nullptr || site->xlink_ != nullptr){
		printf("Error @ site %i_%i: should be unoccupied\n", site->mt_->index_, 
															 site->index_);
		exit(1);
	}
}

void MicrotubuleManagement::UnoccupiedCheck(int i_mt, int i_site){

	if(mt_list_[i_mt].lattice_[i_site].motor_ != nullptr
	|| mt_list_[i_mt].lattice_[i_site].xlink_ != nullptr){
		printf("Error @ site %i_%i: should be unoccupied\n", i_mt, i_site);
		exit(1);
	}
}

void MicrotubuleManagement::OccupiedCheck(Tubulin *site){
	
	if(site->motor_ == nullptr && site->xlink_ == nullptr){
		printf("Error @ site %i_%i: should be occupied\n", site->mt_->index_, 
														   site->index_);
		exit(1);
	}	
}

void MicrotubuleManagement::OccupiedCheck(int i_mt, int i_site){

	if(mt_list_[i_mt].lattice_[i_site].motor_ == nullptr
	&& mt_list_[i_mt].lattice_[i_site].xlink_ == nullptr){
		printf("Error @ site %i_%i: should be occupied\n", i_mt, i_site);
		exit(1);
	}
}

void MicrotubuleManagement::UpdateNumUnoccupied(){
	
	n_unoccupied_ = 0;
	int n_occupied = 0;
	int n_mts = parameters_->n_microtubules;
	int mt_length = parameters_->length_of_microtubule;
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		Microtubule *mt = &mt_list_[i_mt];
		// Exclude boundary sites
		for(int i_site = 1; i_site < mt_length - 1; i_site++){
			Tubulin *site = &mt->lattice_[i_site];
			if(site->occupied_ == false){
				n_unoccupied_++;

			}
			else{
				n_occupied++;
			}
		}
	}	
	// Verify that statistics didn't go all wonky 
	int n_sites_bulk = n_sites_tot_ - 2*n_mts;
	if(n_unoccupied_ + n_occupied != n_sites_bulk){
		printf("something wrong in update_n_unoccupied\n");
		exit(1);
	}
}

void MicrotubuleManagement::UpdateUnoccupiedList(){

	UpdateNumUnoccupied();
	int i_unoccupied = 0;
	int n_mts = parameters_->n_microtubules;
	int mt_length = parameters_->length_of_microtubule;
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		Microtubule *mt = &mt_list_[i_mt];
		// Exclude boundary sites
		for(int i_site = 1; i_site < mt_length - 1; i_site++){
			Tubulin *site = &mt->lattice_[i_site];
			if(site->occupied_ == false){
				unoccupied_list_[i_unoccupied] = site;
				i_unoccupied++;
			}
		}
	}
	if(i_unoccupied != n_unoccupied_){
		printf("something awful in update_unoccupied_list bruh:\n");
		printf("  %i != %i\n", i_unoccupied, n_unoccupied_);
		exit(1);
	}
}

void MicrotubuleManagement::UpdateNumUnoccupiedPairs(){

	n_unoccupied_pairs_ = 0;
	int n_unoccupied_singles = 0;
	int n_occupied_singles = 0;
	int n_mts = parameters_->n_microtubules;
	int mt_length = parameters_->length_of_microtubule;
	bool pair_flag = false;
//	printf("boing 1\n");
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
//		printf("boing 2\n");
		Microtubule *mt = &mt_list_[i_mt];
		// Exclude boundary sites as either member of a pair
		for (int i_site = 1; i_site < mt_length - 1; i_site++){
			Tubulin *site = &mt->lattice_[i_site];
			// If site is unoccupied, check for an adjacent unoccupied site to pair it with
			if(site->occupied_ == false){
				// If pair_flag is active, pair this site with the site behind it
				if(pair_flag == true){
/*					if(properties_->current_step_ > 30000){
						printf("CHURCH @ %i (%i)\n", i_site, properties_->current_step_);
						properties_->wallace.PrintMicrotubules();
					}
*/					n_unoccupied_pairs_++;
					pair_flag = false;
				}
				// If this site is the last one on the MT, do not attempt to pair it
				else if(i_site == mt_length - 2){
/*					if(properties_->current_step_ > 30000){
						printf("the end @ %i (%i)\n", i_site, properties_->current_step_);
						properties_->wallace.PrintMicrotubules();
					}
*/					n_unoccupied_singles++;
				}
				// Otherwise, activate pair_flag for the next site in the FOR loop
				else{
/*					if(properties_->current_step_ > 30000){
						printf("chargin' @ %i (%i)\n", i_site, properties_->current_step_);
						properties_->wallace.PrintMicrotubules();
					}
*/					pair_flag = true;	
				}
			}
			// If this site is occupied, simply update pair_flag and related statistics
			else if(site->occupied_ == true){
				n_occupied_singles++;
				if(pair_flag == true){
					n_unoccupied_singles++;
					pair_flag = false;
				}
			}
		}
	}
	// Verify that statistics didn't go all wonky
	int n_singles = n_occupied_singles + n_unoccupied_singles;
	int n_sites_bulk = n_sites_tot_ - 2*n_mts;
	if(n_singles + 2*n_unoccupied_pairs_ != n_sites_bulk){
		printf("something terrible in update_num_unoccupied_pairs:\n");
		printf("  %i + 2*%i != %i\n", n_singles, n_unoccupied_pairs_, n_sites_bulk);
		exit(1);
	}
}

void MicrotubuleManagement::UpdateUnoccupiedPairsList(){

	UpdateNumUnoccupiedPairs();
	n_pair_entries_ = 0;
	int i_pair = 0;
	int n_mts = parameters_->n_microtubules;
	int mt_length = parameters_->length_of_microtubule;
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		Microtubule *mt = &mt_list_[i_mt];
		// Exclude boundary sites as either member of a pair
		for(int i_site = 1; i_site < mt_length - 2; i_site++){
			Tubulin *first_site = &mt->lattice_[i_site];
			if(first_site->occupied_ == false){
				Tubulin *second_site = &mt->lattice_[i_site + 1];
				if(second_site->occupied_ == false){
					// Only store a pointer to the pair's rear-most site
					if(i_mt%2 == 0){
						unoccupied_pairs_list_[i_pair] = first_site;
						i_pair++;
						n_pair_entries_++;
					}
					else{
						unoccupied_pairs_list_[i_pair] = second_site;
						i_pair++;
						n_pair_entries_++;
					}
				}
			}
		}
	}
	if(n_pair_entries_ > 2*n_unoccupied_pairs_){
		printf("something awful in update_unoccupied_pairs_list bruh:\n");
		printf("  %i > 2*%i\n", i_pair, n_unoccupied_pairs_);
		exit(1);
	}
}

Tubulin* MicrotubuleManagement::GetUnoccupiedSite(){

	UpdateUnoccupiedList();
	// Make sure an unoccupied site exists
	if(n_unoccupied_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_unoccupied_);	
		Tubulin *site = unoccupied_list_[i_entry];
		UnoccupiedCheck(site);
		return site;
	}
	else{
		printf("Error: GetUnoccupiedSite() called, but no unoccupied sites\n");
		exit(1);
	}
}

Tubulin* MicrotubuleManagement::GetUnoccupiedPair_1(){

	// Make sure at least one pair of unoccupied sites exists
	if(n_unoccupied_pairs_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_pair_entries_);
		Tubulin *site = unoccupied_pairs_list_[i_entry];
		UnoccupiedCheck(site);
		return site;
	}
	else{
		printf("Error: GetUnoccupiedPair_1 called, but no unoccupied pairs\n");
		exit(1);
	}
}

Tubulin* MicrotubuleManagement::GetUnoccupiedPair_2(Tubulin* first_site){

	if(first_site->occupied_ == false){
		int i_mt = first_site->mt_->index_;
		int i_site = first_site->index_;
		int delta_x = first_site->mt_->delta_x_;
		Tubulin *site = &mt_list_[i_mt].lattice_[i_site + delta_x];
		UnoccupiedCheck(site);
		if(site->index_ != first_site->index_ + delta_x){
			printf("error in GUP_2\n");
			exit(1);
		}
		return site;
	}
	else{
		printf("Error: GetUnoccupiedPair_2 called, but first site occupied\n");
		exit(1);
	}

}
