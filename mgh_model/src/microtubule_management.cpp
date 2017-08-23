#include "master_header.h"
#include "microtubule_management.h"

MicrotubuleManagement::MicrotubuleManagement(){
}

void MicrotubuleManagement::Initialize(system_parameters *parameters, 
								       system_properties *properties){

	parameters_ = parameters;
	properties_ = properties;

	GenerateActiveMTs();
	GenerateUnoccupiedList();
}

void MicrotubuleManagement::GenerateActiveMTs(){

	int n_mts = parameters_->n_microtubules;
	active_mts.resize(n_mts);
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		active_mts[i_mt].Initialize(parameters_, properties_, i_mt);
	}
}

void MicrotubuleManagement::GenerateUnoccupiedList(){

	int n_mts = parameters_->n_microtubules;
	int n_sites = parameters_->length_of_microtubule;
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		Microtubule *mt = &active_mts[i_mt];
		for(int i_site = 1; i_site < (n_sites - 1); i_site++){
			Tubulin *site = &mt->lattice[i_site];
			UnoccupiedCheck(i_mt, i_site);
			UnoccupiedCheck(site);
			unoccupied_list.push_back(site);
			n_unoccupied_++;
		}
	}
}

void MicrotubuleManagement::UnoccupiedCheck(Tubulin *site){

	if(site->occupant != nullptr){
		printf("Error @ site %i_%i: should be unoccupied\n", site->parent->index_, 
			site->index_);
		exit(1);
	}

}

void MicrotubuleManagement::UnoccupiedCheck(int i_mt, int i_site){

	if(active_mts[i_mt].lattice[i_site].occupant != nullptr){
		printf("Error @ site %i_%i: should be unoccupied\n", i_mt, i_site);
		exit(1);
	}
}

void MicrotubuleManagement::OccupiedCheck(Tubulin *site){
	
	if(site->occupant == nullptr){
		printf("Error @ site %i_%i: should be unoccupied\n", site->parent->index_, 
			site->index_);
		exit(1);
	}	
}

Tubulin* MicrotubuleManagement::GetUnoccupiedSite(){

	if(unoccupied_list.empty() != true){
		int i_entry = properties_->gsl.GetRanInt(unoccupied_list.size());	
		Tubulin *site = unoccupied_list[i_entry];
		UnoccupiedCheck(site->parent->index_, site->index_);
		return site;
	}
	else{
		printf("Error: GetUnoccupiedSite() called, but unoccupied_list is empty\n");
		exit(1);
	}
}
	
void MicrotubuleManagement::RemoveFromUnoccupiedList(Tubulin *site){

	OccupiedCheck(site);
	// Generate iterator pointing to this site's entry in unoccupied_list
	auto site_entry = std::find(unoccupied_list.begin(), unoccupied_list.end(), site);
	// Get numerical index (i.e. distance from vector's start) of the iterator
	int i_entry = std::distance(unoccupied_list.begin(), site_entry);
	// Erase entry in unoccupied_list that corresponds to this site 
	unoccupied_list.erase(unoccupied_list.begin() + i_entry);
	n_unoccupied_--;
}

void MicrotubuleManagement::AddToUnoccupiedList(Tubulin *site){
	
	UnoccupiedCheck(site);
	unoccupied_list.emplace_back(site);
	n_unoccupied_++;
}
