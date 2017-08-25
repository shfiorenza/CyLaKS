#include "master_header.h"
#include "microtubule.h"

Microtubule::Microtubule(){
}

void Microtubule::Initialize(system_parameters *parameters, 
	system_properties *properties, int i_mt){

	parameters_ = parameters;
	properties_ = properties;
	index_ = i_mt;

	SetParameters();
	GenerateLattice();
}

void Microtubule::SetParameters(){

	n_sites_ = parameters_->length_of_microtubule;
	coord_ = 0;
	if(index_%2 == 0){
		polarity_ = 0;
		plus_end_ = n_sites_ - 1;
		minus_end_ = 0;
		delta_x_ = 1;
		mt_index_adj_ = index_ + 1; 	// FIXME
	}
	else if(index_%2 == 1){
		polarity_ = 1;
		plus_end_ = 0;
		minus_end_ = n_sites_ - 1;
		delta_x_ = -1;
		mt_index_adj_ = index_ - 1;		// FIXME 
	}
}

void Microtubule::GenerateLattice(){

	lattice_.resize(n_sites_);
	for(int i_site = 0; i_site < n_sites_; i_site++){
		lattice_[i_site].Initialize(parameters_, properties_, this, i_site);
	}
}
