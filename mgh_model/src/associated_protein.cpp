#include "master_header.h"
#include "associated_protein.h"

AssociatedProtein::AssociatedProtein(){
}

void AssociatedProtein::Initialize(system_parameters *parameters, 
	system_properties *properties, int ID) {

	ID_ = ID;
	parameters_ = parameters;
	properties_ = properties;
	neighbor_list_binding_.resize(50);
	neighbor_list_tethering_.resize(200);
}

void AssociatedProtein::UpdateNeighborList_Tethering(){
	
	int mt_length = parameters_->length_of_microtubule;
	int i_anchor;
	if(heads_active == 1){
		if(site_one_ != nullptr){
			i_anchor = site_one_->index_;
		}
		else if(site_two_ != nullptr){
			i_anchor = site_two_->index_;
		}
		else	
			printf("associated proteins god damnit \n");
	}
	else if(heads_active == 2){
		int i_first = site_one_->index_;
		int i_second = site_one_->index_;
		i_anchor = (i_first  + i_second)/2;
	}	
	int teth_cutoff = 50;
	for(int i_mt = 0; i_mt < n_mts; i_mt++){ //FIXME for 2+ MTs
		// Scan thru 2*teth_cutoff neighboring sites, centered around i_anchor
		for(int dist = -teth_cutoff; dist < teth_cutoff; dist++){
			int i_site = i_anchor + dist;
			// Start scan at first bulk site if i_site is negative or zero
			if(i_site < 1){
				i_site = 1;
				dist = -i_anchor + 1;
			}
			// End scan at last bulk site 
			else if(i_site > mt_length - 2){
				break;
			}
			Microtubule *mt = &properties_->microtubules.mt_list[i_mt];
			Tubulin *site = &mt->lattice_[i_site];
		}
	}
}
