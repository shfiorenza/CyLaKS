#include "master_header.h"
#include "associated_protein.h"

AssociatedProtein::AssociatedProtein(){
}

void AssociatedProtein::Initialize(system_parameters *parameters, 
	system_properties *properties, int ID) {

	ID_ = ID;
	parameters_ = parameters;
	properties_ = properties;
	neighbor_sites_.resize(2*dist_cutoff_ + 1);		//fix cutoff
	binding_lookup_table_.resize(dist_cutoff_ + 1);	//fix cutoff
	PopulateBindingLookupTable();
}

void AssociatedProtein::PopulateBindingLookupTable(){

	double r_y = 35e-9;		//dist between MTs in meters
	for(int i_dist = 0; i_dist <= dist_cutoff_; i_dist++){
		double r_x = (double)i_dist*site_size_;
		double r_squared = r_y*r_y + r_x*r_x;
		double weight = exp(-r_squared*k_spring_/(2*kbT_));
//		printf("r_squared: %g, weight: %g\n", r_squared, weight);
		binding_lookup_table_[i_dist] = weight;
	}
	
}

void AssociatedProtein::UpdateNeighborSites(){
	
	int mt_length = parameters_->length_of_microtubule;
	n_neighbor_sites_ = 0;
	Tubulin *site = GetActiveHeadSite();
	int i_site = site->index_;
	Microtubule *mt = site->mt_;
	Microtubule *adj_mt = mt->neighbor_;
	int offset = adj_mt->coord_ - mt->coord_;
	int cutoff = 12;	//FIXME
	int i_entry = 0;
	// Scan through all potential neighbor sites; only add unoccupied to list 
	for(int distance = -cutoff; distance <= cutoff; distance++){
		int i_neighbor = i_site + distance - offset;
		// Start index at first bulk site (1) if i_neighbor is 0 or  negative
		if(i_neighbor <= 0){
			distance -= i_neighbor;
			i_neighbor = 1;
		}
		// End scan once last bulk site (mt_length - 2) has been checked
		else if(i_neighbor >= mt_length - 1){
			break;
		} 
		else{
			Tubulin *neighbor = &adj_mt->lattice_[i_neighbor];
			if(neighbor->occupied_ == false){
				n_neighbor_sites_++;
				neighbor_sites_[i_entry] = neighbor;
				i_entry++;
//				printf("site: %i, neighb: %i   ", i_site, i_neighbor);
//				printf("boop %i\n", i_entry);
			}
		}
	}
}

void AssociatedProtein::UpdateNeighborMotors(){
	
	int mt_length = parameters_->length_of_microtubule;
	int n_mts = parameters_->n_microtubules;
	int i_anchor;
	if(heads_active_ == 1){
		if(site_one_ != nullptr){
			i_anchor = site_one_->index_;
		}
		else if(site_two_ != nullptr){
			i_anchor = site_two_->index_;
		}
		else	
			printf("associated proteins god damnit \n");
	}
	else if(heads_active_ == 2){
		int i_first = site_one_->index_;
		int i_second = site_one_->index_;
		i_anchor = (i_first  + i_second)/2;
	}	
	// need to calc cutoff
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
			Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
			Tubulin *site = &mt->lattice_[i_site];
		}
	}
}

Tubulin* AssociatedProtein::GetActiveHeadSite(){

	if(heads_active_ == 1){
		if(site_one_ != nullptr)
			return site_one_;
		else if(site_two_ != nullptr)
			return site_two_;
		else
			printf("what uVX\n");
	}
	else{
		printf("not cool bro...not single bound \n");
		exit(1);
	}
}

Tubulin* AssociatedProtein::GetNeighborSite(){

	// Sample my janky-ass normal dist w/ sigma as argument
	int distance = properties_->gsl.SampleNormalDist(dist_cutoff_/2);
	if(distance > dist_cutoff_)
		distance = dist_cutoff_;
	// Attempt to find neighbor w/ that distance
	for(int i_neighb = 0; i_neighb < n_neighbor_sites_; i_neighb++){
		Tubulin *site = GetActiveHeadSite();
		Tubulin *neighbor = neighbor_sites_[i_neighb];
		int offset = neighbor->mt_->coord_ - site->mt_->coord_;
		int i_bound = site->index_;
		int i_target = neighbor->index_;
		int neighb_dist = abs(i_target - i_bound + offset);
//		printf("%i_%i - (%i|%i|%i)\n", distance, neighb_dist, 
//				i_bound, i_target, offset);
		if(neighb_dist == distance)
			return neighbor;
	}
//	printf("%i\n", distance);
	return nullptr;
}

double AssociatedProtein::GetBindingWeight(Tubulin *neighbor){

	Tubulin *site = GetActiveHeadSite();
	Microtubule *mt = site->mt_;
	Microtubule *adj_mt = neighbor->mt_;
	if(adj_mt != mt->neighbor_)
		printf("why the microtubules tho (in assiociated protein GBW)\n");
	int offset = adj_mt->coord_ - mt->coord_;
	int i_bound = site->index_;
	int i_target = neighbor->index_;
	int dist = abs(i_target - i_bound + offset);
	// get weight from lookup table
	double weight = binding_lookup_table_[dist];
	return weight;	
}
