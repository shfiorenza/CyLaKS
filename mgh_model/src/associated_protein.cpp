#include "master_header.h"
#include "associated_protein.h"

AssociatedProtein::AssociatedProtein(){
}

void AssociatedProtein::Initialize(system_parameters *parameters, 
	system_properties *properties, int ID) {

	ID_ = ID;
	parameters_ = parameters;
	properties_ = properties;
	SetParameters();
	CalculateCutoffs();
	InitiateNeighborLists();
	PopulateBindingLookupTable();
	PopulateTethBindingLookupTable();
	PopulateTethBindingIILookupTable();
	PopulateExtensionLookups();
}

void AssociatedProtein::SetParameters(){

	r_0_ = parameters_->xlinks.r_0;
	k_spring_ = parameters_->xlinks.k_spring;
}

void AssociatedProtein::CalculateCutoffs(){

	int site_size = parameters_->microtubules.site_size; 
	double kbT = parameters_->kbT; 
	double r_y = parameters_->microtubules.y_dist;
	/* First, calculate rest_dist_ in number of sites */
	int rough_rest_dist = sqrt(r_0_*r_0_ - r_y*r_y) / site_size; 
	double rest_scan[3]; 
	double scan_force[3]; 
	for(int i_scan = -1; i_scan <= 1; i_scan++){
		rest_scan[i_scan + 1] = rough_rest_dist + (i_scan * 0.5);
		double rest_scan_length = rest_scan[i_scan + 1] * site_size;
		double r_scan = sqrt(r_y*r_y + rest_scan_length*rest_scan_length); 
		scan_force[i_scan + 1] = (r_scan - r_0_) * k_spring_;
	}
	double min_force = 100; 
	for(int i_scan = -1; i_scan <=1; i_scan++){
		double force = fabs(scan_force[i_scan + 1]); 
		if(force < min_force){
			min_force = force;
			rest_dist_ = rest_scan[i_scan + 1]; 
		}	
	}
	/* Finally, calculate extension distance cutoff */
	for(int x_dist = (int) rest_dist_; x_dist < 1000; x_dist++){
		int r_x = x_dist * site_size;
		double r = sqrt(r_y*r_y + r_x*r_x);
		double dr = r - r_0_; 
		double U = (k_spring_/2)*dr*dr;
		double boltzmann_weight = exp(U/(2*kbT));
		if(boltzmann_weight > 100){
			dist_cutoff_ = x_dist;
			break;
		}
	}
}

void AssociatedProtein::InitiateNeighborLists(){

	int n_mts = parameters_->microtubules.count; 
	int teth_cutoff = properties_->kinesin4.dist_cutoff_; 
	neighbor_sites_.resize((n_mts - 1)*(2*dist_cutoff_ + 1));
	teth_neighbor_sites_.resize(n_mts*(2*teth_cutoff + 1));
	teth_neighbor_sites_ii_.resize((n_mts - 1)*(2*dist_cutoff_ + 1));
}

void AssociatedProtein::PopulateBindingLookupTable(){

	double r_y = parameters_->microtubules.y_dist;
	double kbT = parameters_->kbT;
	double site_size = parameters_->microtubules.site_size;
	binding_weight_lookup_.resize(dist_cutoff_ + 1);
	for(int i_dist = 0; i_dist <= dist_cutoff_; i_dist++){
		double r_x = (double)i_dist*site_size;
		double r = sqrt(r_y*r_y + r_x*r_x);
		double dr = r - r_0_;
		double U = (k_spring_/2)*dr*dr;
		double weight = exp(-U/(2*kbT));
		binding_weight_lookup_[i_dist] = weight;
//		printf("r: %g, dr: %g,  weight: %g\n", r, dr, weight);
	}
}

void AssociatedProtein::PopulateTethBindingLookupTable(){

	double r_y = parameters_->microtubules.y_dist / 2;
	double kbT = parameters_->kbT;
	double site_size = parameters_->microtubules.site_size;	
	int teth_cutoff = properties_->kinesin4.dist_cutoff_;
	int comp_cutoff = properties_->kinesin4.comp_cutoff_;
	teth_binding_weight_lookup_.resize(2*teth_cutoff + 1);
	double k_teth = parameters_->motors.k_spring;
	double k_slack = parameters_->motors.k_slack; 
	double r_0_teth = parameters_->motors.r_0; 
	for(int x_dist_dub = 0; x_dist_dub <= 2*teth_cutoff; x_dist_dub++){
		double r_x = (double) x_dist_dub * site_size / 2; 
		double r = sqrt(r_y*r_y + r_x*r_x);
		double dr = r - r_0_teth; 
		double weight = 0; 
		if(x_dist_dub < 2*comp_cutoff){
			weight = 0;
		}
		else if(dr < 0){
			weight = exp(-dr*dr*k_slack/(4*kbT));
		}
		else{
			weight = exp(-dr*dr*k_teth/(4*kbT));
		}
		teth_binding_weight_lookup_[x_dist_dub] = weight;
//		printf("weight for 2x = %i is %g\n", x_dist_dub, weight);
	}
}

// input: CURRENT x_dist_dub of tether and PROPOSED x_dist xlink
void AssociatedProtein::PopulateTethBindingIILookupTable(){

	int teth_dist_cutoff = properties_->kinesin4.dist_cutoff_;
	int comp_cutoff = properties_->kinesin4.comp_cutoff_;
	teth_binding_weight_ii_to_.resize(2*teth_dist_cutoff + 1);
	teth_binding_weight_ii_from_.resize(2*teth_dist_cutoff + 1);
	double kbT = parameters_->kbT; 
	double r_y = parameters_->microtubules.y_dist;
	double k_teth_spring = properties_->kinesin4.motors_[0].k_spring_;
	double k_teth_slack = properties_->kinesin4.motors_[0].k_slack_;
	double r_0_teth = properties_->kinesin4.motors_[0].r_0_;
	double r_y_teth = parameters_->microtubules.y_dist / 2;
	double rest_dist_teth = properties_->kinesin4.motors_[0].rest_dist_;
	double site_size = parameters_->microtubules.site_size;
	for(int x_dub = 0; x_dub <= 2*teth_dist_cutoff; x_dub++){
		teth_binding_weight_ii_to_[x_dub].resize(2*dist_cutoff_ + 1);
		teth_binding_weight_ii_from_[x_dub].resize(2*dist_cutoff_ + 1);
		// Calc x-distances (in nm) for tether
		double r_x_teth = x_dub * site_size / 2;
		// Calc total r values 
		double r_teth = sqrt(r_x_teth*r_x_teth + r_y_teth*r_y_teth);
		// Calc tether exts for current dist and stepping to/from rest
		double dr_teth = r_teth - r_0_teth;	
		for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
			double r_x = (double)x_dist*site_size;
			double r = sqrt(r_y*r_y + r_x*r_x);
			double dr = r - r_0_;
			double U = (k_spring_/2)*dr*dr;
			double weight = exp(-U/(2*kbT));
			// change in teth extension if 2nd xlink head were to bind
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
			int x_from_rest_dub = abs((int)(2*rest_dist_teth) - x_dub);
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
			double tot_weight_to = 0;
			int x_dub_to = 2 * r_x_teth_to / site_size;
			if(x_dub_to < 2*comp_cutoff
			|| x_dub_to > 2*teth_dist_cutoff){
				tot_weight_to = 0;
			}
			else{
				double weight_to_teth = exp(-dU_to_teth/(2*kbT));
				tot_weight_to = weight_to_teth*weight;	
			}
			double tot_weight_from = 0; 
			int x_dub_from = 2 * r_x_teth_from / site_size;
			if(x_dub_from < 2*comp_cutoff
			|| x_dub_from > 2*teth_dist_cutoff){
				tot_weight_from = 0;
			}
			else{
		   		double weight_from_teth = exp(-dU_from_teth/(2*kbT));
				tot_weight_from = weight_from_teth*weight;
			}
			teth_binding_weight_ii_to_[x_dub][x_dist] = tot_weight_to;
			teth_binding_weight_ii_from_[x_dub][x_dist] = tot_weight_from;
/*
			printf("for 2x=%i, x=%i, w_to = %g\n", x_dub, x_dist, 
				teth_binding_weight_ii_to_[x_dub][x_dist]);
			printf("for 2x=%i, x=%i, w_from = %g\n", x_dub, x_dist, 
				teth_binding_weight_ii_from_[x_dub][x_dist]);
*/
		}
	}
}

void AssociatedProtein::PopulateExtensionLookups(){

	double r_y = parameters_->microtubules.y_dist;
	double site_size = parameters_->microtubules.site_size;
	extension_lookup_.resize(dist_cutoff_ + 1);
	cosine_lookup_.resize(dist_cutoff_ + 1);
	for(int x_dist = 0; x_dist <= dist_cutoff_; x_dist++){
		double r_x = site_size * x_dist;
		double r = sqrt(r_x*r_x + r_y*r_y);
		extension_lookup_[x_dist] = r - r_0_;
		cosine_lookup_[x_dist] = r_x / r;
	}
}

AssociatedProtein::Monomer* AssociatedProtein::GetActiveHead(){

	if(heads_active_ == 1){
		if(head_one_.site_ != nullptr)
			return &head_one_;
		else if(head_two_.site_ != nullptr)
			return &head_two_;
		else{
			printf("what in get active head site\n");
			exit(1);
		}
	}
	else{
		printf("not cool bro...not single bound \n");
		exit(1);
	}
}

Tubulin* AssociatedProtein::GetActiveHeadSite(){

	if(heads_active_ == 1){
		if(head_one_.site_ != nullptr)
			return head_one_.site_;
		else if(head_two_.site_ != nullptr)
			return head_two_.site_;
		else{
			printf("what in get active head site\n");
			exit(1);
		}
	}
	else{
		printf("not cool bro...not single bound \n");
		exit(1);
	}
}

double AssociatedProtein::GetAnchorCoordinate(){

	// If single bound, use that head; assume other's diffusion avgs out
	if(heads_active_ == 1){
		Tubulin *site = GetActiveHeadSite();
		int index = site->index_;
		int mt_coord = site->mt_->coord_;
		double site_coord = (double)(mt_coord + index);
		return site_coord;
	}
	// If double bound, use avg of head's indices 
	else if(heads_active_ == 2){
		int index_one = head_one_.site_->index_;
		int mt_coord_one = head_one_.site_->mt_->coord_;
		double coord_one = (double)(mt_coord_one + index_one);
		int index_two = head_two_.site_->index_;
		int mt_coord_two = head_two_.site_->mt_->coord_;
		double coord_two = (double)(mt_coord_two + index_two);
		double avg_coord = (coord_one + coord_two)/2;
		return avg_coord;
	}
	else{
		printf("not NOT cool bro ... cant get anchor index: %i\n", 
				heads_active_);
		exit(1);
	}
}

int AssociatedProtein::GetPRC1NeighbCount(Monomer* head){

	int n_neighbs = 0;
	Tubulin *site = head->site_;
	int i_plus = site->mt_->plus_end_;    
	int i_minus = site->mt_->minus_end_;    
	int dx = site->mt_->delta_x_; 
	if(site->index_ == i_plus){    
		if(site->mt_->lattice_[site->index_-dx].xlink_head_ != nullptr)    
			n_neighbs++;
	}   
	else if(site->index_ == i_minus){    
		if(site->mt_->lattice_[site->index_+dx].xlink_head_ != nullptr)    
			n_neighbs++;
	}   
	else{    
		if(site->mt_->lattice_[site->index_-dx].xlink_head_ != nullptr)    
			n_neighbs++;    
		if(site->mt_->lattice_[site->index_+dx].xlink_head_ != nullptr)    
			n_neighbs++;    
	}   
	return n_neighbs;
}

void AssociatedProtein::UpdateNeighborSites(){
	
	n_neighbor_sites_ = 0;
	int n_mts = parameters_->microtubules.count;
	if(n_mts > 1){
		Tubulin *site = GetActiveHeadSite();
		int i_site = site->index_;
		Microtubule *mt = site->mt_;
		int mt_length = mt->neighbor_->n_sites_; 
		Microtubule *adj_mt = mt->neighbor_;
		int site_coord = mt->coord_ + i_site; 
		// Scan through all potential neighbors; only add unoccupied to list 
		int i_entry = 0;
		for(int dist = -dist_cutoff_; dist <= dist_cutoff_; dist++){
			int i_neighbor = (site_coord - adj_mt->coord_) + dist;
			// Start index at first site (0) if i_neighb is 0 or neg
			if(i_neighbor < 0){
				dist -= (i_neighbor + 1); // - 1;
			}
			// End scan once last site (mt_length - 1) has been checked
			else if(i_neighbor > mt_length - 1){
				break;
			} 
			else{
				Tubulin *neighbor = &adj_mt->lattice_[i_neighbor];
				if(neighbor->occupied_ == false){
					n_neighbor_sites_++;
					neighbor_sites_[i_entry] = neighbor;
					i_entry++;
				}
			}
		}
	}
}

void AssociatedProtein::UpdateTethNeighborSites(){

	n_teth_neighbor_sites_ = 0;
	if(tethered_ == true
	&& heads_active_ == 0){
		int n_mts = parameters_->microtubules.count;
		int teth_cutoff = properties_->kinesin4.dist_cutoff_; 
		int comp_cutoff = properties_->kinesin4.comp_cutoff_;
		double stalk_coord = motor_->GetStalkCoordinate();
//		printf("stalk coord is %g\n\n", anchor_coord);
		// Scan through all potential neighbor sites; add unoccupied to list 
		int i_entry = 0;
		for(int i_mt = 0; i_mt < n_mts; i_mt++){ 	
			Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
			int mt_length = mt->n_sites_;
			double mt_coord = mt->coord_;
			int i_stalk = stalk_coord - mt_coord; 
			for(int x_dist = -teth_cutoff; x_dist <= teth_cutoff; x_dist++){
				int i_site = i_stalk + x_dist; 
//				printf("i_site is %i (x_dist %i)\n", i_site, x_dist);
				// Start index at first site (0) if site index is <= 0
				if(i_site < 0){
					x_dist -= (i_site + 1);
				}
				// End scan at last site (mt_length - 1)
				else if(i_site > mt_length - 1){
					break;
				}
				else{
					Tubulin *neighbor = &mt->lattice_[i_site];
					double site_coord = i_site + neighbor->mt_->coord_; 
					double x_dist = fabs(stalk_coord - site_coord); 
					int x_dist_dub = 2*x_dist;
					if(x_dist_dub >= 2*comp_cutoff
					&& x_dist_dub <= 2*teth_cutoff
					&& neighbor->occupied_ == false){
						n_teth_neighbor_sites_++;
						teth_neighbor_sites_[i_entry] = neighbor;
						i_entry++;
					}
				}
			}   
		}
	}
	else{
		printf("error in XLINK update teth neighbor sites\n");
		exit(1);
	}
}

void AssociatedProtein::UpdateTethNeighborSitesII(){

	n_teth_neighbor_sites_ii_ = 0;
	int n_mts = parameters_->microtubules.count;
	if(n_mts > 1){
		Tubulin *site = GetActiveHeadSite();
		int i_site = site->index_;
		Microtubule *mt = site->mt_;
		int mt_length = mt->neighbor_->n_sites_;
		Microtubule *adj_mt = mt->neighbor_;
		int site_coord = mt->coord_ + i_site; 
		// Scan through all potential neighbors; only add unoccupied to list 
		int i_entry = 0;
		for(int dist = -dist_cutoff_; dist <= dist_cutoff_; dist++){
			int i_neighbor = (site_coord - adj_mt->coord_) + dist;
			// Start index at first site (0) if i_neighb is 0 or neg
			if(i_neighbor < 0){
				dist -= (i_neighbor + 1);
			}
			// End scan once last bulk site (mt_length - 1) has been checked
			else if(i_neighbor > mt_length - 1){
				break;
			} 
			else{
				Tubulin *neighbor = &adj_mt->lattice_[i_neighbor];
				if(neighbor->occupied_ == false){
					n_teth_neighbor_sites_ii_++;
					teth_neighbor_sites_ii_[i_entry] = neighbor;
					i_entry++;
				}
			}
		}
	}
}

void AssociatedProtein::UpdateExtension(){

	if(heads_active_ == 2){
		int x_dist_pre = x_dist_;
		// Calculate first head's coordinate
		int i_head_one = head_one_.site_->index_;
		int mt_coord_one = head_one_.site_->mt_->coord_;
		int coord_one = mt_coord_one + i_head_one;
		// Calculate second head's coordinate
		int i_head_two = head_two_.site_->index_;
		int mt_coord_two = head_two_.site_->mt_->coord_;
		int coord_two = mt_coord_two + i_head_two;
		// Calculate x_distance in # of sites
		int x_dist = abs(coord_one - coord_two);	
		x_dist_ = x_dist; 
		if(x_dist <= dist_cutoff_){
			extension_ = extension_lookup_[x_dist_];
			cosine_ = cosine_lookup_[x_dist_];
		}
		else{
			ForceUnbind(x_dist_pre);
//			printf("forced an unbind event >:O\n");
		}
	}
	else if(heads_active_ == 1){
		x_dist_ = 0;
		extension_ = 0;
		cosine_ = 0;
	}
	else{
		printf("some kinda error in assoc. protein update_extension\n");
		exit(1);
	}
}

void AssociatedProtein::ForceUnbind(int x_dist_pre){

	// different stuff depending on whether or not xlink is tethered
	if(tethered_ == false){
		double ran = properties_->gsl.GetRanProb(); 
		if(ran < 0.5){
			head_one_.site_->xlink_head_ = nullptr;
			head_one_.site_->occupied_ = false;
			head_one_.site_ = nullptr;
		}
		else{
			head_two_.site_->xlink_head_ = nullptr;
			head_two_.site_->occupied_ = false;
			head_two_.site_ = nullptr;	
		}
		heads_active_--;
		x_dist_ = 0;
		extension_ = 0; 
		cosine_ = 0;
	}
	else if(motor_->heads_active_ > 0){
		// Update motor ext (unbinding 2nd head changes anchor coord) 
		Tubulin *site_to = GetSiteFartherFromTethRest(); 
		int neighbs_to = site_to->GetPRC1NeighborCount();
		Tubulin *site_fr = GetSiteCloserToTethRest();
		int neighbs_fr = site_fr->GetPRC1NeighborCount();

		int x_dub_pre = motor_->x_dist_doubled_;
		double p_unbind_to = properties_->prc1.p_unbind_ii_to_teth_
			[neighbs_to][x_dub_pre][x_dist_pre]; 
		double p_unbind_from = properties_->prc1.p_unbind_ii_fr_teth_
			[neighbs_fr][x_dub_pre][x_dist_pre]; 
		double p_tot = p_unbind_to + p_unbind_from;


		double ran = properties_->gsl.GetRanProb();
		if(ran < p_unbind_to / p_tot){
			site_to->xlink_head_ = nullptr;
			site_to->occupied_ = false;
			site_to = nullptr;
		}
		else{
			site_fr->xlink_head_ = nullptr;
			site_fr->occupied_ = false;
			site_fr = nullptr;
		}
		heads_active_--;
		x_dist_ = 0;
		extension_ = 0; 
		cosine_ = 0;
		motor_->UpdateExtension();
	}
	// xlinks tethered to free motors diffuse as if untethered
	else{
		double ran = properties_->gsl.GetRanProb(); 
		if(ran < 0.5){
			head_one_.site_->xlink_head_ = nullptr;
			head_one_.site_->occupied_ = false;
			head_one_.site_ = nullptr;
		}
		else{
			head_two_.site_->xlink_head_ = nullptr;
			head_two_.site_->occupied_ = false;
			head_two_.site_ = nullptr;	
		}
		heads_active_--;
		x_dist_ = 0;
		extension_ = 0; 
		cosine_ = 0;
	}
}

void AssociatedProtein::UntetherSatellite(){

	if(tethered_ == true){
		if(heads_active_ == 0){
			printf("Error in Untether_Satellite() in ASSOC. PROT.\n");
			exit(1);
		}
		// Remove satellite motor from active_ list, replace with last entry
		int i_last = properties_->kinesin4.n_active_ - 1;
		Kinesin *last_entry = properties_->kinesin4.active_[i_last];
		int i_this = motor_->active_index_;
		if(i_this != i_last){
			properties_->kinesin4.active_[i_this] = last_entry;
			last_entry->active_index_ = i_this;
		}
		properties_->kinesin4.n_active_--;
		// Update motor details
		motor_->tethered_ = false;
		motor_->xlink_ = nullptr;
		// Update xlink details
		tethered_ = false;
		motor_ = nullptr;
	}
	else{
		printf("Error in xlink UntethSatellite()\n");
		exit(1);
	}
}

int AssociatedProtein::GetDirectionTowardRest(Tubulin *site){

	if(heads_active_ == 2){
		double anchor_coord = GetAnchorCoordinate();
		int site_coord = site->index_ + site->mt_->coord_;
		if(site_coord == anchor_coord){
			double ran = properties_->gsl.GetRanProb();
			if(ran < 0.5)
				return -1;
			else
				return 1;
		}
		else if(site_coord > anchor_coord)
			return -1;
		else
			return 1;
	}
	else{
		printf("error in get dir. toward rest (xlink)\n");
		exit(1);
	}
}

double AssociatedProtein::GetBindingWeight(Tubulin *neighbor){

	Tubulin *site = GetActiveHeadSite();
	Microtubule *mt = site->mt_;
	Microtubule *adj_mt = neighbor->mt_;
	if(adj_mt != mt->neighbor_){
		printf("adj: %i, neighb: %i\n", 
				adj_mt->index_, mt->neighbor_->index_);
		printf("why the microtubules tho (in assiociated protein GBW)\n");
		exit(1);
	}
	int offset = adj_mt->coord_ - mt->coord_;
	// Calculate distance (in x-dir.) between site and neighbor in # of sites 
	int i_site = site->index_;
	int i_neighbor = neighbor->index_;
	int x_dist = abs(i_neighbor - i_site + offset);
	// Get binding weight that corresponds to this x-distance
	double weight = binding_weight_lookup_[x_dist];
	return weight;	
}

double AssociatedProtein::GetTethBindingWeight(Tubulin *neighbor){

	double stalk_coord = motor_->GetStalkCoordinate(); 
	double site_coord = neighbor->index_ + neighbor->mt_->coord_; 
	double x_dist = fabs(stalk_coord - site_coord);
	int x_dist_dub = 2*x_dist;
	double weight = teth_binding_weight_lookup_[x_dist_dub];
	return weight; 
}

double AssociatedProtein::GetTethBindingWeightII(Tubulin *neighbor){

	// Get x_dist of xlink extension if 2nd head were to bind to this neighb
	Tubulin *site = GetActiveHeadSite();
	double site_coord = site->index_ + site->mt_->coord_;
	double neighb_coord = neighbor->index_ + neighbor->mt_->coord_;
	int x_dist = fabs(site_coord - neighb_coord);
	// Get current x_dist_dub of tether
	double stalk_coord = motor_->GetStalkCoordinate();
	int x_dub = 2*fabs(stalk_coord - site_coord);
	// Get NEW x_dist_dub if 2nd head were to bind to this neighb
	double anchor_coord_post = (neighb_coord + site_coord) / 2;
	int x_dub_post = 2*fabs(stalk_coord - anchor_coord_post);
	int x_dub_rest = 2*properties_->kinesin4.rest_dist_;
	double weight;
	// If extended, ...
	if(x_dub >= x_dub_rest){
		// ... further extension is going away FROM rest
		if(x_dub_post > x_dub)
			weight = teth_binding_weight_ii_from_[x_dub][x_dist];
		else
			weight = teth_binding_weight_ii_to_[x_dub][x_dist];
	}
	// Else, if compressed ...
	else{
		// ... further extension is going TO rest
		if(x_dub_post > x_dub)
			weight = teth_binding_weight_ii_to_[x_dub][x_dist];
		else
			weight = teth_binding_weight_ii_from_[x_dub][x_dist];
	}
	return weight;
}

double AssociatedProtein::GetExtensionForce(Tubulin *site){

	if(heads_active_ == 2){
		UpdateExtension();
		// Make sure we didn't force an unbinding event
		if(heads_active_ == 2){
			double force_mag = extension_ * k_spring_;			// in pN
			double site_coord = site->index_ + site->mt_->coord_;
			double anchor_coord = GetAnchorCoordinate();
			double force;
			if(site_coord < anchor_coord)
				force = force_mag * cosine_;
			else
				force = -1 * force_mag * cosine_;
			return force; 
		}
		else
			return 0;
	}
	else{
		printf("error in get ext force (xlink)\n");
		exit(1);
	}
}

Tubulin* AssociatedProtein::GetSiteCloserToTethRest(){

	if(heads_active_ == 2
	&& tethered_ == true){
		double stalk_coord = motor_->GetStalkCoordinate();
		double site_one_coord 
			= head_one_.site_->index_ + head_one_.site_->mt_->coord_;
		int dx_dub_one = 2 * fabs(stalk_coord - site_one_coord);
		double site_two_coord 
			= head_two_.site_->index_ + head_two_.site_->mt_->coord_;
		int dx_dub_two = 2 * fabs(stalk_coord - site_two_coord);
		double anchor_coord = GetAnchorCoordinate();
		int x_dub = 2 * fabs(anchor_coord - stalk_coord);
		int rest_dist_dub = 2 * properties_->kinesin4.rest_dist_; 
		// twice dx of teth extension if 1 head unbinds
		int dx_teth_dub = x_dist_; 
		// twice the distance teth is drom rest
		int dist_from_rest_dub = abs(x_dub - rest_dist_dub);
		// if unbinding causes teth to go from slack to spring or
		// vise versa, always choose the site that makes it go slack
		if(dx_teth_dub > dist_from_rest_dub){
			if(dx_dub_one < dx_dub_two)
				return head_one_.site_; 
			else
				return head_two_.site_;
		}
		else if(x_dub >= rest_dist_dub){
			if(dx_dub_one < dx_dub_two)
				return head_one_.site_; 
			else
				return head_two_.site_;
		}
		else{
			if(dx_dub_one > dx_dub_two)
				return head_one_.site_;
			else
				return head_two_.site_;
		}
	}
	else{
		printf("Error in get_site_closer_to_teth_rest XLINK\n");
		exit(1);
	}
}

Tubulin* AssociatedProtein::GetSiteFartherFromTethRest(){

	if(heads_active_ == 2
	&& tethered_ == true){
		double stalk_coord = motor_->GetStalkCoordinate();
		double site_one_coord 
			= head_one_.site_->index_ + head_one_.site_->mt_->coord_;
		int dx_dub_one = 2 * fabs(stalk_coord - site_one_coord);
		double site_two_coord 
			= head_two_.site_->index_ + head_two_.site_->mt_->coord_;
		int dx_dub_two = 2 * fabs(stalk_coord - site_two_coord);
		double anchor_coord = GetAnchorCoordinate();
		int x_dub = 2 * fabs(anchor_coord - stalk_coord);
		int rest_dist_dub = 2 * properties_->kinesin4.rest_dist_; 
		// twice dx of teth extension if 1 head unbinds
		int dx_teth_dub = x_dist_; 
		// twice the distance teth is drom rest
		int dist_from_rest_dub = abs(x_dub - rest_dist_dub);
		// if unbinding causes teth to go from slack to spring or
		// vise versa, always choose the site that makes it taut
		if(dx_teth_dub > dist_from_rest_dub){
			if(dx_dub_one > dx_dub_two)
				return head_one_.site_; 
			else
				return head_two_.site_;
		}
		else if(x_dub >= rest_dist_dub){
			if(dx_dub_one > dx_dub_two)
				return head_one_.site_; 
			else
				return head_two_.site_;
		}
		else{
			if(dx_dub_one < dx_dub_two)
				return head_one_.site_;
			else
				return head_two_.site_;
		}
	}
	else{
		printf("Error in get_site_farther_fr_teth_rest XLINK\n");
		exit(1);
	}
}

Tubulin* AssociatedProtein::GetWeightedNeighborSite(){
	
	UpdateNeighborSites();
	double anch_coord = GetAnchorCoordinate(); 
	double p_tot = 0;
	for(int i_site = 0; i_site < n_neighbor_sites_; i_site++){
		Tubulin *site = neighbor_sites_[i_site]; 
		double site_coord = site->index_ + site->mt_->coord_; 
		int x_dist = (int)fabs(anch_coord - site_coord);
//		printf("x_dist is %i\n", x_dist);
		p_tot += binding_weight_lookup_[x_dist]; 
	}
	double ran = properties_->gsl.GetRanProb();
	double p_cum = 0;
	for(int i_site = 0; i_site < n_neighbor_sites_; i_site++){
		Tubulin *site = neighbor_sites_[i_site];
		double site_coord = site->index_ + site->mt_->coord_;
		int x_dist = (int)fabs(anch_coord - site_coord); 
		p_cum += binding_weight_lookup_[x_dist] / p_tot;
		if(ran < p_cum){
			return site;
		}
	}
	return nullptr;
}

Tubulin* AssociatedProtein::GetWeightedTethNeighborSite(){
	
	UpdateTethNeighborSites();
	double stalk_coord = motor_->GetStalkCoordinate();
	double p_tot = 0; 
	for(int i_site = 0; i_site < n_teth_neighbor_sites_; i_site++){
		Tubulin *site = teth_neighbor_sites_[i_site];
		double site_coord = site->index_ + site->mt_->coord_; 
		double x_dist = fabs(stalk_coord - site_coord);
		int x_dist_dub = 2*x_dist;
		p_tot += teth_binding_weight_lookup_[x_dist_dub];
	}
	double ran = properties_->gsl.GetRanProb();
	double p_cum = 0;
	for(int i_site = 0; i_site < n_teth_neighbor_sites_; i_site++){
		Tubulin *site = teth_neighbor_sites_[i_site];
		double site_coord = site->index_ + site->mt_->coord_;
		double x_dist = fabs(stalk_coord - site_coord);
		int x_dist_dub = 2*x_dist;
		p_cum += teth_binding_weight_lookup_[x_dist_dub] / p_tot; 
		if(ran < p_cum){
			return site;
		}
	}
	return nullptr;
}

Tubulin* AssociatedProtein::GetWeightedTethNeighborSiteII(){

	UpdateTethNeighborSitesII();
	double stalk_coord = motor_->GetStalkCoordinate();
	double anchor_coord = GetAnchorCoordinate();
	// All neighbor sites have the same initial x_dub
	int x_dub = 2*fabs(stalk_coord - anchor_coord);
	int x_dub_rest = 2*properties_->kinesin4.rest_dist_; 
	double weight_tot = 0;
	for(int i_site = 0; i_site < n_teth_neighbor_sites_ii_; i_site++){
		Tubulin *site = teth_neighbor_sites_ii_[i_site];
		double site_coord = site->index_ + site->mt_->coord_; 
		int x_dist = (int)fabs(site_coord - anchor_coord);
		double anchor_coord_post = (site_coord + anchor_coord) / 2;
		int x_dub_post = 2*fabs(stalk_coord - anchor_coord_post);
		double weight;
		// If extended, ...
		if(x_dub >= x_dub_rest){
			// ... further extension is going away FROM rest
			if(x_dub_post > x_dub)
				weight = teth_binding_weight_ii_from_[x_dub][x_dist];
			else
				weight = teth_binding_weight_ii_to_[x_dub][x_dist];
		}
		// Else, if compressed ...
		else{
			// ... further extension is going TO rest
			if(x_dub_post > x_dub)
				weight = teth_binding_weight_ii_to_[x_dub][x_dist];
			else
				weight = teth_binding_weight_ii_from_[x_dub][x_dist];
		}
		weight_tot += weight;
	}
	double ran = properties_->gsl.GetRanProb();
	double p_cum = 0;
	for(int i_site = 0; i_site < n_teth_neighbor_sites_ii_; i_site++){
		Tubulin *site = teth_neighbor_sites_ii_[i_site];
		double site_coord = site->index_ + site->mt_->coord_; 
		int x_dist = (int)fabs(site_coord - anchor_coord);
		double anchor_coord_post = (site_coord + anchor_coord) / 2;
		int x_dub_post = 2*fabs(stalk_coord - anchor_coord_post);
		double weight;
		// If extended, ...
		if(x_dub >= x_dub_rest){
			// ... further extension is going away FROM rest
			if(x_dub_post > x_dub)
				weight = teth_binding_weight_ii_from_[x_dub][x_dist];
			else
				weight = teth_binding_weight_ii_to_[x_dub][x_dist];
		}
		// Else, if compressed ...
		else{
			// ... further extension is going TO rest
			if(x_dub_post > x_dub)
				weight = teth_binding_weight_ii_to_[x_dub][x_dist];
			else
				weight = teth_binding_weight_ii_from_[x_dub][x_dist];
		}
		p_cum += weight / weight_tot;
		if(ran < p_cum)
			return site;
	}
	return nullptr;
}
