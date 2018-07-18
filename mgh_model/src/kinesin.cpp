#include "master_header.h"
#include "kinesin.h"

Kinesin::Kinesin(){
}

void Kinesin::Initialize(system_parameters *parameters, 
	system_properties *properties, int ID){

	ID_ = ID; 
	parameters_ = parameters;
	properties_ = properties;
	SetParameters();
	CalculateCutoffs();
	InitiateNeighborLists(); 
	PopulateTetheringLookupTable();
	PopulateBindingLookupTable();
}

void Kinesin::SetParameters(){

	r_0_ = parameters_->motors.r_0;
	k_spring_ = parameters_->motors.k_spring;
	k_slack_ = parameters_->motors.k_slack;
}

void Kinesin::CalculateCutoffs(){

	int site_size = parameters_->microtubules.site_size; 
	double r_y = parameters_->microtubules.y_dist / 2;
	double kbT = parameters_->kbT; 
	/* First, calculate rest_dist_ in number of sites */
	int rough_rest_dist = sqrt(r_0_*r_0_ - r_y*r_y) / site_size; 
	double rest_scan[3]; 
	double scan_force[3]; 
	for(int i_scan = -1; i_scan <= 1; i_scan++){
		rest_scan[i_scan + 1] = rough_rest_dist + (i_scan * 0.5);
		double rest_scan_length = rest_scan[i_scan + 1] * site_size;
		double r_scan = sqrt(r_y*r_y + rest_scan_length*rest_scan_length); 
		if(r_scan >= r_0_)
			scan_force[i_scan + 1] = (r_scan - r_0_) * k_spring_;
		else
			scan_force[i_scan + 1] = (r_scan - r_0_) * k_slack_; 
	}
	double min_force = 100; 
	for(int i_scan = -1; i_scan <=1; i_scan++){
		double force = fabs(scan_force[i_scan + 1]); 
		if(force < min_force){
			min_force = force;
			rest_dist_ = rest_scan[i_scan + 1]; 
		}	
	}
	comp_cutoff_ = 0;
	/* Next, calculate compression distance cutoff */ 
	for(int x_dist = (int)rest_dist_; x_dist >= 0; x_dist--){
		int r_x = x_dist * site_size; 
		double r = sqrt(r_y*r_y + r_x*r_x); 
		double dr = r - r_0_; 
		double boltzmann_weight; 
		if(r < r_0_)
			boltzmann_weight = exp(-dr*dr*k_slack_/(2*kbT));
		else
			boltzmann_weight = exp(-dr*dr*k_spring_/(2*kbT));
		if(boltzmann_weight < 1e-9){
			comp_cutoff_ = x_dist;
			break;
		}
	}
	/* Finally, calculate extension distance cutoff */
	for(int x_dist = (int) rest_dist_; x_dist < 1000; x_dist++){
		int r_x = x_dist * site_size;
		double r = sqrt(r_y*r_y + r_x*r_x);
		double dr = r - r_0_; 
		double boltzmann_weight;
		if(r >= r_0_)
			boltzmann_weight = exp(-dr*dr*k_spring_/(2*kbT));
		else
			boltzmann_weight = exp(-dr*dr*k_slack_/(2*kbT));
		if(boltzmann_weight < 1e-9){
			dist_cutoff_ = x_dist;
			break;
		}
	}
}

void Kinesin::InitiateNeighborLists(){

	int n_mts = parameters_->microtubules.count; 
	// Serialize this bitch so we just roll one random number 
	neighbor_xlinks_.resize(n_mts*(2*dist_cutoff_ + 1));
	neighbor_sites_.resize(n_mts*(2*dist_cutoff_ + 1));
}

/* XXX currently the exact same as binding lookup table...appropriate? XXX */
void Kinesin::PopulateTetheringLookupTable(){

	double r_y = parameters_->microtubules.y_dist / 2;
	double kbT = parameters_->kbT;
	double site_size = parameters_->microtubules.site_size;
	tethering_weight_lookup_.resize(2*dist_cutoff_ + 1);
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		double r_x = x_dist_dub * site_size / 2;
		double r = sqrt(r_x*r_x + r_y*r_y);
		double dr = r - r_0_;
		double weight = 0;
		if(x_dist_dub < 2*comp_cutoff_){
			weight = 0;
		}
		// For compression, assume the tail can bend, lowering its k_spring
		else if(dr < 0){
			weight = exp(-dr*dr*k_slack_/(2*kbT));
		}
		// For extension, treat tail as a spring
		else{
			weight = exp(-dr*dr*k_spring_/(2*kbT));
		}
		tethering_weight_lookup_[x_dist_dub] = weight;
//		printf("TETH_ #%i_ r: %g, dr: %g, weight: %g\n", 
//				x_dist_dub, r, dr, weight);
	}
}

void Kinesin::PopulateBindingLookupTable(){

	double r_y = parameters_->microtubules.y_dist / 2;
	double kbT = parameters_->kbT;
	double site_size = parameters_->microtubules.site_size;
	// Use 2*dist_cutoff because we can have half-integer distances
	binding_weight_lookup_.resize(2*dist_cutoff_ + 1);	
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		double r_x = x_dist_dub * site_size / 2;
		double r = sqrt(r_x*r_x + r_y*r_y);
		double dr = r - r_0_;
		double weight = 0;
		if(x_dist_dub < 2*comp_cutoff_){
			weight = 0;
		}
		// For compression, assume the tail can bend, lowering its k_spring
		else if(dr < 0){
			weight = exp(-dr*dr*k_slack_/(2*kbT));
		}
		// For extension, treat tail as a spring
		else{
			weight = exp(-dr*dr*k_spring_/(2*kbT));
		}
		binding_weight_lookup_[x_dist_dub] = weight;
//		printf("BIND_ #%i_ r: %g, dr: %g, weight: %g\n", 
//				x_dist_dub, r, dr, weight);
	}
}

void Kinesin::UpdateNeighborXlinks(){

	if(tethered_ == false
	&& heads_active_ == 2){
		n_neighbor_xlinks_ = 0;
		int n_mts = parameters_->microtubules.count;
		int mt_length = parameters_->microtubules.length;
		double stalk_coord = GetStalkCoordinate();
		// Scan through all potential neighbor sites; add untethered xlinks
		int i_entry = 0;
		// FIXME this only works for two MTs as of now FIXME
		for(int i_mt = 0; i_mt < n_mts; i_mt++){
			Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
			int i_stalk = stalk_coord - mt_->coord_;
//			printf("i_stalk is %i\n\n", i_stalk);
			for(int x_dist = -dist_cutoff_; x_dist <= dist_cutoff_; x_dist++){
				int i_site = i_stalk + x_dist;
//				printf("i_site is %i\n", i_site);
				// XXX BOUNDARY SITES ACCESSIBLE -- DISABLE FOR ALPHA/BETA
				if(i_site < 0){
					x_dist -= (i_site + 1);
				}
				// End scan once last bulk site (mt_length - 2) is reached
				else if(i_site > mt_length - 1){
					break;
				}
				else{
					Tubulin *site = &mt->lattice_[i_site];
					if(site->xlink_ != nullptr){
						AssociatedProtein *xlink = site->xlink_;
						double anchor_coord = xlink->GetAnchorCoordinate();
						double x_dist = fabs(anchor_coord - stalk_coord);
						int x_dist_dub = 2*x_dist;
						if(x_dist_dub >= 2*comp_cutoff_ 
						&& x_dist_dub <= 2*dist_cutoff_
						&& xlink->tethered_ == false){
							if(xlink->heads_active_ == 1){
								n_neighbor_xlinks_++;
								neighbor_xlinks_[i_entry] = xlink; 
								i_entry++;
							}
							else if(xlink->heads_active_ == 2
							&& i_mt == 0){				// FIXME this is bad
								n_neighbor_xlinks_++;
								neighbor_xlinks_[i_entry] = xlink;
								i_entry++;		
							}
						}
					}
				}
			}
		}
	}
	else{
		printf("error in update neighbor xlinks\n");
		exit(1);
	}
}

void Kinesin::UpdateNeighborSites(){

	if(tethered_ == true
	&& heads_active_ == 0){
		n_neighbor_sites_ = 0;
		int n_mts = parameters_->microtubules.count;
		int mt_length = parameters_->microtubules.length;
		double anchor_coord = xlink_->GetAnchorCoordinate();
//		printf("anchor coord is %g\n\n", anchor_coord);
		// Scan through all potential neighbor sites; add unoccupied to list 
		int i_entry = 0;
		for(int i_mt = 0; i_mt < n_mts; i_mt++){ 	
			Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
			double mt_coord = mt->coord_;
			int i_anchor = anchor_coord - mt_coord; 
			for(int x_dist = -dist_cutoff_; x_dist <= dist_cutoff_; x_dist++){
				int i_site = i_anchor + x_dist; 
//				printf("i_site is %i (x_dist %i)\n", i_site, x_dist);
				// Start index at first bulk site (1) if site index is <= 0
				// XXX BOUNDARY SITES ACCESSIBLE -- DISABLE FOR ALPHA/BETA XXX
				if(i_site < 0){
					x_dist -= (i_site + 1);
				}
				// End scan at last bulk site (mt_length - 2)
				else if(i_site > mt_length - 1){
					break;
				}
				else{
					Tubulin *neighbor = &mt->lattice_[i_site];
					double site_coord = i_site + neighbor->mt_->coord_; 
					double x_dist = fabs(anchor_coord - site_coord); 
					int x_dist_dub = 2*x_dist;
					if(x_dist_dub >= 2*comp_cutoff_
					&& x_dist_dub <= 2*dist_cutoff_
					&& neighbor->occupied_ == false){
						n_neighbor_sites_++;
						neighbor_sites_[i_entry] = neighbor;
						i_entry++;
					}
				}
			}   
		}
	}
	else{
		printf("error in update neighbor sites\n");
		exit(1);
	}
}

void Kinesin::UpdateExtension(){

	double site_size = parameters_->microtubules.site_size;
	if(heads_active_ == 0
	|| tethered_ == false){
		x_dist_doubled_ = 0;
		extension_ = 0;
	}
	else if(heads_active_ == 1){
		int x_dub_pre = x_dist_doubled_;
		Tubulin *site = GetActiveHeadSite();
		double site_coord = site->index_ + mt_->coord_;
		double anchor_coord = xlink_->GetAnchorCoordinate();
		double x_dist = fabs(anchor_coord - site_coord);
		int x_dist_dub = 2*x_dist;
		x_dist_doubled_ = x_dist_dub;
		if(x_dist_doubled_ > 2*dist_cutoff_
		|| x_dist_doubled_ < 2*comp_cutoff_){
			ForceUntether(x_dub_pre);
//			printf("forced SINGLE-HEAD UNTETHER ???\n");
//			fflush(stdout);
		}
		else{
			double r_y = parameters_->microtubules.y_dist / 2;
			double r_x = site_size * x_dist_doubled_ / 2;
			double r = sqrt(r_x*r_x + r_y*r_y);
			double extension = r - r_0_;
			extension_ = extension;
			cosine_ = r_x / r; 
		}
	}
	else if(heads_active_ == 2){
		int x_dub_pre = x_dist_doubled_;
		double stalk_coord = GetStalkCoordinate();
		double anchor_coord = xlink_->GetAnchorCoordinate();
		double x_dist = fabs(anchor_coord - stalk_coord);
		int x_dist_dub = 2*x_dist;
		x_dist_doubled_ = x_dist_dub; 
		if(x_dist_doubled_ > 2*dist_cutoff_
		|| x_dist_doubled_ < 2*comp_cutoff_){
			ForceUntether(x_dub_pre); 
//			printf("forced an untether >:O\n");
//			fflush(stdout);
//			exit(1);
		}
		else{
			// Calculate new extension
			double r_y = parameters_->microtubules.y_dist / 2;
			double r_x = site_size * x_dist_doubled_/2;
			double r = sqrt(r_y*r_y + r_x*r_x);
			double extension = r - r_0_; 
			extension_ = extension;
			cosine_ = r_x / r;
		}
	}
	else{
		printf("ok error in update motor extension\n");
		exit(1);
	}
}

void Kinesin::ForceUntether(int x_dub_pre){
	
	// Update statistics
	if(xlink_->heads_active_ == 1){
		properties_->prc1.n_sites_i_untethered_++;
		properties_->prc1.n_sites_i_tethered_[x_dub_pre]--;
	}
	else if(xlink_->heads_active_ == 2){
		int x_dist = xlink_->x_dist_;
		properties_->prc1.n_sites_ii_untethered_[x_dist] += 2;
		properties_->prc1.n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
	}
	if(heads_active_ == 2){
		properties_->kinesin4.n_bound_untethered_++;
		properties_->kinesin4.n_bound_tethered_tot_--;
		properties_->kinesin4.n_bound_tethered_[x_dub_pre]--;
	}
	properties_->prc1.n_untethered_++;
	// Update motor	
	tethered_ = false;
	x_dist_doubled_ = 0;
	extension_ = 0; 
	// Update xlink
	xlink_->motor_ = nullptr;
	xlink_->tethered_ = false; 
	xlink_ = nullptr;
}	

bool Kinesin::AtCutoff(){

	int dx = mt_->delta_x_; 
	int drest = GetDirectionTowardRest();

	if((x_dist_doubled_  >= (2*dist_cutoff_ - 1) && dx == -drest)
		 ||(x_dist_doubled_ <= (2*comp_cutoff_ + 1) && dx == -drest))
		return true;
	else
		return false;
}

int Kinesin::GetDirectionTowardRest(){

	if(tethered_ == true){
		if(heads_active_ == 1){
			Tubulin *site = GetActiveHeadSite();
			int site_coord = site->index_ + mt_->coord_;
			double rest_coord = GetRestLengthCoordinate();
			if(site_coord > rest_coord)
				return -1;
			else if(site_coord < rest_coord)
				return 1;
			else{
				double anchor_coord = xlink_->GetAnchorCoordinate();
				if(site_coord > anchor_coord)
					return 1;
				else
					return -1;
			}
		}
		else if(heads_active_ == 2){
			double stalk_coord = GetStalkCoordinate();
			double rest_coord = GetRestLengthCoordinate();
			if(stalk_coord > rest_coord)
				return -1;
			else if(stalk_coord < rest_coord)
				return 1;
			else{
				double anchor_coord = xlink_->GetAnchorCoordinate();
				if(stalk_coord > anchor_coord)
					return 1;
				else
					return -1;
			}
		}
		else{
			printf("this is a different type of error in get dir torest\n");
			exit(1);
		}
	}
	else{
		printf("error in get dir. toward xlink\n");
		exit(1);
	}
}

double Kinesin::GetRestLengthCoordinate(){

	if(tethered_ == true){
		if(heads_active_ == 1){
			Tubulin *site = GetActiveHeadSite();
			int site_coord = site->index_ + mt_->coord_;
			double anchor_coord = xlink_->GetAnchorCoordinate();
			if(site_coord > anchor_coord){
				double rest_coord = anchor_coord + rest_dist_;
				return rest_coord;
			}
			else if(site_coord < anchor_coord){
				double rest_coord = anchor_coord - rest_dist_;
				return rest_coord;
			}
			else{
				printf("lameness in get rest length coord kinesin ONE\n");
				exit(1);
				}
			}
			else if(heads_active_ == 2){
				double stalk_coord = GetStalkCoordinate();
				double anchor_coord = xlink_->GetAnchorCoordinate();
				if(stalk_coord > anchor_coord){
					double rest_coord = anchor_coord + rest_dist_;
					return rest_coord;
				}
				else if(stalk_coord < anchor_coord){
					double rest_coord = anchor_coord - rest_dist_;
					return rest_coord;
				}
				else{
					printf("lameness in get rest length coord kinesin TWO\n");
					printf("stalk: %g | anchor: %g\n", stalk_coord, 
							anchor_coord);
					exit(1);
				}
			}
			else{
				printf("god damnit in get rest length coord (kinesin)\n");
				exit(1);
			}
	}
	else{
		printf("error in get_rest_length_coord (kinesin)\n");
		exit(1);
	}
}

double Kinesin::GetStalkCoordinate(){

	if(heads_active_ == 1){
		Tubulin *site = GetActiveHeadSite();
		int i_site = site->index_;
		int mt_coord = site->mt_->coord_;
		double site_coord = i_site + mt_coord;
		return site_coord; 
	}
	else if(heads_active_ == 2){
		int i_one = front_site_->index_;
		int i_two = rear_site_->index_;
		double avg_index = ((double)(i_one + i_two))/2;
		int mt_coord = mt_->coord_;
		double stalk_coord = mt_coord + avg_index;
		return stalk_coord;
	}
	else{
		printf("error: can NOT get this stalk index bro.\n");
		exit(1);
	}
}

double Kinesin::GetTetheringWeight(AssociatedProtein *xlink){

	double stalk_coord = GetStalkCoordinate();
	double anchor_coord = xlink->GetAnchorCoordinate();
	double x_dist = fabs(anchor_coord - stalk_coord);
	// Multiply by 2 to guarentee an integer
	int adj_distance = x_dist*2; 
	// Get tethering weight that corressponds to this adj. x-dist
	double weight = tethering_weight_lookup_[adj_distance];
	return weight; 
}

double Kinesin::GetBindingWeight(Tubulin *site){

	// Since only one foot of the motor attaches initially, 
	// we only need to consider the coordinate of one site
	// (assuming the diffusion of the other foot avgs out) 
	int i_site = site->index_;	
	int mt_coord = site->mt_->coord_; 
	double site_coord = (double)(mt_coord + i_site);
	double anchor_coord = xlink_->GetAnchorCoordinate();
	double x_dist = fabs(anchor_coord - site_coord); 
	// Multiply by 2 to guarentee an integer 
	int adj_distance = x_dist*2; 
	// Get binding weight that corresponds to this adj. x-dist
	double weight =  binding_weight_lookup_[adj_distance]; 
	return weight;
}

double Kinesin::GetTetherForce(Tubulin *site){

	if(tethered_ == true
	&& heads_active_ > 0){
		UpdateExtension();
		// Make sure we didn't force an untether event
		if(tethered_ == true){
			double force_mag;
			if(x_dist_doubled_ <= 2*rest_dist_)
				force_mag = extension_ * k_slack_;
			else
				force_mag = extension_ * k_spring_; 
			double stalk_coord = GetStalkCoordinate();
			double anchor_coord = xlink_->GetAnchorCoordinate();
			int newtons_third;
			if(site->mt_ == mt_)
				newtons_third = -1;
			else
				newtons_third = 1;	
			double force;
			if(anchor_coord < stalk_coord){
				force = newtons_third * force_mag * cosine_;
			}
			else
				force = -newtons_third * force_mag * cosine_; 
			return force;
		}
	}
	else{
		printf("error in get teth force (motor)\n");
		exit(1);
	}
}

Tubulin* Kinesin::GetActiveHeadSite(){

	if(heads_active_ == 1){
		if(front_site_ != nullptr
		&& rear_site_ == nullptr)
			return front_site_;
		else if(rear_site_ != nullptr
		&& front_site_ == nullptr)
			return rear_site_;
		else{
			printf("bad error in kinesin head tracking\n");
			exit(1);
		}
	}
	else{
		printf("how am i supposed to find an ACTIVE KINESIN HEAD?!\n");
		exit(1);
	}
}

Tubulin* Kinesin::GetSiteCloserToRest(){

	if(tethered_ == true
	&& heads_active_ == 2){
		double anchor_coord = xlink_->GetAnchorCoordinate();
		double stalk_coord = GetStalkCoordinate();
		double rest_coord;
		if(stalk_coord > anchor_coord)
			rest_coord = anchor_coord + rest_dist_;
		else
			rest_coord = anchor_coord - rest_dist_;
		int front_coord = front_site_->index_ + mt_->coord_;
		int rear_coord = rear_site_->index_ + mt_->coord_;
		double dist_front = abs(rest_coord - front_coord);
		double dist_rear = abs(rest_coord - rear_coord);
		if(dist_front < dist_rear)
			return front_site_;
		else
			return rear_site_;
	}
	else{
		printf("Error in getting site closest to xlink (motor)\n");
		exit(1);
	}
}

Tubulin* Kinesin::GetSiteFartherFromRest(){

	if(tethered_ == true
	&& heads_active_ == 2){
		double anchor_coord = xlink_->GetAnchorCoordinate();
		double stalk_coord = GetStalkCoordinate();
		double rest_coord;
		if(stalk_coord > anchor_coord)
			rest_coord = anchor_coord + rest_dist_;
		else
			rest_coord = anchor_coord - rest_dist_;
		int front_coord = front_site_->index_ + mt_->coord_;
		int rear_coord = rear_site_->index_ + mt_->coord_;	
		double dist_front = abs(rest_coord - front_coord);
		double dist_rear = abs(rest_coord - rear_coord);
		if(dist_front > dist_rear)
			return front_site_;
		else
			return rear_site_;
	}
	else{
		printf("Error in getting site farthest xlink (motor)\n");
		exit(1);
	}
}

Tubulin* Kinesin::GetWeightedNeighborSite(int binding_affinity){

	UpdateNeighborSites();
	Tubulin* eligible_neighbor_sites[n_neighbor_sites_]; 
	int n_eligible = 0;
	int i_entry = 0;
	// Scan through all neighbor sites to extract eligible ones
	for(int i_site = 0; i_site < n_neighbor_sites_; i_site++){
		Tubulin* site = neighbor_sites_[i_site];
		if(site->binding_affinity_ == binding_affinity){
			eligible_neighbor_sites[i_entry] = site;
			i_entry++;
			n_eligible++;
		}
	}
	double anch_coord = xlink_->GetAnchorCoordinate();
	double p_tot = 0;
	// Get total binding probability of all eligible sites for normalization
	for(int i_site = 0; i_site < n_eligible; i_site++){
		Tubulin* site = eligible_neighbor_sites[i_site];
		double site_coord = site->mt_->coord_ + site->index_;
		double x_dist = fabs(anch_coord - site_coord);
		int x_dist_dub = 2 * x_dist;
		p_tot += binding_weight_lookup_[x_dist_dub];
	}
	// Scan through eligible sites; pick one randomly based on weights
	double ran = properties_->gsl.GetRanProb();
	double p_cum = 0;
	for(int i_site = 0; i_site < n_eligible; i_site++){
		Tubulin* site = eligible_neighbor_sites[i_site];
		double site_coord = site->mt_->coord_ + site->index_;
		double x_dist = fabs(anch_coord - site_coord);
		int x_dist_dub = 2 * x_dist;
		p_cum += binding_weight_lookup_[x_dist_dub] / p_tot; 
		if(ran < p_cum){
			return site; 
		}
	}
	return nullptr;
}

AssociatedProtein* Kinesin::GetWeightedNeighborXlink(){

	UpdateNeighborXlinks();
	double stalk_coord = GetStalkCoordinate();
	double p_tot = 0;
	for(int i_xlink = 0; i_xlink < n_neighbor_xlinks_; i_xlink++){
		AssociatedProtein* xlink = neighbor_xlinks_[i_xlink];
		double anch_coord = xlink->GetAnchorCoordinate();
		double x_dist = fabs(stalk_coord - anch_coord);
		int x_dist_dub = 2 * x_dist;
		p_tot += tethering_weight_lookup_[x_dist_dub];
	}
	double ran = properties_->gsl.GetRanProb();
	double p_cum = 0; 
	for(int i_xlink = 0; i_xlink < n_neighbor_xlinks_; i_xlink++){
		AssociatedProtein* xlink = neighbor_xlinks_[i_xlink];
		double anch_coord = xlink->GetAnchorCoordinate();
		double x_dist = fabs(stalk_coord - anch_coord);
		int x_dist_dub = 2 * x_dist;
		p_cum += tethering_weight_lookup_[x_dist_dub] / p_tot;
		if(ran < p_cum){
			return xlink;
		}
	}
	return nullptr;
}
