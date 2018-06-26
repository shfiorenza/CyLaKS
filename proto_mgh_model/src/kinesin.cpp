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
	InitiateNeighborLists(); 
	PopulateTetheringLookupTable();
	PopulateBindingLookupTable();
}

void Kinesin::SetParameters(){

	kbT_ = parameters_->kbT;
	site_size_ = parameters_->site_size;
	r_0_ = parameters_->r_0_motor;
	k_spring_ = parameters_->k_spring_motor;
	k_eff_slack_ = parameters_->k_slack_motor;
	stall_force_ = parameters_->stall_force;
}

void Kinesin::InitiateNeighborLists(){

	neighbor_xlinks_.resize(2*dist_cutoff_ + 1);
	int n_mts = parameters_->n_microtubules; 
	// Serialize this bitch so we just roll one random number 
	neighbor_sites_.resize(n_mts*(2*dist_cutoff_ + 1));
}

/* XXX currently the exact same as binding lookup table...appropriate? XXX */
void Kinesin::PopulateTetheringLookupTable(){

	double r_y = 17.5;
	tethering_weight_lookup_.resize(2*dist_cutoff_ + 1);
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		double r_x = x_dist_dub * site_size_ / 2;
		double r = sqrt(r_x*r_x + r_y*r_y);
		double dr = r - r_0_;
		double weight = 0;
		if(x_dist_dub < 2*comp_cutoff_){
			weight = 0;
		}
		// For compression, assume the tail can bend, lowering its k_spring
		else if(dr < 0){
			weight = exp(-dr*dr*k_eff_slack_/(2*kbT_));
		}
		// For extension, treat tail as a spring
		else{
			weight = exp(-dr*dr*k_spring_/(2*kbT_));
		}
		tethering_weight_lookup_[x_dist_dub] = weight;
//		printf("TETH_ #%i_ r: %g, dr: %g, weight: %g\n", 
//				x_dist_dub, r, dr, weight);
	}
}

void Kinesin::PopulateBindingLookupTable(){

	double r_y = 17.5;		// dist from either MT to midpoint of prc1
	// Use 2*dist_cutoff because we can have half-integer distances
	binding_weight_lookup_.resize(2*dist_cutoff_ + 1);	
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		double r_x = x_dist_dub * site_size_ / 2;
		double r = sqrt(r_x*r_x + r_y*r_y);
		double dr = r - r_0_;
		double weight = 0;
		if(x_dist_dub < 2*comp_cutoff_){
			weight = 0;
		}
		// For compression, assume the tail can bend, lowering its k_spring
		else if(dr < 0){
			weight = exp(-dr*dr*k_eff_slack_/(2*kbT_));
		}
		// For extension, treat tail as a spring
		else{
			weight = exp(-dr*dr*k_spring_/(2*kbT_));
		}
		binding_weight_lookup_[x_dist_dub] = weight;
//		printf("BIND_ #%i_ r: %g, dr: %g, weight: %g\n", 
//				x_dist_dub, r, dr, weight);
	}
}

void Kinesin::UpdateNeighborXlinks(){

	n_neighbor_xlinks_ = 0;
	int n_mts = parameters_->n_microtubules;
	int mt_length = parameters_->length_of_microtubule;
	double stalk_coord = GetStalkCoordinate();
	// Scan through all potential neighbor sites; add only untethered xlinks
	int i_entry = 0;
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
		int offset = mt->coord_;
//		printf("for mt %i, stalk %g: ", i_mt, stalk_coord); 
		for(int x_dist = -dist_cutoff_; x_dist <= dist_cutoff_; x_dist++){
			int i_site = (int)stalk_coord + x_dist - offset;	// FIXME 
			// Start index at first bulk site (1) if site index is 0 or less
			// XXX BOUNDARY SITES ACCESSIBLE -- DISABLE FOR ALPHA/BETA XXX
			if(i_site < 0){
				x_dist -= (i_site - 1);
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
					double x_dist = abs(anchor_coord - stalk_coord);
					if(x_dist >= comp_cutoff_ && x_dist <= dist_cutoff_){
						if(xlink->heads_active_ == 1
						&& xlink->tethered_ == false){
							n_neighbor_xlinks_++;
							neighbor_xlinks_[i_entry] = xlink; 
							i_entry++;
						}
						else if(xlink->heads_active_ == 2
						&& xlink->tethered_ == false
						&& i_mt == 0){					// FIXME this is bad
//							printf("%iii_", i_site);
							n_neighbor_xlinks_++;
							neighbor_xlinks_[i_entry] = xlink;
							i_entry++;		
						}
					}
				}
			}
		}
		//		printf("\n");
	}
}

bool Kinesin::NeighborXlinkExists(int x_dist_doubled){

	AssociatedProtein *neighbor_xlink = GetNeighborXlink(x_dist_doubled);
	if(neighbor_xlink == nullptr)
		return false;
	else
		return true;
}

void Kinesin::UpdateNeighborSites(){

	if(tethered_ == true
	&& heads_active_ == 0){
		n_neighbor_sites_ = 0;
		int n_mts = parameters_->n_microtubules;
		int mt_length = parameters_->length_of_microtubule;
		double anchor_coord = xlink_->GetAnchorCoordinate();
		// FIXME this only works for two MTs as of now FIXME
		// Scan through all potential neighbor sites; add unoccupied to list 
		int i_entry = 0;
		for(int i_mt = 0; i_mt < n_mts; i_mt++){ 	
			Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
			int offset = mt->coord_;
			for(int x_dist = -dist_cutoff_; x_dist <= dist_cutoff_; x_dist++){
				int i_site = (int)anchor_coord + x_dist - offset;
				// Start index at first bulk site (1) if site index is <= 0
				// XXX BOUNDARY SITES ACCESSIBLE -- DISABLE FOR ALPHA/BETA XXX
				if(i_site < 0){
					x_dist -= (i_site - 1);
				}
				// End scan at last bulk site (mt_length - 2)
				else if(i_site > mt_length - 1){
					break;
				}
				else if(abs(x_dist) >= comp_cutoff_){
					Tubulin *neighbor = &mt->lattice_[i_site];
					if(neighbor->occupied_ == false){
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

bool Kinesin::NeighborSiteExists(int x_dist_doubled){

		Tubulin *neighbor_site = GetNeighborSite(x_dist_doubled);
		if(neighbor_site == nullptr)
			return false;
		else
			return true;
}

void Kinesin::UpdateExtension(){

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
			double r_y = 17.5;
			double r_x = site_size_ * x_dist_doubled_ / 2;
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
			double r_y = 17.5;
			double r_x = site_size_*x_dist_doubled_/2;
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
	properties_->prc1.n_tethered_--;
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

int Kinesin::SampleTailExtensionDoubled(){

	int rest_dist_dub = 2*rest_dist_;
	// Total weight of extension 'side' of energy profile
	double ext_weight = 0;
	// Total weight of slack 'side' of energy profile
	double slack_weight = 0;
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		if(x_dist_dub < rest_dist_dub)
			slack_weight += tethering_weight_lookup_[x_dist_dub];
		else
			ext_weight += tethering_weight_lookup_[x_dist_dub];
	}
	double tot_weight = ext_weight + slack_weight;
	// Roll to determine which side of the profile we sample
	RandomNumberManagement* gsl = &properties_->gsl;
	double ran = properties_->gsl.GetRanProb();
	// Sample a normal distribution around rest length if we get extension
	if(ran < ext_weight/tot_weight){
		double sigma = sqrt(kbT_ / k_spring_) / site_size_;
		int extension = gsl->SampleAbsNormalDist(sigma);
		int x_dist_dub = rest_dist_dub + extension;
		if(x_dist_dub > 2*dist_cutoff_)
			x_dist_dub = 2*dist_cutoff_;
		return x_dist_dub; 
	}
	// Otherwise, select for 'slack' based on different sigma
	else{
		double sigma = sqrt(kbT_ / k_eff_slack_) / site_size_;
		int compression = gsl->SampleAbsNormalDist(sigma);
		int x_dist_dub = rest_dist_dub - compression;
		if(x_dist_dub < 2*comp_cutoff_)
			x_dist_dub = 2*comp_cutoff_;
		return x_dist_dub;
	}
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
		AssociatedProtein *xlink = xlink_; 
		UpdateExtension();
		// Make sure we didn't force an untether event
		if(tethered_ == true){
			double force_mag;
			if(x_dist_doubled_ <= 2*rest_dist_)
				force_mag = extension_ * k_eff_slack_;
			else
				force_mag = extension_ * k_spring_; 
			double stalk_coord = GetStalkCoordinate();
			double anchor_coord = xlink->GetAnchorCoordinate();
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

Tubulin* Kinesin::GetNeighborSite(int x_dist_doubled){

	UpdateNeighborSites();
	double x_dist = ((double)x_dist_doubled)/2;
	if(tethered_ == true
	&& heads_active_ == 0){
		double anchor_coord = xlink_->GetAnchorCoordinate();
		for(int i_entry = 0; i_entry < n_neighbor_sites_; i_entry++){
			Tubulin *neighbor_site = neighbor_sites_[i_entry];
			int neighb_index = neighbor_site->index_;
			int mt_coord = neighbor_site->mt_->coord_; 
			double neighb_coord = mt_coord + neighb_index; 
			double neighb_dist = fabs(anchor_coord - neighb_coord);
			if(neighb_dist == x_dist)
				return neighbor_site;
		}
		return nullptr;
	}
	else{
		printf("error: neighb sites but untethered.\n");
		exit(1);
	}
}

AssociatedProtein* Kinesin::GetNeighborXlink(int x_dist_doubled){

	UpdateNeighborXlinks();
	double x_dist = (double)x_dist_doubled/2; 
	if(heads_active_ == 2){
		double stalk_coord = GetStalkCoordinate();
		for(int i_entry = 0; i_entry < n_neighbor_xlinks_; i_entry++){
			AssociatedProtein *neighbor_xlink = neighbor_xlinks_[i_entry];
			double anchor_coord = neighbor_xlink->GetAnchorCoordinate();
			double neighb_dist = fabs(stalk_coord - anchor_coord);
//			printf("n_d: %g, x_d: %g\n", neighb_dist, x_dist);
			if(neighb_dist == x_dist)
				return neighbor_xlink;
		}
		return nullptr;
	}
	else{
		printf("error: neighb xlinks but not double bound\n");
		exit(1);
	}
}
