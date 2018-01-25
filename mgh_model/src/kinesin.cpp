#include "master_header.h"
#include "kinesin.h"

Kinesin::Kinesin(){
}

void Kinesin::Initialize(system_parameters *parameters, 
	system_properties *properties, int ID){

	ID_ = ID; 
	parameters_ = parameters;
	properties_ = properties;
	InitiateNeighborLists(); 
	PopulateTetheringLookupTable();
	PopulateBindingLookupTable();
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
	for(int i_dist = 0; i_dist <= 2*dist_cutoff_; i_dist++){
		double r_x = (((double)i_dist)/2)*site_size_;
		double r = sqrt(r_x*r_x + r_y*r_y);
		double dr = r - r_0_;
		double weight = 0;
		// For extension, treat tail as a spring
		if(dr >= 0){
			weight = exp(-dr*dr*k_spring_/(2*kbT_));
		}
		// For compression, assume tail is floppy enough so that it just
		// bends slightly and barely affects binding weight vs rest length
		else if(dr >= (-r_0_/6)){
			weight = 0.25;
		}	
		else if(dr >= (-r_0_/3)){
			weight = 0.20;
		}
		else if(dr >= (-r_0_/2)){
			weight = 0.15;
		}
		else if(dr >= (-2*r_0_/3)){
			weight = 0.10;
		}
		else if(dr >= (-5*r_0_/6)){
			weight = 0.05;
		}
		else{
			weight = 0.005;
		}
		tethering_weight_lookup_[i_dist] = weight;
	}
}

void Kinesin::PopulateBindingLookupTable(){

	double r_y = 17.5;		// dist from either MT to midpoint of prc1
	// Use 2*dist_cutoff because we can have half-integer distances
	binding_weight_lookup_.resize(2*dist_cutoff_ + 1);	
	for(int i_dist = 0; i_dist <= 2*dist_cutoff_; i_dist++){
		double r_x = (((double)i_dist)/2)*site_size_;
		double r = sqrt(r_x*r_x + r_y*r_y);
		double dr = r - r_0_;
		double weight = 0;
		// For extension, treat tail as a spring
		if(dr >= 0){
			weight = exp(-dr*dr*k_spring_/(2*kbT_));
		}
		// For compression, assume tail is floppy enough so that it just
		// bends slightly and barely affects binding weight vs rest length
		else if(dr >= (-r_0_/6)){
			weight = 0.5;
		}	
		else if(dr >= (-r_0_/3)){
			weight = 0.25;
		}
		else if(dr >= (-r_0_/2)){
			weight = 0.10;
		}
		else if(dr >= (-2*r_0_/3)){
			weight = 0.05;
		}
		else if(dr >= (-5*r_0_/6)){
			weight = 0.0005;
		}
		else{
			weight = 0.00005;
		}
		binding_weight_lookup_[i_dist] = weight;
//		printf("#%i_ r: %g, dr: %g, weight: %g\n", i_dist, r, dr, weight);
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
			if(i_site <= 0){
				x_dist -= (i_site - 1);
			}
			// End scan once last bulk site (mt_length - 2) is reached
			else if(i_site >= mt_length - 1){
				break;
			}
			else{
				Tubulin *site = &mt->lattice_[i_site];
				if(site->xlink_ != nullptr){
					AssociatedProtein *xlink = site->xlink_;
					if(xlink->heads_active_ == 1
					&& xlink->tethered_ == false){
//						printf("%ii_", i_site);
						n_neighbor_xlinks_++;
						neighbor_xlinks_[i_entry] = xlink; 
						i_entry++;
					}
					else if(xlink->heads_active_ == 2
					&& xlink->tethered_ == false
					&& i_mt == 0){					// FIXME this is bad
						double anchor_coord = xlink->GetAnchorCoordinate();
						double extension = fabs(anchor_coord - stalk_coord);
						if(extension <= dist_cutoff_){
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
				if(i_site <= 0){
					x_dist -= (i_site - 1);
				}
				// End scan at last bulk site (mt_length - 2)
				else if(i_site >= mt_length - 1){
					break;
				}
				else{
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

	int x_dub_pre = x_dist_doubled_;
	if(heads_active_ == 0
	|| tethered_ == false){
		x_dist_doubled_ = 0;
		extension_ = 0;
	}
	else if(heads_active_ == 1){
		printf("ok error in update extension (motor)\n");
		exit(1);
	}
	else if(heads_active_ == 2){
		double stalk_coord = GetStalkCoordinate();
		double anchor_coord = xlink_->GetAnchorCoordinate();
		double x_dist = fabs(anchor_coord - stalk_coord);
		int x_dist_dub = 2*x_dist;
		x_dist_doubled_ = x_dist_dub; 
		if(x_dist_doubled_ > 2*dist_cutoff_){
			ForceUntether(x_dub_pre); 
			printf("forced an untether >:O\n");
//			exit(1);
			fflush(stdout);
		}
		else{
			// Calculate new extension
			double r_y = 17.5;
			double r_x = site_size_*x_dist_doubled_/2;
			double r = sqrt(r_y*r_y + r_x*r_x);
			double extension = r - r_0_; 
			if(extension < 0)
				extension = 0;
			extension_ = extension;
		}
	}
}

void Kinesin::ForceUntether(int x_dub_pre){

	if(xlink_->heads_active_ == 1){
		properties_->prc1.n_sites_single_untethered_++;
		properties_->prc1.n_sites_single_tethered_--;
	}
	else if(xlink_->heads_active_ == 2){
		properties_->prc1.n_sites_double_untethered_ += 2;
		properties_->prc1.n_sites_double_tethered_ -= 2;
	}
	// Update motor	
	tethered_ = false;
	x_dist_doubled_ = 0;
	extension_ = 0; 
	// Update xlink
	xlink_->motor_ = nullptr;
	xlink_->tethered_ = false; 
	xlink_ = nullptr;
	// Update statistics
	properties_->kinesin4.n_bound_untethered_++;
	properties_->kinesin4.n_bound_tethered_tot_--;
	properties_->kinesin4.n_bound_tethered_[x_dub_pre]--;
	properties_->prc1.n_untethered_++;
	properties_->prc1.n_tethered_--;
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
	// Roll a random probability to determine which side of the profile we sample
	double ran = properties_->gsl.GetRanProb();
	// Sample a normal distribution around rest length if we get extension
	if(ran < ext_weight/tot_weight){
		double sigma = ((double)dist_cutoff_)/3; 
		int x_dist_dub = properties_->gsl.SampleNormalDist(sigma, rest_dist_dub);
		if(x_dist_dub > 2*dist_cutoff_)
			x_dist_dub = 2*dist_cutoff_;
		return x_dist_dub; 
	}
	// Otherwise, select for 'slack' based on relative weights
	else{
		double ran2 = properties_->gsl.GetRanProb();
		double cum_weight = 0;
		for(int x_dist_dub = 0; x_dist_dub < rest_dist_dub; x_dist_dub++){
			cum_weight += tethering_weight_lookup_[x_dist_dub] / slack_weight;
			if(cum_weight >= ran2)
				return x_dist_dub; 
		}
	}

}

double Kinesin::GetStalkCoordinate(){

	if(heads_active_ == 2){
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

Tubulin* Kinesin::GetActiveHeadSite(){

	if(heads_active_ == 1){
		if(front_site_ != nullptr)
			return front_site_;
		else if(rear_site_ != nullptr)
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
