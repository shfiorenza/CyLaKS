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
	PopulateBindingLookupTable();
	PopulateTethBindingLookupTable();
}

void AssociatedProtein::SetParameters(){

	r_0_ = parameters_->xlinks.r_0;
	k_spring_ = parameters_->xlinks.k_spring;
	neighbor_sites_.resize(2*dist_cutoff_ + 1);
	int n_mts = parameters_->microtubules.count; 
	int teth_cutoff = properties_->kinesin4.dist_cutoff_; 
	teth_neighbor_sites_.resize(n_mts*(2*teth_cutoff + 1));
}

void AssociatedProtein::PopulateBindingLookupTable(){

	double r_y = 35;		//dist between MTs in nm; static as of now
	double kbT = parameters_->kbT;
	double site_size = parameters_->microtubules.site_size;
	binding_weight_lookup_.resize(dist_cutoff_ + 1);
	for(int i_dist = 0; i_dist <= dist_cutoff_; i_dist++){
		double r_x = (double)i_dist*site_size;
		double r = sqrt(r_y*r_y + r_x*r_x);
		double dr = r - r_0_;
		double weight = exp(-dr*dr*k_spring_/(2*kbT));
		binding_weight_lookup_[i_dist] = weight;
//		printf("r: %g, dr: %g,  weight: %g\n", r, dr, weight);
	}
}

void AssociatedProtein::PopulateTethBindingLookupTable(){

	double r_y = 17.5;
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
			weight = exp(-dr*dr*k_slack/kbT);
		}
		else{
			weight = exp(-dr*dr*k_teth/kbT);
		}
		teth_binding_weight_lookup_[x_dist_dub] = weight;
//		printf("weight for 2x = %i is %g\n", x_dist_dub, weight);
	}
//	properties_->wallace.PauseSim(2);
}

void AssociatedProtein::UpdateNeighborSites(){
	
	n_neighbor_sites_ = 0;
	int n_mts = parameters_->microtubules.count;
	int mt_length = parameters_->microtubules.length;
	if(n_mts > 1){
		Tubulin *site = GetActiveHeadSite();
		int i_site = site->index_;
		Microtubule *mt = site->mt_;
		Microtubule *adj_mt = mt->neighbor_;
		int site_coord = mt->coord_ + i_site; 
		// Scan through all potential neighbors; only add unoccupied to list 
		int i_entry = 0;
		for(int dist = -dist_cutoff_; dist <= dist_cutoff_; dist++){
			int i_neighbor = (site_coord - adj_mt->coord_) + dist;
			// Start index at first bulk site (1) if i_neighb is 0 or neg
			// XXX BOUNDARIES CURRENTLY ACCESSED--DISABLE FOR ALPHA/BETA XXX
			if(i_neighbor < 0){
				dist -= (i_neighbor + 1); // - 1;
			}
			// End scan once last bulk site (mt_length - 2) has been checked
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

	if(tethered_ == true
	&& heads_active_ == 0){
		n_teth_neighbor_sites_ = 0;
		int n_mts = parameters_->microtubules.count;
		int mt_length = parameters_->microtubules.length;
		int teth_cutoff = properties_->kinesin4.dist_cutoff_; 
		int comp_cutoff = properties_->kinesin4.comp_cutoff_;
		double stalk_coord = motor_->GetStalkCoordinate();
//		printf("stalk coord is %g\n\n", anchor_coord);
		// Scan through all potential neighbor sites; add unoccupied to list 
		int i_entry = 0;
		for(int i_mt = 0; i_mt < n_mts; i_mt++){ 	
			Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
			double mt_coord = mt->coord_;
			int i_stalk = stalk_coord - mt_coord; 
			for(int x_dist = -teth_cutoff; x_dist <= teth_cutoff; x_dist++){
				int i_site = i_stalk + x_dist; 
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

void AssociatedProtein::UpdateExtension(){

	double site_size = parameters_->microtubules.site_size;
	if(heads_active_ == 2){
		int x_dist_pre = x_dist_;
		// Calculate first head's coordinate
		int i_head_one = site_one_->index_;
		int mt_coord_one = site_one_->mt_->coord_;
		int coord_one = mt_coord_one + i_head_one;
		// Calculate second head's coordinate
		int i_head_two = site_two_->index_;
		int mt_coord_two = site_two_->mt_->coord_;
		int coord_two = mt_coord_two + i_head_two;
		// Calculate x_distance in # of sites
		int x_dist = abs(coord_one - coord_two);	
		if(x_dist <= dist_cutoff_){
			x_dist_ = x_dist; 
			double r_y = 35;					// static as of now
			double r_x = site_size*x_dist_;
			double r = sqrt(r_y*r_y + r_x*r_x);
			double extension = r - r_0_; 
			extension_ = extension;
			cosine_ = r_x / r;
		}
		else{
			ForceUnbind(x_dist_pre);
//			printf("forced an unbind event >:O\n");
//			fflush(stdout); 
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

	double ran = properties_->gsl.GetRanProb(); 
	if(ran < 0.5){
		// Remove xlink head from site
		site_one_->xlink_ = nullptr;
		site_one_->occupied_ = false;
		// Update xlink details
		site_one_ = nullptr;
		heads_active_--;
		x_dist_ = 0;
		extension_ = 0; 
		cosine_ = 0;
	}
	else{
		// Remove xlink head from site
		site_two_->xlink_ = nullptr;
		site_two_->occupied_ = false;
		// Update xlink details	
		site_two_ = nullptr;	
		heads_active_--;
		x_dist_ = 0;
		extension_ = 0;
		cosine_ = 0;
	}
	// Update statistics
	if(tethered_ == true){
		if(motor_->heads_active_ > 0){
			// Update motor ext (unbinding 2nd head changes anchor coord) 
			int x_dub_pre = motor_->x_dist_doubled_;
			properties_->prc1.n_sites_ii_tethered_
				[x_dub_pre][x_dist_pre] -= 2;
			motor_->UpdateExtension();
			// If untether event didn't occur, update tethered and motor stats
			if(tethered_ == true){
				int x_dub_post = motor_->x_dist_doubled_;
				properties_->prc1.n_sites_i_tethered_[x_dub_post]++;
				// Update kinesin stats
				KinesinManagement *kinesin4 = &properties_->kinesin4;
				if(motor_->heads_active_ == 2){
					kinesin4->n_bound_tethered_[x_dub_pre]--;
					kinesin4->n_bound_tethered_[x_dub_post]++;
				}
			}
			// If untether event DOES occur, counteract force_untether stats
			else{
				// NOT a typo; see kinesin ForceUntether
				properties_->prc1.n_sites_i_tethered_[x_dub_pre]++;
			}
		}
		else{
			properties_->prc1.n_sites_ii_untethered_[x_dist_pre] -= 2;
			properties_->prc1.n_sites_i_untethered_++;
		}
	}
	else{
		properties_->prc1.n_sites_ii_untethered_[x_dist_pre] -= 2;
		properties_->prc1.n_sites_i_untethered_++;
	}
	properties_->prc1.n_double_bound_[x_dist_pre]--;
	properties_->prc1.n_single_bound_++;
}

int AssociatedProtein::GetDirectionTowardRest(Tubulin *site){

	Microtubule *mt = site->mt_;
	if(heads_active_ == 2){
		double anchor_coord = GetAnchorCoordinate();
		int i_site = site->index_;
		int site_coord = mt->coord_ + i_site;
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

int AssociatedProtein::SampleSpringExtension(){

	// Scale sigma so that the avg for binomial is greater than 1
	double kbT = parameters_->kbT;
	double site_size = parameters_->microtubules.site_size;
	double sigma = sqrt(kbT / k_spring_) / site_size * 100;
	int x_dist = properties_->gsl.SampleNormalDist(sigma);
	x_dist = x_dist / 100;
	if(x_dist > dist_cutoff_)
		x_dist = dist_cutoff_;
	else if(x_dist < -dist_cutoff_)
		x_dist = -dist_cutoff_;
	return x_dist;
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
		int index_one = site_one_->index_;
		int mt_coord_one = site_one_->mt_->coord_;
		double coord_one = (double)(mt_coord_one + index_one);
		int index_two = site_two_->index_;
		int mt_coord_two = site_two_->mt_->coord_;
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

double AssociatedProtein::GetBindingWeight(Tubulin *neighbor){

	Tubulin *site = GetActiveHeadSite();
	Microtubule *mt = site->mt_;
	Microtubule *adj_mt = neighbor->mt_;
	if(adj_mt != mt->neighbor_){
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
	}
	else{
		printf("error in get ext force (xlink)\n");
		exit(1);
	}
}

Tubulin* AssociatedProtein::GetActiveHeadSite(){

	if(heads_active_ == 1){
		if(site_one_ != nullptr)
			return site_one_;
		else if(site_two_ != nullptr)
			return site_two_;
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

Tubulin* AssociatedProtein::GetWeightedNeighborSite(){
	
	UpdateNeighborSites();
	double anch_coord = GetAnchorCoordinate(); 
	double p_tot = 0;
	for(int i_site = 0; i_site < n_neighbor_sites_; i_site++){
		Tubulin *site = neighbor_sites_[i_site]; 
		double site_coord = site->index_ + site->mt_->coord_; 
		int x_dist = abs(anch_coord - site_coord);
//		printf("x_dist is %i\n", x_dist);
		p_tot += binding_weight_lookup_[x_dist]; 
	}
	double ran = properties_->gsl.GetRanProb();
	double p_cum = 0;
	for(int i_site = 0; i_site < n_neighbor_sites_; i_site++){
		Tubulin *site = neighbor_sites_[i_site];
		double site_coord = site->index_ + site->mt_->coord_;
		int x_dist = abs(anch_coord - site_coord); 
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
	for(int i_site = 0; i_site < n_neighbor_sites_; i_site++){
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
