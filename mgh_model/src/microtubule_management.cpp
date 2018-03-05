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
	UpdateNeighbors();
	UpdateUnoccupiedList();
}

void MicrotubuleManagement::SetParameters(){

	int n_mts = parameters_->n_microtubules;
	int mt_length = parameters_->length_of_microtubule;
	n_sites_tot_ = n_mts*mt_length;
	// XXX BOUNDARY SITES ACCESSIBLE -- DISABLE FOR ALPHA/BETA XXX
	int n_sites_bulk = n_sites_tot_; // - 2*n_mts;
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

void MicrotubuleManagement::UpdateNeighbors(){

	int n_mts = parameters_->n_microtubules;
	for(int i_mt = 0; i_mt < n_mts; i_mt+=2){
		Microtubule *mt = &mt_list_[i_mt];
		Microtubule *mt_adj = &mt_list_[i_mt+1];
		mt->neighbor_ = mt_adj;
		mt_adj->neighbor_ = mt;
	}

}

void MicrotubuleManagement::UpdateNumUnoccupied(){
	
	n_unoccupied_ = 0;
	int n_occupied = 0;
	int n_mts = parameters_->n_microtubules;
	int mt_length = parameters_->length_of_microtubule;
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		Microtubule *mt = &mt_list_[i_mt];
		// XXX BOUNDARY SITES INCLUDED - DISABLE FOR ALPHA/BETA !! XXX 
		// Exclude boundary sites
		for(int i_site = 0; i_site <= mt_length - 1; i_site++){
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
	int n_sites_bulk = n_sites_tot_; // - 2*n_mts;
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
		// XXX BOUNDARY SITES INCLUDED - DISABLE FOR ALPHA/BETA !! XXX 
		// Exclude boundary sites
		for(int i_site = 0; i_site <= mt_length - 1; i_site++){
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
		printf("Error: GetUnoccupiedSite called, but no unoccupied sites\n");
		exit(1);
	}
}

void MicrotubuleManagement::RunDiffusion(){

	int n_mts = parameters_->n_microtubules;
	int n_sites = parameters_->length_of_microtubule;
	double forces_summed[n_mts];
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		forces_summed[i_mt] = 0;
		for(int i_site = 0; i_site < n_sites; i_site++){
			Tubulin *site = &mt_list_[i_mt].lattice_[i_site];
			// Check if site is occupied by a motor
			if(site->motor_ != nullptr){
				Kinesin *motor = site->motor_;
				// Motors can only exert forces if they're tethered
				if(motor->tethered_ == true){
					AssociatedProtein *xlink = motor->xlink_;
					// If xlink is single bound, only add force if
					// it's on a different microtubule
					if(xlink->heads_active_ == 1){
						Tubulin *xlink_site = xlink->GetActiveHeadSite();
						if(site->mt_ != xlink_site->mt_){
							if(motor->heads_active_ == 1){
								forces_summed[i_mt] += 
									motor->GetTetherForce(site);
							}
							else if(site == motor->front_site_){
								forces_summed[i_mt] += 
									motor->GetTetherForce(site); 
//								printf("%g from motor\n", motor->
//										GetTetherForce(site));
							}
						}
					}
					// If xlink is double bound, add force
					else if(xlink->heads_active_ == 2){
						if(motor->heads_active_ == 1){
							forces_summed[i_mt] += 
										motor->GetTetherForce(site);
						}
						// Only count force from 1st foot (no double counting)
						else if(site == motor->front_site_){
							forces_summed[i_mt] += 
										motor->GetTetherForce(site); 
						}
					}
				}
			}
			// Otherwise, check if site is occupied by an xlink
			else if(site->xlink_ != nullptr){
				AssociatedProtein *xlink = site->xlink_;
				// Xlinks can only exert forces if they're double bound
				if(xlink->heads_active_ == 2){
					forces_summed[i_mt] += xlink->GetExtensionForce(site);
				}
				// Motors tethered to this xlink can also exert forces
				if(xlink->tethered_ == true){
					Kinesin *motor = xlink->motor_;
					// To avoid double counting, make sure xlink's motor
					// is on a different microtubule when double bound
					if(site->mt_ != motor->mt_
					&& motor->heads_active_ > 0){
						forces_summed[i_mt] += 
							motor->GetTetherForce(site);
					}
				}
			}
		}
	}
	// Check for symmetry
	double tolerance = 0.0001; 
	for(int i_mt = 0; i_mt < n_mts; i_mt += 2){
		double delta = abs(forces_summed[i_mt] + forces_summed[i_mt + 1]);
		if(delta > tolerance){
				printf("aw man in MT diffusion\n");
				printf("for mt %i: %g, for mt %i: %g\n", 
						i_mt, forces_summed[i_mt], 
						i_mt + 1, - forces_summed[i_mt + 1]);
				exit(1);
		}
	}
	// Calculate MT displacements for this timestep 
	double kbT = mt_list_[0].kbT_;
	double site_size = mt_list_[0].site_size_;
	double delta_t = parameters_->delta_t;
	double gamma = mt_list_[0].gamma_;
	// variance of the gaussian to be sampled below
	double sigma = sqrt(2 * kbT * delta_t / gamma);
	double displacement[n_mts];
//	printf("amp: %g\n", amplitude);
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		double velocity = forces_summed[i_mt] / gamma;
		// gaussian noise is added into the calculated displacement
		double noise = properties_->gsl.GetGaussianNoise(sigma);
		double raw_displacement = velocity * delta_t + noise;
//		printf("dx: %g (%g noise) nm / s\n", raw_displacement, 
//				noise);
//		printf("gamma be %g\n", gamma);
		double site_displacement = (raw_displacement) / site_size;
//		printf("%g sites\n", site_displacement);
		// Get number of sites MT is expected to move
		int n_sites = (int) site_displacement;
		// Use leftover as a probability to roll for another step
		double leftover = abs(site_displacement - n_sites);
		double ran = properties_->gsl.GetRanProb();
		if(ran < leftover
		&& n_sites > 0)
			n_sites++;
		else if(ran < leftover
		&& n_sites < 0)
			n_sites--;
		else if(ran < leftover
		&& n_sites == 0
		&& site_displacement > 0)
			n_sites++;
		else if(ran < leftover
		&& n_sites == 0
		&& site_displacement < 0)
			n_sites--;
		// Store value in array
		displacement[i_mt] = n_sites;
//		if(n_sites > 0)
//			printf("MT diffusion for #%i: %i\n", i_mt, n_sites);
	}
	// Run through MT list and update displacements
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		Microtubule *mt = &mt_list_[i_mt];
		int n_steps = displacement[i_mt];
		int dx = 0;
		if(n_steps > 0)
			dx = 1;
		else if(n_steps < 0){
			dx = -1;
			n_steps = abs(n_steps);
		}
		for(int i_step = 0; i_step < n_steps; i_step++){
			mt->coord_ += dx; 
			mt->UpdateExtensions();
		}
	}
}

