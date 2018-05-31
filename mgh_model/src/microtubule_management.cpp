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
/*	Sum up all forces exerted on the MTs 	*/
	double forces_summed[n_mts];
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		Microtubule *mt = &mt_list_[i_mt]; 
		forces_summed[i_mt] = mt->GetNetForce();
	}
	// Check for symmetry
	double tolerance = 0.000001; 
	for(int i_mt = 0; i_mt < n_mts; i_mt += 2){
//		if(forces_summed[i_mt] != 0){
//		printf("for mt # %i: %g\n", i_mt, forces_summed[i_mt]);
//		printf("for mt # %i: %g\n", i_mt + 1, forces_summed[i_mt + 1]);
//		properties_->wallace.PrintMicrotubules(0.1);
//		}
		double delta = abs(forces_summed[i_mt] + forces_summed[i_mt + 1]);
		if(delta > tolerance){
				printf("aw man in MT diffusion\n");
				printf("for mt %i: %g, for mt %i: %g\n", 
						i_mt, forces_summed[i_mt], 
						i_mt + 1, - forces_summed[i_mt + 1]);
				exit(1);
		}
	}
/*	Calculate MT displacements for this timestep */
	double kbT = mt_list_[0].kbT_;
	double site_size = mt_list_[0].site_size_;
	double delta_t = parameters_->delta_t;
	double gamma = mt_list_[0].gamma_;
	// variance of the gaussian to be sampled below
	double sigma = sqrt(2 * kbT * delta_t / gamma);
	//	Imposed velocity on top MT (being pulled against dir. it slides)
	double imposed_vel = parameters_->top_mt_imposed_velocity;
	// step in KMC sim at which top MT will be free to slide
	int step_unpin = parameters_->top_mt_pinned_until / delta_t; 
	int cur_step = properties_->current_step_;
	int displacement[n_mts];
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		double velocity = forces_summed[i_mt] / gamma;
		// gaussian noise is added into the calculated displacement
		double noise = properties_->gsl.GetGaussianNoise(sigma);
		double raw_displacement = velocity * delta_t + noise;
//		printf("dx: %g (%g noise) sites for mt #%i\n", 
//				raw_displacement / site_size, noise / site_size, i_mt);
		double site_displacement = (raw_displacement) / site_size;
		// Get number of sites MT is expected to move
		int n_steps = (int) site_displacement;
		// Use leftover as a probability to roll for another step
		double leftover = fabs(site_displacement - n_steps);
		double ran = properties_->gsl.GetRanProb();
		if(ran < leftover
		&& site_displacement > 0)
			n_steps++;
		else if(ran < leftover
		&& site_displacement < 0)
			n_steps--;
		// Store value in array
		displacement[i_mt] = n_steps;
		// Add imposed velocity to top MT only
		if(i_mt == 1){
			if(imposed_vel != 0
			&& cur_step >= step_unpin){
				double imposed_disp = delta_t * imposed_vel / site_size; 
				int step_dir = mt_list_[i_mt].delta_x_; 
				double ran = properties_->gsl.GetRanProb();
				if(ran < imposed_disp)
					displacement[i_mt] += step_dir;
			}
			if(cur_step < step_unpin){
				displacement[i_mt] = 0;
			}
		}
	}
/*  Run through MT list and update displacementsi */
	int i_mt_start;
	if(parameters_->bot_mt_pinned == true){
		i_mt_start = 1;
	}	
	else{
		i_mt_start = 0;
	}
	for(int i_mt = i_mt_start; i_mt < n_mts; i_mt++){ 
		Microtubule *mt = &mt_list_[i_mt];
		Microtubule *neighb_mt; 
		if(n_mts > 1){
			int i_neighb_mt = 1 - i_mt; 
			neighb_mt = &mt_list_[i_neighb_mt];
		}
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
			neighb_mt->UpdateExtensions();
		}
	}
}

