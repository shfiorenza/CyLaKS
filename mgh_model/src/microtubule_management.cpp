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

	int n_mts = parameters_->microtubules.count;
	int mt_length = parameters_->microtubules.length;
	n_sites_tot_ = n_mts*mt_length;
	// int n_sites_bulk = n_sites_tot_ - 2*n_mts;	
	unoccupied_list_.resize(n_sites_tot_); 
}

void MicrotubuleManagement::GenerateMicrotubules(){

	int n_mts = parameters_->microtubules.count;
	mt_list_.resize(n_mts);
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		mt_list_[i_mt].Initialize(parameters_, properties_, i_mt);
	}
}

void MicrotubuleManagement::UnoccupiedCheck(Tubulin *site){

	if(site->motor_ != nullptr || site->xlink_ != nullptr){
		printf("Error @ site %i_%i: should be unoccupied\n", 
				site->mt_->index_, site->index_);
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

	int n_mts = parameters_->microtubules.count; 
	for(int i_mt = 0; i_mt < n_mts; i_mt+=2){
		Microtubule *mt = &mt_list_[i_mt];
		Microtubule *mt_adj = &mt_list_[i_mt+1];
		mt->neighbor_ = mt_adj;
		mt_adj->neighbor_ = mt;
	}

}

void MicrotubuleManagement::UpdateUnoccupiedList(){

	int n_mts = parameters_->microtubules.count;
	int mt_length = parameters_->microtubules.length;
	n_unoccupied_ = 0; 
	int i_unoccupied = 0; 
	int n_occupied = 0;
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		Microtubule *mt = &mt_list_[i_mt];
		// XXX BOUNDARY SITES INCLUDED - DISABLE FOR ALPHA/BETA !! XXX 
		for(int i_site = 0; i_site <= mt_length - 1; i_site++){
			Tubulin *site = &mt->lattice_[i_site];
			if(site->occupied_ == false){
				unoccupied_list_[i_unoccupied] = site; 
				i_unoccupied++;
				n_unoccupied_++;
			}
			else{
				n_occupied++; 
			}
		}
	}
	if(n_unoccupied_ + n_occupied != n_sites_tot_){
		printf("something awful in update_unoccupied_list bruh: \n");
		printf("%i != %i + %i\n", n_sites_tot_, 
				n_unoccupied_, n_occupied);
		exit(1);
	}
}

Tubulin* MicrotubuleManagement::GetUnoccupiedSite(){

	UpdateUnoccupiedList();
	int n_unoccupied = n_unoccupied_;
	// Make sure an unoccupied site exists
	if(n_unoccupied > 0){
		int i_entry = properties_->gsl.GetRanInt(n_unoccupied);	
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


	int n_mts = parameters_->microtubules.count;
	int current_step = properties_->current_step_; 
	double delta_t = parameters_->delta_t;
	double site_size = parameters_->microtubules.site_size;
	double kbT = parameters_->kbT;
	double gamma = mt_list_[0].gamma_;
	long unpin_step[n_mts]; 
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		double unpin_time = parameters_->microtubules.immobile_until[i_mt]; 
		unpin_step[i_mt] = unpin_time / delta_t; 
	}
	/* To ensure smooth diffusion without having to use timesteps of 
	  ~1e-7 seconds, divide all forces by N and run loop N times */
	int n_iterations = 100;
	double delta_t_eff = delta_t / n_iterations; 
	double sigma = sqrt(2 * kbT * delta_t_eff / gamma); 
	for(int i_itr = 0; i_itr < n_iterations; i_itr++){
		/*	Sum up all forces exerted on the MTs 	*/
		double forces_summed[n_mts];
		for(int i_mt = 0; i_mt < n_mts; i_mt++){
			Microtubule *mt = &mt_list_[i_mt]; 
			if(i_itr == 0){		// For checking symmetry
				mt->UpdateExtensions(); 
				forces_summed[i_mt] = mt->GetNetForce();
			}
			else if(current_step >= unpin_step[i_mt]){
				mt->UpdateExtensions();
				forces_summed[i_mt] = mt->GetNetForce();
			}
		}
		// Check for symmetry on first iteration only
		if(i_itr == 0){
			double tolerance = 0.000001; 
			for(int i_mt = 0; i_mt < n_mts; i_mt += 2){
				double delta = 
					fabs(forces_summed[i_mt] + forces_summed[i_mt+1]);
				if(delta > tolerance){
					printf("aw man in MT diffusion\n");
					printf("for mt %i: %g, for mt %i: %g\n", 
							i_mt, forces_summed[i_mt], 
							i_mt + 1, - forces_summed[i_mt + 1]);
					exit(1);
				}
			}
		}
		/*	Calculate MT displacements for this timestep */
		int displacement[n_mts];
		for(int i_mt = 0; i_mt < n_mts; i_mt++){
			if(current_step >= unpin_step[i_mt]){
				double velocity = forces_summed[i_mt] / gamma; 
				// gaussian noise is added into the calculated displacement
				double noise = properties_->gsl.GetGaussianNoise(sigma); 
				double raw_displacement = velocity*delta_t_eff + noise;
//				printf("dx: %g (%g noise) sites for mt #%i\n", 
//					raw_displacement/site_size, noise_eff/site_size, i_mt);
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
			}
			else{
				displacement[i_mt] = 0;	
			}
		}
		/*  Run through MT list and update displacementsi */
		for(int i_mt = 0; i_mt < n_mts; i_mt++){ 
			Microtubule *mt = &mt_list_[i_mt];
			Microtubule *neighb = &mt_list_[mt->mt_index_adj_];
			int n_steps = displacement[i_mt];
			int dx = 0;
			if(n_steps > 0){
				dx = 1;
			}
			else if(n_steps < 0){
				dx = -1;
				n_steps = abs(n_steps);
			}
			for(int i_step = 0; i_step < n_steps; i_step++){
				mt->coord_ += dx; 
				mt->UpdateExtensions();
				neighb->UpdateExtensions();
			}
		}
	}
}
