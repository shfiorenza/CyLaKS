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
//	KinesinManagement *kinesin4 = &properties_->kinesin4; 
//	int n_binding_affinities = kinesin4->motor_list_[0].n_binding_affinities_;
	// XXX BOUNDARY SITES ACCESSIBLE -- DISABLE FOR ALPHA/BETA XXX
	// int n_sites_bulk = n_sites_tot_ - 2*n_mts;	
	unoccupied_list_.resize(n_sites_tot_); 	
//	for(int i_aff = 0; i_aff < n_binding_affinities; i_aff++){
//		unoccupied_list_[i_aff].resize(n_sites_tot_);
//	}
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

	n_unoccupied_ = 0;
	int n_occupied = 0;
	int i_unoccupied = 0;
	int n_mts = parameters_->microtubules.count;
	int mt_length = parameters_->microtubules.length;
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

	int n_mts = parameters_->microtubules.count;
	/*	Sum up all forces exerted on the MTs 	*/
	double forces_summed[n_mts];
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		forces_summed[i_mt] = 0;
		Microtubule *mt = &mt_list_[i_mt]; 
		forces_summed[i_mt] = mt->GetNetForce();
	}
	// Check for symmetry
	double tolerance = 0.000001; 
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
	/*	Calculate MT displacements for this timestep */
	int displacement[n_mts];
	double site_size = parameters_->microtubules.site_size;
	double gamma = mt_list_[0].gamma_;
	double kbT = parameters_->kbT;
	double delta_t = parameters_->delta_t;
	// standard deviation of the gaussian to be sampled below
	double sigma = sqrt(2 * kbT * delta_t / gamma);
	// Imposed velocity on top MT (being pulled against dir. it slides)
	std::vector<double> imp_vel = parameters_->microtubules.imposed_velocity;
	std::vector<double> immobile_until = 
		parameters_->microtubules.immobile_until;
	int current_step = properties_->current_step_; 
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		long unpin_step = immobile_until[i_mt] / delta_t; 
		int sites_per_sec = imp_vel[i_mt] / site_size; 
		double sites_per_timestep = sites_per_sec * delta_t; 
		int timesteps_per_site;
		if(imp_vel[i_mt] != 0)
			timesteps_per_site = 1 / sites_per_timestep;
		else
			timesteps_per_site = 2 * parameters_->n_steps; //locally infinite
		if(current_step >= unpin_step){
			double velocity = forces_summed[i_mt] / gamma; 
			// gaussian noise is added into the calculated displacement
			double noise = properties_->gsl.GetGaussianNoise(sigma);
			double raw_displacement = velocity * delta_t + noise;
//			printf("dx: %g (%g noise) sites for mt #%i\n", 
//					raw_displacement / site_size, noise / site_size, i_mt);
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
			if(current_step % timesteps_per_site == 0
			&& imp_vel[i_mt] > 0)
				n_steps++;
			else if(current_step % timesteps_per_site == 0
			&& imp_vel[i_mt] < 0)
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
