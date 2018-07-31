#include "master_header.h"
#include "microtubule.h"

Microtubule::Microtubule(){
}

void Microtubule::Initialize(system_parameters *parameters, 
	system_properties *properties, int i_mt){

	parameters_ = parameters;
	properties_ = properties;
	index_ = i_mt;

	SetParameters();
	GenerateLattice();
}

void Microtubule::SetParameters(){

	n_sites_ = parameters_->microtubules.length;
	coord_ = 0;
	if(index_%2 == 0){
		polarity_ = 0;
		plus_end_ = 0;
		minus_end_ = n_sites_ - 1;
		delta_x_ = -1;
		mt_index_adj_ = index_ + 1; 	// FIXME
		coord_ = parameters_->microtubules.start_coord[0];
	}
	else if(index_%2 == 1){
		polarity_ = 1;
		plus_end_ = n_sites_ - 1;
		minus_end_ = 0;
		delta_x_ = 1;
		mt_index_adj_ = index_ - 1;		// FIXME 
		coord_ = parameters_->microtubules.start_coord[1]; 
	}
	int mt_length = parameters_->microtubules.length;
	double site_size = parameters_->microtubules.site_size;
	double big_l = mt_length * site_size;
	double radius = parameters_->microtubules.radius;
	double height = parameters_->microtubules.elevation;
	double eta = parameters_->eta;
	// see radhika sliding paper for any of this to make sense
	double numerator = 2 * 3.14159 * big_l * (eta / 1000000);;
	double denom = log(2 * height / radius);
	gamma_ = (numerator / denom);
}

void Microtubule::GenerateLattice(){

	lattice_.resize(n_sites_);
	for(int i_site = 0; i_site < n_sites_; i_site++){
		lattice_[i_site].Initialize(parameters_, properties_, this, i_site);
	}
}

void Microtubule::UpdateExtensions(){

	KinesinManagement *kinesin4 = &properties_->kinesin4;
	AssociatedProteinManagement *prc1 = &properties_->prc1;
	int n_sites = n_sites_;
	// Run through all sites on MT
	for(int i_site = 0; i_site < n_sites; i_site++){
		Tubulin *site = &lattice_[i_site]; 
		// Check if site is occupied by a motor
		if(site->motor_ != nullptr){
			Kinesin *motor = site->motor_;
			// Update tether extension if motor is tethered
			if(motor->tethered_ == true){
				int x_dub_pre = motor->x_dist_doubled_;
				if(motor->xlink_->heads_active_ > 0){
					motor->UpdateExtension();
					// Make sure we didn't force an untether event
					if(motor->tethered_ == true){
						int x_dub_post = motor->x_dist_doubled_;
						if(x_dub_pre != x_dub_post){
							// Only update kinesin stats if double-bound
							if(motor->heads_active_ == 2){
								kinesin4->n_bound_ii_tethered_[x_dub_pre]--;
								kinesin4->n_bound_ii_tethered_[x_dub_post]++;
							}
							// Update site statistics for prc1
							AssociatedProtein *xlink = motor->xlink_;
							if(xlink->heads_active_ == 1){
								prc1->n_sites_i_tethered_[x_dub_pre]--;
								prc1->n_sites_i_tethered_[x_dub_post]++;
							}
							else if(xlink->heads_active_ == 2){
								int x_dist = xlink->x_dist_;
								prc1->n_sites_ii_tethered_[x_dub_pre][x_dist] 
									-= 2;
								prc1->n_sites_ii_tethered_[x_dub_post][x_dist]
									+= 2;
							}
						}
					}
				}
			}
		}
		// Check if site is occupied by an xlink
		if(site->xlink_ != nullptr){
			AssociatedProtein *xlink = site->xlink_;
			if(xlink->heads_active_ == 2){
				int x_pre = xlink->x_dist_;
				xlink->UpdateExtension();
				// Make sure we didn't force an unbind event
				if(xlink->heads_active_ == 2){
					int x_post = xlink->x_dist_;
					if(x_pre != x_post){
						prc1->n_double_bound_[x_pre]--;
						prc1->n_double_bound_[x_post]++;
						if(xlink->tethered_ == false){
							prc1->n_sites_ii_untethered_[x_pre] -= 2;
							prc1->n_sites_ii_untethered_[x_post] += 2;
						}
						else{
							Kinesin *motor = xlink->motor_;
							int x_dub_pre = motor->x_dist_doubled_;
							motor->UpdateExtension();
							// Make sure we didn't force an untether event
							if(motor->tethered_ == true){
								// Only bound motors contribute to teth stats
								if(motor->heads_active_ > 0){
									int x_dub_post = motor->x_dist_doubled_;
									if(x_dub_pre != x_dub_post){
										prc1->n_sites_ii_tethered_
											[x_dub_pre][x_pre] -= 2;
										prc1->n_sites_ii_tethered_
											[x_dub_post][x_post] += 2;
										// Update kinesin statistics
										if(motor->heads_active_ == 2){
											kinesin4->n_bound_ii_tethered_
												[x_dub_pre]--;
											kinesin4->n_bound_ii_tethered_
												[x_dub_post]++;
										}
									}
									else{
										prc1->n_sites_ii_tethered_
											[x_dub_pre][x_pre] -= 2;
										prc1->n_sites_ii_tethered_
											[x_dub_pre][x_post] += 2;
									}
								}
								// Xlinks attached to free motors diffuse 
								// as though they are untethered
								else{
									prc1->n_sites_ii_untethered_[x_pre] -= 2;
									prc1->n_sites_ii_untethered_[x_post] += 2;
								}
							}
							// If untether event was forced, correct stats
							else{
								prc1->n_sites_ii_tethered_
									[x_dub_pre][x_pre] -= 2;	
								prc1->n_sites_ii_tethered_
									[x_dub_pre][x_post] += 2;
							}
						}
					}
				}
			}
			// If xlink is single bound, check if motor is on another MT
			else if(xlink->tethered_ == true){
				Kinesin *motor = xlink->motor_;
				if(motor->heads_active_ > 0
				&& motor->mt_ != this){
					int x_dub_pre = motor->x_dist_doubled_;
					motor->UpdateExtension();
					// Make sure we didn't force an unbinding event
					if(motor->tethered_ == true){
						int x_dub_post = motor->x_dist_doubled_;
						if(x_dub_pre != x_dub_post){
							prc1->n_sites_i_tethered_[x_dub_pre]--;
							prc1->n_sites_i_tethered_[x_dub_post]++;
							// Update kinesin statistics
							if(motor->heads_active_ == 2){
								kinesin4->n_bound_ii_tethered_[x_dub_pre]--;
								kinesin4->n_bound_ii_tethered_[x_dub_post]++;
							}
						}
					}
				}
			}
		}
	}
}

double Microtubule::GetNetForce(){

	double forces_summed = 0;
	for(int i_site = 0; i_site < n_sites_; i_site++){
		Tubulin *site = &lattice_[i_site];
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
							forces_summed += motor->GetTetherForce(site);
						}
						// Only count force from 1st foot
						else if(site == motor->front_site_){
							forces_summed += motor->GetTetherForce(site);
						}
					}
				}
				// If xlink is double bound, add force regardless
				else if(xlink->heads_active_ == 2){
					if(motor->heads_active_ == 1){
						forces_summed += motor->GetTetherForce(site);
					}
					// Only count force from 1st foot (no double counting)
					else if(site == motor->front_site_){
						forces_summed += motor->GetTetherForce(site);
					}
				}
			}
		}
		// Otherwise, check if site is occupied by an xlink
		else if(site->xlink_ != nullptr){
			AssociatedProtein *xlink = site->xlink_;
			// Xlinks can only exert forces if they're double bound
			if(xlink->heads_active_ == 2){
				forces_summed += xlink->GetExtensionForce(site);
			}
			// Motors tethered to this xlink can also exert forces
			if(xlink->tethered_ == true){
				Kinesin *motor = xlink->motor_;
				// To avoid double counting, make sure xlink's motor
				// is on a different microtubule when double bound
				if(site->mt_ != motor->mt_ 
				&& motor->heads_active_ > 0){
					forces_summed += motor->GetTetherForce(site);
				}
			}
		}
	}
	return forces_summed; 
}

double Microtubule::GetNetForce_Motors(){

	double forces_summed = 0;
	for(int i_site = 0; i_site < n_sites_; i_site++){
		Tubulin *site = &lattice_[i_site];
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
							forces_summed += motor->GetTetherForce(site);
						}
						// Only count force from 1st foot (no double counting)
						else if(site == motor->front_site_){
							forces_summed += motor->GetTetherForce(site);
						}
					}
				}
				// If xlink is double bound, add force regardless
				else if(xlink->heads_active_ == 2){
					if(motor->heads_active_ == 1){
						forces_summed += motor->GetTetherForce(site);
					}
					// Only count force from 1st foot (no double counting)
					else if(site == motor->front_site_){
						forces_summed += motor->GetTetherForce(site);
					}
				}
			}
		}
		// Otherwise, check if site is occupied by an xlink
		else if(site->xlink_ != nullptr){
			AssociatedProtein *xlink = site->xlink_;
			// Motors tethered to this xlink can also exert forces
			if(xlink->tethered_ == true){
				Kinesin *motor = xlink->motor_;
				// To avoid double counting, make sure xlink's motor
				// is on a different microtubule when double bound
				if(site->mt_ != motor->mt_ 
				&& motor->heads_active_ > 0){
					forces_summed += motor->GetTetherForce(site);
				}
			}
		}
	}
	return forces_summed; 
}

double Microtubule::GetNetForce_Xlinks(){

	double forces_summed = 0;
	for(int i_site = 0; i_site < n_sites_; i_site++){
		Tubulin *site = &lattice_[i_site];
		// Check if site is occupied by an xlink
		if(site->xlink_ != nullptr){
			AssociatedProtein *xlink = site->xlink_;
			// Xlinks can only exert forces if they're double bound
			if(xlink->heads_active_ == 2){
				forces_summed += xlink->GetExtensionForce(site);
			}
		}
	}
	return forces_summed; 
}
