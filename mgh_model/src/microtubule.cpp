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

	n_sites_ = parameters_->length_of_microtubule;
	kbT_ = parameters_->kbT;
	radius_ = parameters_->mt_radius;
	height_ = parameters_->mt_height;
	eta_inverse_ = parameters_->eta_inverse;
	site_size_ = parameters_->site_size;
	coord_ = 0;
	if(index_%2 == 0){
		polarity_ = 0;
		plus_end_ = n_sites_ - 1;
		minus_end_ = 0;
		delta_x_ = 1;
		mt_index_adj_ = index_ + 1; 	// FIXME
	}
	else if(index_%2 == 1){
		polarity_ = 1;
		plus_end_ = 0;
		minus_end_ = n_sites_ - 1;
		delta_x_ = -1;
		mt_index_adj_ = index_ - 1;		// FIXME 
	}
	int mt_length = parameters_->length_of_microtubule;
	big_l_ = mt_length * site_size_;
	// see radhika paper for any of this to make sense
	double delta_t = parameters_->delta_t;
	double numerator = 2 * 3.14159 * big_l_ / (eta_inverse_ * 1000000);;
	double denom = log(2 * height_ / radius_);
	// XXX divide by delta_t to make units work??
	gamma_ = (numerator / denom);
	printf("gamma: %g = (%g / %g) / %g\n", gamma_, numerator, denom, delta_t);
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
	int n_sites = parameters_->length_of_microtubule;
	for(int i_site = 0; i_site < n_sites; i_site++){
		Tubulin *site = &lattice_[i_site]; 
		if(site->motor_ != nullptr){
			Kinesin *motor = site->motor_;
			if(motor->tethered_ == true){
				int x_dub_pre = motor->x_dist_doubled_;
				motor->UpdateExtension();
				// Make sure we didn't force an untether event
				if(motor->tethered_ == true){
					int x_dub_post = motor->x_dist_doubled_;
					if(x_dub_pre != x_dub_post){
						// Only update kinesin stats if double-bound
						if(motor->heads_active_ == 2){
							kinesin4->n_bound_tethered_[x_dub_pre]--;
							kinesin4->n_bound_tethered_[x_dub_post]++;
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
		else if(site->xlink_ != nullptr){
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
											kinesin4->n_bound_tethered_
												[x_dub_pre]--;
											kinesin4->n_bound_tethered_
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
								kinesin4->n_bound_tethered_[x_dub_pre]--;
								kinesin4->n_bound_tethered_[x_dub_post]++;
							}
						}
					}
				}
			}
		}
	}
}
