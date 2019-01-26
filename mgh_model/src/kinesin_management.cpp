#include "master_header.h"
#include "kinesin_management.h"

KinesinManagement::KinesinManagement(){
}

void KinesinManagement::Initialize(system_parameters *parameters, 
        system_properties *properties){

    parameters_ = parameters;
    properties_ = properties;

    GenerateMotors();	
    SetParameters();	
	InitializeLists();
	InitializeEvents();
}

void KinesinManagement::GenerateMotors(){

    int n_mts = parameters_->microtubules.count;
    int n_sites = parameters_->microtubules.length;
    // Since only one head has to be bound, the most that will ever
    // be needed (all single-bound) is the total number of sites 
    n_motors_ = n_mts*n_sites;
    motors_.resize(n_motors_);
    for(int ID = 0; ID < n_motors_; ID++){
        motors_[ID].Initialize(parameters_, properties_, ID);
    }
}

void KinesinManagement::SetParameters(){

	// Simulation constants
    double delta_t = parameters_->delta_t;
	double site_size = parameters_->microtubules.site_size;
	double kbT = parameters_->kbT;
	double r_0 = motors_[0].r_0_;
	double k_spring = motors_[0].k_spring_;
	double k_slack = motors_[0].k_slack_;
	double r_y = parameters_->microtubules.y_dist / 2;

	// Non-tethered statistics
    double k_on = parameters_->motors.k_on;
    double c_motor = parameters_->motors.c_bulk;
	p_bind_i_ = k_on * c_motor * delta_t;
	double k_on_ATP = parameters_->motors.k_on_ATP;
	double c_ATP = parameters_->motors.c_ATP;
	p_bind_ATP_ = k_on_ATP * c_ATP * delta_t;
	double k_hydrolyze = parameters_->motors.k_hydrolyze;
	p_hydrolyze_ = k_hydrolyze * delta_t;
	double c_eff = parameters_->motors.c_eff_bind;
	p_bind_ii_ = k_on * c_eff * delta_t; 
	double k_off_ii = parameters_->motors.k_off_ii;
	p_unbind_ii_ = k_off_ii * delta_t; 
	double k_off_i = parameters_->motors.k_off_i;
	p_unbind_i_ = k_off_i * delta_t; 
	// Sound the alarm if our timestep is too large
	if(p_bind_i_ > 1) 
		printf("WARNING: p_bind_i=%g for motors\n", p_bind_i_);
	if(p_bind_ATP_ > 1)
		printf("WARNING: p_bind_ATP=%g for motors\n", p_bind_ATP_);	
	if(p_hydrolyze_ > 1)
		printf("WARNING: p_phos=%g for motors\n", p_hydrolyze_);	
	if(p_bind_ii_ > 1)
		printf("WARNING: p_bind_ii=%g for motors\n", p_bind_ii_);
	if(p_unbind_ii_ > 1)
		printf("WARNING: p_unbind_ii=%g for motors\n", p_unbind_ii_);
	if(p_unbind_i_ > 1)
		printf("WARNING: p_unbind_i=%g for motors\n", p_unbind_i_);

	// Get tether information from motors
	rest_dist_ = motors_[0].rest_dist_;
	comp_cutoff_ = motors_[0].comp_cutoff_;
	dist_cutoff_ = motors_[0].dist_cutoff_;
	if(parameters_->motors.tethers_active){
		printf("For motors:\n");
		printf("  rest_dist is %g\n", rest_dist_);
		printf("  comp_cutoff is %i\n", comp_cutoff_);
		printf("  dist_cutoff is %i\n", dist_cutoff_);
	}
	// Tethered statistics
	double k_tether = parameters_->motors.k_tether;
	double k_untether= parameters_->motors.k_untether;
	double c_eff_teth = parameters_->motors.c_eff_tether;
	if(!parameters_->motors.tethers_active){
		k_tether = 0;
		k_untether = 0;
		c_eff_teth = 0;
	}
	p_tether_free_ = k_tether * c_motor * delta_t;
	p_untether_free_ = k_untether * delta_t;
	// Sound the alarm if our timestep is too large
	if(p_tether_free_ > 1) 
		printf("WARNING: p_teth_free=%g for motors\n", p_tether_free_);
	if(p_untether_free_ > 1)
		printf("WARNING: p_unteth_free=%g for motors\n", p_untether_free_);	
	// For events that depend on tether stretch, each different extension 
	// has its own rate; "base" refers to when the tether is unstretched 
	double p_bind_base = k_on * c_eff_teth * delta_t; 
	double p_unbind_base = p_unbind_i_;
	double p_teth_base = k_tether * c_eff_teth * delta_t;
	double p_unteth_base = p_untether_free_;
	p_bind_i_tethered_.resize(2*dist_cutoff_ + 1);
	p_unbind_i_tethered_.resize(2*dist_cutoff_ + 1);
	p_tether_bound_.resize(2*dist_cutoff_ + 1);
	p_untether_bound_.resize(2*dist_cutoff_ + 1);
	for(int x_dub = 0; x_dub <= 2*dist_cutoff_; x_dub++){
		// Calculate tether length for this x_dist 
		double r_x = (double)x_dub * site_size / 2;
		double r = sqrt(r_y*r_y + r_x*r_x);
		// Calculate extension of tether for given x_dub
		double dr = r - r_0; 
		// Calculate potential energy of this tether extension
		double U_teth;
		if(dr > 0)
			U_teth = k_spring * dr * dr / 2;
		else
			U_teth = k_spring * dr * dr / 2;
		// Weight of creating this configuration
		double weight_create = exp(-U_teth/(2*kbT));
		// Weight of annihilating this configuration
		double weight_annihilate = exp(U_teth/(2*kbT));
		// If tethering is disabled, all weights are automatically zero
		if(!parameters_->motors.tethers_active){
			weight_create = 0;
			weight_annihilate = 0; 
		}
		// Calculate appropriately-weighted probabilities
		p_bind_i_tethered_[x_dub] = weight_create * p_bind_base;
		p_unbind_i_tethered_[x_dub] = weight_annihilate * p_unbind_base;
		p_tether_bound_[x_dub] = weight_create * p_teth_base;
		p_untether_bound_[x_dub] = weight_annihilate * p_unteth_base;
		// Sound the alarm if our timestep is too large
		if(p_bind_i_tethered_[x_dub] > 1)
			printf("WARNING: p_bind_i_teth=%g for 2x=%i\n", 
					p_bind_i_tethered_[x_dub], x_dub);
		if(p_unbind_i_tethered_[x_dub] > 1)
			printf("WARNING: p_unbind_i_teth=%g for 2x=%i\n", 
					p_unbind_i_tethered_[x_dub], x_dub);
		if(p_tether_bound_[x_dub] > 1)
			printf("WARNING: p_teth_bound=%g for 2x=%i\n", 
					p_tether_bound_[x_dub], x_dub);
		if(p_untether_bound_[x_dub] > 1)
			printf("WARNING: p_unteth_bound=%g for 2x=%i\n", 
					p_untether_bound_[x_dub], x_dub);
	}
}

void KinesinManagement::InitializeLists(){

	// One dimensional stuff
	active_.resize(n_motors_);
	docked_.resize(n_motors_); 
	bound_NULL_.resize(n_motors_);
	bound_ATP_.resize(n_motors_);
	bound_ADPP_i_.resize(n_motors_);
	bound_ADPP_ii_.resize(n_motors_);
	free_tethered_.resize(n_motors_);
	bound_untethered_.resize(n_motors_);
	// Two dimensional stuff
	n_bound_tethered_.resize(2*dist_cutoff_ + 1);
	bound_tethered_.resize(2*dist_cutoff_ + 1);
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		bound_tethered_[x_dist_dub].resize(n_motors_); 
		n_bound_tethered_[x_dist_dub] = 0;
	}
}

void KinesinManagement::InitializeEvents(){

	int tot_size = 6;
/*
	if(parameters_->motors.tethers_active)
		tot_size = 16 + 16*dist_cutoff_ + 1;
*/
	serial_events_.resize(tot_size);

	// create modular binomial function to use in event data structures
	auto binomial = [&](double p, int n, int ID){
		if(n > 0){
			return properties_->gsl.SampleBinomialDist_Kinesin(p, n, ID);
		}
		else return 0;
	};

	/*
	event bind_i = {"bind_i", "unocc", 0, binomial_ptr, 
		&properties_->microtubules.n_unoccupied_, &p_bind_i_, 
			0, 
			std::bind(&KinesinManagement::KMC_Bind_I, this, 0)};
	serial_events_[bind_i.index_] = bind_i;

	event bind_ATP = {"bind_ATP", "bound_NULL", &n_bound_NULL_, 
			20, 1, 0, binomial, p_bind_ATP_, 0};
	serial_events_[bind_ATP.index_] = bind_ATP;

	event hydrolyze = {"hydrolyze", "bound_ATP", &n_bound_ATP_,
			30, 2, 0, binomial, p_hydrolyze_, 0};
	serial_events_[hydrolyze.index_] = hydrolyze;

	event bind_ii = {"bind_ii", "docked", &n_docked_, 
			40, 3, 0, binomial, p_bind_ii_, 0};
	serial_events_[bind_ii.index_] = bind_ii;

	event unbind_ii = {"unbind_ii", "bound_ADPP_ii", &n_bound_ADPP_ii_, 
			50, 4, 0, binomial, p_unbind_ii_, 0};
	serial_events_[unbind_ii.index_] = unbind_ii;

	event unbind_i = {"unbind_i", "bound_ADPP_i", &n_bound_ADPP_i_, 
			60, 5, 0, binomial, p_unbind_i_, 0};
	serial_events_[unbind_i.index_] = unbind_i;

	*/
	/*
	// If tethers ARE active, serialize all extension-based events
	if(parameters_->motors.tethers_active){
		event bind_i_tethered = {&n_free_tethered_, 0, 
			"bind_i_tethered", 1, 0};
		serial_events_[1] = bind_i_tethered;

		event bind_ii = {&n_bound_i_bindable_, 0, "bind_ii", 2, 0};
		serial_events_[2] = bind_ii;

		for(int x_dub(0); x_dub <= 2*dist_cutoff_; x_dub++){
			event bind_ii_to_teth = {&n_bound_i_bindable_to_teth_[x_dub], 
				x_dub, "bind_ii_to_teth", 3 + x_dub, 0};
			serial_events_[3 + x_dub] = bind_ii_to_teth;

			event bind_ii_fr_teth = {&n_bound_i_bindable_fr_teth_[x_dub], 
				x_dub, "bind_ii_fr_teth", 4 + 2*dist_cutoff_ + x_dub, 0};
			serial_events_[4 + 2*dist_cutoff_ + x_dub] = bind_ii_fr_teth; 
		}
		int offset = 4 + 4*dist_cutoff_;

		event unbind_i = {&n_bound_i_, 0, "unbind_i", offset + 1, 0};
		serial_events_[offset + 1] = unbind_i;

		for(int x_dub(0); x_dub <= 2*dist_cutoff_; x_dub++){
			event unbind_i_tethered = {&n_bound_i_tethered_[x_dub], 
				x_dub, "unbind_i_tethered", offset + 2 + x_dub, 0};
			serial_events_[offset + 2 + x_dub] = unbind_i_tethered;
		}
		int offset2 = 6 + 6*dist_cutoff_; 

		event unbind_ii = {&n_bound_ii_, 0, "unbind_ii", offset2 + 1, 0};
		serial_events_[offset2 + 1] = unbind_ii;

		for(int x_dub(0); x_dub <= 2*dist_cutoff_; x_dub++){
			event unbind_to_teth = {&n_bound_ii_tethered_[x_dub],
			   	x_dub, "unbind_ii_to_teth", offset2 + 2 + x_dub, 0};
			serial_events_[offset2 + 2 + x_dub] = unbind_to_teth;

			event unbind_fr_teth = {&n_bound_ii_tethered_[x_dub], 
				x_dub, "unbind_ii_fr_teth", 
				offset2 + 3 + 2*dist_cutoff_ + x_dub, 0};
			serial_events_[offset2 + 3 + 2*dist_cutoff_ + x_dub] = 
				unbind_fr_teth;
		}
		int offset3 = 9 + 10*dist_cutoff_;

		event tether_free = {&properties_->prc1.n_bound_untethered_, 
			0,  "tether_free", offset3 + 1, 0};
		serial_events_[offset3 + 1] = tether_free;

		event tether_bound = {&n_bound_untethered_, 0, "tether_bound", 
			offset3 + 2, 0};
		serial_events_[offset3 + 2] = tether_bound;

		event untether_free = {&n_free_tethered_, 0, "untether_free",
			offset3 + 3, 0};
		serial_events_[offset3 + 3] = untether_free;

		for(int x_dub(0); x_dub <= 2*dist_cutoff_; x_dub++){
			event untether_bound = {&n_bound_tethered_[x_dub], 
				x_dub, "untether_bound", offset3 + 4 + x_dub, 0};
			serial_events_[offset3 + 4 + x_dub] = untether_bound; 
		}
		int offset4 = 13 + 12*dist_cutoff_; 

		event step = {&n_stepable_, 0, "step", offset4 + 1, 0};
		serial_events_[offset4 + 1] = step;

		for(int x_dub(0); x_dub <= 2*dist_cutoff_; x_dub++){
			event step_to_teth = {&n_stepable_to_teth_[x_dub], 
				x_dub, "step_to_teth", offset4 + 2 + x_dub, 0};
			serial_events_[offset4 + 2 + x_dub] = step_to_teth;

			event step_fr_teth = {&n_stepable_fr_teth_[x_dub], 
				x_dub, "step_fr_teth",
				offset4 + 3 + 2*dist_cutoff_ + x_dub, 0};
			serial_events_[offset4 + 3 + 2*dist_cutoff_ + x_dub] = 
				step_fr_teth; 
		}
	}
	*/
}

/*
void KinesinManagement::InitializeSamplingFunctions(){

	// Use lambda expressions to create functions that sample statistical
	// distributions with the appropriate p & n values for each population,
	// then bind them to a string key via std::make_pair and store in map

	// Function that gets num to bind_i
	auto bind_i = [&](int x_dub, int index){
		if(properties_->microtubules.n_unoccupied_ > 0){
			return properties_->gsl.SampleBinomialDist_Kinesin(
					p_bind_i_, 
					properties_->microtubules.n_unoccupied_, 
					index);
		}
		else return 0;
	};
	sampling_functs_.push_back(std::make_pair("bind_i", bind_i));

	// If tethering isn't enabled, population types are simple
	if(!parameters_->motors.tethers_active){
		// Function that gets num to bind_ii
		auto bind_ii = [&](int x_dub, int index){
			if(n_bound_i_bindable_ > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_bind_ii_, 
						n_bound_i_bindable_, 
						index);
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("bind_ii", bind_ii));

		// Function that gets num to unbind_i
		auto unbind_i = [&](int x_dub, int index){
			if(n_bound_i_ > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_unbind_i_, 
						n_bound_i_, 
						index);
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("unbind_i", unbind_i));
		// Function that gets num to unbind_ii
		auto unbind_ii = [&](int x_dub, int index){
			if(n_bound_ii_ > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_unbind_ii_,
						n_bound_ii_,
						index);
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("unbind_ii", unbind_ii));
		// Function that gets num to step
		auto step = [&](int x_dub, int index){
			if(n_stepable_ > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_step_, 
						n_stepable_,
						index);
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("step", step));
	}
	// Otherwise, serialize all populations of different tether extension
	else{
		// Function that gets num to bind_i_tethered
		auto bind_i_tethered = [&](int x_dub, int index){
			if(n_free_tethered_ > 0){
				double weight = GetWeightBindITethered(); 
				if(weight > 0){
					return properties_->gsl.SamplePoissonDist_Kinesin(
							p_bind_i_tethered_ * weight,
							index);
				}
				else return 0;
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("bind_i_tethered", 
					bind_i_tethered));
		// Function that gets num to bind_ii
		auto bind_ii = [&](int x_dub, int index){
			if(n_bound_i_bindable_ > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_bind_ii_, 
						n_bound_i_bindable_, 
						index);
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("bind_ii", bind_ii));
		// Functions that gets num to bind_ii_to/fr_teth
		auto bind_ii_to_teth = [&](int x_dub, int index){
			if(n_bound_i_bindable_to_teth_[x_dub] > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_bind_ii_to_teth_[x_dub], 
						n_bound_i_bindable_to_teth_[x_dub], 
						index);
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("bind_ii_to_teth", 
					bind_ii_to_teth));
		auto bind_ii_fr_teth = [&](int x_dub, int index){
			if(n_bound_i_bindable_fr_teth_[x_dub] > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_bind_ii_fr_teth_[x_dub], 
						n_bound_i_bindable_fr_teth_[x_dub],
						index);
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("bind_ii_fr_teth", 
					bind_ii_fr_teth));
		// Function that gets num to unbind_i
		auto unbind_i = [&](int x_dub, int index){
			if(n_bound_i_ > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_unbind_i_, 
						n_bound_i_, 
						index);
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("unbind_i", unbind_i));
		// Function that gets num to unbind_i_tethered
		auto unbind_i_tethered = [&](int x_dub, int index){
			if(n_bound_i_tethered_[x_dub] > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_unbind_i_tethered_[x_dub],
						n_bound_i_tethered_[x_dub], 
						index);
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("unbind_i_tethered", 
					unbind_i_tethered));
		// Function that gets num to unbind_ii
		auto unbind_ii = [&](int x_dub, int index){
			if(n_bound_ii_ > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_unbind_ii_,
						n_bound_ii_,
						index);
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("unbind_ii", unbind_ii));
		// Functions that get num to unbind_ii_to/fr_teth
		auto unbind_ii_to_teth = [&](int x_dub, int index){
			if(n_bound_ii_tethered_[x_dub] > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_unbind_ii_to_teth_[x_dub], 
						n_bound_ii_tethered_[x_dub], 
						index);
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("unbind_ii_to_teth", 
					bind_ii_to_teth));
		auto unbind_ii_fr_teth = [&](int x_dub, int index){
			if(n_bound_ii_tethered_[x_dub] > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_unbind_ii_fr_teth_[x_dub], 
						n_bound_ii_tethered_[x_dub], 
						index);
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("unbind_ii_fr_teth", 
					bind_ii_fr_teth));
		// Function that gets num tether_free
		auto tether_free = [&](int x_dub, int index){
			if(properties_->prc1.n_bound_untethered_ > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_tether_free_, 
						properties_->prc1.n_bound_untethered_, 
						index);
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("tether_free", 
					tether_free));
		// Function that gets num tether_bound
		auto tether_bound = [&](int x_dub, int index){
			if(n_bound_untethered_ > 0){
				double weight = GetWeightTetherBound();
				if(weight > 0){
					return properties_->gsl.SamplePoissonDist_Kinesin(
							p_tether_bound_ * weight, 
							index);
				}
				else return 0;
			}
			else return 0;	
		};
		sampling_functs_.push_back(std::make_pair("tether_bound", 
					tether_bound));
		// Function that gets num untether_free
		auto untether_free = [&](int x_dub, int index){
			if(n_free_tethered_ > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_untether_free_, 
						n_free_tethered_, 
						index);
			}
			else return 0;	
		};
		sampling_functs_.push_back(std::make_pair("untether_free", 
					untether_free));
		// Function that gets num untether_bound
		auto untether_bound = [&](int x_dub, int index){
			if(n_bound_tethered_[x_dub] > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_untether_bound_[x_dub], 
						n_bound_tethered_[x_dub], 
						index);
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("untether_bound", 
					untether_bound));
		// Function that gets num to step
		auto step = [&](int x_dub, int index){
			if(n_stepable_ > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_step_, 
						n_stepable_,
						index);
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("step", step));
		// Functions that get num to step_to/fr_teth
		auto step_to_teth = [&](int x_dub, int index){
			if(n_stepable_to_teth_[x_dub] > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_step_to_teth_[x_dub], 
						n_stepable_to_teth_[x_dub], 
						index);
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("step_to_teth", 
				step_to_teth));	
		auto step_fr_teth = [&](int x_dub, int index){
			if(n_stepable_fr_teth_[x_dub] > 0){
				return properties_->gsl.SampleBinomialDist_Kinesin(
						p_step_fr_teth_[x_dub], 
						n_stepable_fr_teth_[x_dub], 
						index);
			}
			else return 0;
		};
		sampling_functs_.push_back(std::make_pair("step_fr_teth", 
					step_fr_teth));
	}
}
*/

Kinesin* KinesinManagement::GetFreeMotor(){

	// Randomly pick a motor from the reservoir
	int i_motor = properties_->gsl.GetRanInt(n_motors_);
	Kinesin *motor = &motors_[i_motor];
	int attempts = 0;
	while(motor->heads_active_ > 0 
	|| motor->tethered_ == true){
		i_motor++;
		if(i_motor == n_motors_) i_motor = 0;
		motor = &motors_[i_motor];
		attempts++;
		if(attempts > n_motors_){
			printf("error in get_free_motor\n");
			exit(1);
		}
	}
	return motor;
}

Kinesin* KinesinManagement::GetBoundUntetheredMotor(){

	UpdateBoundUntethered();
	int i_motor = properties_->gsl.GetRanInt(n_bound_untethered_); 
	return bound_untethered_[i_motor];
}

int KinesinManagement::GetNumBoundUntethered(){

	UpdateBoundUntethered();
	return n_bound_untethered_;
}

void KinesinManagement::UpdateAllLists(){

	properties_->microtubules.UpdateUnoccupied();
	UpdateDocked();
	UpdateBoundNULL();
	UpdateBoundATP();
	UpdateBoundADPP_I();
	UpdateBoundADPP_II();
	if(parameters_->motors.tethers_active){
		properties_->prc1.UpdateUntethered();
		UpdateFreeTethered();
		UpdateBoundADPP_I_Tethered();
		UpdateBoundUntethered();
		UpdateBoundTethered();
	}
}

void KinesinManagement::UpdateFreeTethered(){

	n_free_tethered_ = 0;
	for(int i_motor = 0; i_motor < n_active_; i_motor++){
		Kinesin *motor = active_[i_motor];
		if(motor->heads_active_ == 0
		&& motor->tethered_ == true){
			free_tethered_[n_free_tethered_] = motor;
			n_free_tethered_++;
		}
	}
}

void KinesinManagement::UpdateDocked(){

	n_docked_ = 0;
	for(int i_motor = 0; i_motor < n_active_; i_motor++){
		Kinesin *motor = active_[i_motor];
		if(motor->heads_active_ == 1){
			if(motor->GetActiveHead()->ligand_ == "ADPP"){
				double site_coord = motor->GetDockedCoordinate();
				int i_site = site_coord - motor->mt_->coord_;
				if(i_site >= 0 
				&& i_site <= parameters_->microtubules.length - 1){
					// Ensure site isn't occupied (otherwise, cannot bind_ii)
					if(!motor->mt_->lattice_[i_site].occupied_){
						docked_[n_docked_] = motor->GetDockedHead(); 
						n_docked_++;
					}
				}
			}
		}
	}
}

void KinesinManagement::UpdateBoundNULL(){

	n_bound_NULL_ = 0;
	for(int i_motor = 0; i_motor < n_active_; i_motor++){
		Kinesin *motor = active_[i_motor];
		if(motor->heads_active_ == 1){
			if(motor->head_one_.site_ != nullptr
			&& motor->head_one_.ligand_ == "NULL"){
				bound_NULL_[n_bound_NULL_] = &motor->head_one_;
				n_bound_NULL_++;
			}
			if(motor->head_two_.site_ != nullptr
			&& motor->head_two_.ligand_ == "NULL"){
				bound_NULL_[n_bound_NULL_] = &motor->head_two_;
				n_bound_NULL_++;
			}
		}
		else if(motor->heads_active_ == 2){
			if(motor->head_one_.site_ != nullptr
			&& motor->head_one_.ligand_ == "NULL"
			&& motor->head_one_.trailing_){
				bound_NULL_[n_bound_NULL_] = &motor->head_one_;
				n_bound_NULL_++;
			}
			if(motor->head_two_.site_ != nullptr
			&& motor->head_two_.ligand_ == "NULL"
			&& motor->head_two_.trailing_){
				bound_NULL_[n_bound_NULL_] = &motor->head_two_;
				n_bound_NULL_++;
			}
		}
	}
}

void KinesinManagement::UpdateBoundATP(){

	n_bound_ATP_ = 0;
	for(int i_motor = 0; i_motor < n_active_; i_motor++){
		Kinesin *motor = active_[i_motor];
		if(motor->head_one_.site_ != nullptr
		&& motor->head_one_.ligand_ == "ATP"){
			bound_ATP_[n_bound_ATP_] = &motor->head_one_;
			n_bound_ATP_++;
		}
		if(motor->head_two_.site_ != nullptr
		&& motor->head_two_.ligand_ == "ATP"){
			bound_ATP_[n_bound_ATP_] = &motor->head_two_;
			n_bound_ATP_++;
		}
	}
}

void KinesinManagement::UpdateBoundADPP_I(){

	n_bound_ADPP_i_ = 0;
	for(int i_motor = 0; i_motor < n_active_; i_motor++){
		Kinesin *motor = active_[i_motor];
		if(motor->heads_active_ == 1){
			if(motor->head_one_.site_ != nullptr
			&& motor->head_one_.ligand_ == "ADPP"){
				bound_ADPP_i_[n_bound_ADPP_i_] = &motor->head_one_;
				n_bound_ADPP_i_++;
			}
			if(motor->head_two_.site_ != nullptr
			&& motor->head_two_.ligand_ == "ADPP"){
				bound_ADPP_i_[n_bound_ADPP_i_] = &motor->head_two_;
				n_bound_ADPP_i_++;
			}
		}
	}
}

void KinesinManagement::UpdateBoundADPP_I_Tethered(){

}

void KinesinManagement::UpdateBoundADPP_II(){

	n_bound_ADPP_ii_ = 0;
	for(int i_motor = 0; i_motor < n_active_; i_motor++){
		Kinesin *motor = active_[i_motor];
		if(motor->heads_active_ == 2){
			if(motor->head_one_.ligand_ == "ADPP"){
				bound_ADPP_ii_[n_bound_ADPP_ii_] = &motor->head_one_;
				n_bound_ADPP_ii_++;
			}
			if(motor->head_two_.ligand_ == "ADPP"){
				bound_ADPP_ii_[n_bound_ADPP_ii_] = &motor->head_two_;
				n_bound_ADPP_ii_++;
			}
		}
	}
}

void KinesinManagement::UpdateBoundUntethered(){

	n_bound_untethered_ = 0;
	for(int i_motor = 0; i_motor < n_active_; i_motor++){
		Kinesin *motor = active_[i_motor];
		if(motor->tethered_ == false
		&& motor->heads_active_ > 0){
			bound_untethered_[n_bound_untethered_] = motor;
			n_bound_untethered_++; 
		}
	}
}

void KinesinManagement::UpdateBoundTethered(){

	for(int x_dub = 0; x_dub <= 2*dist_cutoff_; x_dub++){
		n_bound_tethered_[x_dub] = 0;
	}
	for(int i_motor = 0; i_motor < n_active_; i_motor++){
		Kinesin *motor = active_[i_motor];
		if(motor->heads_active_ > 0
		&& motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				motor->UpdateExtension(); 
				if(motor->tethered_ == true){
					int x_dub = motor->x_dist_doubled_; 
					int index = n_bound_tethered_[x_dub];
					bound_tethered_[x_dub][index] = motor;
					n_bound_tethered_[x_dub]++; 
				}
			}
		}
	}
}

void KinesinManagement::GenerateKMCList(){

	sys_time start = sys_clock::now();

	UpdateAllLists();
	UpdateEvents();

	sys_time finish = sys_clock::now();
	auto elapsed = std::chrono::duration_cast<t_microsec>(finish - start);
	properties_->wallace.t_motors_[1] += elapsed.count();

	start = sys_clock::now();

	// currently hardcoded ... FIXME eventually 
	// Stat correcton; ensure there aren't more events than population size
	while(serial_events_[3].n_events_ + serial_events_[5].n_events_ 
	> n_bound_ADPP_i_){
		double ran = properties_->gsl.GetRanProb();
		double p_tot = serial_events_[3].p_event_ 
					 + serial_events_[5].p_event_;
		if(ran < serial_events_[3].p_event_/p_tot
		&& serial_events_[3].n_events_ > 0)
			serial_events_[3].n_events_--;
		else if(serial_events_[5].n_events_ > 0) 
			serial_events_[5].n_events_--;
	}

	int n_events = 0;
	event pre_array[5*serial_events_.size()];
	// Scan over all events; record those with >0 expected in this timestep
	for(int i_entry = 0; i_entry < serial_events_.size(); i_entry++){
		for(int i = 0; i < serial_events_[i_entry].n_events_; i++){
			pre_array[n_events] = serial_events_[i_entry];
			n_events++; 
			// Make sure we don't bamboozle ourselves here
			if(n_events > 5 * serial_events_.size()){
				printf("Error in GenerateKMCList for motors!!\n");
				exit(1);
			}
		}
	}
	// If total expected events is greater than 0, construct kmc_list_
	if(n_events > 0){
		// Trim array to appropriate size/data
		event reduced_array[n_events];
		for(int i_entry = 0; i_entry < n_events; i_entry++){
			reduced_array[i_entry] = pre_array[i_entry];
		}
		// If there are more than 1 expected events, shuffle their order
		if(n_events > 1) gsl_ran_shuffle(properties_->gsl.rng_, 
										 reduced_array, n_events, 
										 sizeof(event));
		// Transfer shuffled array into kmc_list_ structure
        kmc_list_.resize(n_events);
        for(int i_entry = 0; i_entry < n_events; i_entry++){
            kmc_list_[i_entry] = reduced_array[i_entry];
        }
	}
	// Otherwise, simply clear kmc_list_
    else{
        kmc_list_.clear();
    }
	// Conserve memory
	kmc_list_.shrink_to_fit();

	finish = sys_clock::now();
	elapsed = std::chrono::duration_cast<t_microsec>(finish - start);
	properties_->wallace.t_motors_[2] += elapsed.count();
	/*

	int n_bind_i,
		n_bind_ii, 
		n_unbind_i, 
		n_unbind_ii, 
		n_step;
	
	int n_bind_i_tethered, 
		n_bind_ii_to_teth[2*dist_cutoff_ + 1],
		n_bind_ii_fr_teth[2*dist_cutoff_ + 1],
		n_unbind_i_tethered[2*dist_cutoff_ + 1],
		n_unbind_ii_to_teth[2*dist_cutoff_ + 1],
		n_unbind_ii_fr_teth[2*dist_cutoff_ + 1],
		n_tether_free, 
		n_tether_bound,
		n_untether_free,
		n_untether_bound[2*dist_cutoff_ + 1], 
		n_step_to_teth[2*dist_cutoff_ + 1], 
		n_step_fr_teth[2*dist_cutoff_ + 1];

	n_bind_i = serial_events_[0].n_events_; 
	n_events += n_bind_i;

	// Without tethers, only 5 events to manage
	if(!parameters_->motors.tethers_active){
		n_bind_ii = serial_events_[1].n_events_;
		n_unbind_i = serial_events_[2].n_events_;
		// Ensure there aren't more KMC events than singly-bound motors
		while(n_unbind_i + n_bind_ii > n_bound_i_){
			double p_tot = p_bind_ii_ + p_unbind_i_;
			double ran = properties_->gsl.GetRanProb();
			if(ran < p_bind_ii_/p_tot
			&& n_bind_ii > 0)
				n_bind_ii--;
			else if(n_unbind_i > 0)
				n_unbind_i--;
			else if(n_bind_ii > 0)
				n_bind_ii--;
		}
		n_events += n_bind_ii;
		n_events += n_unbind_i;

		n_unbind_ii = serial_events_[3].n_events_;
		n_step = serial_events_[4].n_events_;
		// Ensure there aren't more KMC events than doubly-bound motors
		while(n_unbind_ii + n_step > n_stepable_){
			double p_tot = p_step_ + p_unbind_ii_;
			double ran = properties_->gsl.GetRanProb();
			if(ran < p_step_/p_tot
			&& n_step > 0)
				n_step--;
			else if(n_unbind_ii > 0)
				n_unbind_ii--;
			else if(n_step > 0)
				n_step--;
		}
		n_events += n_unbind_ii;
		n_events += n_step; 
	}
	// If tethers are active, things are a bit more complicated
	else{
		n_bind_ii = serial_events_[2].n_events_;
		n_unbind_i = serial_events_[5 + 4*dist_cutoff_].n_events_;
		// Ensure there aren't more KMC events than singly-bound motors
		while(n_unbind_i + n_bind_ii > n_bound_i_){
			double p_tot = p_bind_ii_ + p_unbind_i_;
			double ran = properties_->gsl.GetRanProb();
			if(ran < p_bind_ii_/p_tot
			&& n_bind_ii > 0)
				n_bind_ii--;
			else if(n_unbind_i > 0)
				n_unbind_i--;
			else if(n_bind_ii > 0)
				n_bind_ii--;
		}
		n_events += n_bind_ii;
		n_events += n_unbind_i;

		n_unbind_ii = serial_events_[7 + 6*dist_cutoff_].n_events_;
		n_step = serial_events_[14 + 12*dist_cutoff_].n_events_;
		// Ensure there aren't more KMC events than doubly-bound motors
		while(n_unbind_ii + n_step > n_stepable_){
			double p_tot = p_step_ + p_unbind_ii_;
			double ran = properties_->gsl.GetRanProb();
			if(ran < p_step_/p_tot
			&& n_step > 0)
				n_step--;
			else if(n_unbind_ii > 0)
				n_unbind_ii--;
			else if(n_step > 0)
				n_step--;
		}
		n_events += n_unbind_ii;
		n_events += n_step;

		n_bind_i_tethered = serial_events_[1].n_events_;
		n_untether_free = serial_events_[12 + 10*dist_cutoff_].n_events_;
		// Ensure there aren't more KMC events than free tethered motors
		while(n_untether_free + n_bind_i_tethered > n_free_tethered_){
			double p_unteth = p_untether_free_;
			double p_bind = p_bind_i_tethered_;
			double p_tot = p_untether_free_ + p_bind_i_tethered_; 
			double ran = properties_->gsl.GetRanProb();
			if(ran < p_untether_free_/p_tot
			&& n_bind_i_tethered > 0)
				n_bind_i_tethered--;
			else if(n_untether_free > 0)
				n_untether_free--;
			else if(n_bind_i_tethered > 0)
				n_bind_i_tethered--;
		}
		n_events += n_bind_i_tethered;
		n_events += n_untether_free;

		n_tether_free = serial_events_[10 + 10*dist_cutoff_].n_events_;
		n_events += n_tether_free;
		// FIXME correct tether_bound stat in loop??
		n_tether_bound = serial_events_[11 + 10*dist_cutoff_].n_events_;
		n_events += n_tether_bound;

		// Handle the statistics of differnt tether extensions separately 
		for(int x_dub = 2*comp_cutoff_; x_dub <= 2*dist_cutoff_; x_dub++){

			n_bind_ii_to_teth[x_dub] = 
				serial_events_[3 + x_dub].n_events_;
			n_bind_ii_fr_teth[x_dub] = 
				serial_events_[4 + 2*dist_cutoff_ + x_dub].n_events_;
			n_unbind_i_tethered[x_dub] = 
				serial_events_[6 + 4*dist_cutoff_ + x_dub].n_events_;
			n_unbind_ii_to_teth[x_dub] = 
				serial_events_[8 + 6*dist_cutoff_ + x_dub].n_events_;
			n_unbind_ii_fr_teth[x_dub] = 
				serial_events_[9 + 8*dist_cutoff_ + x_dub].n_events_;
			n_untether_bound[x_dub] = 
				serial_events_[13 + 10*dist_cutoff_ + x_dub].n_events_;
			n_step_to_teth[x_dub] = 
				serial_events_[15 + 12*dist_cutoff_ + x_dub].n_events_;
			n_step_fr_teth[x_dub] = 
				serial_events_[16 + 14*dist_cutoff_ + x_dub].n_events_;

			int n_bind_to = n_bind_ii_to_teth[x_dub];
			int n_bind_from = n_bind_ii_fr_teth[x_dub];
			int n_unbind_teth = n_unbind_i_tethered[x_dub];
			int n_unteth = n_untether_bound[x_dub];
			int n_bound_i_teth = n_bound_i_tethered_[x_dub]; 
			// Make sure there aren't too many events for singly-bound teth'd
			while(n_bind_to + n_bind_from + 
			n_unbind_teth + n_unteth > n_bound_i_teth){
				printf("%i avail\n", n_bound_i_teth);
				printf("%i, %i, %i, %i\n", n_bind_to, n_bind_from, 
						n_unbind_teth, n_unteth);
				printf("x_dub is %i\n", x_dub);
				double p_to = p_bind_ii_to_teth_[x_dub];
				double p_from = p_bind_ii_fr_teth_[x_dub];
				double p_unbind = p_unbind_i_tethered_[x_dub];
				double p_unteth = p_untether_bound_[x_dub];
				double p_tot = p_to + p_from + p_unbind + p_unteth; 
				double ran = properties_->gsl.GetRanProb();
				if(ran < p_to / p_tot
				&& n_bind_to > 0){
					n_bind_ii_to_teth[x_dub]--;
					n_bind_to = n_bind_ii_to_teth[x_dub]; 
				}
				else if(ran < (p_to + p_from) / p_tot
				&& n_bind_from > 0){
					n_bind_ii_fr_teth[x_dub]--;
					n_bind_from = n_bind_ii_fr_teth[x_dub];
				}
				else if(ran < (p_to + p_from + p_unbind) / p_tot
				&& n_unbind_teth > 0){
					n_unbind_i_tethered[x_dub]--;
					n_unbind_teth = n_unbind_i_tethered[x_dub];
				}
				else if(n_unteth > 0){
					n_untether_bound[x_dub]--;
					n_unteth = n_untether_bound[x_dub];
				}

			}
			int n_unbind_to = n_unbind_ii_to_teth[x_dub];
			int n_unbind_from = n_unbind_ii_fr_teth[x_dub];
			int n_step_to = n_step_to_teth[x_dub];
			int n_step_from = n_step_fr_teth[x_dub];
			int n_tethered = n_bound_ii_tethered_[x_dub];
			// Prevent too many steps from occuring (imagine if 2 motors
			// are bound and we roll for 1 unbind but 2 steps)
			while(n_unteth + n_step_to + n_step_from +
			n_unbind_to + n_unbind_from > n_tethered){
				// First, remove stepping events before any else
				double p_to = p_step_to_teth_[x_dub];
				double p_from = p_step_fr_teth_[x_dub];
				double p_tot = p_to + p_from;
				double ran = properties_->gsl.GetRanProb();
				if(ran < p_to / p_tot
				&& n_step_to > 0){
					n_step_to_teth[x_dub]--;
					n_step_to = n_step_to_teth[x_dub];
				}
				else if(n_step_from > 0){
					n_step_fr_teth[x_dub]--;
					n_step_from = n_step_fr_teth[x_dub];
				}
				else if(n_step_to > 0){
					n_step_to_teth[x_dub]--;
					n_step_to = n_step_to_teth[x_dub];
				}
				// If removing stepping events resolved issue, break
				if(n_unteth + n_step_to + n_step_from + 
				n_unbind_to + n_unbind_from <= n_tethered)
					break;
				// Otherwise, remove unbind/unteth events appropriately
				double p_unbind_to = p_unbind_ii_to_teth_[x_dub];
				double p_unbind_from = p_unbind_ii_fr_teth_[x_dub];
				double p_unteth = p_untether_bound_[x_dub];
				double p_tot2 = p_unbind_to + p_unbind_from + p_unteth;
				double ran2 = properties_->gsl.GetRanProb();
				if(ran2 < p_unbind_to / p_tot2
				&& n_unbind_to > 0){
					n_unbind_ii_to_teth[x_dub]--;
					n_unbind_to = n_unbind_ii_to_teth[x_dub];
				}
				else if(ran2 < (p_unbind_to + p_unbind_from) / p_tot2
				&& n_unbind_from > 0){
					n_unbind_ii_fr_teth[x_dub]--;
					n_unbind_from = n_unbind_ii_fr_teth[x_dub];
				}
				else if (n_unteth > 0){
					n_untether_bound[x_dub]--;
					n_unteth = n_untether_bound[x_dub];
				}
			} 
			n_events += n_bind_ii_to_teth[x_dub];
			n_events += n_bind_ii_fr_teth[x_dub];
			n_events += n_unbind_i_tethered[x_dub];
			n_events += n_unbind_ii_to_teth[x_dub];
			n_events += n_unbind_ii_fr_teth[x_dub];
			n_events += n_untether_bound[x_dub];
			n_events += n_step_to_teth[x_dub];
			n_events += n_step_fr_teth[x_dub];
		}
	}
	if(n_events > 0){
//		printf("n_events: %i \n", n_events);
		int pre_list[n_events];
		int kmc_index = 0;
		for(int i_event = 0; i_event < n_bind_i; i_event++){
//			printf("10\n");
			pre_list[kmc_index] = 10;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_bind_ii; i_event++){
//			printf("12\n");
			pre_list[kmc_index] = 12;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_unbind_i; i_event++){
//			printf("20\n");
			pre_list[kmc_index] = 20;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_unbind_ii; i_event++){
//			printf("22\n");
			pre_list[kmc_index] = 22;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_step; i_event++){
//			printf("60\n");
			pre_list[kmc_index] = 60;
			kmc_index++;
		}
		if(parameters_->motors.tethers_active){
			for(int i_event = 0; i_event < n_bind_i_tethered; i_event++){
//				printf("11\n");
				pre_list[kmc_index] = 11;
				kmc_index++;
			}
			for(int i_event = 0; i_event < n_tether_free; i_event++){
//				printf("30\n");
				pre_list[kmc_index] = 30;
				kmc_index++;
			}
			for(int i_event = 0; i_event < n_tether_bound; i_event++){
//				printf("31\n");
				pre_list[kmc_index] = 31;
				kmc_index++;
			}
			for(int i_event = 0; i_event < n_untether_free; i_event++){
//				printf("50\n");
				pre_list[kmc_index] = 50;
				kmc_index++;
			}
			for(int x_dub = 2*comp_cutoff_; x_dub <= 2*dist_cutoff_;x_dub++){
				int n_bind_to = n_bind_ii_to_teth[x_dub];
				for(int i_event = 0; i_event < n_bind_to; i_event++){
//					printf("%i\n", 1300 + x_dub);
					pre_list[kmc_index] = 1300 + x_dub; 
					kmc_index++;
				}
				int n_bind_from = n_bind_ii_fr_teth[x_dub];
				for(int i_event = 0; i_event < n_bind_from; i_event++){
//					printf("%i\n", 1400 + x_dub);
					pre_list[kmc_index] = 1400 + x_dub; 
					kmc_index++;
				}
				int n_unbind_teth = n_unbind_i_tethered[x_dub];
				for(int i_event = 0; i_event < n_unbind_teth; i_event++){
//					printf("%i\n", 2100 + x_dub);
					pre_list[kmc_index] = 2100 + x_dub; 
					kmc_index++;
				}
				int n_unbind_to = n_unbind_ii_to_teth[x_dub];
				for(int i_event = 0; i_event < n_unbind_to; i_event++){
//					printf("%i\n", 2300 + x_dub);
					pre_list[kmc_index] = 2300 + x_dub; 
					kmc_index++;
				}
				int n_unbind_from = n_unbind_ii_fr_teth[x_dub];
				for(int i_event = 0; i_event < n_unbind_from; i_event++){
//					printf("%i\n", 2400 + x_dub);
					pre_list[kmc_index] = 2400 + x_dub; 
					kmc_index++;
				}
				int n_unteth_b = n_untether_bound[x_dub];
				for(int i_event = 0; i_event < n_unteth_b; i_event++){
					pre_list[kmc_index] = 500 + x_dub; 	
//					printf("%i\n", 500 + x_dub);
					kmc_index++;
				}
				int n_step_to = n_step_to_teth[x_dub];
				for(int i_event = 0; i_event < n_step_to; i_event++){
					pre_list[kmc_index] = 600 + x_dub;
//					printf("%i\n", 600 + x_dub);
					kmc_index++;
				}
				int n_step_from = n_step_fr_teth[x_dub];
				for(int i_event = 0; i_event < n_step_from; i_event++){
					pre_list[kmc_index] = 700 + x_dub;
//					printf("%i\n", 700 + x_dub);
					kmc_index++;
				}
			}
		}
		if(n_events > 1){
			gsl_ran_shuffle(properties_->gsl.rng_, pre_list, n_events, 
					sizeof(int));
		}
		// Transfer pre_list to class' kmc_list_ structure
        kmc_list_.resize(n_events);
        for(int i_event = 0; i_event < n_events; i_event++){
            kmc_list_[i_event] = pre_list[i_event];
        }
    }
    else{
        kmc_list_.clear();
    }
	*/
}

void KinesinManagement::UpdateEvents(){

	for(int i_entry = 0; i_entry < serial_events_.size(); i_entry++){
		event *entry = &serial_events_[i_entry];
		/*
		printf("event \"");
		std::cout << entry->label_;
		printf("\" targets population \"");
		std::cout << entry->target_pop_;
		printf("\", which has a current size of %i\n", *entry->pop_ptr_);
		*/
		/*
		if(*entry->pop_ptr_ > 0){
			entry->n_events_ = entry->sample_stats_(entry->p_event_, 
													 *entry->pop_ptr_, 
													  entry->index_);
													  */
//		}
		// If population is 0, no. of kmc events is automatically 0
//		else entry->n_events_ = 0;
	}
}

/*
double KinesinManagement::GetWeightBindITethered(){

	double weights_summed = 0;
	// Sum over all tethered but unbound motors
	for(int i_motor = 0; i_motor < n_free_tethered_; i_motor++){
		Kinesin *motor = free_tethered_[i_motor];
		motor->UpdateNeighborSites();
		// Get weight of all neighbor sites
		int n_neighbs = motor->n_neighbor_sites_;
		for(int i_neighb = 0; i_neighb < n_neighbs; i_neighb++){
			Tubulin *site = motor->neighbor_sites_[i_neighb];
			double weight = motor->GetBindingWeight(site);
			weights_summed += weight;
		}
    }
	return weights_summed; 
}

double KinesinManagement::GetWeightTetherBound(){

	double weights_summed = 0;
	// Sum over all bound (both singly- and doubly-) but untethered motors
	for(int i_motor = 0; i_motor < n_bound_untethered_; i_motor++){
		Kinesin *motor = bound_untethered_[i_motor];
		motor->UpdateNeighborXlinks();
		// Get weight of all neighbor xlinks
		int n_neighbs = motor->n_neighbor_xlinks_; 
		for(int i_neighb = 0; i_neighb < n_neighbs; i_neighb++){
			AssociatedProtein *xlink = motor->neighbor_xlinks_[i_neighb];
			double weight = motor->GetTetheringWeight(xlink);
			weights_summed += weight;
		}
	}
	return weights_summed;
}
*/

void KinesinManagement::RunKMC(){

	sys_time start1 = sys_clock::now();
    GenerateKMCList();

	sys_time start2 = sys_clock::now();

	if(!kmc_list_.empty()){
//		printf("\nStart of Kinesin KMC cycle\n");
		for(int i_entry = 0; i_entry < kmc_list_.size(); i_entry++){
			int event_code_ = kmc_list_[i_entry].num_code_;
			switch(event_code_){
				case 10:
//					printf("Bind_I\n");
					KMC_Bind_I();
					break;
				case 11:
					printf("Bind_I_Tethered\n");
	//				KMC_Bind_I_Tethered();
					break;
				case 20:
//					printf("Bind_ATP\n");
					KMC_Bind_ATP();
					break;
				case 30:
//					printf("Hydrolyze\n");
					KMC_Hydrolyze();
					break;
				case 40:
//					printf("Bind_II\n");
					KMC_Bind_II();
					break;
				case 50:
//					printf("Unbind_II\n");
//					properties_->wallace.PrintMicrotubules();
					KMC_Unbind_II();
					break;
				case 60:
//					printf("Unbind_I\n");
//					properties_->wallace.PrintMicrotubules();
					KMC_Unbind_I();
					break;
			}
		}
	}
	sys_time finish = sys_clock::now(); 
	auto elapsed = std::chrono::duration_cast<t_microsec>(finish - start1);
	properties_->wallace.t_motors_[0] += elapsed.count();
	elapsed = std::chrono::duration_cast<t_microsec>(finish - start2);
	properties_->wallace.t_motors_[3] += elapsed.count();
/*
    if(kmc_list_.empty() == false){
        int n_events = kmc_list_.size();
//		printf("%i MOTOR KMC EVENTS\n", n_events);
		int x_dist_doubled = 0;
        for(int i_event = 0; i_event < n_events; i_event++){
            int kmc_event = kmc_list_[i_event];

			if(kmc_event >= 1300 && kmc_event < 1400){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 13; 
			}
			if(kmc_event >= 1400 && kmc_event < 1500){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 14;
			}
			if(kmc_event >= 2100 && kmc_event < 2200){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 21;
			}
			if(kmc_event >= 2300 && kmc_event < 2400){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 23;
			}
			if(kmc_event >= 2400 && kmc_event < 2500){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 24;
			}
			if(kmc_event >= 500 && kmc_event < 600){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 51;
			}
			if(kmc_event >= 600 && kmc_event < 700){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 62; 
			}	
			if(kmc_event >= 700 && kmc_event < 800){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 63;
			}

            switch(kmc_event){
				case 10:
//						printf("free motor pseudo-bound\n");
						KMC_Bind_I();
                        break;
                case 11:
//						printf("tethered motor pseudo-bound\n");
						KMC_Bind_I_Tethered();
                        break;
				case 12:
//						printf("pseudo motor bound\n");
						KMC_Bind_II();
						break;
				case 13:
//						printf("bind ii to teth\n");
						KMC_Bind_II_To_Teth_Rest(x_dist_doubled); 
						break;
				case 14: 
//						printf("bind ii from teth\n");
						KMC_Bind_II_From_Teth_Rest(x_dist_doubled);
						break;
				case 20:
//						printf("pseudo-bound motor unbound\n");
						KMC_Unbind_I();
						break;
				case 21:
//						printf("unbind i tethered\n");
						KMC_Unbind_I_Tethered(x_dist_doubled);
						break;
				case 22:
//						printf("untethered stepable motor unbound\n");
						KMC_Unbind_II();
                        break;
				case 23:
///						printf("unbind ii to teth\n");
						KMC_Unbind_II_To_Teth_Rest(x_dist_doubled);
						break;
				case 24:
//						printf("unbind ii from teth\n");
						KMC_Unbind_II_From_Teth_Rest(x_dist_doubled);
						break;
                case 30:
//						printf("free motor tethered\n");
						KMC_Tether_Free();
                        break;
				case 31:
//						printf("bound motor tethered\n");
						KMC_Tether_Bound();
						break;
				case 50:
//						printf("free motor untethered\n");
						KMC_Untether_Free();
						break;
				case 51:
//						printf("bound motor (ext %i) untethered\n", 
//								x_dist_doubled);
						KMC_Untether_Bound(x_dist_doubled);
						break;
				case 60:
//						printf("untethered motor stepped\n");
						KMC_Step();
						break;
						
				case 62:
//						printf("tethered motor (ext %i) stepped\n", 
//								x_dist_doubled);
						KMC_Step_To_Teth_Rest(x_dist_doubled);
						break;
				case 63:
//						printf("tethered motor (ext %i) stepped\n", 
//								x_dist_doubled);
						KMC_Step_From_Teth_Rest(x_dist_doubled);
						break;
            }
        }
    }
	*/
}

void KinesinManagement::KMC_Bind_I(){

	properties_->microtubules.UpdateUnoccupied();
	if(properties_->microtubules.n_unoccupied_ > 0){
		// Get random free motor
		Kinesin *motor = GetFreeMotor();
		// Get random unoccupied site
		Tubulin *site = properties_->microtubules.GetUnoccupiedSite();
		// Update site details
		site->motor_head_ = &motor->head_one_;
		site->occupied_ = true;
		// Update motor details
		motor->mt_ = site->mt_;
		motor->head_one_.site_ = site;
		motor->head_one_.ligand_ = "NULL";
		motor->head_one_.trailing_ = false;
		motor->head_two_.trailing_ = true;
		motor->heads_active_++;
		// Update active_ list
		active_[n_active_] = motor;
		motor->active_index_ = n_active_;
		n_active_++;

	}
	else{
		printf("Failed to Bind_I (motors): no unoccupied sites.\n");
	}
}

void KinesinManagement::KMC_Bind_ATP(){

	UpdateBoundNULL();
	if(n_bound_NULL_ > 0){
		// Get a random bound_NULL motor
		int i_entry = properties_->gsl.GetRanInt(n_bound_NULL_);
		Kinesin::head *head = bound_NULL_[i_entry];
		// Verify correctness
		if(head->ligand_ != "NULL"
		|| head->site_ == nullptr){
			printf("Error in KMC_Bind_ATP()\n");
			exit(1);
		}
		// Update motor head
		head->ligand_ = "ATP";
		Microtubule *mt = head->site_->mt_;
//		if(!head->trailing_){
			int i_site = head->site_->index_;
			int dx = mt->delta_x_;
			int mt_length = parameters_->microtubules.length - 1;
			if(parameters_->motors.endpausing_active){
				if(!(i_site == 0 && dx == -1) 
				&& !(i_site == mt_length && dx == 1)){
					if(!mt->lattice_[i_site + dx].occupied_){
						head->motor_->ChangeConformation();
					}
				}
			}
			else head->motor_->ChangeConformation(); 
			/*	
		 	else {
				if((i_site == 0 && dx == -1)
				|| (i_site == mt_length && dx == 1)){
					head->motor_->ChangeConformation();
				}
				else if(!mt->lattice_[i_site + dx].occupied_){
					head->motor_->ChangeConformation();
				}
			}
			*/
//		}
//		else printf("atp bound to trailing head\n");
	}
	else{
		printf("Failed to Bind_ATP: no bound_NULL motors.\n");
	}
}

void KinesinManagement::KMC_Hydrolyze(){

	UpdateBoundATP();
	if(n_bound_ATP_ > 0){
		// Get a random bound_ATP motor
		int i_entry = properties_->gsl.GetRanInt(n_bound_ATP_);
		Kinesin::head *head = bound_ATP_[i_entry];
		// Verify correctness
		if(head->ligand_ != "ATP"
		|| head->site_ == nullptr){
			printf("Error in KMC_Hydrolyze()\n");
			exit(1);
		}
		head->ligand_ = "ADPP";
	}
	else{
		printf("Failed to Hydrolyze: no bound_ATP motors.\n");
	}
}

void KinesinManagement::KMC_Bind_II(){

	UpdateDocked();
	if(n_docked_ > 0){
		// Get a random docked motor
		int i_entry = properties_->gsl.GetRanInt(n_docked_);
		Kinesin::head *docked_head = docked_[i_entry];
		Kinesin *motor = docked_head->motor_; 
		// Verify correctness
		if(docked_head->ligand_ != "ADP"
		|| motor->GetActiveHead()->ligand_ != "ADPP"
		|| motor->frustrated_){
			printf("Error in KMC_Bind_II()\n");
			exit(1);
		}
		// Verify that proposed site is unoccupied
		int i_dock = motor->GetDockedCoordinate() - motor->mt_->coord_;
		Tubulin *dock_site = &motor->mt_->lattice_[i_dock];
		if(dock_site->occupied_){
			printf("Error in KMC_Bind_II(): Dock is occupied!!\n");
			exit(1);
		}	
		// Update site
		dock_site->motor_head_ = docked_head; 
		dock_site->occupied_ = true;
		// Update motor
		docked_head->site_ = dock_site;
		docked_head->ligand_ = "NULL";
		motor->heads_active_++;
	}
	else{
		printf("Failed to Bind_II: no docked motors.\n");
		properties_->wallace.PrintMicrotubules(2);
	}
}

void KinesinManagement::KMC_Unbind_II(){

	UpdateBoundADPP_II();
	if(n_bound_ADPP_ii_ > 0){
		// Get a random bound_ADPP motor head
		int i_entry = properties_->gsl.GetRanInt(n_bound_ADPP_ii_);
		Kinesin::head *head = bound_ADPP_ii_[i_entry];
		if(head->site_ == nullptr
		|| head->ligand_ != "ADPP"){
			printf("Error in KMC_Unbind (motors): ");
			std::cout << head->ligand_;
			printf(" bound to head\n");
			exit(1);
		}
		if(head->motor_->heads_active_ == 1){
			printf("Error TWO in KMC_Unbind (motors).\n");
			exit(1);
		}
		// Update site
		head->site_->occupied_ = false;
		head->site_->motor_head_ = nullptr;
		// Update motor
		head->site_ = nullptr;
		head->ligand_ = "ADP";
		head->motor_->heads_active_--;
		if(head->motor_->frustrated_){
			if(head->trailing_){
//				printf("STEPPED TRAILING HEAD:\n");
				head->motor_->ChangeConformation();
				if(head->motor_->frustrated_){
					printf("Error THREE in KMC_Unbind().\n");
					exit(1);
				}
			}
			else{
//				printf("STEPPED LEADING HEAD:\n");
				head->motor_->frustrated_ = false; 
			}
		}
	}
	else{
		printf("Failed to KMC_Unbind_II: no bound_ADPP motors.\n");
	}
}

void KinesinManagement::KMC_Unbind_I(){

	UpdateBoundADPP_I();
	if(n_bound_ADPP_i_ > 0){
		// Get a random bound_ADPP motor head
		int i_entry = properties_->gsl.GetRanInt(n_bound_ADPP_i_);
		Kinesin::head *head = bound_ADPP_i_[i_entry];
		if(head->site_ == nullptr
		|| head->ligand_ != "ADPP"){
			printf("Error in KMC_Unbind (motors): ");
			std::cout << head->ligand_;
			printf(" bound to head\n");
			exit(1);
		}
		if(head->motor_->heads_active_ == 2){
			printf("Error TWO in KMC_Unbind (motors).\n");
			exit(1);
		}
		// Update site
		head->site_->occupied_ = false;
		head->site_->motor_head_ = nullptr;
		// Update motor
		head->site_ = nullptr;
		head->ligand_ = "ADP";
		head->motor_->heads_active_--;
		head->motor_->mt_ = nullptr; 
		// Remove this motor from active_, replace with last entry
		Kinesin *last_entry = active_[n_active_ - 1];
		int this_index = head->motor_->active_index_; 
		if(this_index != n_active_ - 1){
			active_[this_index] = last_entry; 
			last_entry->active_index_ = this_index; 
		}
		n_active_--;
	}
	else{
		printf("Failed to KMC_Unbind_I: no bound_ADPP motors.\n");
	}
}

/*
void KinesinManagement::KMC_Bind_I(){

    // Make sure that at least one unbound motor exists
	properties_->microtubules.UpdateUnoccupied();
	if(properties_->microtubules.n_unoccupied_ > 0){	
        Kinesin *motor = GetFreeMotor();
		MicrotubuleManagement *mts = &properties_->microtubules;
		Tubulin *site = mts->GetUnoccupiedSite();
		Microtubule *mt = site->mt_;
		// Update site details
		site->motor_ = motor;
		site->occupied_ = true;
		// Update motor details
		motor->mt_ = mt;
		motor->front_site_ = site;
		motor->heads_active_++;
		// Update active_ list
		active_[n_active_] = motor;
		motor->active_index_ = n_active_;
	//	printf("added motor %i to active_[%i]\n", motor->ID_, n_active_);
		n_active_++;
	}
	else{
		printf("Error in Bind_I: no unoccupied sites.\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Bind_II(){

	UpdateBoundIBindable();
	if(n_bound_i_bindable_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound_i_bindable_);
		Kinesin *motor = bound_i_bindable_[i_entry];
		if(motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				printf("tried to bind_ii motor w/ an ACTIVE tether??\n");
				exit(1);
			}
		}
		Tubulin *bound_site = motor->GetActiveHeadSite();
		Microtubule *mt = motor->mt_;
		int i_site = bound_site->index_;
		int dx = motor->mt_->delta_x_;
		int i_plus = motor->mt_->plus_end_;
		int i_minus = motor->mt_->minus_end_;
		Tubulin *front_site;
		Tubulin	*rear_site; 
		// Choose where to bind second head
		if(i_site == i_plus){
			if(mt->lattice_[i_site - dx].occupied_ == false){
				front_site = bound_site;
				rear_site = &mt->lattice_[i_site - dx];	
			}
			else{
				printf("something wrong in bind_ii ONE\n");
				exit(1);
			}
		}
		else if(i_site == i_minus){
			if(mt->lattice_[i_site + dx].occupied_ == false){
				front_site = &mt->lattice_[i_site + dx];
				rear_site = bound_site;
			}
			else{
				printf("something wrong in bind_ii TWO\n");
				exit(1);
			}
		}
		else if(mt->lattice_[i_site + dx].occupied_ == false
		&& mt->lattice_[i_site - dx].occupied_ == false){
			double ran = properties_->gsl.GetRanProb();
			if(ran < 0.5){
				front_site = &mt->lattice_[i_site + dx];
				rear_site = bound_site;
			}
			else{
				front_site = bound_site;
				rear_site = &mt->lattice_[i_site - dx];
			}
		}
		else if(mt->lattice_[i_site + dx].occupied_ == false){
			front_site = &mt->lattice_[i_site + dx];
			rear_site = bound_site;
		}
		else if(mt->lattice_[i_site - dx].occupied_ == false){
			front_site = bound_site;
			rear_site = &mt->lattice_[i_site - dx];
		}
		else{
			printf("something wrong in bind_ii THREE\n");
			exit(1); 
		}
		// Update site details
		front_site->motor_ = motor;
		front_site->occupied_ = true;
		rear_site->motor_ = motor;
		rear_site->occupied_ = true;
		// Update motor details
		motor->front_site_ = front_site;
		motor->rear_site_ = rear_site;
		motor->heads_active_++;
	}
	else{
//		printf("Error in Bind_II: no singly-bound motors. \n");
//		exit(1);
	}
}	

void KinesinManagement::KMC_Unbind_I(){

	UpdateBoundI();
	if(n_bound_i_ > 0){
		// Randomly pick a motor
		int i_entry = properties_->gsl.GetRanInt(n_bound_i_);
		Kinesin *motor = bound_i_[i_entry];
		Tubulin *site = motor->GetActiveHeadSite();
		// Update site details
		site->motor_ = nullptr;
		site->occupied_ = false;
		// Update motor details
		motor->mt_ = nullptr;
		motor->front_site_ = nullptr;
		motor->rear_site_ = nullptr;
		motor->heads_active_--;
		// Check to see if we had a satellite xlink; if so, untether it
		if(motor->tethered_ == true){
			if(motor->xlink_->heads_active_ == 0)
				motor->UntetherSatellite();
			else{
				printf("ummm whywouldudothis??? mot_unbind_i\n");
				exit(1);
			}
		}
		// Remove this motor from active_, replace with last entry
		Kinesin *last_entry = active_[n_active_ - 1];
		int this_index = motor->active_index_; 
		if(this_index != n_active_ - 1){
			active_[this_index] = last_entry; 
			last_entry->active_index_ = this_index; 
		}
		n_active_--;
	}
	else{
		printf("Error in Unbind_I: no pseudo bound motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Unbind_II(){

	UpdateBoundII();
    if(n_bound_ii_ > 0){
		// Randomly pick a motor
		int i_entry = properties_->gsl.GetRanInt(n_bound_ii_);
		Kinesin *motor = bound_ii_[i_entry];
		// Roll and randomly pick a head to unbind
		double ran = properties_->gsl.GetRanProb();
		// Update site details
		if(ran < 0.5){
			motor->front_site_->motor_ = nullptr;
			motor->front_site_->occupied_ = false;
			motor->front_site_ = nullptr;
		}
		else{
			motor->rear_site_->motor_ = nullptr;
			motor->rear_site_->occupied_ = false;
			motor->rear_site_ = nullptr;
		}
		// Update motor statistics
		motor->heads_active_--;
	}
	else{
		printf("Error in Unbind_II: no bound untethered motors!\n");
//      exit(1);
    }
}

void KinesinManagement::KMC_Step(){

    // Make sure there is at least one stepable untethered motor
	UpdateStepable();
	if(n_stepable_ > 0){
		// Randomly pick a motor
		int i_entry = properties_->gsl.GetRanInt(n_stepable_);
		Kinesin *motor = stepable_[i_entry];
		Tubulin *old_front_site = motor->front_site_; 
		Tubulin *old_rear_site = motor->rear_site_;
		int i_old_front = old_front_site->index_;
		int dx = motor->mt_->delta_x_;
		Tubulin *new_front_site = &motor->mt_->lattice_[i_old_front + dx];
		Tubulin *new_rear_site = old_front_site; 
		// Update front head of motor
		motor->front_site_ = new_front_site;
		new_front_site->motor_ = motor;
		new_front_site->occupied_ = true;
		// Update rear head of motor 
		motor->rear_site_ = new_rear_site;
		old_rear_site->motor_ = nullptr;
		old_rear_site->occupied_ = false;
	}
	else{
		printf("Error in Step: no stepable untethered motors\n");
//		exit(1);
    }
}

void KinesinManagement::KMC_Bind_I_Tethered(){
	
	UpdateFreeTethered();
	properties_->microtubules.UpdateUnoccupied();
	if(n_free_tethered_ > 0
	&& properties_->microtubules.n_unoccupied_ > 0){
		// Pick a random tethered_free motor to bind
		int i_motor = properties_->gsl.GetRanInt(n_free_tethered_);
		Kinesin *motor = free_tethered_[i_motor];
		Tubulin *site = motor->GetWeightedNeighborSite();
		int attempts = 0;
		while(site == nullptr){
			if(attempts > 10*n_free_tethered_){
				break;
			}
			i_motor = properties_->gsl.GetRanInt(n_free_tethered_);
			motor = free_tethered_[i_motor];
			site = motor->GetWeightedNeighborSite();
			attempts++;	
		}
		if(site != nullptr){
			// Update site details
			site->motor_ = motor;
			site->occupied_ = true;
			// Update motor details
			motor->mt_ = site->mt_;
			motor->front_site_ = site;
			motor->heads_active_++;
			motor->UpdateExtension();
		}
		else{
			printf("Failed to bind_i_teth\n");
			exit(1);
		}
	}
	else if(properties_->microtubules.n_unoccupied_ == 0){
		printf("Error in Bind_I_Tethered: no unoccupied sites\n");
//		exit(1);
	}
	else{
		printf("Error in Bind_I_Tethered: no tethered free motors\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Bind_II_To_Teth_Rest(int x_dist_doubled){

	UpdateBindableToTeth();
	int x_dub = x_dist_doubled; 
	int n_bound = n_bound_i_bindable_to_teth_[x_dub];
	if(n_bound > 0){
		// Randomly pick a motor
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_i_bindable_to_teth_[x_dub][i_entry];
		if(motor->x_dist_doubled_ != x_dist_doubled){
			printf("noplz in bind_ii_to_teth_rest (motors)\n");
			exit(1);
		}
		Tubulin *bound_site = motor->GetActiveHeadSite();
		Microtubule *mt = motor->mt_;
		int i_site = bound_site->index_;
		int mt_dx = bound_site->mt_->delta_x_;
		int dx = motor->GetDirectionTowardRest();
		int i_plus = motor->mt_->plus_end_;
		int i_minus = motor->mt_->minus_end_;
		Tubulin *front_site;
		Tubulin *rear_site; 
		// Choose where to bind second head
		if(i_site == i_plus){
			if(dx != mt_dx){
				if(mt->lattice_[i_site + dx].occupied_ == false){
					front_site = bound_site;
					rear_site = &mt->lattice_[i_site + dx];	
				}
			}
			else{
				printf("something wrong in bind_ii_to ONE\n");
				exit(1);
			}
		}
		else if(i_site == i_minus){
			if(dx == mt_dx){
				if(mt->lattice_[i_site + dx].occupied_ == false){
					front_site = &mt->lattice_[i_site + dx];
					rear_site = bound_site;
				}
			}
			else{
				printf("something wrong in bind_ii_to TWO\n");
				exit(1);
			}
		}
		else if(mt->lattice_[i_site + dx].occupied_ == false){
			if(dx == mt_dx){
				front_site = &mt->lattice_[i_site + dx];
				rear_site = bound_site;
			}
			else{
				rear_site = &mt->lattice_[i_site + dx];
				front_site = bound_site;
			}
		}
		else{
			printf("something wrong in bind_ii_to THREE\n");
			exit(1);
		}
		// Update site details
		front_site->motor_ = motor;
		front_site->occupied_ = true;
		rear_site->motor_ = motor;
		rear_site->occupied_ = true;
		// Update motor details
		motor->front_site_ = front_site;
		motor->rear_site_ = rear_site;
		motor->heads_active_++;
		motor->UpdateExtension();
	}
	else{
		printf("Error in Bind_II_To_Teth: no singly-bound motors. \n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Bind_II_From_Teth_Rest(int x_dist_doubled){

	UpdateBindableFromTeth();
	int x_dub = x_dist_doubled; 
	int n_bound = n_bound_i_bindable_fr_teth_[x_dub];
	if(n_bound > 0){
		// Randomly pick a motor
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_i_bindable_fr_teth_[x_dub][i_entry];
		if(motor->x_dist_doubled_ != x_dist_doubled){
			printf("noplz in bind_ii_fr_teth_rest (motors)\n");
			exit(1);
		}
		Tubulin *bound_site = motor->GetActiveHeadSite();
		Microtubule *mt = motor->mt_;
		int i_site = bound_site->index_;
		int mt_dx = bound_site->mt_->delta_x_;
		int dx = -1 * motor->GetDirectionTowardRest();
		int i_plus = motor->mt_->plus_end_;
		int i_minus = motor->mt_->minus_end_;
		Tubulin *front_site;
		Tubulin *rear_site; 
		// Choose where to bind second head
		if(i_site == i_plus){
			if(dx != mt_dx){
				if(mt->lattice_[i_site + dx].occupied_ == false){
					front_site = bound_site;
					rear_site = &mt->lattice_[i_site + dx];	
				}
			}
			else{
				printf("something wrong in bind_ii_fr ONE\n");
				exit(1);
			}
		}
		else if(i_site == i_minus){
			if(dx == mt_dx){
				if(mt->lattice_[i_site + dx].occupied_ == false){
					rear_site = bound_site;
					front_site = &mt->lattice_[i_site + dx];
				}
			}
			else{
				printf("something wrong in bind_ii_fr TWO\n");
				exit(1);
			}
		}
		else if(mt->lattice_[i_site + dx].occupied_ == false){
			if(dx == mt_dx){
				front_site = &mt->lattice_[i_site + dx];
				rear_site = bound_site;
			}
			else{
				rear_site = &mt->lattice_[i_site + dx];
				front_site = bound_site;
			}
		}
		else{
			printf("something wrong in bind_ii_fr THREE\n");
			exit(1);
		}
		// Update site details
		front_site->motor_ = motor;
		front_site->occupied_ = true;
		rear_site->motor_ = motor;
		rear_site->occupied_ = true;
		// Update motor details
		motor->front_site_ = front_site;
		motor->rear_site_ = rear_site;
		motor->heads_active_++;
		motor->UpdateExtension();
	}
	else{
		printf("Error in Bind_II_From_Teth: no singly-bound motors. \n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Unbind_I_Tethered(int x_dist_doubled){

	UpdateBoundITethered();
	int n_bound = n_bound_i_tethered_[x_dist_doubled]; 
	if(n_bound > 0){
		// Randomly pick a motor
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_i_tethered_[x_dist_doubled][i_entry];
		if(motor->x_dist_doubled_ != x_dist_doubled){
			printf("\"not a fan\" -- unbind_i_teth (motors)\n");
			exit(1);
		}
		Tubulin *site = motor->GetActiveHeadSite();
		// Update site details
		site->motor_ = nullptr;
		site->occupied_ = false;
		// Update motor details
		motor->mt_ = nullptr;
		motor->front_site_ = nullptr;
		motor->rear_site_ = nullptr;
		motor->heads_active_--;
	}
	else{
		printf("Error in Unbind_I_Tethered: no pseudo bound motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Unbind_II_To_Teth_Rest(int x_dist_doubled){

	UpdateBoundIITethered();
	int x_dub = x_dist_doubled;
	int n_bound = n_bound_ii_tethered_[x_dub];
    if(n_bound > 0){
		// Randomly pick a motor
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_ii_tethered_[x_dub][i_entry];
		if(motor->x_dist_doubled_ != x_dist_doubled){
			printf("pretty whack in unbind_ii_to_teth (motor)\n");
			exit(1);
		}
		// Unbind site farther from rest to bring stalk closer to rest
		Tubulin *unbind_site = motor->GetSiteFartherFromRest();
		// Update site details
		if(unbind_site == motor->front_site_){
			motor->front_site_->motor_ = nullptr;
			motor->front_site_->occupied_ = false;
			motor->front_site_ = nullptr;
		}
		else{
			motor->rear_site_->motor_ = nullptr;
			motor->rear_site_->occupied_ = false;
			motor->rear_site_ = nullptr;
		}
		motor->heads_active_--;
		motor->UpdateExtension();
	}
	else{
		printf("Error in Unbind_II_To_Teth: no teth bound motors\n");
		exit(1);
	}
}

void KinesinManagement::KMC_Unbind_II_From_Teth_Rest(int x_dist_doubled){

	UpdateBoundIITethered();
	int x_dub = x_dist_doubled;
	int n_bound = n_bound_ii_tethered_[x_dub];
    if(n_bound > 0){
		// Randomly pick a motor
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_ii_tethered_[x_dub][i_entry];
		if(motor->x_dist_doubled_ != x_dist_doubled){
			printf("pretty whack in unbind_ii_fr_teth (motor)\n");
			exit(1);
		}
		// Unbind site closer to rest to bring stalk further from rest
		Tubulin *unbind_site = motor->GetSiteCloserToRest();
		// Update site details
		if(unbind_site == motor->front_site_){
			motor->front_site_->motor_ = nullptr;
			motor->front_site_->occupied_ = false;
			motor->front_site_ = nullptr;
		}
		else{
			motor->rear_site_->motor_ = nullptr;
			motor->rear_site_->occupied_ = false;
			motor->rear_site_ = nullptr;
		}
		// Update motor statistics
		motor->heads_active_--;
		motor->UpdateExtension();
	}
	else{
		printf("Error in Unbind_II_From_Teth: no teth bound motors\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Tether_Free(){

	if(properties_->prc1.n_bound_untethered_ > 0){
		Kinesin *motor = GetFreeMotor();
		// Randomly pick an xlink
		AssociatedProtein *xlink = properties_->prc1.GetUntetheredXlink();
		// Update motor and xlink details
		motor->xlink_ = xlink;
		motor->tethered_ = true;
		xlink->tethered_ = true;
		xlink->motor_ = motor; 	
		// Update active_ list
		active_[n_active_] = motor;
		motor->active_index_ = n_active_;
		n_active_++;
	}
	else{
		printf("Error in Tether_Free: no bound untethered xlinks\n");
//		exit(1);		
	}
}

void KinesinManagement::KMC_Tether_Bound(){

	UpdateBoundUntethered();
	if(n_bound_untethered_ > 0){
		// Randomly pick a motor
		int i_motor = properties_->gsl.GetRanInt(n_bound_untethered_);
		Kinesin *motor = bound_untethered_[i_motor];
		AssociatedProtein* xlink = motor->GetWeightedNeighborXlink();
		// it's possible to pick a motor w/ no neighbs, 
		// so make sure we reroll if that is the case
		int attempts = 0; 
		while(xlink == nullptr){
			if(attempts > 10*n_bound_untethered_){
				break;
			}
			i_motor = properties_->gsl.GetRanInt(n_bound_untethered_);
			motor = bound_untethered_[i_motor];
			xlink = motor->GetWeightedNeighborXlink();
			attempts++;
		}
		if(xlink != nullptr){
			// Update motor and xlink details
			motor->xlink_ = xlink;
			motor->tethered_ = true;
			xlink->tethered_ = true;
			xlink->motor_ = motor;
			motor->UpdateExtension();
			if(motor->tethered_ == false){
				printf("what in tether_bound nation\n");
				exit(1);
			}
		}
		else{
			printf("Failed to tether bound motor\n");
		}
	}
	else{
		printf("Error in Tether_Bound: no bound untethered motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Untether_Bound(int x_dist_doubled){

	UpdateBoundTethered();
	int n_bound_tethered = n_bound_tethered_[x_dist_doubled];
	if(n_bound_tethered  > 0){
		// Randomly pick a motor
		int i_entry = properties_->gsl.GetRanInt(n_bound_tethered);
		Kinesin *motor = bound_tethered_[x_dist_doubled][i_entry];
		if(motor->x_dist_doubled_ != x_dist_doubled){
			printf("error in Untether_Bound (motor) x_dub: ");
			printf("%i received; %i in motor\n", x_dist_doubled, 
					motor->x_dist_doubled_);
			exit(1);
		}
		AssociatedProtein *xlink = motor->xlink_;
		// Update motor and xlink details
		xlink->motor_ = nullptr;
		xlink->tethered_ = false;
		motor->xlink_ = nullptr;
		motor->tethered_ = false; 
		// Update statistics
		motor->UpdateExtension();
	}
	else{
		printf("Error in Untether_Bound: no bound tethered motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Untether_Free(){

	UpdateFreeTethered();
	if(n_free_tethered_ > 0){
		// Randomly pick a motor
		int i_entry = properties_->gsl.GetRanInt(n_free_tethered_);
		Kinesin *motor = free_tethered_[i_entry];
		AssociatedProtein *xlink = motor->xlink_;
		// Update motor and xlink detail
		xlink->motor_ = nullptr; 
		xlink->tethered_ = false; 
		motor->xlink_ = nullptr;
		motor->tethered_ = false;
		// Remove this motor from active_, replace with last entry
		// (only if there are more than 1 active motor)
		if(n_active_ > 1){
			int this_index = motor->active_index_; 
			Kinesin *last_entry = active_[n_active_ - 1];
			active_[n_active_ - 1] = nullptr; 
			active_[this_index] = last_entry; 
		}
		n_active_--;
	}
	else{
		printf("Error in Untether_Free: no free tethered motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Step_To_Teth_Rest(int x_dist_doubled){

	UpdateStepableTethered();
	int n_stepable = n_stepable_to_teth_[x_dist_doubled]; 
	if(n_stepable > 0){
		// Randomly pick a motor
		int i_entry = properties_->gsl.GetRanInt(n_stepable);
		Kinesin *motor = stepable_to_teth_[x_dist_doubled][i_entry];
		if(motor->x_dist_doubled_ != x_dist_doubled){
			printf("error in Step_To_Teth_Rest (motor) x_dub: ");
			printf("%i received; %i in motor\n", x_dist_doubled, 
					motor->x_dist_doubled_);
			exit(1);
		}
		if(motor->mt_->delta_x_ != motor->GetDirectionTowardRest()){
			printf("wat in step_to_rest motor KMC\n");
			exit(1);
		}
		Tubulin *old_front_site = motor->front_site_;
		Tubulin *old_rear_site = motor->rear_site_;
		int i_old_front = old_front_site->index_;
		int dx = motor->mt_->delta_x_;
		Tubulin *new_front_site = &motor->mt_->lattice_[i_old_front + dx];
		Tubulin *new_rear_site = old_front_site; 
		// Update front head of motor
		motor->front_site_ = new_front_site;
		new_front_site->motor_ = motor;
		new_front_site->occupied_ = true;
		// Update rear head of motor 
		motor->rear_site_ = new_rear_site;
		old_rear_site->motor_ = nullptr;
		old_rear_site->occupied_ = false;
		// Update statistics
		motor->UpdateExtension();
		if(motor->x_dist_doubled_ == x_dist_doubled){
			printf("error in step_to_teth (motor)\n");
			exit(1);
		}
	}
	else{
		printf("Error in Step_To_Teth_Rest: no stepable motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Step_From_Teth_Rest(int x_dist_doubled){

	UpdateStepableTethered();
	int n_stepable = n_stepable_fr_teth_[x_dist_doubled]; 
	if(n_stepable > 0){
		// Randomly pick a motor
		int i_entry = properties_->gsl.GetRanInt(n_stepable);
		Kinesin *motor = stepable_fr_teth_[x_dist_doubled][i_entry];
		if(motor->x_dist_doubled_ != x_dist_doubled){
			printf("error in Step_From_Teth (motor) x_dub: ");
			printf("%i received; %i in motor\n", x_dist_doubled,
					motor->x_dist_doubled_);
			exit(1);
		}
		if(motor->mt_->delta_x_ != -1 * motor->GetDirectionTowardRest()){
			printf("wat in step_fr_rest motor KMC\n");
			exit(1);
		}
		Tubulin *old_front_site = motor->front_site_;
		Tubulin *old_rear_site = motor->rear_site_;
		int i_old_front = old_front_site->index_;
		int dx = motor->mt_->delta_x_;
		Tubulin *new_front_site = &motor->mt_->lattice_[i_old_front + dx];
		Tubulin *new_rear_site = old_front_site; 
		// Update front head of motor
		motor->front_site_ = new_front_site;
		new_front_site->motor_ = motor;
		new_front_site->occupied_ = true;
		// Update rear head of motor 
		motor->rear_site_ = new_rear_site;
		old_rear_site->motor_ = nullptr;
		old_rear_site->occupied_ = false;
		// Update statistics
		motor->UpdateExtension();
		if(motor->x_dist_doubled_ == x_dist_doubled){
			printf("error in step_fr_teth (motor)\n");
			exit(1);
		}
	}
	else{
		printf("Error in Step_From_Teth_Rest: no stepable motors!\n");
//		exit(1);
	}
}
*/
