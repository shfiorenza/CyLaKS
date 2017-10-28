#include "master_header.h"
#include "curator.h"

Curator::Curator(){
}

void Curator::Initialize(system_parameters *parameters, 
						 system_properties *properties){

	parameters_ = parameters;
	properties_ = properties;
	start_ = clock();
	SetParameters();
	OutputSimDetails();
}

void Curator::SetParameters(){
	
	int n_steps = parameters_->n_steps;
	int n_datapoints = parameters_->n_datapoints;
	data_threshold_ = parameters_->data_threshold;

	range_of_data_ = n_steps - data_threshold_;
	n_pickup_ = range_of_data_/n_datapoints;
	equil_milestone_ = data_threshold_/10;
	data_milestone_ = range_of_data_/10;
}

void Curator::OutputSimDetails(){

	int n_steps = parameters_->n_steps;
	double k_on = parameters_->k_on;
	double c_motor = parameters_->c_motor;
	double k_off = parameters_->k_off;
	double switch_rate = parameters_->switch_rate;
	double motor_speed = parameters_->motor_speed;
	double delta_t = parameters_->delta_t;
	double p_bind = k_on*c_motor*delta_t;
	double p_unbind = k_off*delta_t;
	double p_switch = switch_rate*delta_t;
	double p_step = motor_speed*delta_t;
    printf("Total simulation duration: %g seconds\n", delta_t*n_steps);
    printf("Timestep duration: %g seconds\n\n", delta_t);
    printf("Motor binding frequency: %g per timestep\n", p_bind);
    printf("Motor unbinding frequency: %g per timestep\n", p_unbind);
    printf("Motor switching frequency: %g per timestep\n", p_switch);
    printf("Motor stepping frequency: %g sites per timestep\n\n", p_step);
    fflush(stdout);
}

void Curator::PrintMicrotubules(){

    int n_mts = parameters_->n_microtubules;
    int mt_length = parameters_->length_of_microtubule;
    for(int i_mt = 0; i_mt < n_mts; i_mt++){
		Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
        for(int i_site = 0; i_site < mt_length; i_site++){
            if(mt->lattice_[i_site].occupied_ == false)
                printf("=");
			else if(mt->lattice_[i_site].xlink_ != nullptr)
				printf("X");
            else if(mt->lattice_[i_site].motor_ != nullptr){
				int i_front = mt->lattice_[i_site].motor_->front_site_->index_;
				int i_rear = mt->lattice_[i_site].motor_->rear_site_->index_;
				if(i_front > i_rear){
					if(i_site == i_rear)
						printf("(");
					if(i_site == i_front)
						printf(")");
				} 
				if(i_front < i_rear){
					if(i_site == i_front)
						printf("(");
					if(i_site == i_rear)
						printf(")");
				}
           	}
        }
		printf(" %i\n", mt->polarity_);
    }   
	printf("\n");
}

void Curator::PrintMicrotubules(double pause_duration){

	PrintMicrotubules();
	PauseSim(pause_duration);

}

void Curator::OutputData(){

	FILE *occupancy_file = properties_->occupancy_file_;
	FILE *motor_ID_file = properties_->motor_ID_file_;
	FILE *xlink_ID_file = properties_->xlink_ID_file_;
	int n_mts = parameters_->n_microtubules;
	int mt_length = parameters_->length_of_microtubule;
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
		int occupancy_array[mt_length], 
			motor_ID_array[mt_length],
			xlink_ID_array[mt_length];
		int *occupancy_ptr = occupancy_array, 
			*motor_ID_ptr = motor_ID_array, 
			*xlink_ID_ptr = xlink_ID_array;
		for(int i_site = 0; i_site < mt_length; i_site++){
			Tubulin *site = &mt->lattice_[i_site];
			// If unoccupied, store the speciesID of tubulin to occupancy file
			// and an ID of -1 (null) to motor/xlink ID files 
			if(site->occupied_ == false){
				occupancy_array[i_site] = site->speciesID_;
				motor_ID_array[i_site] = -1;
				xlink_ID_array[i_site] = -1;
			}
			// If occupied by xlink, store its species ID to occupancy_file,
			// its unique ID to the xlink ID file, and -1 to motor ID file
			else if(site->xlink_ != nullptr){
				occupancy_array[i_site] = site->xlink_->speciesID_;
				motor_ID_array[i_site] = -1;
				xlink_ID_array[i_site] = site->xlink_->ID_;
			}
			// If occupied by motor, store its species ID to occupancy_file, 
			// its unique ID to the motor ID file, and -1 to xlink ID file
			else if(site->motor_ != nullptr){
				occupancy_array[i_site] = site->motor_->speciesID_;
				motor_ID_array[i_site] = site->motor_->ID_;
				xlink_ID_array[i_site] = -1;
			}
		}
		// Write the data to respective files one microtubule at a time
		fwrite(occupancy_ptr, sizeof(int), mt_length, occupancy_file);
		fwrite(motor_ID_ptr, sizeof(int), mt_length, motor_ID_file);
		fwrite(xlink_ID_array, sizeof(int), mt_length, xlink_ID_file);
	}	
}

void Curator::UpdateTimestep(int i_step){

	properties_->current_step_ = i_step;
	// Give updates on equilibrium process (every 10%)
	if(i_step < data_threshold_ && i_step%equil_milestone_ == 0){
		printf("Equilibration is %i percent complete (step number %i)\n", 
							  (int)(i_step/equil_milestone_)*10, i_step);
	}
	// Start data collection at appropriate step threshold
	else if(i_step >= data_threshold_){
		int delta = i_step - data_threshold_;
		// Collect data every n_pickup timesteps
		if(delta%n_pickup_ == 0){
			OutputData();
		}
		if(delta%data_milestone_ == 0){
			printf("Data collection is %i percent complete (step number %i)\n",
					    		      (int)(delta/data_milestone_)*10, i_step);
		}
		else if(delta == range_of_data_ - 1){
			printf("Done!");
		}
	}
}

void Curator::PauseSim(double duration){
	
	// Duration should be input in seconds
	pause_dur_.tv_sec = (int) duration;	 
	pause_dur_.tv_nsec = (duration - (int)duration)*100000000;
	nanosleep(&pause_dur_, NULL);
}

void Curator::OutputSimDuration(){

	finish_ = clock();
	sim_duration_ = (double)(finish_ - start_)/CLOCKS_PER_SEC;
	stream_ = fopen("sim_duration.dat", "w");
	fprintf(stream_, "Time to execute sim: %f seconds.\n", sim_duration_);
	fclose(stream_);
	printf(" Time to execute: %f seconds.\n\n", sim_duration_);
}

void Curator::CleanUp(){

	fclose(properties_->occupancy_file_);
	fclose(properties_->motor_ID_file_);
	fclose(properties_->xlink_ID_file_);
}