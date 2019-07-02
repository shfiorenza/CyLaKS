#include "master_header.h"
#include "tubulin.h"

Tubulin::Tubulin(){
}

void Tubulin::Initialize(system_parameters *parameters, 
	system_properties *properties, Microtubule *mt, int i_site){

	parameters_ = parameters;
	properties_ = properties;
	index_ = i_site;
	coord_ = i_site;
	mt_ = mt;
}

bool Tubulin::EquilibriumInSameDirection(){

	if(occupied_ == true
	&& xlink_ != nullptr){
		if(xlink_->heads_active_ == 2
		&& xlink_->tethered_ == true){
			Kinesin* motor = xlink_->motor_;
			double site_coord = index_ + mt_->coord_;
			double xlink_rest = xlink_->GetAnchorCoordinate();
			double motor_rest = motor->GetRestLengthCoordinate();
			// When tether is extended past rest length, ...
			if(motor->x_dist_doubled_ > 2*motor->rest_dist_){
				// ... equils are in same dir. if coords are on same side
				if((xlink_rest > site_coord && motor_rest > site_coord)
				|| (xlink_rest < site_coord && motor_rest < site_coord)){
					return true;
				}
				// ... equils are in oppo dir. if coords are on oppo side
				else if((xlink_rest > site_coord && motor_rest < site_coord)
				|| (xlink_rest < site_coord && motor_rest > site_coord)){
					return false;
				}
				else{
					printf("error in Tubulin::SpringEquilOnSameSide\n");
					exit(1);
				}
			}
			// When tether at or compressed below rest length, ...
			else{
				// ... equils are in same dir. if coords are on oppo side
				if((xlink_rest > site_coord && motor_rest < site_coord)
				|| (xlink_rest < site_coord && motor_rest > site_coord)){
					return true;
				}
				// ... equils are in oppo dir. if coords are on same side
				else if((xlink_rest > site_coord && motor_rest > site_coord)
		 		|| (xlink_rest < site_coord && motor_rest < site_coord)){
					return false;
				}
				else{
					printf("error in Tubulin::SpringEquilOnSameSide TWO\n");
					exit(1);
				}
			}
		}
		else{
			printf("error NO. 2 in tubulin: spring equil on same side ??\n");
			exit(1);
		}
	}
	else{
		printf("error in tubulin:spring equil on same side\n");
		exit(1);
	}
}

int Tubulin::GetPRC1NeighborCount(){

	int n_neighbs = 0;
	int i_plus = mt_->plus_end_;                                    
	int i_minus = mt_->minus_end_;                                  
	int dx = mt_->delta_x_; 
	if(index_ == i_plus){                                                 
		if(mt_->lattice_[index_-dx].xlink_ != nullptr)              
			n_neighbs++;
	}
	else if(index_ == i_minus){                                           
		if(mt_->lattice_[index_+dx].xlink_ != nullptr)              
			n_neighbs++;
	}
	else{                                                                 
		if(mt_->lattice_[index_-dx].xlink_ != nullptr)              
			n_neighbs++;                                                  
		if(mt_->lattice_[index_+dx].xlink_ != nullptr)              
			n_neighbs++;                                                  
	}
	return n_neighbs;
}
