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

bool Tubulin::SpringEquilOnSameSide(){

	if(occupied_ == true
	&& xlink_ != nullptr){
		if(xlink_->heads_active_ == 2
		&& xlink_->tethered_ == true){
			Kinesin *motor = xlink_->motor_;
			double site_coord = index_ + mt_->coord_;
			double xlink_rest = xlink_->GetAnchorCoordinate();
			double motor_rest = motor->GetRestLengthCoordinate();
			if((xlink_rest > site_coord && motor_rest > site_coord)
			|| (xlink_rest < site_coord && motor_rest < site_coord)){
				return true;
			}
			else if((xlink_rest > site_coord && motor_rest < site_coord)
				 || (xlink_rest < site_coord && motor_rest > site_coord)){
				return false;
			}
			else{
				printf("BACK THE F UP!! in spring equil on same (tubby)\n");
				exit(1);
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
