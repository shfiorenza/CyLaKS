#include "master_header.h"
#include "kinesin.h"

Kinesin::Kinesin(){
}

void Kinesin::Initialize(system_parameters *parameters, 
	system_properties *properties, int ID){

	ID_ = ID; 
	parameters_ = parameters;
	properties_ = properties;
}
