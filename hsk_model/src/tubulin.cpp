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
	parent = mt;
}
