#ifndef _SYSTEM_PROPERTIES
#define _SYSTEM_PROPERTIES
#include "curator.h"
#include "rng_management.h"
#include "microtubule_management.h"
#include "kinesin_management.h"

typedef struct system_properties{

	Curator wallace;

	RandomNumberManagement gsl; 

	MicrotubuleManagement microtubules;

	KinesinManagement kinesin4; 

} system_properties;
#endif

