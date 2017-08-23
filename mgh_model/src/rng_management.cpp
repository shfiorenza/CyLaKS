#include "rng_management.h"

RandomNumberManagement::RandomNumberManagement(){
}

void RandomNumberManagement::Initialize(int seed){

	// Mersenne Twister up in this bitch 
	generator_type_ = gsl_rng_mt19937;
	rng = gsl_rng_alloc(generator_type_);
	gsl_rng_set(rng, seed);
}

/* NOTE: These functions could be wrapped to simply return 0 if
   given 0 as an input, but this can mask more fundamental errors
   in the simulation, so the segfaults from GSL should be observed */

int RandomNumberManagement::GetRanInt(int n){
	return gsl_rng_uniform_int(rng, n);
}

double RandomNumberManagement::GetRanProb(){
	return gsl_rng_uniform(rng);
}

int RandomNumberManagement::SampleBinomialDist(double p, int n){
	return gsl_ran_binomial(rng, p, n);
}

int RandomNumberManagement::SamplePoissonDist(double n_avg){
	return gsl_ran_poisson(rng, n_avg);
}
