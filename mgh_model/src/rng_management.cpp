#include "master_header.h"
#include "rng_management.h"

RandomNumberManagement::RandomNumberManagement(){
}

void RandomNumberManagement::Initialize(system_parameters *parameters, 
		system_properties *properties){

	parameters_ = parameters;
	properties_ = properties;
	long seed = parameters->seed; 
	// Mersenne Twister up in this bitch 
	generator_type_ = gsl_rng_mt19937;
	rng = gsl_rng_alloc(generator_type_);
	gsl_rng_set(rng, seed);

	int n_types = properties->kinesin4.serial_pop_.size();
	rngs.resize(n_types); 
	for(int i_type(0); i_type < n_types; i_type++){
		rngs[i_type] = gsl_rng_alloc(generator_type_);
		long local_seed = seed + i_type*1000;
		gsl_rng_set(rngs[i_type], local_seed);
	}
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

double RandomNumberManagement::GetGaussianNoise(double sigma){

	double noise = gsl_ran_gaussian(rng, sigma);
	return noise; 
}

int RandomNumberManagement::SampleNormalDist(double sigma){

	double p = 0.00001;
	double n = sigma*sigma/(p*(1 - p));
	double avg = p*n;
	int sample = SampleBinomialDist(p, n);
	int result = sample - avg;
	return result;
}

int RandomNumberManagement::SampleAbsNormalDist(double sigma){

	/*  recall, for a binomial distribution:

		mean = n * p, N is no. of trials, p is probability for success
		sigma^2 = n * p * (1 - p) 
	
	This is intended to sample the normal distribution around 0 */

	double p = 0.0001;
	double n = sigma*sigma/(p*(1 - p));
	double avg = p*n;
	int sample = SampleBinomialDist(p, n);
	int result;
	if(sample > avg)
		result = sample - avg;
	else
		result = avg - sample;
	return result;
}

int RandomNumberManagement::SampleBinomialDist(double p, int n){
	return gsl_ran_binomial(rng, p, n);
}

int RandomNumberManagement::SampleBinomialDist(double p, int n, int i_event){
	return gsl_ran_binomial(rngs[i_event], p, n);
}

int RandomNumberManagement::SamplePoissonDist(double n_avg){
	return gsl_ran_poisson(rng, n_avg);
}
