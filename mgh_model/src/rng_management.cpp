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
//	printf("sample: %i, avg: %g\n", sample, avg);
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

/*
int RandomNumberManagement::SampleNormalDist(double sigma, int center){

	// Same as above, but centered around some input value 

	double p = 0.0001;
	double n = sigma*sigma/(p*(1-p));
	double avg = p*n;
	int sample = SampleBinomialDist(p, n);
	int result;
	if(sample > avg)
		result = sample - avg + center;
	else
		result = avg - sample + center;
	return result;
}
*/

int RandomNumberManagement::SampleBinomialDist(double p, int n){
	return gsl_ran_binomial(rng, p, n);
}

int RandomNumberManagement::SamplePoissonDist(double n_avg){
	return gsl_ran_poisson(rng, n_avg);
}
