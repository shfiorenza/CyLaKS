#include "master_header.h"
#include "rng_management.h"

RandomNumberManagement::RandomNumberManagement(){
}

void RandomNumberManagement::Initialize(system_parameters *parameters, 
		system_properties *properties){

	int world_rank; 
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// initialize "root" RNG
	long seed = parameters->seed; 
	rng_ = gsl_rng_alloc(generator_type_);
	gsl_rng_set(rng_, seed);

	// initialize RNG for each KMC event on each MPI node
	int n_kin_events = properties->kinesin4.serial_kmc_.size();
	kinesin_rngs_.resize(n_kin_events); 
	for(int i_event(0); i_event < n_kin_events; i_event++){
		long event_seed = seed + (world_rank+1)*10000 + (i_event+1)*100;
//		printf("seed for kin_rng #%i is %lu\n", i_event, event_seed);
		kinesin_rngs_[i_event] = gsl_rng_alloc(generator_type_);
		gsl_rng_set(kinesin_rngs_[i_event], event_seed);
	}
	printf("\n");
	parameters_ = parameters;
	properties_ = properties;
}

/* NOTE: These functions could be wrapped to simply return 0 if
   given 0 as an input, but this can mask more fundamental errors
   in the simulation, so the segfaults from GSL should be observed */

int RandomNumberManagement::GetRanInt(int n){
	return gsl_rng_uniform_int(rng_, n);
}

double RandomNumberManagement::GetRanProb(){
	return gsl_rng_uniform(rng_);
}

double RandomNumberManagement::GetGaussianNoise(double sigma){

	double noise = gsl_ran_gaussian(rng_, sigma);
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
	return gsl_ran_binomial(rng_, p, n);
}

int RandomNumberManagement::SampleBinomialDist_Kinesin(double p, 
		int n, int i_rng){
	return gsl_ran_binomial(kinesin_rngs_[i_rng], p, n);
}

int RandomNumberManagement::SampleBinomialDist_Xlink(double p, 
		int n, int i_rng){
	return gsl_ran_binomial(xlink_rngs_[i_rng], p, n);
}

int RandomNumberManagement::SamplePoissonDist(double n_avg){
	return gsl_ran_poisson(rng_, n_avg);
}

int RandomNumberManagement::SamplePoissonDist_Kinesin(double n_avg, 
		int i_rng){
	return gsl_ran_poisson(kinesin_rngs_[i_rng], n_avg);
}
