#include "master_header.h"
#include "rng_management.h"

RandomNumberManagement::RandomNumberManagement(){
}

void RandomNumberManagement::Initialize(system_parameters *parameters, 
		system_properties *properties){

	parameters_ = parameters;
	properties_ = properties;

	int world_rank; 
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// initialize "root" RNG
	long seed = parameters->seed; 
	rng_ = gsl_rng_alloc(generator_type_);
	gsl_rng_set(rng_, seed);

	// initialize RNG for each kinesin population on each MPI node
	int n_kin_pop = properties->kinesin4.serial_pop_.size();
	kinesin_rngs_.resize(n_kin_pop); 
	long last_seed = seed; 
	for(int i_pop(0); i_pop < n_kin_pop; i_pop++){
		long pop_seed = seed + 1 + world_rank*(n_kin_pop + 1) + i_pop;
		kinesin_rngs_[i_pop] = gsl_rng_alloc(generator_type_);
		gsl_rng_set(kinesin_rngs_[i_pop], pop_seed);
		// Ensure seeds aren't duplicated amongst threads
		if(pop_seed == last_seed){
			printf("Error in seeding kinsin RNGs: duplicate seeds!\n");
			exit(1);
		}
		last_seed = pop_seed; 
	}

	// initialize RNG for each crosslinker KMC event on each MPI node
	// (there will always be more diffusion events than KMC)
	int n_xl_pop = properties->prc1.serial_dif_.size();
	crosslinker_rngs_.resize(n_xl_pop);
	long seed_offset = seed + 1 + world_size*(n_kin_pop + 1) + n_kin_pop;
	for(int i_pop(0); i_pop < n_xl_pop; i_pop++){
		long pop_seed = seed_offset + 1 + world_rank*(n_xl_pop + 1) + i_pop;
		crosslinker_rngs_[i_pop] = gsl_rng_alloc(generator_type_);
		gsl_rng_set(crosslinker_rngs_[i_pop], pop_seed);
		if(pop_seed == last_seed){
			printf("Error in seeding crossinker RNGs: duplicate seeds!\n");
			exit(1);
		}
		last_seed = pop_seed;
	}
}

void RandomNumberManagement::CleanUp(){

	gsl_rng_free(rng_);
	for(int i_rng(0); i_rng < kinesin_rngs_.size(); i_rng++){
		gsl_rng_free(kinesin_rngs_[i_rng]);
	}
	for(int i_rng(0); i_rng < crosslinker_rngs_.size(); i_rng++){
		gsl_rng_free(crosslinker_rngs_[i_rng]);
	}
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

int RandomNumberManagement::SampleBinomialDist_Crosslinker(double p, 
		int n, int i_rng){
	return gsl_ran_binomial(crosslinker_rngs_[i_rng], p, n);
}

int RandomNumberManagement::SamplePoissonDist(double n_avg){
	return gsl_ran_poisson(rng_, n_avg);
}

int RandomNumberManagement::SamplePoissonDist_Kinesin(double n_avg, 
		int i_rng){
	return gsl_ran_poisson(kinesin_rngs_[i_rng], n_avg);
}

int RandomNumberManagement::SamplePoissonDist_Crosslinker(double n_avg,
		int i_rng){
	return gsl_ran_poisson(crosslinker_rngs_[i_rng], n_avg);
}
