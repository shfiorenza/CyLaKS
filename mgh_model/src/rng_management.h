#ifndef _RNG_MANAGEMENT_H
#define _RNG_MANAGEMENT_H
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

struct system_parameters;
struct system_properties;

class RandomNumberManagement{
	private:

	public:
		const gsl_rng_type *generator_type_;
		gsl_rng *rng;		
		std::vector<gsl_rng*> rngs; 

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
	private:
	
	public:
		RandomNumberManagement();

		void Initialize(system_parameters *parameters, 
				system_properties *properties);

		int GetRanInt(int n);
		double GetRanProb();
		double GetGaussianNoise(double sigma);

		int SampleNormalDist(double sigma);
		int SampleAbsNormalDist(double sigma);
		int SampleBinomialDist(double p, int n);
		int SampleBinomialDist(double p, int n, int i_event);
		int SamplePoissonDist(double n_avg);
};
#endif
