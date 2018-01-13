#ifndef _RNG_MANAGEMENT_H
#define _RNG_MANAGEMENT_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class RandomNumberManagement{
	private:

	public:
		const gsl_rng_type *generator_type_;
		gsl_rng *rng;		
	private:
	
	public:
		RandomNumberManagement();

		void Initialize(int seed);

		int GetRanInt(int n);
		double GetRanProb();

		int SampleNormalDist(double sigma);
		int SampleNormalDist(double sigma, int center); 
		int SampleBinomialDist(double p, int n);
		int SamplePoissonDist(double n_avg);
};
#endif
