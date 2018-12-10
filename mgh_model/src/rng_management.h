#ifndef _RNG_MANAGEMENT_H
#define _RNG_MANAGEMENT_H
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

struct system_parameters;
struct system_properties;

class RandomNumberManagement{
	private:
		const gsl_rng_type *generator_type_ = gsl_rng_mt19937;
	public:
		gsl_rng *rng_;		
		std::vector<gsl_rng*> kinesin_rngs_; 
		std::vector<gsl_rng*> xlink_dif_rngs_;
		std::vector<gsl_rng*> xlink_kmc_rngs_;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
	private:

	public:
		RandomNumberManagement();
		void Initialize(system_parameters *parameters, 
				system_properties *properties);
		void CleanUp();

		int GetRanInt(int n);
		double GetRanProb();
		double GetGaussianNoise(double sigma);

		int SampleNormalDist(double sigma);
		int SampleAbsNormalDist(double sigma);
		int SampleBinomialDist(double p, int n);
		int SampleBinomialDist_Kinesin(double p, int n, int i_rng);
		int SampleBinomialDist_XlinkDif(double p, int n, int i_rng); 
		int SampleBinomialDist_XlinkKMC(double p, int n, int i_rng); 
		int SamplePoissonDist(double n_avg);
		int SamplePoissonDist_Kinesin(double n_avg, int i_rng);
		int SamplePoissonDist_Xlink(double n_avg, int i_rng);
};
#endif
