#include "rng_management.h"
#include "master_header.h"

RandomNumberManagement::RandomNumberManagement() {}

void RandomNumberManagement::Initialize(system_parameters *parameters,
                                        system_properties *properties) {

  parameters_ = parameters;
  properties_ = properties;

  // initialize "root" RNG
  long seed = parameters->seed;
  rng_ = gsl_rng_alloc(generator_type_);
  gsl_rng_set(rng_, seed);
}

void RandomNumberManagement::CleanUp() { gsl_rng_free(rng_); }

/* NOTE: These functions could be wrapped to simply return 0 if
   given 0 as an input, but this can mask more fundamental errors
   in the simulation, so the segfaults from GSL should be observed */

int RandomNumberManagement::GetRanInt(int n) {

  return gsl_rng_uniform_int(rng_, n);
}

double RandomNumberManagement::GetRanProb() { return gsl_rng_uniform(rng_); }

double RandomNumberManagement::GetGaussianNoise(double sigma) {

  return gsl_ran_gaussian(rng_, sigma);
}

double RandomNumberManagement::GetGaussianPDF(double x, double sigma) {

  return gsl_ran_gaussian_pdf(x, sigma);
}

int RandomNumberManagement::SampleNormalDist(double sigma) {

  double p = 0.00001;
  double n = sigma * sigma / (p * (1 - p));
  double avg = p * n;
  int sample = SampleBinomialDist(p, n);
  int result = sample - avg;
  return result;
}

int RandomNumberManagement::SampleAbsNormalDist(double sigma) {

  /*  recall, for a binomial distribution:

          mean = n * p, N is no. of trials, p is probability for success
          sigma^2 = n * p * (1 - p)

  This is intended to sample the normal distribution around 0 */

  double p = 0.0001;
  double n = sigma * sigma / (p * (1 - p));
  double avg = p * n;
  int sample = SampleBinomialDist(p, n);
  int result;
  if (sample > avg)
    result = sample - avg;
  else
    result = avg - sample;
  return result;
}

int RandomNumberManagement::SampleBinomialDist(double p, int n) {
  return gsl_ran_binomial(rng_, p, n);
}

int RandomNumberManagement::SamplePoissonDist(double n_avg) {
  return gsl_ran_poisson(rng_, n_avg);
}
