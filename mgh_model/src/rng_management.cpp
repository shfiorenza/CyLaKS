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

void RandomNumberManagement::SetRanIndices(int indices[], int n, int m) {

  int integer_pool[m];
  for (int i{0}; i < m; i++) {
    integer_pool[i] = i;
  }
  gsl_ran_choose(rng_, indices, n, integer_pool, m, sizeof(int));
}

double RandomNumberManagement::GetRanProb() { return gsl_rng_uniform(rng_); }

double RandomNumberManagement::GetGaussianNoise(double sigma) {

  double noise = gsl_ran_gaussian(rng_, sigma);
  return noise;
}

void RandomNumberManagement::Shuffle(void *array, size_t length, size_t size) {

  gsl_ran_shuffle(properties_->gsl.rng_, array, length, size);
}

int RandomNumberManagement::SampleBinomialDist(double p, int n) {
  return gsl_ran_binomial(rng_, p, n);
}

int RandomNumberManagement::SamplePoissonDist(double n_avg) {
  return gsl_ran_poisson(rng_, n_avg);
}