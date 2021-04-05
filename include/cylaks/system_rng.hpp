#ifndef _CYLAKS_SYSTEM_RNG_HPP_
#define _CYLAKS_SYSTEM_RNG_HPP_
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

// SysRNG: A wrapper for GSL so that we dont have to keep passing the RNG object
struct SysRNG {
private:
  inline static const gsl_rng_type *generator_type_{gsl_rng_mt19937};
  inline static gsl_rng *rng_;

public:
  SysRNG() {}
  ~SysRNG() { gsl_rng_free(rng_); }
  static void Initialize(int seed) {
    rng_ = gsl_rng_alloc(generator_type_);
    gsl_rng_set(rng_, seed);
  }
  // NOTE: These functions could be wrapped to simply return 0 if
  // given 0 as an input, but this can mask more fundamental errors
  // in the simulation, so the seg-faults from GSL should be observed
  static int SampleBinomial(double p, int n) {
    return gsl_ran_binomial(rng_, p, n);
  }
  static int SamplePoisson(double n_avg) {
    return gsl_ran_poisson(rng_, n_avg);
  }
  static int GetRanInt(int n) { return gsl_rng_uniform_int(rng_, n); }
  static double GetRanProb() { return gsl_rng_uniform(rng_); }
  static double GetGaussianPDF(double x, double sigma) {
    return gsl_ran_gaussian_pdf(x, sigma);
  }
  static double GetGaussianNoise(double sigma) {
    return gsl_ran_gaussian(rng_, sigma);
  }
  static void Shuffle(void *array, int length, int element_size) {
    gsl_ran_shuffle(rng_, array, length, element_size);
  }
  static void SetRanIndices(int indices[], int n, int m) {
    int integer_pool[m];
    for (int i{0}; i < m; i++) {
      integer_pool[i] = i;
    }
    gsl_ran_choose(rng_, indices, n, integer_pool, m, sizeof(int));
  }
};
#endif
