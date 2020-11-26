#ifndef _CYLAKS_RNG_MANAGER_HPP_
#define _CYLAKS_RNG_MANAGER_HPP_
#include "definitions.hpp"

class RandomNumberManager {
private:
  const gsl_rng_type *generator_type_{gsl_rng_mt19937};
  gsl_rng *rng_;

public:
private:
public:
  RandomNumberManagement();
  void Initialize(system_parameters *parameters);
  void CleanUp();

  int GetRanInt(int n);
  double GetRanProb();
  double GetGaussianPDF(double x, double sigma);

  void SetRanIndices(int indices[], int n, int m);
  void Shuffle(void *array, size_t length, size_t element_size);

  int SampleBinomialDist(double p, int n);
  int SamplePoissonDist(double n_avg);
};
#endif
