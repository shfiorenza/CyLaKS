#ifndef _RNG_MANAGEMENT_H
#define _RNG_MANAGEMENT_H
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <vector>

struct system_parameters;
struct system_properties;

class RandomNumberManagement {
private:
  const gsl_rng_type *generator_type_ = gsl_rng_mt19937;

public:
  gsl_rng *rng_;
  std::vector<gsl_rng *> kinesin_rngs_;
  std::vector<gsl_rng *> crosslinker_rngs_;

  system_parameters *parameters_ = nullptr;
  system_properties *properties_ = nullptr;

private:
public:
  RandomNumberManagement();
  void Initialize(system_parameters *parameters, system_properties *properties);
  void CleanUp();

  int GetRanInt(int n);
  void SetRanIndices(int indices[], int n, int m);
  double GetRanProb();
  double GetGaussianNoise(double sigma);

  void Shuffle(void *array, size_t length, size_t element_size);

  int SampleBinomialDist(double p, int n);
  int SamplePoissonDist(double n_avg);
};
#endif
