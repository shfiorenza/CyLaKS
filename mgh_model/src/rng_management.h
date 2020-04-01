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
  double GetGaussianNoise(double sigma);

  void SetRanIndices(int indices[], int n, int m);
  void Shuffle(void *array, size_t length, size_t element_size);

  int SampleBinomialDist(double p, int n);
  int SamplePoissonDist(double n_avg);
};
#endif
