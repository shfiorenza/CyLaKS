#ifndef _MICROTUBULE_PARAMETERS_H
#define _MICROTUBULE_PARAMETERS_H
#include <vector>

struct microtubule_parameters {
private:
  using vec_t = std::vector<double>;

public:
  int count;               // Number of MTs in simulation
  std::vector<int> length; // Length of each MT in sites (tubulin dimers)
  double y_dist;           // Vertical distance between each MT
  double site_size;        // Length of tubulin dimer; nm
  double radius;           // Outer radius of MT barrel; nm
  double elevation;        // Distance above glass slide; nm
  std::vector<double> start_coord;    // Coords of each MT at sim start; nm
  std::vector<double> immobile_until; // MTs immobilzed until this time; s
  int n_iterations;     // Number of iterations for each diffusion step
  double applied_force; // Perpetually applied force on each MT (pN)
  bool printout_on;     // Whether or not ASCII printout is enabled
  bool diffusion_on;    // Whether or not MT diffusion is enabled
};
#endif
