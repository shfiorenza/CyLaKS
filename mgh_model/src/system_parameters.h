#ifndef _PARAMETERS_H
#define _PARAMETERS_H
#include "associated_protein_parameters.h"
#include "kinesin_parameters.h"
#include "microtubule_parameters.h"

struct system_parameters {

  // Simulation stuff
  long seed;          /* Random number seed */
  int n_steps;        /* Number of timesteps in simulation */
  int n_datapoints;   /* No. data points to sample during sim run */
  int data_threshold; /* Step after which to begin data recording */
  double delta_t;     /* How much time passes in one time step (s) */

  // Physical constants
  double kbT; /* boltzmann constant * temp, in pN * nm */
  double eta; /* Viscosity of liquid; in (pN/um^2)*s */

  // Parameter structures for various sim objects
  microtubule_parameters microtubules;
  kinesin_parameters motors;
  associated_protein_parameters xlinks;
};
#endif
