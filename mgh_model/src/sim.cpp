#include "master_header.h"

int main(int argc, char *argv[]) {

  system_parameters parameters;
  system_properties properties;
  // Initialize sim: parse params & make objects (MTs, kinesin, MAPs, etc)
  properties.wallace.InitializeSimulation(argv, &properties, &parameters);
  // Main KMC loop
  for (unsigned long i_step{0}; i_step < parameters.n_steps; i_step++) {
    properties.wallace.UpdateTimestep(i_step);
    // properties.kinesin4.RunKMC();
    properties.prc1.RunKMC();
    properties.microtubules.RunDiffusion();
  }
  // properties.kinesin4.ReportProbabilities();
  // properties.prc1.ReportProbabilities();
  // Clean up simulation before exiting
  properties.wallace.CleanUp();
  return 0;
}