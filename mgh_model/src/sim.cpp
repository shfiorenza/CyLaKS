#include "master_header.h"

int main(int argc, char *argv[])
{

  system_parameters parameters;
  system_properties properties;

  // Initialize sim: parse params & make objects (MTs, kinesin, MAPs, etc)
  properties.wallace.InitializeSimulation(argv[0], argv[1], argv[2], argc,
                                          &properties, &parameters);

  // Main KMC loop
  for (int i_step = 0; i_step < parameters.n_steps; i_step++)
  {
    properties.wallace.UpdateTimestep(i_step);
    properties.kinesin4.Run_KMC();
    properties.prc1.Run_KMC();
    properties.microtubules.RunDiffusion();
  }

  properties.wallace.CleanUp();
  return 0;
}
