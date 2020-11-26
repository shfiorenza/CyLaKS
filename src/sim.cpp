#include "curator.h"

int main(int argc, char *argv[]) {

  Curator wallace(argv);
  while (wallace.sim_running_) {
    wallace.EvolveSimulation();
  }
  return 0;
}