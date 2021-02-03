#include "curator.hpp"

int main(int argc, char *argv[]) {

  Curator wallace(argc, argv);
  while (Sys::running_) {
    wallace.EvolveSimulation();
  }
  return 0;
}