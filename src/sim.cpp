#include "curator.hpp"

int main(int argc, char *argv[]) {

  for (int i{0}; i < 9; i++) {
    Sys::i_picked_[i] = 0;
  }

  Curator wallace(argc, argv);
  while (Sys::running_) {
    wallace.EvolveSimulation();
  }
  for (int i{0}; i < 9; i++) {
    printf("%i picked %i times\n", i, Sys::i_picked_[i]);
  }
  return 0;
}