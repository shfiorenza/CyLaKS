#include "cylaks/filament_tester.hpp"
#include "cylaks/protein_tester.hpp"

void FilamentTester::Initialize(ProteinTester *proteins) {

  recorded_force_.resize(Sys::n_steps_run_);
  proteins_ = proteins;
  FilamentManager::proteins_ = dynamic_cast<ProteinManager *>(proteins_);
  FilamentManager::SetParameters();
  FilamentManager::GenerateFilaments();
}

void FilamentTester::UpdateForces() {

  FilamentManager::UpdateForces();
  // Add necessary force to get desired sliding velocity
  // (only applies to x dimension, so i = 0 index)
  double f_required{Sys::slide_velocity_ * protofilaments_[1].gamma_[0]};
  if (Sys::constant_velocity_) {
    double f_applied{f_required - protofilaments_[1].force_[0]};
    if (Sys::i_step_ < Sys::i_pause_ or Sys::i_step_ >= Sys::i_resume_) {
      protofilaments_[1].force_[0] = f_required;
      // Record applied force
      if (recorded_force_[Sys::i_step_] == 0) {
        recorded_force_[Sys::i_step_] = f_applied;
      }
    }
  } else {
    protofilaments_[1].force_[0] += f_required;
    f_required = Sys::slide_velocity_ * protofilaments_[0].gamma_[0];
    protofilaments_[0].force_[0] -= f_required;
  }
  // else {
  //   printf("i_step is %zu\n", Sys::i_step_);
  // }
}
