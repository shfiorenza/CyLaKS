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
      protofilaments_[1].force_[0] = -f_required;
      // Record applied force
      if (recorded_force_[Sys::i_step_] == 0) {
        recorded_force_[Sys::i_step_] = -f_applied;
      }
    }
  } else {
    /* FORCE-DEP UNBINDING */
    double f_stall{6.0}; // pN
    // Add motor force to 'bottom' microtubule
    size_t n_motors_bot{
        (size_t)std::round(double(protofilaments_[0].sites_.size()) / 20.0)};
    double f_bot{protofilaments_[0].force_[0]};
    double f_per_mot_bot{f_bot / n_motors_bot};
    double v_mot_bot{Sys::slide_velocity_ *
                     (1.0 - std::fabs(f_per_mot_bot / f_stall))};
    protofilaments_[0].force_[0] += v_mot_bot * protofilaments_[0].gamma_[0];
    // Add motor force to 'top' microtubule
    size_t n_motors_top{
        (size_t)std::round(double(protofilaments_[1].sites_.size()) / 20.0)};
    double f_top{protofilaments_[1].force_[0]};
    double f_per_mot_top{f_bot / n_motors_top};
    double v_mot_top{Sys::slide_velocity_ *
                     (1.0 - std::fabs(f_per_mot_top / f_stall))};
    protofilaments_[1].force_[0] -= v_mot_top * protofilaments_[1].gamma_[0];
    /* END FORCE-DEP UNBINDING*/

    // printf("f_per_mot: %g // %g\n", f_per_mot_bot, f_per_mot_top);
    // printf("v_exp: %g\n\n", v_mot_top);

    /* OLD PRIMITIVE FORCE (NON) DEPENDENCE BELOW */
    // protofilaments_[1].force_[0] -= f_required;
    // f_required = Sys::slide_velocity_ * protofilaments_[0].gamma_[0];
    // protofilaments_[0].force_[0] += f_required;
    /* END OLD BOI */
  }
  // else {
  //   printf("i_step is %zu\n", Sys::i_step_);
  // }
}
