#ifndef _SYSTEM_PROPERTIES
#define _SYSTEM_PROPERTIES
#include "associated_protein_management.h"
#include "curator.h"
#include "kinesin_management.h"
#include "microtubule_management.h"
#include "rng_management.h"

struct system_properties {

  Curator wallace;
  RandomNumberManagement gsl;
  MicrotubuleManagement microtubules;
  KinesinManagement kinesin4;
  AssociatedProteinManagement prc1;

  size_t current_step_{0};

  bool sim_equilibrating_{true};
  bool sim_running_{true};

  FILE *log_file_, *occupancy_file_, *motor_ID_file_, *xlink_ID_file_,
      *tether_coord_file_, *mt_coord_file_, *motor_extension_file_,
      *xlink_extension_file_, *motor_force_file_, *xlink_force_file_,
      *total_force_file_, *motor_head_status_file_;
};
#endif
