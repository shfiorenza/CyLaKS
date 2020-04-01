#ifndef _SYSTEM_PROPERTIES
#define _SYSTEM_PROPERTIES
#include "associated_protein_management.h"
#include "curator.h"
#include "kinesin_management.h"
#include "microtubule_management.h"
#include "rng_management.h"

struct system_properties {

  struct system_file {
    std::string name_{"nope"};
    FILE *file_ptr_;
    template <typename DATA_T>
    void WriteData(DATA_T *array, unsigned int count) {
      int n_chars_written{fwrite(array, sizeof(DATA_T), count, file_ptr_)};
    }
    void OpenFile(char *sim_name) {
      std::string file_name{*sim_name + "_" + name_ + ".file"};
      file_ptr_ = fopen(name_.c_str(), "w");
      if (file_ptr_ == nullptr) {
        printf("Cannot open %s\n", name_.c_str());
        exit(1);
      }
    }
    void CloseFile() { fclose(file_ptr_); }
  };

  Curator wallace;
  RandomNumberManagement gsl;
  MicrotubuleManagement microtubules;
  KinesinManagement kinesin4;
  AssociatedProteinManagement prc1;

  unsigned long current_step_{0};

  std::unordered_map<std::string, system_file> files_; // ADD THIS BRUHH

  FILE *occupancy_file_, *motor_ID_file_, *xlink_ID_file_, *tether_coord_file_,
      *mt_coord_file_, *motor_extension_file_, *xlink_extension_file_,
      *motor_force_file_, *xlink_force_file_, *total_force_file_,
      *motor_head_status_file_;
};
#endif
