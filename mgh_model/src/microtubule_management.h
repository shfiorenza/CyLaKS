#ifndef _MICROTUBULE_MANAGEMENT_
#define _MICROTUBULE_MANAGEMENT_
#include "entry.hpp"
#include "microtubule.h" // Includes <vector> lib as well
class Tubulin;
struct system_parameters;
struct system_properties;

class MicrotubuleManagement {
private:
  template <class DATA_T> using Vec = std::vector<DATA_T>;
  int n_sites_tot_{0};
  int max_neighbs_xlink_{2};
  int max_neighbs_motor_{0}; // {2};
  int n_affinities_{1};      // {11};
  bool lists_up_to_date_{false};
  system_parameters *parameters_{nullptr};
  system_properties *properties_{nullptr};

public:
  Vec<int> n_unocc_xlink_;
  Vec<Vec<int>> n_unocc_motor_;
  Vec<Vec<ENTRY_T>> unocc_xlink_;
  Vec<Vec<Vec<ENTRY_T>>> unocc_motor_;
  Vec<Microtubule> mt_list_;

private:
  void SetParameters();
  void GenerateMicrotubules();

public:
  MicrotubuleManagement();
  void Initialize(system_parameters *parameters, system_properties *properties);

  void FlagForUpdate();
  void UpdateNeighbors();
  void UpdateUnoccupied();

  void RunDiffusion();
};
#endif
