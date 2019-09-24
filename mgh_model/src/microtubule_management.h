#ifndef _MICROTUBULE_MANAGEMENT_H
#define _MICROTUBULE_MANAGEMENT_H
#include "entry.h"
#include "microtubule.h" // Includes <vector> lib as well
class Tubulin;
struct system_parameters;
struct system_properties;

class MicrotubuleManagement {
private:
  template <class DATA_T> using Vec = std::vector<DATA_T>;

public:
  int n_sites_tot_ = 0;
  int n_unoccupied_ = 0;
  std::vector<int> n_unoccupied_xl_;

  system_parameters *parameters_ = nullptr;
  system_properties *properties_ = nullptr;

  std::vector<Microtubule> mt_list_;
  Vec<ENTRY_T> unoccupied_list_;
  Vec<Vec<ENTRY_T>> unoccupied_list_xl_;

private:
public:
  MicrotubuleManagement();

  void Initialize(system_parameters *parameters, system_properties *properties);

  void SetParameters();
  void GenerateMicrotubules();

  void UnoccupiedCheck(Tubulin *site);
  void UnoccupiedCheck(int i_mt, int i_site);
  void OccupiedCheck(Tubulin *site);
  void OccupiedCheck(int i_mt, int i_site);

  void UpdateNeighbors();
  void UpdateUnoccupied();

  Tubulin *GetUnoccupiedSite();
  Tubulin *GetUnoccupiedSite(int n_neighbs);

  void RunDiffusion();
};
#endif
