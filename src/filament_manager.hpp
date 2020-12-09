#ifndef _CYLAKS_FILAMENT_MANAGER_HPP_
#define _CYLAKS_FILAMENT_MANAGER_HPP_
#include "population.hpp"
#include "protofilament.hpp"

class BindingSite;
class Curator;
struct SysParams;
struct SysRNG;

class FilamentManager {
private:
  bool up_to_date_{false};

  bool immobile_{true};

  Vec<BindingSite *> sites_;

  Curator *wally_{nullptr};
  SysParams *params_{nullptr};
  SysRNG *gsl_{nullptr};

public:
  bool mobile_{false};
  Vec<Protofilament> list_;

  UMap<Str, Population<BindingSite>> unocc_;

private:
  void GenerateFilaments();

public:
  FilamentManager() {}
  void Initialize(Curator *wally, SysParams *params) {
    wally_ = wally;
    params_ = params;
    GenerateFilaments();
  }

  void FlagForUpdate() { up_to_date_ = false; }
  void UpdateUnoccupied();

  void RunBD();
};
#endif