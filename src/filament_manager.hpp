#ifndef _CYLAKS_FILAMENT_MANAGER_HPP_
#define _CYLAKS_FILAMENT_MANAGER_HPP_
#include "protofilament.hpp"
#include "reservoir.hpp"

class BindingSite;
class Curator;
struct SysParameters;

class FilamentManager {
private:
  friend Reservoir<BindingSite>;
  bool up_to_date_{false};

  Vec<BindingSite *> sites_;

  Curator *wally_;
  SysParameters *params_;

public:
  bool mobile_{false};
  Vec<Protofilament> list_;

  UMap<Str, Reservoir<BindingSite>::PopEntry> unocc_;

private:
  void GenerateFilaments();

public:
  FilamentManager();
  void Initialize(Curator *wally, SysParameters *params) {
    wally_ = wally;
    params_ = params;
    GenerateFilaments();
  }

  void FlagForUpdate() { up_to_date_ = false; }
  void UpdateUnoccupied();

  void RunBD();
};
#endif