#ifndef _CYLAKS_FILAMENT_MANAGER_HPP_
#define _CYLAKS_FILAMENT_MANAGER_HPP_
#include "curator.hpp"
#include "protofilament.hpp"
#include "reservoir.hpp"

class BindingSite;
struct SysParameters;

class FilamentManager {
private:
  friend Reservoir<BindingSite>;
  bool up_to_date_{false};

  Vec<BindingSite *> sites_;

  Curator *wally_;
  SysParameters *params_;

public:
  UMap<Str, Reservoir<BindingSite>::PopEntry> unocc_;

private:
  void GenerateFilaments();

public:
  FilamentManager();
  void Initialize(Curator *wally) {
    wally_ = wally;
    params_ = &wally_->params_;
    GenerateFilaments();
  }

  void FlagForUpdate() { up_to_date_ = false; }
  void UpdateUnoccupied();

  void RunBD();
};
#endif