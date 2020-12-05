#ifndef _CYLAKS_PROTOFILAMENT_HPP_
#define _CYLAKS_PROTOFILAMENT_HPP_
#include "binding_site.hpp"
#include "rigid_rod.hpp"

class Curator;
struct SysParams;

class Protofilament : public RigidRod {
protected:
  Curator *wally_{nullptr};
  SysParams *params_{nullptr};

public:
  int dx_{0}; // Towards plus end
  size_t n_sites_;
  std::vector<BindingSite> sites_;

  BindingSite *plus_end_{nullptr};
  BindingSite *minus_end_{nullptr};
  Protofilament *neighbor_{nullptr};

private:
  void SetParameters();
  void GenerateSites();

public:
  Protofilament(Curator *wally, SysParams *params, size_t sid, size_t id,
                double length)
      : wally_{wally}, params_{params}, RigidRod(sid, id, length) {
    GenerateSites();
  }
};
#endif