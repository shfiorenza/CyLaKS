#ifndef _CYLAKS_BINDING_SITE_H_
#define _CYLAKS_BINDING_SITE_H_
#include "sphere.hpp"

class BindingHead;
class Protofilament;

class BindingSite : public Sphere {
protected:
  Vec<BindingSite *> neighbors_;
  double weight_bind_{0.0};
  double weight_unbind_{0.0};

public:
  size_t index_{0};
  BindingHead *occupant_{nullptr};
  Protofilament *filament_{nullptr};

private:
public:
  BindingSite() {}
  void Initialize(size_t sid, size_t id, double radius, size_t index,
                  Protofilament *filament) {
    Sphere::Initialize(sid, id, radius);
    index_ = index;
    filament_ = filament;
  }

  void AddNeighbor(BindingSite *site) { neighbors_.emplace_back(site); }

  bool IsOccupied() {
    if (occupant_ == nullptr) {
      return false;
    }
    return true;
  }

  double GetWeight_Bind() { return weight_bind_; }
  double GetWeight_Unbind() { return weight_unbind_; }

  int GetNumNeighborsOccupied() {
    if (_n_neighbs_max == 0) {
      return 0;
    }
    int n_neighbs{0};
    for (auto const &site : neighbors_) {
      if (site->occupant_ != nullptr) {
        n_neighbs++;
      }
    }
    return n_neighbs;
  }
  BindingSite *GetNeighbor(int dir);
};
#endif