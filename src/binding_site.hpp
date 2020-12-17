#ifndef _CYLAKS_BINDING_SITE_H_
#define _CYLAKS_BINDING_SITE_H_
#include "sphere.hpp"

class BindingHead;
class Protofilament;

class BindingSite : public Sphere {
protected:
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

  bool IsOccupied() {
    if (occupant_ == nullptr) {
      return false;
    }
    return true;
  }

  double GetWeight_Bind() { return weight_bind_; }
  double GetWeight_Unbind() { return weight_unbind_; }

  int GetNeighborCount() { return 0; }
};
#endif