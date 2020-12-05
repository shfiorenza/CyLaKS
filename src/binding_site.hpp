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
  bool occupied_{false};
  BindingHead *occupant_{nullptr};
  Protofilament *filament_{nullptr};

private:
public:
  BindingSite(size_t sid, size_t id, double radius, Protofilament *filament)
      : Sphere(sid, id, radius), filament_{filament} {}

  double GetWeight_Bind() { return weight_bind_; }
  double GetWeight_Unbind() { return weight_unbind_; }

  bool HeadTrailing();

  size_t GetNeighborCount() { return 0; }
};
#endif