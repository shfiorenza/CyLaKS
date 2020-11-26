#ifndef _CYLAKS_BINDING_SITE_H_
#define _CYLAKS_BINDING_SITE_H_
#include "sphere.hpp"

class BindingSite : public Sphere {
protected:
  double weight_bind_{0.0};
  double weight_unbind_{0.0};

public:
  bool occupied_{false};
  Object *occupant_{nullptr};
  Object *filament_{nullptr};

private:
public:
  BindingSite(size_t sid, size_t id, double radius, Object *filament)
      : filament_{filament}, Sphere(sid, id, radius) {}

  double GetWeight_Bind() { return weight_bind_; }
  double GetWeight_Unbind() { return weight_unbind_; }

  size_t GetNeighborCount() {}
};
#endif