#ifndef _CYLAKS_BINDING_SITE_H_
#define _CYLAKS_BINDING_SITE_H_
#include "sphere.hpp"

class BindingHead;
class CatalyticHead;
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
  BindingSite(size_t sid, size_t id, double radius, Object *filament)
      : filament_{filament}, Sphere(sid, id, radius) {}

  double GetWeight_Bind() { return weight_bind_; }
  double GetWeight_Unbind() { return weight_unbind_; }

  bool HeadTrailing() {
    if (occupant_ == nullptr) {
      return false;
    }
    if (occupant_->GetSpeciesID() != _id_motor) {
      return false;
    }
    return dynamic_cast<CatalyticHead *>(occupant_)->trailing_;
  }

  size_t GetNeighborCount() {}
};
#endif