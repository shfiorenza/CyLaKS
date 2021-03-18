#ifndef _CYLAKS_SPHERE_HPP_
#define _CYLAKS_SPHERE_HPP_
#include "object.hpp"

class Sphere : public Object {
private:
  double radius_{0.0};

public:
  Sphere() {}
  void Initialize(size_t sid, size_t id, double radius) {
    Object::Initialize(sid, id);
    radius_ = radius;
  }
};
#endif