#ifndef _CYLAKS_SPHERE_H_
#define _CYLAKS_SPHERE_H_
#include "object.hpp"

class Sphere : public Object {
private:
  double radius_{0.0};

public:
  Sphere(Object *parent, double radius)
      : radius_{radius}, Object(parent->GetSpeciesID(), parent->GetID()) {}
};
#endif