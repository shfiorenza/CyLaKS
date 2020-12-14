#ifndef _CYLAKS_LINEAR_SPRING_HPP_
#define _CYLAKS_LINEAR_SPRING_HPP_
#include "object.hpp"

class LinearSpring : public Object {
protected:
  double r_cutoff_{0.0};
  double k_spring_{0.0};
  Object *parent_;

public:
  Vec<double> force_;

private:
public:
  LinearSpring() {}
  void Initialize(size_t sid, size_t id, Object *parent) {
    Object::Initialize(sid, id);
    parent_ = parent;
  }
  void SetSpringConstant(double k_spring) { k_spring_ = k_spring; }
  void SetForce();
};

#endif