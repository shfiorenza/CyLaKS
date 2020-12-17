#ifndef _CYLAKS_BINDING_HEAD_HPP_
#define _CYLAKS_BINDING_HEAD_HPP_
#include "sphere.hpp"

class BindingSite;
class Protein;

class BindingHead : public Sphere {
protected:
public:
  BindingSite *site_{nullptr};

  Protein *parent_{nullptr};
  BindingHead *other_head_{nullptr};

private:
public:
  BindingHead() {}
  void Initialize(size_t sid, size_t id, double radius, Protein *parent_ptr,
                  BindingHead *other_head_ptr) {
    Sphere::Initialize(sid, id, radius);
    parent_ = parent_ptr;
    other_head_ = other_head_ptr;
  }
  void Initialize(size_t sid, size_t id, double radius) {
    Sphere::Initialize(sid, id, radius);
  }

  virtual bool Trailing() { return false; }
  int GetNeighborCount() {
    printf("head\n");
    return 0;
  }

  bool Unbind();

  Protein *GetParent() { return parent_; }
  BindingHead *GetOtherHead() { return other_head_; }
  BindingSite *GetSite() { return site_; }
};
#endif