#ifndef _CYLAKS_BINDING_HEAD_HPP_
#define _CYLAKS_BINDING_HEAD_HPP_
#include "sphere.hpp"

class BindingHead : public Sphere {
protected:
public:
  BindingSite *site_{nullptr};

  Protein *parent_{nullptr};
  BindingHead *other_head_{nullptr};

private:
public:
  BindingHead(Object *parent, BindingHead *partner, double radius)
      : other_head_{partner}, Sphere(parent, radius) {}
  virtual ~BindingHead();

  Protein *GetParent() { return parent_; }
  BindingHead *GetOtherHead() { return other_head_; }
  BindingSite *GetSite() { return site_; }
};
#endif