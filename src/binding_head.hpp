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
  // FIXME sid & id
  BindingHead(Protein *parent, BindingHead *other_head, double radius)
      : Sphere(0, 0, radius), other_head_{other_head} {}
  BindingHead(size_t sid, size_t id, double r) : Sphere(sid, id, r) {}

  virtual ~BindingHead();
  virtual bool Trailing() { return false; }

  Protein *GetParent() { return parent_; }
  BindingHead *GetOtherHead() { return other_head_; }
  BindingSite *GetSite() { return site_; }
};
#endif