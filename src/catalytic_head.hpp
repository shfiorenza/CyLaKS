#ifndef _CYLAKS_CATALYTIC_HEAD_HPP_
#define _CYLAKS_CATALYTIC_HEAD_HPP_
#include "binding_head.hpp"

class BindingSite;
class Motor;

class CatalyticHead : public BindingHead {
private:
public:
  enum Ligand { NONE, ATP, ADPP, ADP };
  Ligand ligand_{ADP};

  bool trailing_{false};

  Motor *parent_{nullptr};
  CatalyticHead *other_head_{nullptr};

private:
public:
  CatalyticHead() {}
  void Initialize(size_t sid, size_t id, double radius, Motor *parent_ptr,
                  CatalyticHead *other_head_ptr) {
    BindingHead::Initialize(sid, id, radius);
    parent_ = parent_ptr;
    other_head_ = other_head_ptr;
  }

  Ligand GetLigand() { return ligand_; }
  CatalyticHead *GetOtherHead() { return other_head_; }

  double GetWeight_Unbind_II();
  bool Unbind();

  bool Trailing() { return trailing_; }
};
#endif