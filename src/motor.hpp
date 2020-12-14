#ifndef _CYLAKS_MOTOR_HPP_
#define _CYLAKS_MOTOR_HPP_
#include "catalytic_head.hpp"
#include "linear_spring.hpp"
#include "protein.hpp"

class BindingSite;

class Motor : public Protein {
protected:
  Str ligands_{"yuhh yuh"};

public:
  CatalyticHead head_one_, head_two_;
  LinearSpring tether_;

private:
public:
  Motor() {}
  void Initialize(size_t sid, size_t id) {
    Object::Initialize(sid, id);
    head_one_.Initialize(sid, id, _r_motor_head, this, &head_two_);
    head_two_.Initialize(sid, id, _r_motor_head, this, &head_one_);
    tether_.Initialize(sid, id, this);
  }

  BindingSite *GetDockSite();
  CatalyticHead *GetActiveHead();
  void ChangeConformation();
  bool Bind(BindingSite *site, CatalyticHead *head);
  bool Unbind(CatalyticHead *head);
};

#endif